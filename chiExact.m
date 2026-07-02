% function x = chiExact(h, lambda, n)
%     arguments
%         h      (1,1) double
%         lambda (1,1) double {mustBePositive}
%         n      (1,1) double {mustBePositive}
%     end
% 
%     % Numerically stable log(cosh(x))
%     % Avoids overflow of cosh(x) when |x| is large.
%     log2 = log(2);
%     logcosh = @(x) stable_logcosh(x, log2);
% 
%     % Log of (unnormalized) density: g(psi)
%     % g(psi) = -n * ( psi^2 / (2*lambda) - logcosh(h + psi) )
%     g = @(psi) -n .* ( psi.^2 ./ (2*lambda) - logcosh(h + psi) );
% 
%     % Find approximate mode of g(psi) to shift the exponent
%     % Width of the peak is O(sqrt(lambda/n)), use that to set a search range.
%     R   = max(10, 10*sqrt(lambda/n) + abs(h));  % fairly generous
%     obj = @(x) -g(x);                           % minimize -g -> maximize g
%     [psi0, objval] = fminbnd(obj, -R, R);
%     g0  = -objval;                              % maximum value of g
% 
%     % Shifted log-density so that max is at 0
%     raw_logf = @(psi) g(psi) - g0;
% 
%     % Clamp tiny positive roundoff to 0 to *guarantee* no overflow in exp
%     logf = @(psi) min(raw_logf(psi), 0);
% 
%     f = @(psi) exp(logf(psi));                  % now always in (0, 1]
% 
%     tanh_ = @(psi) tanh(h + psi);
%     sech2 = @(psi) 1 ./ cosh(h + psi).^2;
% 
%     integ = @(fun) integral( ...
%         fun, -Inf, Inf, ...
%         'RelTol',     1e-10, ...
%         'AbsTol',     1e-12, ...
%         'ArrayValued', true);
% 
%     % Normalization factor; exp(g0) cancels algebraically.
%     Zcore = integ(f);
%     Etanh = integ(@(psi) f(psi) .* tanh_(psi)) / Zcore;
% 
%     Etanh   = integ(@(psi) f(psi).*tanh_(psi)) / Zcore;
%     Etanh2  = integ(@(psi) f(psi).*(tanh_(psi).^2)) / Zcore;
%     Esech2  = integ(@(psi) f(psi).*sech2(psi)) / Zcore;
% 
%     x = -(- Esech2 - n*(Etanh2 - Etanh.^2));
% end
% 
% % --- helper (subfunction) ---
% 
% function y = stable_logcosh(x, log2)
%     % Compute log(cosh(x)) without overflowing cosh(x)
%     %
%     % For |x| small, just use log(cosh(x)).
%     % For |x| large, use:
%     %   log(cosh(x)) = |x| + log( (1 + exp(-2|x|)) / 2 )
%     %                = |x| + log1p(exp(-2|x|)) - log(2)
%     ax = abs(x);
%     y  = zeros(size(x));
% 
%     % Threshold can be tuned; 20 is a very safe choice
%     small = ax < 20;
% 
%     if any(small, 'all')
%         xs       = x(small);
%         y(small) = log(cosh(xs));
%     end
% 
%     if any(~small, 'all')
%         xl        = ax(~small);
%         y(~small) = xl + log1p(exp(-2*xl)) - log2;
%     end
% end


function x = chiExact(h, lambda, n)

    arguments
        h      (1,1) double
        lambda (1,1) double
        n      (1,1) double {mustBePositive}
    end

    if lambda == 0
        error("lambda must be nonzero.");
    end

    log2 = log(2);

    % Stable log(cosh(z)) for complex z = a + i b.
    % Uses complex log-sum-exp:
    % log(cosh z) = log(exp(z) + exp(-z)) - log(2)
    %             = m + log(exp(z-m) + exp(-z-m)) - log(2),
    % where m = max(real(z), real(-z)) = abs(real(z)).
    logcosh = @(z) stable_logcosh_complex(z, log2);

    % Choose the contour / integrand based on sign of lambda
    if lambda > 0
        zarg  = @(psi) h + psi;          % real line
        quad  = @(psi) -(psi.^2) ./ (2*lambda);
        tanh_ = @(psi) tanh(h + psi);
        sech2 = @(psi) 1 ./ cosh(h + psi).^2;
    else
        zarg  = @(psi) h + 1i*psi;       % imaginary direction
        quad  = @(psi)  (psi.^2) ./ (2*lambda);  % lambda<0 => negative real quadratic
        tanh_ = @(psi) tanh(h + 1i*psi);
        sech2 = @(psi) 1 ./ cosh(h + 1i*psi).^2;
    end

    % Exponent (generally complex if lambda<0)
    g = @(psi) n .* ( quad(psi) + logcosh(zarg(psi)) );

    % Find a shift g0 so that real(g - g0) <= 0 (prevents overflow, keeps |exp|<=1).
    % The integrand is essentially Gaussian in psi with width ~ sqrt(|lambda|/n).
    R   = max(10, 10*sqrt(abs(lambda)/n) + abs(h));
    obj = @(x) -real(g(x));                      % maximize real(g)
    [~, objval] = fminbnd(obj, -R, R);
    g0 = -objval;                                % approx max of real(g)

    raw_logf = @(psi) g(psi) - g0;
    % Clamp any tiny positive real roundoff so exp never exceeds magnitude 1
    logf = @(psi) raw_logf(psi) - max(real(raw_logf(psi)), 0);

    f = @(psi) exp(logf(psi));                   % complex-valued allowed

    integ = @(fun) integral(fun, -Inf, Inf, ...
        'RelTol',      1e-14, ...
        'AbsTol',      1e-14, ...
        'ArrayValued', true);

    Zcore = integ(f);
    Etanh = integ(@(psi) f(psi) .* tanh_(psi)) / Zcore;
    Etanh2  = integ(@(psi) f(psi).*(tanh_(psi).^2)) / Zcore;
    Esech2  = integ(@(psi) f(psi).*sech2(psi)) / Zcore;

    x = -(- Esech2 - n*(Etanh2 - Etanh.^2));  % may be complex for lambda<0; caller can take real(u) if desired
end

% ---- helper ----
function y = stable_logcosh_complex(z, log2)
    m = abs(real(z));  % max(real(z), real(-z))
    y = m + log( exp(z - m) + exp(-z - m) ) - log2;
end
