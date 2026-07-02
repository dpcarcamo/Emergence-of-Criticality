function x = dchidlam(h, lambda, n)
    arguments
        h      (1,1) double
        lambda (1,1) double {mustBePositive}
        n      (1,1) double {mustBePositive}
    end

    % Numerically stable log(cosh(x))
    % Avoids overflow of cosh(x) when |x| is large.
    log2 = log(2);
    logcosh = @(x) stable_logcosh(x, log2);

    % Log of (unnormalized) density: g(psi)
    % g(psi) = -n * ( psi^2 / (2*lambda) - logcosh(h + psi) )
    g = @(psi) -n .* ( psi.^2 ./ (2*lambda) - logcosh(h + psi) );

    % Find approximate mode of g(psi) to shift the exponent
    % Width of the peak is O(sqrt(lambda/n)), use that to set a search range.
    R   = max(10, 10*sqrt(lambda/n) + abs(h));  % fairly generous
    obj = @(x) -g(x);                           % minimize -g -> maximize g
    [psi0, objval] = fminbnd(obj, -R, R);
    g0  = -objval;                              % maximum value of g

    % Shifted log-density so that max is at 0
    raw_logf = @(psi) g(psi) - g0;

    % Clamp tiny positive roundoff to 0 to *guarantee* no overflow in exp
    logf = @(psi) min(raw_logf(psi), 0);

    f = @(psi) exp(logf(psi));                  % now always in (0, 1]

    tanh_ = @(psi) tanh(h + psi);
    sech2 = @(psi) 1 ./ cosh(h + psi).^2;

    integ = @(fun) integral( ...
        fun, -Inf, Inf, ...
        'RelTol',     1e-10, ...
        'AbsTol',     1e-12, ...
        'ArrayValued', true);

    % Normalization factor; exp(g0) cancels algebraically.
    Zcore = integ(f);
    Etanh = integ(@(psi) f(psi) .* tanh_(psi)) / Zcore;

    Etanh   = integ(@(psi) f(psi).*tanh_(psi)) / Zcore;
    Etanh2  = integ(@(psi) f(psi).*(tanh_(psi).^2)) / Zcore;
    Etanh3  = integ(@(psi) f(psi).*(tanh_(psi).^3)) / Zcore;
    Etanh4  = integ(@(psi) f(psi).*(tanh_(psi).^4)) / Zcore;


    x = (n*(3*n-2) + 2.*n.*(n-1).*(3.*n-4).*Etanh2 + n*(n-1).*(n-2).*(n-3).*Etanh4  - (n.*(n-1).*Etanh2 + n).^2 - 2.*(n.*Etanh).*(n.*(n-1).*(n-2).*Etanh3 + 3.*n.*(n-1).*Etanh + n.*Etanh) + 2.*(n.*(n-1).*Etanh2 + n).*(n.*Etanh).^2)./(2.*n.^2);
end

% --- helper (subfunction) ---

function y = stable_logcosh(x, log2)
    % Compute log(cosh(x)) without overflowing cosh(x)
    %
    % For |x| small, just use log(cosh(x)).
    % For |x| large, use:
    %   log(cosh(x)) = |x| + log( (1 + exp(-2|x|)) / 2 )
    %                = |x| + log1p(exp(-2|x|)) - log(2)
    ax = abs(x);
    y  = zeros(size(x));

    % Threshold can be tuned; 20 is a very safe choice
    small = ax < 20;

    if any(small, 'all')
        xs       = x(small);
        y(small) = log(cosh(xs));
    end

    if any(~small, 'all')
        xl        = ax(~small);
        y(~small) = xl + log1p(exp(-2*xl)) - log2;
    end
end
