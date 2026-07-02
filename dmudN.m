function u = dmudN(h, lambda, n)

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
    else
        zarg  = @(psi) h + 1i*psi;       % imaginary direction
        quad  = @(psi)  (psi.^2) ./ (2*lambda);  % lambda<0 => negative real quadratic
        tanh_ = @(psi) tanh(h + 1i*psi);
    end

    % Exponent (generally complex if lambda<0)
    g = @(psi) n .* ( quad(psi) + logcosh(zarg(psi)) );

    F = @(psi) -( quad(psi) + logcosh(zarg(psi)) );

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
    Ef = integ(@(psi) f(psi) .* F(psi))/ Zcore;
    Eftanh = integ(@(psi) f(psi) .* F(psi) .* tanh_(psi)) / Zcore;

    u = Ef.*Etanh - Eftanh;  % may be complex for lambda<0; caller can take real(u) if desired
end

% ---- helper ----
function y = stable_logcosh_complex(z, log2)
    m = abs(real(z));  % max(real(z), real(-z))
    y = m + log( exp(z - m) + exp(-z - m) ) - log2;
end