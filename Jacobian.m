function [Jac, J11, J12, J21, J22, dmdn, dchidn] = Jacobian(h, lambda, n)
%MUCHIEXACT2SPIN  Exact finite-N magnetization and susceptibility for the
% fully-connected 2-spin (p=2) mean-field Ising model.
%
% Model (with beta absorbed into h, lambda):
%   H(s) = - n * [ h * m + (lambda/2) * m^2 ],
%   m = (1/n) * sum_i s_i,   s_i in {+1,-1}.
%
% Partition function written exactly as a sum over k = #(+1) spins:
%   m_k = (2k - n)/n
%   weight_k = C(n,k) * exp( n * [ h m_k + (lambda/2) m_k^2 ] )
%
% Outputs:
%   mu  = <m>
%   chi = d< m >/dh = n * Var(m)    (beta=1 convention / beta absorbed)

    arguments
        h      (1,1) double
        lambda (1,1) double
        n      (1,1) double {mustBePositive}
    end

    % Enumerate k = number of +1 spins
    k  = (0:n).';
    mk = (2*k - n) ./ n;  % magnetization values

    % log binomial coefficient: log C(n,k)
    logC = gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1);

    % "energy" exponent (beta absorbed): n*(h m + (lambda/6)m^3)
    % (Use elementwise ops; mk is (n+1)x1)
    logw = logC + n .* (h .* mk + (lambda/2) .* (mk.^2));

    % logZ via log-sum-exp for stability
    logZ = logsumexp(logw);

    % normalized weights
    wk = exp(logw - logZ);

    F = -h.*mk - lambda.*mk.^2 ./2 - logC./n;
    dmkdn = -(mk + 1)./n;

    dfdn = (gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1))./n.^2 + (-psi(n+1) + psi(k+1).*((1+mk)/2 + n.*dmkdn/2)+ psi(n-k+1).*((1-mk)/2 -n.*dmkdn/2))/n - h.*dmkdn - lambda.*mk.*dmkdn;

    % moments
    mu  = sum(wk .* mk);
    m2  = sum(wk .* (mk.^2));
    m3  = sum(wk .* (mk.^3));
    m4  = sum(wk .* (mk.^4));
    f = sum(wk .* F);
    mf = sum(wk .* F .* mk);
    df = sum(wk.*dfdn);
    mdf = sum(wk .* mk.*dfdn);
    m2f = sum(wk .* mk.^2 .*F);
    m2df = sum(wk .* mk.^2.*dfdn);

    % susceptibility with respect to h (beta absorbed): chi = n * Var(m)
    J11 = n.*(m2-mu.^2);
    J12 = n.*(0.5*m3-0.5*mu*m2);
    J21 = n.^2 .*(m3 - 3*mu*m2 + 2*mu^3 );
    J22 = n.^2 .*(0.5*m4 - 0.5*m2.^2 - mu.*m3 + mu.^2 .* m2);

    dmdn = - n.*mdf + n.*mu.*df - mf + mu.*f;

    dchidn = J11/n + n*(m2*f +  n*m2*df - m2f- n*m2df - 2*mu^2 * f - 2*mu^2 * n*df + 2*mu*mf + 2*mu*n*mdf);

    Jac = J11*J22 - J12*J21;
end

% ---- helper: stable log(sum(exp(a))) for real a ----
function s = logsumexp(a)
    amax = max(a);
    s = amax + log(sum(exp(a - amax)));
end

% ---- validators ----
function mustBeInteger(x)
    if any(abs(x - round(x)) > 0)
        error("n must be an integer.");
    end
end

% function [Jac, J11, J12, J21, J22, dmdn] = Jacobian(h, lambda, n)
% 
%     arguments
%         h      (1,1) double
%         lambda (1,1) double
%         n      (1,1) double {mustBeInteger, mustBePositive}
%     end
% 
%     % Enumerate k = number of +1 spins
%     k  = (0:n).';
% 
%     % Magnetization values
%     mk = (2*k - n) ./ n;
% 
%     % Log binomial coefficient
%     logC = gammaln(n+1) ...
%          - gammaln(k+1) ...
%          - gammaln(n-k+1);
% 
%     % Model weight
%     logw = logC + n .* (h .* mk + (lambda/2) .* mk.^2);
% 
%     % Normalize weights stably
%     logZ = logsumexp(logw);
%     wk = exp(logw - logZ);
% 
%     % Moments
%     mu = sum(wk .* mk);
%     m2 = sum(wk .* mk.^2);
%     m3 = sum(wk .* mk.^3);
%     m4 = sum(wk .* mk.^4);
% 
%     % Jacobian terms
%     J11 = n .* (m2 - mu.^2);
% 
%     J12 = n .* (0.5*m3 - 0.5*mu*m2);
% 
%     J21 = n.^2 .* (m3 - 3*mu*m2 + 2*mu.^3);
% 
%     J22 = n.^2 .* (0.5*m4 - 0.5*m2.^2 - mu.*m3 + mu.^2.*m2);
% 
%     % ------------------------------------------------------------
%     % Partial derivative d mu / d n
%     % treating k as fixed in the finite sum.
%     %
%     % mk = (2k - n)/n = 2k/n - 1
%     % dmkdn = -2k/n^2 = -(mk + 1)/n
%     % ------------------------------------------------------------
% 
%     dmkdn = -(mk + 1) ./ n;
% 
%     dlogCdn = psi(n+1) - psi(n-k+1);
% 
%     dlogwdn = dlogCdn ...
%             + (h .* mk + (lambda/2) .* mk.^2) ...
%             + n .* (h + lambda .* mk) .* dmkdn;
% 
%     dmdn = sum(wk .* dmkdn) ...
%          + sum(wk .* mk .* dlogwdn) ...
%          - mu .* sum(wk .* dlogwdn);
% 
%     % Determinant of Jacobian
%     Jac = J11*J22 - J12*J21;
% 
% end
% 
% 
% % ---- helper: stable log(sum(exp(a))) for real a ----
% function s = logsumexp(a)
%     amax = max(a);
%     s = amax + log(sum(exp(a - amax)));
% end
% 
% 
% % ---- validators ----
% function mustBeInteger(x)
%     if any(abs(x - round(x)) > 0)
%         error("n must be an integer.");
%     end
% end