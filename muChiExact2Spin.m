function [mu, chi, dudl, ddzdldl, C, wk] = muChiExact2Spin(h, lambda, n)
%MUCHIEXACT3SPIN  Exact finite-N magnetization and susceptibility for the
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

    % moments
    mu  = sum(wk .* mk);
    m2  = sum(wk .* (mk.^2));
    m3  = sum(wk .* (mk.^3));
    m4  = sum(wk .* (mk.^4));
    E = -n .* (h .* mk + (lambda/2) .* (mk.^2));

    % susceptibility with respect to h (beta absorbed): chi = n * Var(m)
    chi = n .* (m2 - mu.^2);

    dudl = n .* (m3 -mu.*m2)/2; 
    ddzdldl = (n / 4)*(m4 - m2.^2);
    %C = (n^2 * lambda^2 / 4) *(m4 - m2^2) + h*lambda*n^2 *(m3-m2*mu) + h^2*n^2 *(m2 - mu^2 ); 
    C = sum(wk.*E.^2) - sum(wk.*E).^2;
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