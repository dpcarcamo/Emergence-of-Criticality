%% KEEP PLOTS
% Used for creating the cartoon double well plot

l = 1.5;
h = 0.01;
N = 20;

figure
hold on

x = linspace(-1,1);


Ns = [20,50, 100, 200, 1000];
labels = strings(size(Ns));  % store legend entries

for k = 1:length(Ns)
    N = Ns(k);

    [h, l] = hlambda(-0.962,0.24,N);
    
    y = (N*h*x + N*l*x.^2/2 ...
        + logGammaLanczos(N+1) ...
        - logGammaLanczos(N*(1+x)/2 + 1) ...
        - logGammaLanczos(N*(1-x)/2 + 1)) / N;

    % Find peaks and take top two
    y_padded = [min(y)-1, y, min(y)-1];
    [pks, locs_padded] = findpeaks(y_padded);
    %[pks, ~] = findpeaks(y);
    top = maxk(pks, 2);

delta = diff(top); %#ok<NASGU> % retained for optional spacing diagnostics

    plot(x, -(y - top(end)), 'LineWidth', 1.5)

    % Store label
    labels(k) = sprintf('N = %d', N);
end

plot(x, x.*0, 'k--', 'DisplayName', 'y = 0')
box on
ax = gca;
ax.FontSize = 20;      % tick labels

ylabel('f(\mu)', 'FontSize', 18)
xlabel('\mu', 'FontSize', 18)
axis square

legend(labels, 'Location', 'best')

%%






function y = logGammaLanczos(x)
%LOGGAMMALANCZOS Vectorized log(Gamma(x)) using Lanczos approximation
% Works for real x (avoids poles at non-positive integers)

    % Coefficients (g=7, n=9)
    p = [ ...
        0.99999999999980993
        676.5203681218851
       -1259.1392167224028
        771.32342877765313
       -176.61502916214059
        12.507343278686905
       -0.13857109526572012
        9.9843695780195716e-6
        1.5056327351493116e-7];

    g = 7;

    % Initialize output
    y = zeros(size(x));

    % ---- Reflection region: x < 0.5 ----
    mask_reflect = (x < 0.5);
    xr = x(mask_reflect);

    if any(mask_reflect)
        y(mask_reflect) = log(pi) ...
            - log(abs(sin(pi * xr))) ...
            - logGammaLanczos(1 - xr);
    end

    % ---- Main Lanczos region: x >= 0.5 ----
    mask_main = ~mask_reflect;
    xm = x(mask_main);

    if any(mask_main)
        xm = xm - 1;

        % Compute a = p(1) + sum p(i)/(x+i-1)
        a = p(1) * ones(size(xm));
        for i = 2:length(p)
            a = a + p(i) ./ (xm + i - 1);
        end

        t = xm + g + 0.5;

        y(mask_main) = 0.5*log(2*pi) ...
            + (xm + 0.5).*log(t) ...
            - t ...
            + log(a);
    end
end
