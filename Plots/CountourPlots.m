%% KEEP PLOTS

%% Contour plots for Jeffreys prior and Jacobian

clear; clc;

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repoRoot);

%% Colormap

colors = [
    [255, 216, 118] / 256
    [235, 172, 71] / 256
    [195, 91, 29] / 256
    [180, 20, 0] / 256
];
cmap = interpolateColors(colors, 256);

%% Grid

h_vals = -logspace(0, -3, 1000);
lambda_vals = linspace(0.00001, 3, 1000);
plot_lambda_vals = lambda_vals(lambda_vals > 0.01);

Ns = [20, 100, 1000];

F_jeffreys = zeros(length(lambda_vals), length(h_vals), length(Ns));
F_jacobian = zeros(length(lambda_vals), length(h_vals), length(Ns));

%% Compute contour values

for count = 1:length(Ns)
    N = Ns(count);
    disp(N)

    for i = 1:length(lambda_vals)
        l = lambda_vals(i);
        disp(l)

        for j = 1:length(h_vals)
            h = h_vals(j);

            [~, m11, m12, m22] = muChiExact2Spin(h, l, N);

            F_jeffreys(i, j, count) = N^2 * (m11 * m22 - m12^2);
            F_jacobian(i, j, count) = Jacobian(h, l, N);
        end
    end
end

%% Plot Jeffreys prior

plotContourSet( ...
    F_jeffreys, h_vals, lambda_vals, plot_lambda_vals, Ns, cmap, ...
    [-7 8], 'Contour of log Jefferies Prior (h,\lambda)');

%% Plot Jacobian

plotContourSet( ...
    F_jacobian, h_vals, lambda_vals, plot_lambda_vals, Ns, cmap, ...
    [-6 5], 'Contour of log Jacobian (h,\lambda)');

function plotContourSet(F, h_vals, lambda_vals, plot_lambda_vals, Ns, cmap, c_limits, plot_title)
    x = linspace(1, 3);

    for count = 1:length(Ns)
        figure
        hold on
        contourf(h_vals, plot_lambda_vals, log10(abs(F(lambda_vals > 0.01, :, count))), 15)
        colormap(cmap)
        clim(c_limits)
        colorbar('FontSize', 26)
        xlabel('h', 'FontSize', 30)
        ylabel('\lambda', 'FontSize', 30)
        title(sprintf('%s, N = %d', plot_title, Ns(count)))
        axis square
        scatter(0, 1, 'or', 'filled')
        plot(0 * x, x, 'r')

        xlim([-1, -0.0005])
        xscale log
        xticks([-1, -0.1, -0.01, -0.001, -0.0005])
        xticklabels({'-10^0', '-10^{-1}', '-10^{-2}', '-10^{-3}', '0'})

        set(gca, 'TickDir', 'both')
        ax = gca;
        ax.FontSize = 26;
        box on
    end
end

function cmap = interpolateColors(colors, nColors)
    x = linspace(0, 1, size(colors, 1));
    xq = linspace(0, 1, nColors);
    cmap = zeros(nColors, 3);

    for channel = 1:3
        cmap(:, channel) = interp1(x, colors(:, channel), xq, 'linear');
    end
end
