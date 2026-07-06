%% KEEP PLOTS

%% Contour plots for Jeffreys prior and Jacobian
% Figure mapping: the Jeffreys-prior contours are Fig. 3a. The Jacobian
% contours are an alternate diagnostic using the same grid.

clear; clc;

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repoRoot);
plotsDir = fileparts(mfilename('fullpath'));

%% Colormap

colors = [
    [10, 10, 10] / 256
    [180, 20, 0] / 256
    [195, 91, 29] / 256
    [235, 172, 71] / 256
    [245, 245, 245] / 256
];
cmap = interpolateColors(colors, 256);
axisLabelFontSize = 30;
tickLabelFontSize = 26;
textFontName = 'Calibri';
labelFont = ['\fontname{' textFontName '}'];

%% Grid

h_vals = -logspace(0, -3, 1000);
lambda_vals = linspace(0.00001, 3, 1000);
plot_lambda_vals = lambda_vals(lambda_vals > 0.01);

Ns = [20, 100, 1000];

%% Load or compute contour values

forceRecompute = false;
cacheFile = fullfile(plotsDir, 'CountourPlots_cache.mat');
cacheVersion = 1;

useCache = false;

if ~forceRecompute && exist(cacheFile, 'file') == 2
    cached = load(cacheFile);
    useCache = isContourCacheValid(cached, h_vals, lambda_vals, plot_lambda_vals, Ns, cacheVersion);

    if useCache
        disp(['Loading contour values from ', cacheFile])
        F_jeffreys = cached.F_jeffreys;
        F_jacobian = cached.F_jacobian;
    else
        disp('Contour cache does not match current grid. Recomputing.')
    end
end

if ~useCache
    [F_jeffreys, F_jacobian] = computeContourValues(h_vals, lambda_vals, Ns);
    save(cacheFile, 'F_jeffreys', 'F_jacobian', ...
        'h_vals', 'lambda_vals', 'plot_lambda_vals', 'Ns', 'cacheVersion', '-v7.3')
    disp(['Saved contour values to ', cacheFile])
end

%% Plot Jeffreys prior

plotContourSet( ...
    F_jeffreys, h_vals, lambda_vals, plot_lambda_vals, Ns, cmap, ...
    [-7 8], 'Contour of log Jefferies Prior (h,J)', ...
    labelFont, textFontName, axisLabelFontSize, tickLabelFontSize);

%% Plot Jacobian

plotContourSet( ...
    F_jacobian, h_vals, lambda_vals, plot_lambda_vals, Ns, cmap, ...
    [-9 5], 'Contour of log Jacobian (h,J)', ...
    labelFont, textFontName, axisLabelFontSize, tickLabelFontSize);

function plotContourSet(F, h_vals, lambda_vals, plot_lambda_vals, Ns, cmap, c_limits, plot_title, labelFont, textFontName, axisLabelFontSize, tickLabelFontSize)
    x = linspace(1, 3);

    for count = 1:length(Ns)
        figure
        hold on
        contourf(h_vals, plot_lambda_vals, log10(abs(F(lambda_vals > 0.01, :, count))), 15)
        colormap(cmap)
        clim(c_limits)
        colorbar('FontSize', tickLabelFontSize, 'FontName', textFontName)
        xlabel([labelFont 'External field {\ith}'], ...
            'Interpreter', 'tex', ...
            'FontSize', axisLabelFontSize)
        ylabel([labelFont 'Interaction strength {\itJ}'], ...
            'Interpreter', 'tex', ...
            'FontSize', axisLabelFontSize)
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
        ax.FontName = textFontName;
        ax.FontSize = tickLabelFontSize;
        box on
    end
end

function [F_jeffreys, F_jacobian] = computeContourValues(h_vals, lambda_vals, Ns)
    F_jeffreys = zeros(length(lambda_vals), length(h_vals), length(Ns));
    F_jacobian = zeros(length(lambda_vals), length(h_vals), length(Ns));

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
end

function tf = isContourCacheValid(cached, h_vals, lambda_vals, plot_lambda_vals, Ns, cacheVersion)
    requiredFields = {'F_jeffreys', 'F_jacobian', 'h_vals', ...
        'lambda_vals', 'plot_lambda_vals', 'Ns', 'cacheVersion'};

    tf = all(isfield(cached, requiredFields));

    if ~tf
        return
    end

    tf = cached.cacheVersion == cacheVersion ...
        && isequal(cached.h_vals, h_vals) ...
        && isequal(cached.lambda_vals, lambda_vals) ...
        && isequal(cached.plot_lambda_vals, plot_lambda_vals) ...
        && isequal(cached.Ns, Ns);
end

function cmap = interpolateColors(colors, nColors)
    x = linspace(0, 1, size(colors, 1));
    xq = linspace(0, 1, nColors);
    cmap = zeros(nColors, 3);

    for channel = 1:3
        cmap(:, channel) = interp1(x, colors(:, channel), xq, 'linear');
    end
end
