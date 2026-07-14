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
    [0.25 0.25 0.25]
    [0.9 0 0]
];
cmap = interpolateColors(colors, 256);
axisLabelFontSize = 30;
tickLabelFontSize = 26;
textFontName = 'Helvetica';
labelFont = ['\fontname{' textFontName '}'];
plotPanelsSeparately = false;
criticalHPlot = -0.0005;
criticalLineWidth = 2;
criticalPointSize = 85;
criticalColor = 'r';

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
    labelFont, textFontName, axisLabelFontSize, tickLabelFontSize, ...
    plotPanelsSeparately, criticalHPlot, criticalLineWidth, criticalPointSize, criticalColor);
plotStandaloneLogColorbar(cmap, [-7 8], ...
    'Jeffreys prior', textFontName, tickLabelFontSize);

%% Plot Jacobian

plotContourSet( ...
    F_jacobian, h_vals, lambda_vals, plot_lambda_vals, Ns, cmap, ...
    [-9 6], 'Contour of log Jacobian (h,J)', ...
    labelFont, textFontName, axisLabelFontSize, tickLabelFontSize, ...
    plotPanelsSeparately, criticalHPlot, criticalLineWidth, criticalPointSize, criticalColor);
plotStandaloneLogColorbar(cmap, [-9 6], ...
    'Jacobian determinant', textFontName, tickLabelFontSize);

function plotContourSet(F, h_vals, lambda_vals, plot_lambda_vals, Ns, cmap, c_limits, plot_title, labelFont, textFontName, axisLabelFontSize, tickLabelFontSize, plotPanelsSeparately, criticalHPlot, criticalLineWidth, criticalPointSize, criticalColor)
    if ~plotPanelsSeparately
        figure('Name', plot_title, 'Color', 'w')
        tiledlayout(1, length(Ns), ...
            'TileSpacing', 'compact', ...
            'Padding', 'compact');
    end

    for count = 1:length(Ns)
        if plotPanelsSeparately
            figure('Name', sprintf('%s, N = %d', plot_title, Ns(count)), 'Color', 'w')
        else
            nexttile
        end

        hold on
        contourf(h_vals, plot_lambda_vals, log10(abs(F(lambda_vals > 0.01, :, count))), 15)
        colormap(cmap)
        clim(c_limits)
        %cb = colorbar('FontSize', tickLabelFontSize, 'FontName', textFontName);
        %formatLogColorbar(cb, c_limits, textFontName, tickLabelFontSize)
        xlabel([labelFont 'External field {\ith}'], ...
            'Interpreter', 'tex', ...
            'FontSize', axisLabelFontSize)
        if plotPanelsSeparately || count == 1
            ylabel([labelFont 'Interaction strength {\itJ}'], ...
                'Interpreter', 'tex', ...
                'FontSize', axisLabelFontSize)
        else
            ylabel('')
        end
        %title(sprintf('%s, N = %d', plot_title, Ns(count)))
        axis square

        xlim([-1, criticalHPlot])
        xscale log
        xticks([-1, -0.1, -0.01, -0.001, criticalHPlot])
        xticklabels({'-10^0', '-10^{-1}', '-10^{-2}', '-10^{-3}', '0'})

        set(gca, 'TickDir', 'both')
        ax = gca;
        ax.FontName = textFontName;
        ax.FontSize = tickLabelFontSize;
        box on

        plot([criticalHPlot, criticalHPlot], [1, max(plot_lambda_vals)], ...
            'Color', criticalColor, ...
            'LineWidth', criticalLineWidth, ...
            'Clipping', 'off')
        scatter(criticalHPlot, 1, criticalPointSize, ...
            'MarkerFaceColor', criticalColor, ...
            'MarkerEdgeColor', criticalColor, ...
            'Clipping', 'off')
    end
end

function plotStandaloneLogColorbar(cmap, c_limits, colorbarLabel, textFontName, tickLabelFontSize)
    figure('Name', [colorbarLabel, ' colorbar'], 'Color', 'w')
    ax = axes;
    colormap(ax, cmap)
    clim(ax, c_limits)
    axis(ax, 'off')

    cb = colorbar(ax);
    cb.FontSize = tickLabelFontSize;
    cb.FontName = textFontName;
    cb.Label.String = colorbarLabel;
    cb.Label.FontSize = tickLabelFontSize;
    cb.Label.FontName = textFontName;

    formatLogColorbar(cb, c_limits, textFontName, tickLabelFontSize)
end

function formatLogColorbar(cb, c_limits, textFontName, tickLabelFontSize)
    ticks = logColorbarTicks(c_limits);

    cb.Ticks = ticks;
    cb.TickLabels = arrayfun(@(x) sprintf('10^{%g}', x), ticks, 'UniformOutput', false);
    cb.TickLabelInterpreter = 'tex';
    cb.FontName = textFontName;
    cb.FontSize = tickLabelFontSize;
end

function ticks = logColorbarTicks(c_limits)
    minTick = ceil(c_limits(1));
    maxTick = floor(c_limits(2));
    allTicks = minTick:maxTick;

    if numel(allTicks) <= 6
        ticks = allTicks;
    else
        step = ceil((maxTick - minTick) / 5);
        ticks = minTick:step:maxTick;
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
