%% KEEP PLOTS


%% RG vector fields for fixed N values
% Standalone fixed-N vector-field plot.
% Figure mapping: this script generates the Fig. 3b vector-field panels.

clearvars
clc

currentPath = pwd;
folders = split(currentPath, '\');
newPath = join(folders(1:length(folders)-1),"\");
addpath(newPath{1});
plotsDir = fileparts(mfilename('fullpath'));

%% Plot settings

NValues = [20 100 1000];

hlim = [-0.9 -0.001];
lambdalim = [0.001 3];

Nh = 15;
Nlambda = 15;

arrowScale = 0.21;
arrowLengthMode = 'logMagnitude';  % Options: 'fixed' or 'logMagnitude'
arrowHeadSize = 2.8;
arrowHeadBaseLength = 0.018;
arrowHeadBaseWidth = 0.014;
plotPanelsSeparately = false;
logMagnitudeReference = 0;
logMagnitudeFloor = -8;
axisLabelFontSize = 30;
tickLabelFontSize = 26;
colorbarFontSize = 30;
textFontName = 'Helvetica';
labelFont = ['\fontname{' textFontName '}'];
criticalLineWidth = 2;
criticalPointSize = 85;
criticalColor = 'r';
forceRecompute = false;
cacheFile = fullfile(plotsDir, 'RGvectorFieldFixedN_cache.mat');
cacheVersion = 1;

colors = [
    [0.25 0.25 0.25]
    [0.9 0 0]
];
cmap = interpolateColors(colors, 256);

%% Grid

hValues = -logspace(log10(-hlim(2)), log10(-hlim(1)), Nh);
lambdaValues = linspace(lambdalim(1), lambdalim(2), Nlambda);
[hGrid, lambdaGrid] = meshgrid(hValues, lambdaValues);

plotH = -log10(-hGrid);
xMax = max(plotH, [], 'all') - log10(0.5);
criticalXPlot = xMax;
xPlotMax = criticalXPlot;

%% Load or compute vector fields

useCache = false;

if ~forceRecompute && exist(cacheFile, 'file') == 2
    cached = load(cacheFile);
    useCache = isVectorFieldCacheValid(cached, hGrid, lambdaGrid, NValues, cacheVersion);

    if useCache
        disp(['Loading vector fields from ', cacheFile])
        dhdNFields = cached.dhdNFields;
        dlambdadNFields = cached.dlambdadNFields;
    else
        disp('Vector-field cache does not match current grid. Recomputing.')
    end
end

if ~useCache
    [dhdNFields, dlambdadNFields] = computeVectorFields(hGrid, lambdaGrid, NValues);
    save(cacheFile, 'dhdNFields', 'dlambdadNFields', ...
        'hGrid', 'lambdaGrid', 'NValues', 'cacheVersion', '-v7.3')
    disp(['Saved vector fields to ', cacheFile])
end

%% Plot vector fields

if ~plotPanelsSeparately
    figure('Name', 'RG vector fields for fixed N', 'Color', 'w');
    tiledlayout(1, numel(NValues), ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');
end

for idx = 1:numel(NValues)
    N = NValues(idx);

    dhdN = dhdNFields(:, :, idx);
    dlambdadN = dlambdadNFields(:, :, idx);

    plotDh = dhdN./(abs(hGrid)*log(10));
    plotDlambda = dlambdadN;

    if plotPanelsSeparately
        figure('Name', sprintf('RG vector field N = %g', N), 'Color', 'w');
    else
        nexttile
    end

    coloredVectorField( ...
        plotH, lambdaGrid, plotDh, plotDlambda, ...
        'NormalizeArrows', true, ...
        'MagnitudeMode', 'original', ...
        'ArrowLengthMode', arrowLengthMode, ...
        'ArrowScale', arrowScale, ...
        'LogMagnitudeReference', logMagnitudeReference, ...
        'LogMagnitudeFloor', logMagnitudeFloor, ...
        'ColorMap', cmap, ...
        'MaxHeadSize', arrowHeadSize, ...
        'UseCustomArrowHeads', true, ...
        'ArrowHeadSize', arrowHeadSize, ...
        'ArrowHeadBaseLength', arrowHeadBaseLength, ...
        'ArrowHeadBaseWidth', arrowHeadBaseWidth);

    xlim([0, xPlotMax])
    ylim(lambdalim)
    xticks([0 1 2 3 criticalXPlot])
    xticklabels({'10^0', '10^{-1}', '10^{-2}', '10^{-3}', '0'})

    xlabel([labelFont 'External field {\ith}'], ...
        'Interpreter', 'tex', ...
        'FontSize', axisLabelFontSize)
    if plotPanelsSeparately || idx == 1
        ylabel([labelFont 'Interaction strength {\itJ}'], ...
            'Interpreter', 'tex', ...
            'FontSize', axisLabelFontSize)
    else
        ylabel('')
    end
    %title(sprintf('N = %g', N))

    ax = gca;
    ax.FontName = textFontName;
    ax.FontSize = tickLabelFontSize;
    ax.LineWidth = 1.3;
    ax.TickDir = 'out';
    ax.Layer = 'bottom';
    box on
    axis square

    hold on
    plot([criticalXPlot, criticalXPlot], [1, lambdalim(2)], ...
        'Color', criticalColor, ...
        'LineWidth', criticalLineWidth, ...
        'Clipping', 'off')
    scatter(criticalXPlot, 1, criticalPointSize, ...
        'MarkerFaceColor', criticalColor, ...
        'MarkerEdgeColor', criticalColor, ...
        'Clipping', 'off')
end

colormap(gcf, cmap)
plotStandaloneLogColorbar(cmap, [logMagnitudeFloor logMagnitudeReference], ...
    'Vector magnitude', colorbarFontSize);

%% Local functions

function [dhdNFields, dlambdadNFields] = computeVectorFields(hGrid, lambdaGrid, NValues)
    dhdNFields = zeros([size(hGrid), numel(NValues)]);
    dlambdadNFields = zeros([size(lambdaGrid), numel(NValues)]);

    for idx = 1:numel(NValues)
        N = NValues(idx);
        fprintf('Evaluating vector field for N = %g\n', N);

        dhdN = zeros(size(hGrid));
        dlambdadN = zeros(size(lambdaGrid));

        for ii = 1:numel(hGrid)
            h = hGrid(ii);
            lambda = lambdaGrid(ii);

            dhdN(ii) = real(DhDNtot(h, lambda, N));
            dlambdadN(ii) = real(DlamDNtot(h, lambda, N));
        end

        dhdNFields(:, :, idx) = dhdN;
        dlambdadNFields(:, :, idx) = dlambdadN;
    end
end

function tf = isVectorFieldCacheValid(cached, hGrid, lambdaGrid, NValues, cacheVersion)
    requiredFields = {'dhdNFields', 'dlambdadNFields', ...
        'hGrid', 'lambdaGrid', 'NValues', 'cacheVersion'};

    tf = all(isfield(cached, requiredFields));

    if ~tf
        return
    end

    tf = cached.cacheVersion == cacheVersion ...
        && isequal(cached.hGrid, hGrid) ...
        && isequal(cached.lambdaGrid, lambdaGrid) ...
        && isequal(cached.NValues, NValues);
end

function hq = coloredVectorField(X, Y, U, V, varargin)
    p = inputParser;

    addParameter(p, 'ColorByMagnitude', true);
    addParameter(p, 'ColorMap', turbo(256));
    addParameter(p, 'NormalizeArrows', false);
    addParameter(p, 'ArrowScale', 0.6);
    addParameter(p, 'LineWidth', 2);
    addParameter(p, 'MaxHeadSize', 0.7);
    addParameter(p, 'MagnitudeMode', 'original');
    addParameter(p, 'ArrowLengthMode', 'fixed');
    addParameter(p, 'LogMagnitudeReference', 0);
    addParameter(p, 'LogMagnitudeFloor', -8);
    addParameter(p, 'UseCustomArrowHeads', false);
    addParameter(p, 'ArrowHeadSize', 1);
    addParameter(p, 'ArrowHeadBaseLength', 0.018);
    addParameter(p, 'ArrowHeadBaseWidth', 0.014);
    addParameter(p, 'Skip', 1);

    parse(p, varargin{:});
    opts = p.Results;

    skip = opts.Skip;
    X = X(1:skip:end, 1:skip:end);
    Y = Y(1:skip:end, 1:skip:end);
    U = U(1:skip:end, 1:skip:end);
    V = V(1:skip:end, 1:skip:end);

    Uorig = U;
    Vorig = V;

    switch lower(opts.MagnitudeMode)
        case 'original'
            M = sqrt(Uorig.^2 + Vorig.^2);
        case 'plotted'
            M = sqrt(U.^2 + V.^2);
        otherwise
            error('MagnitudeMode must be either original or plotted.')
    end

    logM = log10(M);
    quiverScale = opts.ArrowScale;

    switch lower(opts.ArrowLengthMode)
        case 'fixed'
            if opts.NormalizeArrows
                [U, V] = normalizeVectors(U, V);
            end

        case 'logmagnitude'
            [U, V] = normalizeVectors(U, V);

            if opts.LogMagnitudeReference <= opts.LogMagnitudeFloor
                error('LogMagnitudeReference must be larger than LogMagnitudeFloor.')
            end

            clippedLogM = min(max(logM, opts.LogMagnitudeFloor), opts.LogMagnitudeReference);
            relativeLength = (clippedLogM - opts.LogMagnitudeFloor) ./ ...
                (opts.LogMagnitudeReference - opts.LogMagnitudeFloor);

            U = U .* opts.ArrowScale .* relativeLength;
            V = V .* opts.ArrowScale .* relativeLength;
            quiverScale = 0;

        otherwise
            error('ArrowLengthMode must be either fixed or logMagnitude.')
    end

    if opts.UseCustomArrowHeads && quiverScale ~= 0
        U = U .* quiverScale;
        V = V .* quiverScale;
        quiverScale = 0;
    end

    Xf = X(:);
    Yf = Y(:);
    Uf = U(:);
    Vf = V(:);
    Mf = logM(:);

    valid = ~(isnan(Xf) | isnan(Yf) | isnan(Uf) | isnan(Vf) | isnan(Mf) | ...
              isinf(Xf) | isinf(Yf) | isinf(Uf) | isinf(Vf) | isinf(Mf));

    Xf = Xf(valid);
    Yf = Yf(valid);
    Uf = Uf(valid);
    Vf = Vf(valid);
    Mf = Mf(valid);

    if opts.LogMagnitudeReference <= opts.LogMagnitudeFloor
        error('LogMagnitudeReference must be larger than LogMagnitudeFloor.')
    end

    cmap = opts.ColorMap;
    nColors = size(cmap, 1);

    clippedMf = min(max(Mf, opts.LogMagnitudeFloor), opts.LogMagnitudeReference);
    colorInds = round(1 + (nColors - 1) * ...
        (clippedMf - opts.LogMagnitudeFloor) ./ ...
        (opts.LogMagnitudeReference - opts.LogMagnitudeFloor));
    colorInds = max(1, min(nColors, colorInds));

    holdState = ishold;
    hold on

    hq = gobjects(length(Xf), 1);

    for k = 1:length(Xf)
        if opts.ColorByMagnitude
            arrowColor = cmap(colorInds(k), :);
        else
            arrowColor = [0 0 0];
        end

        if opts.UseCustomArrowHeads
            hq(k) = line([Xf(k), Xf(k) + Uf(k)], ...
                [Yf(k), Yf(k) + Vf(k)], ...
                'Color', arrowColor, ...
                'LineWidth', opts.LineWidth);
            drawArrowHead(Xf(k), Yf(k), Uf(k), Vf(k), arrowColor, ...
                opts.LineWidth, opts.ArrowHeadSize, ...
                opts.ArrowHeadBaseLength, opts.ArrowHeadBaseWidth);
        else
            hq(k) = quiver(Xf(k), Yf(k), Uf(k), Vf(k), quiverScale, ...
                'Color', arrowColor, ...
                'LineWidth', opts.LineWidth, ...
                'MaxHeadSize', opts.MaxHeadSize);
        end
    end

    if opts.ColorByMagnitude
        colormap(cmap)
        clim([opts.LogMagnitudeFloor opts.LogMagnitudeReference])
    end

    axis tight
    box on

    if ~holdState
        hold off
    end
end

function drawArrowHead(x, y, u, v, arrowColor, lineWidth, arrowHeadSize, baseLength, baseWidth)
    arrowLength = sqrt(u^2 + v^2);

    if arrowLength <= 0
        return
    end

    direction = [u, v] ./ arrowLength;
    normal = [-direction(2), direction(1)];

    headLength = min(baseLength * arrowHeadSize, 0.8 * arrowLength);
    headWidth = baseWidth * arrowHeadSize;

    tip = [x + u, y + v];
    baseCenter = tip - headLength * direction;

    xPatch = [
        tip(1)
        baseCenter(1) + 0.5 * headWidth * normal(1)
        baseCenter(1) - 0.5 * headWidth * normal(1)
    ];
    yPatch = [
        tip(2)
        baseCenter(2) + 0.5 * headWidth * normal(2)
        baseCenter(2) - 0.5 * headWidth * normal(2)
    ];

    patch(xPatch, yPatch, arrowColor, ...
        'EdgeColor', arrowColor, ...
        'LineWidth', 0.5 * lineWidth, ...
        'Clipping', 'on');
end

function plotStandaloneLogColorbar(cmap, c_limits, colorbarLabel, tickLabelFontSize)
    figure('Name', [colorbarLabel, ' colorbar'], 'Color', 'w')
    ax = axes;
    colormap(ax, cmap)
    clim(ax, c_limits)
    axis(ax, 'off')

    cb = colorbar(ax);
    cb.FontSize = tickLabelFontSize;
    cb.Label.String = colorbarLabel;
    cb.Label.FontSize = tickLabelFontSize;

    formatLogColorbar(cb, c_limits, tickLabelFontSize)
end

function formatLogColorbar(cb, c_limits, tickLabelFontSize)
    ticks = logColorbarTicks(c_limits);

    cb.Ticks = ticks;
    cb.TickLabels = arrayfun(@(x) sprintf('10^{%g}', x), ticks, 'UniformOutput', false);
    cb.TickLabelInterpreter = 'tex';
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

function [Uout, Vout] = normalizeVectors(U, V)
    L = sqrt(U.^2 + V.^2);
    mask = L > 0;

    Uout = zeros(size(U));
    Vout = zeros(size(V));

    Uout(mask) = U(mask)./L(mask);
    Vout(mask) = V(mask)./L(mask);
end

function cmap = interpolateColors(colors, nColors)
    nAnchors = size(colors, 1);
    xAnchors = linspace(0, 1, nAnchors);
    xInterp = linspace(0, 1, nColors);
    cmap = interp1(xAnchors, colors, xInterp, 'linear');
    cmap = max(min(cmap, 1), 0);
end
