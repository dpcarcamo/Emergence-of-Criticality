%% KEEP PLOTS

%% Combined contour, vector-field, and dot-plot panels
% Row 1: Jacobian contours in (h,lambda)
% Row 2: RG vector fields in (h,lambda)
% Row 3: mapped dot plots in (h,lambda)

clearvars
clc

plotsDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(plotsDir);
addpath(repoRoot);

%% Shared plot settings

Ns = [20, 100, 1000];

hlim = [-0.9 -0.001];
lambdalim = [0.001 3];
lambdaPlotLim = [0 3];
lambdaTicks = 0:1:3;

axisLabelFontSize = 15;%30;
tickLabelFontSize = 10;%26;
colorbarFontSize = 30;
textFontName = 'Helvetica';
labelFont = ['\fontname{' textFontName '}'];

criticalHPlot = -0.0005;
criticalLineWidth = 2;%4;
criticalPointSize = 75;%105;
criticalColor = 'r';

vectorLineWidth = 1;%2;
dotSize = 15;%45;
edgeWidth = 3;%5;
insetFontScale = 1.8;
insetAxisLabelFontSize = round(axisLabelFontSize * insetFontScale);
insetTickLabelFontSize = round(tickLabelFontSize * insetFontScale);

colors = [
    1.0000    1.0000    1.0000
    0.8800    0.8400    1.0000
    0.4500    0.5200    1.0000
    0.0400    0.1500    0.6500
         0         0         0
];

xAnchors = linspace(0, 1, size(colors, 1));
xInterp = linspace(0, 1, 256);
cmap = interp1(xAnchors, colors, xInterp, 'linear');
cmap = max(min(cmap, 1), 0);

%% Contour grid and cache

h_vals = -logspace(0, -3, 1000);
lambda_vals = linspace(0.00001, 3, 1000);
plot_lambda_vals = lambda_vals(lambda_vals > 0.01);

contourCLimits = [-9 6];
forceRecomputeContour = false;
contourCacheVersion = 1;
contourCacheFile = fullfile(plotsDir, 'Figure3CombinedPanels_contour_cache.mat');
originalContourCacheFile = fullfile(plotsDir, 'CountourPlots_cache.mat');

useContourCache = false;
contourCacheCandidates = {contourCacheFile, originalContourCacheFile};

for cacheIdx = 1:numel(contourCacheCandidates)
    currentContourCacheFile = contourCacheCandidates{cacheIdx};

    if forceRecomputeContour || exist(currentContourCacheFile, 'file') ~= 2
        continue
    end

    cached = load(currentContourCacheFile);
    requiredContourFields = {'F_jacobian', 'h_vals', 'lambda_vals', ...
        'plot_lambda_vals', 'Ns', 'cacheVersion'};

    if all(isfield(cached, requiredContourFields)) ...
            && cached.cacheVersion == contourCacheVersion ...
            && isequal(cached.h_vals, h_vals) ...
            && isequal(cached.lambda_vals, lambda_vals) ...
            && isequal(cached.plot_lambda_vals, plot_lambda_vals) ...
            && isequal(cached.Ns, Ns)
        disp(['Loading Jacobian contour values from ', currentContourCacheFile])
        F_jacobian = cached.F_jacobian;
        useContourCache = true;
        break
    end
end

if ~useContourCache
    F_jacobian = zeros(length(lambda_vals), length(h_vals), length(Ns));

    for count = 1:length(Ns)
        N = Ns(count);
        disp(['Computing Jacobian contours for N = ', num2str(N)])

        for i = 1:length(lambda_vals)
            l = lambda_vals(i);
            disp(l)

            for j = 1:length(h_vals)
                h = h_vals(j);
                F_jacobian(i, j, count) = Jacobian(h, l, N);
            end
        end
    end

    cacheVersion = contourCacheVersion;
    save(contourCacheFile, 'F_jacobian', ...
        'h_vals', 'lambda_vals', 'plot_lambda_vals', 'Ns', 'cacheVersion', '-v7.3')
    disp(['Saved Jacobian contour values to ', contourCacheFile])
end

%% Vector-field grid and cache

NValues = Ns;
Nh = 15;
Nlambda = 15;

arrowScale = 0.21;
arrowLengthMode = 'logMagnitude';  % Options: 'fixed' or 'logMagnitude'
arrowHeadSize = 2.8;
arrowHeadBaseLength = 0.018;
arrowHeadBaseWidth = 0.014;
logMagnitudeReference = 0;
logMagnitudeFloor = -5;

hValues = -logspace(log10(-hlim(2)), log10(-hlim(1)), Nh);
lambdaValues = linspace(lambdalim(1), lambdalim(2), Nlambda);
[hGrid, lambdaGrid] = meshgrid(hValues, lambdaValues);

plotH = -log10(-hGrid);
xMax = max(plotH, [], 'all') - log10(0.5);
criticalXPlot = xMax;
xPlotMax = criticalXPlot;

forceRecomputeVectorField = false;
vectorFieldCacheVersion = 1;
vectorFieldCacheFile = fullfile(plotsDir, 'RGvectorFieldFixedN_cache.mat');

useVectorFieldCache = false;

if ~forceRecomputeVectorField && exist(vectorFieldCacheFile, 'file') == 2
    cached = load(vectorFieldCacheFile);
    requiredVectorFields = {'dhdNFields', 'dlambdadNFields', ...
        'hGrid', 'lambdaGrid', 'NValues', 'cacheVersion'};

    useVectorFieldCache = all(isfield(cached, requiredVectorFields)) ...
        && cached.cacheVersion == vectorFieldCacheVersion ...
        && isequal(cached.hGrid, hGrid) ...
        && isequal(cached.lambdaGrid, lambdaGrid) ...
        && isequal(cached.NValues, NValues);

    if useVectorFieldCache
        disp(['Loading vector fields from ', vectorFieldCacheFile])
        dhdNFields = cached.dhdNFields;
        dlambdadNFields = cached.dlambdadNFields;
    else
        disp('Vector-field cache does not match current grid. Recomputing.')
    end
end

if ~useVectorFieldCache
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

    cacheVersion = vectorFieldCacheVersion;
    save(vectorFieldCacheFile, 'dhdNFields', 'dlambdadNFields', ...
        'hGrid', 'lambdaGrid', 'NValues', 'cacheVersion', '-v7.3')
    disp(['Saved vector fields to ', vectorFieldCacheFile])
end

%% Dot-plot grid and cache

mus  = linspace(-0.9,-0.02,24);
chis = linspace(1,3.5,24);

edgePts = 50;

muMin = min(mus);
muMax = max(mus);
chiMin = min(chis);
chiMax = max(chis);

muBottom  = linspace(muMin, muMax, edgePts);
chiBottom = chiMin * ones(size(muBottom));

chiRight = linspace(chiMin, chiMax, edgePts);
muRight  = muMax * ones(size(chiRight));

muTop  = linspace(muMax, muMin, edgePts);
chiTop = chiMax * ones(size(muTop));

chiLeft = linspace(chiMax, chiMin, edgePts);
muLeft  = muMin * ones(size(chiLeft));

muEdges  = {muBottom,  muRight,  muTop,  muLeft};
chiEdges = {chiBottom, chiRight, chiTop, chiLeft};

edgeColors = {
    [0.00 0.55 0.60 0.75], ... % bottom edge
    [0.85 0.05 0.45 0.75], ... % right edge
    [0.00 0.45 0.08 0.75], ... % top edge
    [0.00 0.10 0.70 0.75]  ... % left edge
};

forceRecomputeDotplot = false;
dotplotCacheVersion = 1;
dotplotCacheFile = fullfile(plotsDir, 'dotplots_cache.mat');

useDotplotCache = false;

if ~forceRecomputeDotplot && exist(dotplotCacheFile, 'file') == 2
    cached = load(dotplotCacheFile);
    requiredDotplotFields = {'hs', 'ls', 'hEdges', 'lEdges', 'Ns', 'mus', ...
        'chis', 'edgePts', 'muEdges', 'chiEdges', 'cacheVersion'};

    useDotplotCache = all(isfield(cached, requiredDotplotFields)) ...
        && cached.cacheVersion == dotplotCacheVersion ...
        && isequal(cached.Ns, Ns) ...
        && isequal(cached.mus, mus) ...
        && isequal(cached.chis, chis) ...
        && isequal(cached.edgePts, edgePts) ...
        && isequal(cached.muEdges, muEdges) ...
        && isequal(cached.chiEdges, chiEdges);

    if useDotplotCache
        disp(['Loading dotplot data from ', dotplotCacheFile])
        hs = cached.hs;
        ls = cached.ls;
        hEdges = cached.hEdges;
        lEdges = cached.lEdges;
    else
        disp('Dotplot cache does not match current grid. Recomputing.')
    end
end

if ~useDotplotCache
    hs = zeros(length(mus), length(chis), length(Ns));
    ls = zeros(length(mus), length(chis), length(Ns));
    hEdges = cell(length(Ns),4);
    lEdges = cell(length(Ns),4);

    for a = 1:length(Ns)
        N = Ns(a);
        fprintf('Computing dotplot map for N = %d\n', N)

        for j = 1:length(mus)
            for k = 1:length(chis)
                fprintf('N index %d, mu index %d, chi index %d\n', a, j, k)

                mu = mus(j);
                chi = chis(k);

                [h, lambda] = hlambda(mu, chi, N);

                hs(j,k,a) = h;
                ls(j,k,a) = lambda;
            end
        end

        for e = 1:4
            hEdge = zeros(size(muEdges{e}));
            lEdge = zeros(size(muEdges{e}));

            for q = 1:length(muEdges{e})
                fprintf('N index %d, edge %d, point %d\n', a, e, q)
                [hTemp, lTemp] = hlambda(muEdges{e}(q), chiEdges{e}(q), N);
                hEdge(q) = hTemp;
                lEdge(q) = lTemp;
            end

            hEdges{a,e} = hEdge;
            lEdges{a,e} = lEdge;
        end
    end

    cacheVersion = dotplotCacheVersion;
    save(dotplotCacheFile, 'hs', 'ls', 'hEdges', 'lEdges', ...
        'Ns', 'mus', 'chis', 'edgePts', 'muEdges', 'chiEdges', ...
        'cacheVersion', '-v7.3')
    disp(['Saved dotplot data to ', dotplotCacheFile])
end

%% Helper figure: original box in (m,chi)

figure('Name', 'Original box in statistics space', 'Color', 'w')
hold on
set(gca, 'TickDir', 'both')
box on

for e = 1:4
    plot(muEdges{e}, chiEdges{e}, ...
        'Color', edgeColors{e}, ...
        'LineWidth', edgeWidth);
end

mus_in  = mus(2:end-1);
chis_in = chis(2:end-1);

[MU, CHI] = meshgrid(mus_in, chis_in);
scatter(MU(:), CHI(:), dotSize, 'k', 'filled');

xlabel([labelFont 'Average activity {\itm}'], ...
    'Interpreter', 'tex', ...
    'FontSize', insetAxisLabelFontSize)
ylabel([labelFont 'Correlation \chi'], ...
    'Interpreter', 'tex', ...
    'FontSize', insetAxisLabelFontSize)

ax = gca;
ax.FontName = textFontName;
ax.FontSize = insetTickLabelFontSize;

axis square
xlim([-1, 0])
xticks([-1, -0.5, 0])
xticklabels({'-1', '-0.5', '0'})
ylim([chiMin-0.1, chiMax+0.1])

%% Combined 3x3 figure

figure('Name', 'Combined contour, vector-field, and dot plots', 'Color', 'w')
panelLayout = tiledlayout(3, length(Ns), ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

%% Row 1: Jacobian contours

for count = 1:length(Ns)
    nexttile(panelLayout, count)
    hold on

    [~, contourHandle] = contourf(h_vals, plot_lambda_vals, ...
        log10(abs(F_jacobian(lambda_vals > 0.01, :, count))), 15);
    contourHandle.LineColor = 'none';
    colormap(cmap)
    clim(contourCLimits)

    xlabel([labelFont 'External field {\ith}'], ...
        'Interpreter', 'tex', ...
        'FontSize', axisLabelFontSize)
    if count == 1
        ylabel([labelFont 'Interaction strength \lambda'], ...
            'Interpreter', 'tex', ...
            'FontSize', axisLabelFontSize)
    else
        ylabel('')
    end

    axis square

    xlim([-1, criticalHPlot])
    ylim(lambdaPlotLim)
    xscale log
    xticks([-1, -0.1, -0.01, -0.001, criticalHPlot])
    xticklabels({'-10^0', '-10^{-1}', '-10^{-2}', '-10^{-3}', ' 0'})
    yticks(lambdaTicks)
    yticklabels(arrayfun(@num2str, lambdaTicks, 'UniformOutput', false))

    set(gca, 'TickDir', 'both')
    ax = gca;
    ax.FontName = textFontName;
    ax.FontSize = tickLabelFontSize;
    box on

    plot([criticalHPlot, criticalHPlot], [1, lambdaPlotLim(2)], ...
        'Color', criticalColor, ...
        'LineWidth', criticalLineWidth, ...
        'Clipping', 'off')
    scatter(criticalHPlot, 1, criticalPointSize, ...
        'MarkerFaceColor', criticalColor, ...
        'MarkerEdgeColor', criticalColor, ...
        'Clipping', 'off')
end

%% Row 2: RG vector fields

for idx = 1:numel(NValues)
    N = NValues(idx);

    dhdN = N*dhdNFields(:, :, idx);
    dlambdadN = N*dlambdadNFields(:, :, idx);

    plotDh = dhdN./(abs(hGrid)*log(10));
    plotDlambda = dlambdadN;

    nexttile(panelLayout, length(Ns) + idx)

    coloredVectorField( ...
        plotH, lambdaGrid, plotDh, plotDlambda, ...
        'NormalizeArrows', true, ...
        'MagnitudeMode', 'original', ...
        'ArrowLengthMode', arrowLengthMode, ...
        'ArrowScale', arrowScale, ...
        'LineWidth', vectorLineWidth, ...
        'LogMagnitudeReference', logMagnitudeReference, ...
        'LogMagnitudeFloor', logMagnitudeFloor, ...
        'ColorMap', cmap, ...
        'MaxHeadSize', arrowHeadSize, ...
        'UseCustomArrowHeads', true, ...
        'ArrowHeadSize', arrowHeadSize, ...
        'ArrowHeadBaseLength', arrowHeadBaseLength, ...
        'ArrowHeadBaseWidth', arrowHeadBaseWidth);

    xlim([0, xPlotMax])
    ylim(lambdaPlotLim)
    xticks([0 1 2 3 criticalXPlot])
    xticklabels({'10^0', '10^{-1}', '10^{-2}', '10^{-3}', '0'})
    yticks(lambdaTicks)
    yticklabels(arrayfun(@num2str, lambdaTicks, 'UniformOutput', false))

    xlabel([labelFont 'External field {\ith}'], ...
        'Interpreter', 'tex', ...
        'FontSize', axisLabelFontSize)
    if idx == 1
        ylabel([labelFont 'Interaction strength \lambda'], ...
            'Interpreter', 'tex', ...
            'FontSize', axisLabelFontSize)
    else
        ylabel('')
    end

    ax = gca;
    ax.FontName = textFontName;
    ax.FontSize = tickLabelFontSize;
    ax.LineWidth = 1.3;
    ax.TickDir = 'out';
    ax.Layer = 'bottom';
    box on
    axis square

    hold on
    plot([criticalXPlot, criticalXPlot], [1, lambdaPlotLim(2)], ...
        'Color', criticalColor, ...
        'LineWidth', criticalLineWidth, ...
        'Clipping', 'off')
    scatter(criticalXPlot, 1, criticalPointSize, ...
        'MarkerFaceColor', criticalColor, ...
        'MarkerEdgeColor', criticalColor, ...
        'Clipping', 'off')
end

%% Row 3: Mapped dot plots

for a = 1:length(Ns)
    N = Ns(a);

    nexttile(panelLayout, 2*length(Ns) + a)
    hold on

    hs_in = hs(2:end-1, 2:end-1, a);
    ls_in = ls(2:end-1, 2:end-1, a);

    scatter(hs_in(:), ls_in(:), dotSize, 'k', 'filled', MarkerFaceAlpha=0.6);

    for e = 1:4
        plot(hEdges{a,e}, lEdges{a,e}, ...
            'Color', edgeColors{e}, ...
            'LineWidth', edgeWidth);
    end

    plot([0, 0], [1, 3], ...
        'Color', criticalColor, ...
        'LineWidth', criticalLineWidth)
    scatter(0, 1, criticalPointSize, ...
        'MarkerFaceColor', criticalColor, ...
        'MarkerEdgeColor', criticalColor)

    xlabel([labelFont 'External field {\ith}'], ...
        'Interpreter', 'tex', ...
        'FontSize', axisLabelFontSize)
    if a == 1
        ylabel([labelFont 'Interaction strength \lambda'], ...
            'Interpreter', 'tex', ...
            'FontSize', axisLabelFontSize)
    else
        ylabel('')
    end

    axis square

    ylim(lambdaPlotLim)
    yticks(lambdaTicks)
    yticklabels(arrayfun(@num2str, lambdaTicks, 'UniformOutput', false))
    set(gca, 'TickDir', 'both')
    ax = gca;
    ax.FontName = textFontName;
    ax.FontSize = tickLabelFontSize;
    box on
end

set(gcf, 'Renderer', 'painters')

%% Helper figure: Jacobian colorbar

figure('Name', 'Jacobian |J| colorbar', 'Color', 'w')
ax = axes;
colormap(ax, cmap)
clim(ax, contourCLimits)
axis(ax, 'off')

cb = colorbar(ax);
cb.FontSize = tickLabelFontSize + 4;
cb.FontName = textFontName;
cb.Label.String = 'Jacobian |J|';
cb.Label.FontSize = tickLabelFontSize + 4;
cb.Label.FontName = textFontName;
cb.Label.Interpreter = 'tex';

minTick = ceil(contourCLimits(1));
maxTick = floor(contourCLimits(2));
allTicks = minTick:maxTick;

if numel(allTicks) <= 6
    colorbarTicks = allTicks;
else
    tickStep = ceil((maxTick - minTick) / 5);
    colorbarTicks = minTick:tickStep:maxTick;
end

cb.Ticks = colorbarTicks;
cb.TickLabels = arrayfun(@(x) sprintf('10^{%g}', x), colorbarTicks, 'UniformOutput', false);
cb.TickLabelInterpreter = 'tex';

%% Helper figure: flow magnitude colorbar

figure('Name', 'Flow magnitude colorbar', 'Color', 'w')
ax = axes;
colormap(ax, cmap)
clim(ax, [logMagnitudeFloor logMagnitudeReference])
axis(ax, 'off')

cb = colorbar(ax);
cb.FontSize = colorbarFontSize;
cb.Label.String = 'Flow magnitude';
cb.Label.FontSize = colorbarFontSize;

minTick = ceil(logMagnitudeFloor);
maxTick = floor(logMagnitudeReference);
allTicks = minTick:maxTick;

if numel(allTicks) <= 6
    colorbarTicks = allTicks;
else
    tickStep = ceil((maxTick - minTick) / 5);
    colorbarTicks = minTick:tickStep:maxTick;
end

cb.Ticks = colorbarTicks;
cb.TickLabels = arrayfun(@(x) sprintf('10^{%g}', x), colorbarTicks, 'UniformOutput', false);
cb.TickLabelInterpreter = 'tex';
cb.FontSize = colorbarFontSize;

%% Local functions

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

function [Uout, Vout] = normalizeVectors(U, V)
    L = sqrt(U.^2 + V.^2);
    mask = L > 0;

    Uout = zeros(size(U));
    Vout = zeros(size(V));

    Uout(mask) = U(mask)./L(mask);
    Vout(mask) = V(mask)./L(mask);
end
