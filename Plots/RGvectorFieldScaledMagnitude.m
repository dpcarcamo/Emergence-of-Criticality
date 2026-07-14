%% KEEP PLOTS


%% RG vector fields with magnitude-scaled arrows
% Standalone fixed-N vector-field plot. Arrow lengths are scaled by the
% same vector magnitude used for coloring, with the largest arrow across all
% N panels assigned the reference length from arrowScale.

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

Nh = 20;
Nlambda = 20;

arrowScale = 0.15;
axisLabelFontSize = 26;
tickLabelFontSize = 18;
forceRecompute = false;
cacheFile = fullfile(plotsDir, 'RGvectorFieldScaledMagnitude_cache.mat');
cacheVersion = 1;
referenceMagnitudeMode = 'last';  % Options: 'all', 'middle', or 'last'

colors = [
    [10, 10, 10] / 256
    [180, 20, 0] / 256
    [195, 91, 29] / 256
    [235, 172, 71] / 256
    [245, 245, 245] / 256
];
cmap = interpolateColors(colors, 256);

%% Grid

hValues = -logspace(log10(-hlim(2)), log10(-hlim(1)), Nh);
lambdaValues = linspace(lambdalim(1), lambdalim(2), Nlambda);
[hGrid, lambdaGrid] = meshgrid(hValues, lambdaValues);

plotH = -log10(-hGrid);
xMax = max(plotH, [], 'all') - log10(0.5);

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

%% Magnitude scale

[plotDhFields, plotDlambdaFields, magnitudeFields] = ...
    preparePlottedFields(dhdNFields, dlambdadNFields, hGrid);

maxMagnitude = referenceMagnitude(magnitudeFields, referenceMagnitudeMode);
fprintf('Using %s reference magnitude: %.6g\n', referenceMagnitudeMode, maxMagnitude);

%% Plot vector fields

figure('Name', 'RG vector fields with magnitude-scaled arrows', 'Color', 'w');
tiledlayout(1, numel(NValues), ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

for idx = 1:numel(NValues)
    N = NValues(idx);

    nexttile
    coloredVectorField( ...
        plotH, lambdaGrid, plotDhFields(:, :, idx), plotDlambdaFields(:, :, idx), ...
        magnitudeFields(:, :, idx), maxMagnitude, ...
        'ArrowScale', arrowScale, ...
        'ColorMap', cmap, ...
        'MaxHeadSize', 1);

    xlim([0, xMax])
    ylim(lambdalim)
    xticks([0 1 2 3 xMax])
    xticklabels({'10^0', '10^{-1}', '10^{-2}', '10^{-3}', '0'})

    xlabel('-h', 'FontSize', axisLabelFontSize)
    ylabel('\lambda', 'FontSize', axisLabelFontSize)
    title(sprintf('N = %g', N))

    ax = gca;
    ax.FontSize = tickLabelFontSize;
    ax.LineWidth = 1.3;
    ax.TickDir = 'out';
    ax.Layer = 'top';
    box on
    axis square
end

colormap(gcf, cmap)

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

function [plotDhFields, plotDlambdaFields, magnitudeFields] = preparePlottedFields(dhdNFields, dlambdadNFields, hGrid)
    plotDhFields = zeros(size(dhdNFields));
    plotDlambdaFields = dlambdadNFields;
    magnitudeFields = zeros(size(dhdNFields));

    for idx = 1:size(dhdNFields, 3)
        plotDhFields(:, :, idx) = dhdNFields(:, :, idx)./(abs(hGrid)*log(10));
        magnitudeFields(:, :, idx) = sqrt( ...
            plotDhFields(:, :, idx).^2 + plotDlambdaFields(:, :, idx).^2);
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

function maxMagnitude = referenceMagnitude(magnitudeFields, referenceMagnitudeMode)
    switch lower(referenceMagnitudeMode)
        case 'all'
            referenceMagnitudes = magnitudeFields;
        case 'middle'
            middleIndex = ceil(size(magnitudeFields, 3) / 2);
            referenceMagnitudes = magnitudeFields(:, :, middleIndex);
        case 'last'
            referenceMagnitudes = magnitudeFields(:, :, end);
        otherwise
            error('referenceMagnitudeMode must be ''all'', ''middle'', or ''last''.')
    end

    referenceMagnitudes = referenceMagnitudes(isfinite(referenceMagnitudes));

    if isempty(referenceMagnitudes)
        error('No finite vector magnitudes are available for scaling.')
    end

    maxMagnitude = max(referenceMagnitudes);
end

function hq = coloredVectorField(X, Y, U, V, M, maxMagnitude, varargin)
    p = inputParser;

    addParameter(p, 'ColorByMagnitude', true);
    addParameter(p, 'ColorMap', turbo(256));
    addParameter(p, 'ArrowScale', 0.6);
    addParameter(p, 'LineWidth', 2);
    addParameter(p, 'MaxHeadSize', 0.7);
    addParameter(p, 'Skip', 1);

    parse(p, varargin{:});
    opts = p.Results;

    skip = opts.Skip;
    X = X(1:skip:end, 1:skip:end);
    Y = Y(1:skip:end, 1:skip:end);
    U = U(1:skip:end, 1:skip:end);
    V = V(1:skip:end, 1:skip:end);
    M = M(1:skip:end, 1:skip:end);

    L = sqrt(U.^2 + V.^2);
    mask = L > 0 & M > 0 & isfinite(M) & isfinite(maxMagnitude) & maxMagnitude > 0;

    Uscaled = zeros(size(U));
    Vscaled = zeros(size(V));
    scaleByMagnitude = opts.ArrowScale * M(mask)./maxMagnitude;

    Uscaled(mask) = U(mask)./L(mask) .* scaleByMagnitude;
    Vscaled(mask) = V(mask)./L(mask) .* scaleByMagnitude;

    Xf = X(:);
    Yf = Y(:);
    Uf = Uscaled(:);
    Vf = Vscaled(:);
    Mf = log10(M(:));

    valid = ~(isnan(Xf) | isnan(Yf) | isnan(Uf) | isnan(Vf) | isnan(Mf) | ...
              isinf(Xf) | isinf(Yf) | isinf(Uf) | isinf(Vf) | isinf(Mf));

    Xf = Xf(valid);
    Yf = Yf(valid);
    Uf = Uf(valid);
    Vf = Vf(valid);
    Mf = Mf(valid);

    cmap = opts.ColorMap;
    nColors = size(cmap, 1);

    Mmin = min(Mf);
    Mmax = max(Mf);

    if Mmax == Mmin
        colorInds = ones(size(Mf));
    else
        colorInds = round(1 + (nColors - 1)*(Mf - Mmin)/(Mmax - Mmin));
    end

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

        hq(k) = quiver(Xf(k), Yf(k), Uf(k), Vf(k), 0, ...
            'Color', arrowColor, ...
            'LineWidth', opts.LineWidth, ...
            'MaxHeadSize', opts.MaxHeadSize);
    end

    if opts.ColorByMagnitude
        colormap(cmap)
        clim([-7 0])
    end

    axis tight
    box on

    if ~holdState
        hold off
    end
end

function cmap = interpolateColors(colors, nColors)
    nAnchors = size(colors, 1);
    xAnchors = linspace(0, 1, nAnchors);
    xInterp = linspace(0, 1, nColors);
    cmap = interp1(xAnchors, colors, xInterp, 'linear');
    cmap = max(min(cmap, 1), 0);
end
