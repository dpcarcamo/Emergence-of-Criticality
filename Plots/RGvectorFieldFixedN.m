%% KEEP PLOTS


%% RG vector fields for fixed N values
% Standalone fixed-N vector-field plot.

clearvars
clc

currentPath = pwd;
folders = split(currentPath, '\');
newPath = join(folders(1:length(folders)-1),"\");
addpath(newPath{1});

%% Plot settings

NValues = [20 100 1000];

hlim = [-0.9 -0.001];
lambdalim = [0.001 3];

Nh = 20;
Nlambda = 20;

arrowScale = 0.15;
axisLabelFontSize = 26;
tickLabelFontSize = 18;

colors = [
    [255, 216, 118]/256
    [235, 172, 71]/256
    [195, 91, 29]/256
    [180, 20, 0]/256
];
cmap = interpolateColors(colors, 256);

%% Grid

hValues = -logspace(log10(-hlim(2)), log10(-hlim(1)), Nh);
lambdaValues = linspace(lambdalim(1), lambdalim(2), Nlambda);
[hGrid, lambdaGrid] = meshgrid(hValues, lambdaValues);

plotH = -log10(-hGrid);
xMax = max(plotH, [], 'all') - log10(0.5);

%% Plot vector fields

figure('Name', 'RG vector fields for fixed N', 'Color', 'w');
tiledlayout(1, numel(NValues), ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

for idx = 1:numel(NValues)
    N = NValues(idx);
    fprintf('Evaluating vector field for N = %g\n', N);

    [dhdN, dlambdadN] = evalFieldOnGrid(hGrid, lambdaGrid, N);

    plotDh = dhdN./(abs(hGrid)*log(10));
    plotDlambda = dlambdadN;

    nexttile
    coloredVectorField( ...
        plotH, lambdaGrid, plotDh, plotDlambda, ...
        'NormalizeArrows', true, ...
        'MagnitudeMode', 'original', ...
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

function [U, V] = evalFieldOnGrid(hGrid, lambdaGrid, N)
    U = zeros(size(hGrid));
    V = zeros(size(lambdaGrid));

    for ii = 1:numel(hGrid)
        h = hGrid(ii);
        lambda = lambdaGrid(ii);

        U(ii) = real(DhDNtot(h, lambda, N));
        V(ii) = real(DlamDNtot(h, lambda, N));
    end
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

    if opts.NormalizeArrows
        L = sqrt(U.^2 + V.^2);
        mask = L > 0;

        Utemp = zeros(size(U));
        Vtemp = zeros(size(V));

        Utemp(mask) = U(mask)./L(mask);
        Vtemp(mask) = V(mask)./L(mask);

        U = Utemp;
        V = Vtemp;
    end

    Xf = X(:);
    Yf = Y(:);
    Uf = U(:);
    Vf = V(:);
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

        hq(k) = quiver(Xf(k), Yf(k), Uf(k), Vf(k), opts.ArrowScale, ...
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
