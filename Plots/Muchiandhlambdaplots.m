%% KEEP PLOTS

% Plots data points in (mu, chi) and stored data points in (h, lambda).
% Figure mapping: Fig. 1 uses the full data points in statistics space and
% model space. Fig. 4 uses the N = 20, 100, 1000, and full-data points.

clear
clc

plotsDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(plotsDir);
addpath(repoRoot)
addpath(fullfile(repoRoot, 'Allen'))
addpath(fullfile(repoRoot, 'Hippo'))
addpath(fullfile(repoRoot, 'Stringer'))

files = ["Allenhldata.mat", "hippomuchidata.mat", "stringerhldata.mat"];
dataMarker = 'o';
colors = ["#2676ad", "#4f4cc4", "#2f682c"];
referenceLineColor = "#939598";
markerSize = 70;
markerEdgeColor = [0 0 0];
markerLineWidth = 0.5;
markerFaceAlpha = 0.75;
fullMarkerSizeRange = [45 220];
fullNRange = fullDataNRange(files);
axisLabelFontSize = 30;
tickLabelFontSize = 18;
textFontName = 'Helvetica';
labelFont = ['\fontname{' textFontName '}'];
hZeroPlot = -0.00005;
criticalPointSize = 45;
modelYLim = [0, 2.5];

% -1 denotes the final/full data point.
figureSets = {[-1]};
%figureSets = {[20 100 1000 -1]};

%% Plot mu and chi from data

for f = 1:length(figureSets)
    figure
    hold on

    currentSet = figureSets{f};
    scaleFullDataMarkers = isFullDataOnlySet(currentSet);

    for i = 1:length(files)
        [Data, ~, ~] = loadDataFile(files(i));
        [Nss, mus, chis] = summarizeData(Data);

        baseColor = hex2rgb(colors(i));
        [c20, c100, c1000, cFull] = dataAmountColors(baseColor);
        dim = datasetSearchDim(i);

        for val = currentSet
            if val == -1
                if scaleFullDataMarkers
                    fullNs = lastofarray(Nss);
                    fullMarkerSizes = markerSizeFromN(fullNs, fullNRange, fullMarkerSizeRange);
                    scatter(lastofarray(mus), lastofarray(chis), fullMarkerSizes, ...
                        'Marker', dataMarker, ...
                        'MarkerEdgeColor', markerEdgeColor, ...
                        'MarkerFaceColor', cFull, ...
                        'MarkerFaceAlpha', markerFaceAlpha, ...
                        'LineWidth', markerLineWidth)
                else
                    scatter(lastofarray(mus), lastofarray(chis), markerSize, ...
                        'Marker', dataMarker, ...
                        'MarkerEdgeColor', markerEdgeColor, ...
                        'MarkerFaceColor', cFull, ...
                        'MarkerFaceAlpha', markerFaceAlpha, ...
                        'LineWidth', markerLineWidth)
                end
            elseif val == 0
                continue
            else
                [~, I] = min(abs(Nss - val), [], dim, "linear");
                pointColor = colorForN(val, c20, c100, c1000, cFull);
                scatter(mus(I), chis(I), markerSize, ...
                    'Marker', dataMarker, ...
                    'MarkerEdgeColor', markerEdgeColor, ...
                    'MarkerFaceColor', pointColor, ...
                    'MarkerFaceAlpha', markerFaceAlpha, ...
                    'LineWidth', markerLineWidth)
            end
        end
    end

    x = linspace(-1, 0, 1000);
    plot(x, 1 - x.^2, 'LineWidth', 2, 'Color', referenceLineColor)

    ylim([0.05, 70])
    xlim([-1, -0.8])
    axis square
    xlabel([labelFont 'Average activity \langle\mu\rangle'], ...
        'Interpreter', 'tex', ...
        'FontSize', axisLabelFontSize)
    ylabel([labelFont 'Susceptibility \chi'], ...
        'Interpreter', 'tex', ...
        'FontSize', axisLabelFontSize)
    set(gca, 'YScale', 'log')
    set(gca, 'TickDir', 'both')
    ax = gca;
    ax.FontName = textFontName;
    ax.FontSize = tickLabelFontSize;
    box on
end

%% Plot h and lambda from data

for f = 1:length(figureSets)
    figure
    hold on

    plot(hZeroPlot + 0 * linspace(1, 2), linspace(1, 2.5), 'r', 'LineWidth', 2)
    plot(linspace(-1.4, hZeroPlot), 0 * linspace(1, 2.5), 'Color', referenceLineColor, 'LineWidth', 2)

    currentSet = figureSets{f};
    scaleFullDataMarkers = isFullDataOnlySet(currentSet);

    for i = 1:length(files)
        [Data, hs, ls] = loadDataFile(files(i));
        [Nss, ~, ~] = summarizeData(Data);

        baseColor = hex2rgb(colors(i));
        [c20, c100, c1000, cFull] = dataAmountColors(baseColor);
        dim = datasetSearchDim(i);

        for val = currentSet
            if val == -1
                if scaleFullDataMarkers
                    fullNs = lastofarray(Nss);
                    fullMarkerSizes = markerSizeFromN(fullNs, fullNRange, fullMarkerSizeRange);
                    scatter(lastofarray(hs), lastofarray(ls), fullMarkerSizes, ...
                        'Marker', dataMarker, ...
                        'MarkerEdgeColor', markerEdgeColor, ...
                        'MarkerFaceColor', cFull, ...
                        'MarkerFaceAlpha', markerFaceAlpha, ...
                        'LineWidth', markerLineWidth)
                else
                    scatter(lastofarray(hs), lastofarray(ls), markerSize, ...
                        'Marker', dataMarker, ...
                        'MarkerEdgeColor', markerEdgeColor, ...
                        'MarkerFaceColor', cFull, ...
                        'MarkerFaceAlpha', markerFaceAlpha, ...
                        'LineWidth', markerLineWidth)
                end
            elseif val == 0
                continue
            else
                [~, I] = min(abs(Nss - val), [], dim, "linear");
                pointColor = colorForN(val, c20, c100, c1000, cFull);
                scatter(hs(I), ls(I), markerSize, ...
                    'Marker', dataMarker, ...
                    'MarkerEdgeColor', markerEdgeColor, ...
                    'MarkerFaceColor', pointColor, ...
                    'MarkerFaceAlpha', markerFaceAlpha, ...
                    'LineWidth', markerLineWidth)
            end
        end
    end

    scatter(hZeroPlot, 1, criticalPointSize, ...
        'Marker', dataMarker, ...
        'MarkerFaceColor', 'r', ...
        'MarkerEdgeColor', 'r')

    axis square
    xlabel([labelFont 'External field {\ith}'], ...
        'Interpreter', 'tex', ...
        'FontSize', axisLabelFontSize)
    ylabel([labelFont 'Interaction strength {\itJ}'], ...
        'Interpreter', 'tex', ...
        'FontSize', axisLabelFontSize)
    xscale log
    box on
    set(gca, 'TickDir', 'both')
    xlim([-1, hZeroPlot])
    xticks([-1, -0.1, -0.01, -0.001, -0.0001, hZeroPlot])
    xticklabels({'-10^{0}', '-10^{-1}', '-10^{-2}', '-10^{-3}', ...
        '-10^{-4}', '0'})
    ylim(modelYLim)
    ax = gca;
    ax.FontName = textFontName;
    ax.FontSize = tickLabelFontSize;
    ax.TickLabelInterpreter = 'tex';
end

%% Plot full data in h and J with a linear h axis

figure
hold on

plot([0, 0], [1,modelYLim(2)], 'r', 'LineWidth', 2)
plot([-1, 1], [0, 0], 'Color', referenceLineColor, 'LineWidth', 2)

for i = 1:length(files)
    [Data, hs, ls] = loadDataFile(files(i));
    [Nss, ~, ~] = summarizeData(Data);

    baseColor = hex2rgb(colors(i));
    [~, ~, ~, cFull] = dataAmountColors(baseColor);

    fullNs = lastofarray(Nss);
    fullMarkerSizes = markerSizeFromN(fullNs, fullNRange, fullMarkerSizeRange);

    scatter(lastofarray(hs), lastofarray(ls), fullMarkerSizes, ...
        'Marker', dataMarker, ...
        'MarkerEdgeColor', markerEdgeColor, ...
        'MarkerFaceColor', cFull, ...
        'MarkerFaceAlpha', markerFaceAlpha, ...
        'LineWidth', markerLineWidth)
end

scatter(0, 1, criticalPointSize, ...
    'Marker', dataMarker, ...
    'MarkerFaceColor', 'r', ...
    'MarkerEdgeColor', 'r')

axis square
xlim([-1, 1])
ylim(modelYLim)
xlabel([labelFont 'External field {\ith}'], ...
    'Interpreter', 'tex', ...
    'FontSize', axisLabelFontSize)
ylabel([labelFont 'Interaction strength {\itJ}'], ...
    'Interpreter', 'tex', ...
    'FontSize', axisLabelFontSize)
box on
set(gca, 'TickDir', 'both')
ax = gca;
ax.FontName = textFontName;
ax.FontSize = tickLabelFontSize;

function [Data, hs, ls] = loadDataFile(fileName)
    obj = load(fileName);

    Data = obj.Data;
    hs = obj.hs;
    ls = obj.ls;

    hs(hs == 0) = nan;
    ls(ls == 0) = nan;
    Data(Data == 0) = nan;
end

function [Nss, mus, chis] = summarizeData(Data)
    Nss = squeeze(mean(squeeze(Data(3, :, :, :)), 1, 'omitnan'));
    mus = squeeze(mean(squeeze(Data(1, :, :, :)), 1, 'omitnan'));
    chis = squeeze(mean(squeeze(Data(2, :, :, :)), 1, 'omitnan'));
end

function dim = datasetSearchDim(datasetIndex)
    dim = 1;
    if datasetIndex == 2
        dim = 2;
    end
end

function [c20, c100, c1000, cFull] = dataAmountColors(baseColor)
    c20 = lightenColor(baseColor, 0.60);
    c100 = lightenColor(baseColor, 0.45);
    c1000 = lightenColor(baseColor, 0.19);
    cFull = baseColor;
end

function tf = isFullDataOnlySet(currentSet)
    tf = numel(currentSet) == 1 && currentSet == -1;
end

function nRange = fullDataNRange(files)
    allFullNs = [];

    for i = 1:length(files)
        [Data, ~, ~] = loadDataFile(files(i));
        [Nss, ~, ~] = summarizeData(Data);
        allFullNs = [allFullNs, lastofarray(Nss)];
    end

    allFullNs = allFullNs(isfinite(allFullNs));
    nRange = [min(allFullNs), max(allFullNs)];
end

function markerSizes = markerSizeFromN(Ns, nRange, sizeRange)
    if nRange(1) == nRange(2)
        markerSizes = mean(sizeRange) * ones(size(Ns));
        return
    end

    scaledN = (Ns - nRange(1)) ./ diff(nRange);
    markerSizes = sizeRange(1) + scaledN .* diff(sizeRange);
end

function c = colorForN(N, c20, c100, c1000, cFull)
    if N == 20
        c = c20;
    elseif N == 100
        c = c100;
    elseif N == 1000
        c = c1000;
    else
        c = cFull;
    end
end

function cOut = lightenColor(cIn, amt)
    cOut = cIn + (1 - cIn) * amt;
end

function rgb = hex2rgb(hex)
    hex = char(erase(string(hex), "#"));
    rgb = sscanf(hex, '%2x%2x%2x', [1 3]) / 255;
end

function out = lastofarray(A)
    if size(A, 1) == 1
        A = A.';
    end

    A(:, isnan(A(1, :))) = [];
    valid = ~isnan(A);
    indices = arrayfun(@(x) find(valid(:, x), 1, 'last'), 1:size(A, 2));
    out = arrayfun(@(x, y) A(x, y), indices, 1:size(A, 2));
end
