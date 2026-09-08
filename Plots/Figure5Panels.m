%% KEEP PLOTS

% Standalone script for the Figure 5 layout. The top row contains the
% full-data statistics panel and the inferred free-energy panel. The bottom
% left region contains the exact and double-well parameter curves. The bottom
% right panel shows double-well parameter trajectories for fixed full-data
% statistics.

clear
clc

plotsDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(plotsDir);
addpath(repoRoot)
addpath(fullfile(repoRoot, 'Allen'))
addpath(fullfile(repoRoot, 'Hippo'))
addpath(fullfile(repoRoot, 'Stringer'))

%% Figure settings

files = ["Allenhldata.mat", "stringerhldata2.mat", "hippomuchidata.mat"];
datasetNames = ["Allen", "Stringer", "Hippocampus"];
datasetColors = ["#2676ad", "#2f682c", "#4f4cc4"];

textFontName = 'Helvetica';
labelFont = ['\fontname{' textFontName '}'];
axisLabelFontSize = 22;
subAxisLabelFontSize = 16;
tickLabelFontSize = 14;
legendFontSize = 12;

referenceLineColor = "#414042";
markerSize = 130;
markerEdgeColor = [0 0 0];
markerLineWidth = 0.7;
markerFaceAlpha = 0.45;
fullMarkerSizeRange = [85 240];
fullDataMarkerSizeScale = 1.15;
nLegendValues = [1000 10000];
statisticsXLim = [-1, -0.8];
statisticsYLim = [0.05, 70];
upperBoundShadeColor = [1, 0, 0];
upperBoundShadeAlpha = 0.15;
upperBoundLineWidth = 2.6;
criticalLabelColor = [0 0 0];

freeEnergyNs = [20, 50, 100, 200, 1000];
freeEnergyColors = ["#f9cd6b", "#edb34f", "#c2591c", "#ba340d", "#b41a02"];
freeEnergyInfinityColor = [0 0 0];
freeEnergyLineWidth = 3.3;
% Statistics used for the inferred free-energy panel.
fixedStatisticsMu = -0.862;
fixedStatisticsChi = 0.74;
muGrid = linspace(-0.999, 0.999, 1000);

theoryLineWidth = 3.5;
infinityLineWidth = 3.3;
plotExternalFieldLogY = true;
panelCYLimits = [
    -10, -1e-3
    -10, -1e-4
    -10, -1e-3
];
panelCYTicks = {
    [-10, -1e-1, -1e-3]
    [-10, -1e-2, -1e-4]
    [-10, -1e-1, -1e-3]
};
panelCYTickLabels = {
    {'-10^{1}', '-10^{-1}', '-10^{-3}'}
    {'-10^{1}', '-10^{-2}', '-10^{-4}'}
    {'-10^{1}', '-10^{-1}', '-10^{-3}'}
};
panelDYLimits = [0, 2.5];
panelDYTicks = 0:floor(panelDYLimits(2));
parameterAxisLabelFontSize = 24;
parameterTickLabelFontSize = 20;
parameterMarkerSize = 400;
parameterPopulationMode = "single";  % "single" or "all"
% In "single" mode, entries are [Allen, Stringer, Hippocampus].
% Use nan to randomly choose one population among those with the most valid N values.
selectedPopulationByDataset = [nan, nan, nan];
autoPopulationSelection = "randomMaxN";
xTickValues = 10.^(0:4);
xTickLabels = arrayfun(@(x) sprintf('10^{%d}', x), 0:4, 'UniformOutput', false);

trajectoryNMin = 1e2;
trajectoryNMax = 1e6;
trajectoryNCount = 50;
trajectoryNValues = unique(round(logspace(log10(trajectoryNMin), ...
    log10(trajectoryNMax), trajectoryNCount)));
trajectoryLineWidth = theoryLineWidth;
trajectoryLowNLightenAmount = 0.60;
trajectoryColorMidpoint = 0.55;
trajectoryHighNColor = [0.02 0.02 0.02];
trajectoryNLegendLineWidth = 4;
trajectoryLegendSegmentCount = 32;
trajectoryLegendX = -10.^linspace(-1.35, -2.85, trajectoryLegendSegmentCount + 1);
trajectoryLegendY = [2.73, 2.52, 2.31];
trajectoryReferenceLineColor = "#939598";
trajectoryHZeroPlot = -10^-6.5;
trajectoryHXLim = [-0.1, trajectoryHZeroPlot];
trajectoryHTickExponents = -1:-1:-6;
trajectoryHTickValues = [-10.^trajectoryHTickExponents, trajectoryHZeroPlot];
trajectoryHTickLabels = arrayfun(@(x) sprintf('-10^{%d}', x), ...
    trajectoryHTickExponents, 'UniformOutput', false);
trajectoryHTickLabels{end + 1} = '    0';
linearTrajectoryHXLim = [-0.05, 0];
linearTrajectoryHTicks = -0.05:0.01:0;
trajectoryLambdaYLim = [0, 3];
linearTrajectoryLambdaYLim = [0, trajectoryLambdaYLim(2)];
trajectoryLambdaTicks = 0:1:3;
trajectoryCriticalLineWidth = 2;
trajectoryCriticalPointSize = 85;
combinedParameterAxisLabelFontSize = 28;
combinedParameterTickLabelFontSize = 22;
combinedParameterMarkerSize = 260;

%% Cached computed data

forceRecompute = false;
cacheVersion = 1;
cacheFile = fullfile(plotsDir, 'Figure5Panels_cache.mat');

cacheSettings = struct();
cacheSettings.cacheVersion = cacheVersion;
cacheSettings.files = files;
cacheSettings.datasetColors = datasetColors;
cacheSettings.statisticsXLim = statisticsXLim;
cacheSettings.statisticsYLim = statisticsYLim;
cacheSettings.fullMarkerSizeRange = fullMarkerSizeRange;
cacheSettings.parameterPopulationMode = parameterPopulationMode;
cacheSettings.selectedPopulationByDataset = selectedPopulationByDataset;
cacheSettings.autoPopulationSelection = autoPopulationSelection;

fprintf('Figure 5: checking cached computed data.\n')
useCache = false;
cacheSettingsMatch = false;
cacheData = struct();

if ~forceRecompute && exist(cacheFile, 'file') == 2
    cached = load(cacheFile, 'cacheData', 'cacheSettings', 'cacheVersion');

    if isfield(cached, 'cacheData') ...
            && isfield(cached, 'cacheSettings') ...
            && isfield(cached, 'cacheVersion') ...
            && cached.cacheVersion == cacheVersion ...
            && isFigure5CacheCompatible(cached.cacheSettings, cacheSettings)
        cacheData = cached.cacheData;
        cacheSettingsMatch = true;

        if isfield(cacheData, 'isComplete') && cacheData.isComplete
            useCache = true;
            fprintf('Figure 5: loaded complete cached data from %s.\n', cacheFile)
        else
            fprintf('Figure 5: found partial cached data. Resuming missing steps.\n')
        end
    else
        fprintf('Figure 5: cache exists but settings changed. Recomputing.\n')
    end
end

if ~useCache
    if ~cacheSettingsMatch
        cacheData = struct();
    end

    if ~isfield(cacheData, 'panelB')
        fprintf('Figure 5: computing panel a full-data points.\n')
        cacheData.panelB = computePanelBData(files, datasetColors, fullMarkerSizeRange);
        saveFigure5Cache(cacheFile, cacheData, cacheSettings, cacheVersion, false)
    else
        fprintf('Figure 5: using cached panel a full-data points.\n')
    end

    fprintf('Figure 5: computing panels c and d parameter curves.\n')
    if ~isfield(cacheData, 'parameterCurves')
        cacheData.parameterCurves = struct('datasetName', {}, 'populationCurves', {});
    end

    cacheData.parameterCurves = computeAllParameterCurveData( ...
        files, datasetNames, datasetColors, ...
        parameterPopulationMode, selectedPopulationByDataset, ...
        cacheFile, cacheSettings, cacheVersion, cacheData);

    saveFigure5Cache(cacheFile, cacheData, cacheSettings, cacheVersion, true)
    fprintf('Figure 5: saved complete computed data to %s.\n', cacheFile)
end

%% Layout

mainFigure = figure('Name', 'Figure 5 panels', ...
    'Color', 'w');

panelAPosition = [0.08 0.56 0.34 0.34];
panelBPosition = [0.56 0.56 0.34 0.34];
panelEPosition = [0.56 0.10 0.34 0.34];

parameterGridLeft = 0.10;
parameterGridBottom = 0.10;
parameterGridSize = 0.80;
parameterGridColumnGap = 0.07;
parameterGridRowGap = 0.055;
parameterSubWidth = (parameterGridSize - parameterGridColumnGap) / 2;
parameterSubHeight = (parameterGridSize - 2 * parameterGridRowGap) / 3;
panelCPositions = zeros(length(files), 4);
panelDPositions = zeros(length(files), 4);

for datasetIdx = 1:length(files)
    rowBottom = parameterGridBottom ...
        + (length(files) - datasetIdx) * (parameterSubHeight + parameterGridRowGap);
    panelCPositions(datasetIdx, :) = [parameterGridLeft, rowBottom, ...
        parameterSubWidth, parameterSubHeight];
    panelDPositions(datasetIdx, :) = [parameterGridLeft + parameterSubWidth + parameterGridColumnGap, ...
        rowBottom, parameterSubWidth, parameterSubHeight];
end

% Panel map:
% a: full neural data relative to the critical region
% b: inferred free energy for fixed statistics
% c: inferred external field versus N
% d: inferred interaction strength versus N
% e: double-well parameter trajectories for fixed full-data statistics
% Positions are [left bottom width height]. Panels a, b, and e are on the
% main figure. The c/d block is plotted in a separate square figure.

%% a. Full neural data relative to the single-minimum region

axA = axes('Parent', mainFigure, 'Position', panelAPosition);
hold(axA, 'on')
fprintf('Figure 5: plotting panel a.\n')

x = linspace(statisticsXLim(1), statisticsXLim(2), 1000);
independentBound = 1 - x.^2;
upperBound = x .* (1 - x.^2) ./ (x - atanh(x) .* (1 - x.^2));
validShade = isfinite(upperBound) & upperBound > 0 & upperBound < statisticsYLim(2);

fill(axA, [x(validShade), fliplr(x(validShade))], ...
    [upperBound(validShade), statisticsYLim(2) * ones(1, nnz(validShade))], ...
    upperBoundShadeColor, ...
    'FaceAlpha', upperBoundShadeAlpha, ...
    'EdgeColor', 'none')
plot(axA, x, independentBound, ...
    'LineWidth', upperBoundLineWidth, ...
    'Color', referenceLineColor)
plot(axA, x(validShade), upperBound(validShade), ...
    'LineWidth', upperBoundLineWidth, ...
    'Color', upperBoundShadeColor)

for i = 1:numel(cacheData.panelB)
    scatter(axA, cacheData.panelB(i).mus, cacheData.panelB(i).chis, ...
        cacheData.panelB(i).markerSizes * fullDataMarkerSizeScale, ...
        'Marker', 'o', ...
        'MarkerEdgeColor', markerEdgeColor, ...
        'MarkerFaceColor', cacheData.panelB(i).color, ...
        'MarkerFaceAlpha', markerFaceAlpha, ...
        'LineWidth', markerLineWidth)
end

text(axA, -0.895, 25, 'Critical', ...
    'Color', criticalLabelColor, ...
    'FontName', textFontName, ...
    'FontSize', legendFontSize, ...
    'HorizontalAlignment', 'center')

fullNsForLegend = [];
for i = 1:numel(cacheData.panelB)
    if isfield(cacheData.panelB, 'fullNs')
        fullNsForLegend = [fullNsForLegend, cacheData.panelB(i).fullNs]; %#ok<AGROW>
    end
end

fullNsForLegend = fullNsForLegend(isfinite(fullNsForLegend));
if isempty(fullNsForLegend)
    fullNRangeForLegend = fullDataNRange(files);
else
    fullNRangeForLegend = [min(fullNsForLegend), max(fullNsForLegend)];
end

nLegendHandles = gobjects(1, numel(nLegendValues));
nLegendLabels = strings(1, numel(nLegendValues));
for k = 1:numel(nLegendValues)
    nLegendHandles(k) = scatter(axA, nan, nan, ...
        markerSizeFromN(nLegendValues(k), fullNRangeForLegend, fullMarkerSizeRange) * fullDataMarkerSizeScale, ...
        'Marker', 'o', ...
        'MarkerEdgeColor', markerEdgeColor, ...
        'MarkerFaceColor', [0.65 0.65 0.65], ...
        'LineWidth', markerLineWidth);
    nLegendLabels(k) = sprintf('N = %d', nLegendValues(k));
end

legend(axA, nLegendHandles, nLegendLabels, ...
    'Interpreter', 'tex', ...
    'FontName', textFontName, ...
    'FontSize', legendFontSize, ...
    'Location', 'southwest')

xlim(axA, statisticsXLim)
ylim(axA, statisticsYLim)
set(axA, 'YScale', 'log')
xlabel(axA, [labelFont 'Average activity {\itm}'], ...
    'Interpreter', 'tex', ...
    'FontSize', axisLabelFontSize)
ylabel(axA, [labelFont 'Correlation \chi'], ...
    'Interpreter', 'tex', ...
    'FontSize', axisLabelFontSize)
formatSquareAxis(axA, textFontName, tickLabelFontSize)

%% b. Inferred free energy at fixed statistics

axB = axes('Parent', mainFigure, 'Position', panelBPosition);
hold(axB, 'on')
fprintf('Figure 5: computing and plotting panel b.\n')

freeEnergyLabels = strings(1, numel(freeEnergyNs) + 1);

for k = 1:numel(freeEnergyNs)
    N = freeEnergyNs(k);
    [h, lambda] = hlambda(fixedStatisticsMu, fixedStatisticsChi, N);

    score = (N * h * muGrid + N * lambda * muGrid.^2 / 2 ...
        + logGammaLanczos(N + 1) ...
        - logGammaLanczos(N * (1 + muGrid) / 2 + 1) ...
        - logGammaLanczos(N * (1 - muGrid) / 2 + 1)) / N;

    peaks = findpeaks(score);
    if ~isempty(peaks)
        referenceValue = peaks(1);
    else
        referenceValue = max(score);
    end
    y = -(score - referenceValue);

    plot(axB, muGrid, y, ...
        'Color', freeEnergyColors(k), ...
        'LineWidth', freeEnergyLineWidth)

    freeEnergyLabels(k) = sprintf('N = %d', N);
end

lambdaInf = atanh(abs(fixedStatisticsMu)) / max(abs(fixedStatisticsMu), eps);
pPlus = (1 + muGrid) / 2;
pMinus = (1 - muGrid) / 2;
entropyTerm = pPlus .* log(pPlus) + pMinus .* log(pMinus);
scoreInf = lambdaInf * muGrid.^2 / 2 - entropyTerm;
peaks = findpeaks(scoreInf);
if ~isempty(peaks)
    referenceValue = peaks(1);
else
    referenceValue = max(scoreInf);
end
yInf = -(scoreInf - referenceValue);

plot(axB, muGrid, yInf, ...
    'Color', freeEnergyInfinityColor, ...
    'LineWidth', freeEnergyLineWidth)
freeEnergyLabels(end) = 'N = \infty';

legend(axB, freeEnergyLabels, ...
    'Interpreter', 'tex', ...
    'FontName', textFontName, ...
    'FontSize', 9, ...
    'Location', 'best')

xlabel(axB, [labelFont 'Activity \mu'], ...
    'Interpreter', 'tex', ...
    'FontSize', subAxisLabelFontSize)
ylabel(axB, [labelFont 'Free energy {\itf}(\mu)'], ...
    'Interpreter', 'tex', ...
    'FontSize', subAxisLabelFontSize)
formatSubAxis(axB, textFontName, tickLabelFontSize)
axis(axB, 'square')

%% Separate c/d figure

sharedParameterNValues = [];

for datasetIdx = 1:length(files)
    datasetCurveData = cacheData.parameterCurves(datasetIdx);

    for popIdx = 1:numel(datasetCurveData.populationCurves)
        sharedParameterNValues = [sharedParameterNValues; ...
            datasetCurveData.populationCurves(popIdx).Nplot(:)]; %#ok<AGROW>
    end
end

sharedParameterNValues = sharedParameterNValues( ...
    isfinite(sharedParameterNValues) & sharedParameterNValues > 0);

if isempty(sharedParameterNValues)
    parameterXLimits = [min(xTickValues), max(xTickValues)];
else
    parameterXLimits = [min(sharedParameterNValues), max(sharedParameterNValues)];
end

parameterFigure = figure('Name', 'Figure 5 panels c and d', ...
    'Color', 'w');

%% c. Exact inferred field and double-well approximation

for datasetIdx = 1:length(files)
    fprintf('Figure 5: plotting panel c for %s.\n', datasetNames(datasetIdx))
    axC = axes('Parent', parameterFigure, 'Position', panelCPositions(datasetIdx, :));
    hold(axC, 'on')

    datasetCurveData = cacheData.parameterCurves(datasetIdx);

    if isfield(datasetCurveData, 'selectedPopulation') ...
            && ~isempty(datasetCurveData.selectedPopulation) ...
            && isfinite(datasetCurveData.selectedPopulation)
        fprintf('  plotting population %d of %s.\n', ...
            datasetCurveData.selectedPopulation, datasetNames(datasetIdx))
    end

    for popIdx = 1:numel(datasetCurveData.populationCurves)
        curve = datasetCurveData.populationCurves(popIdx);
        popColor = curve.popColor;

        if isfield(curve, 'sourcePopulation') && ~isempty(curve.sourcePopulation)
            fprintf('  panel c curve %d uses source population %d.\n', ...
                popIdx, curve.sourcePopulation)
        end

        % Circles are exact inversions. Solid lines are the double-well approximation.
        scatter(axC, curve.Nplot, curve.hFit, parameterMarkerSize, ...
            'Marker', 'o', ...
            'MarkerEdgeColor', markerEdgeColor, ...
            'MarkerFaceColor', popColor, ...
            'MarkerFaceAlpha', markerFaceAlpha, ...
            'LineWidth', markerLineWidth)
        plot(axC, curve.Nplot, curve.hTheory, ...
            'Color', popColor, ...
            'LineWidth', theoryLineWidth)

        if ~plotExternalFieldLogY
            plot(axC, [min(curve.Nplot), max(curve.Nplot)], curve.hInfinity * [1, 1], ...
                '--', ...
                'Color', popColor, ...
                'LineWidth', infinityLineWidth)
        end
    end

    ylabel(axC, [labelFont 'External field {\ith}'], ...
        'Interpreter', 'tex', ...
        'FontSize', parameterAxisLabelFontSize)

    formatParameterAxis(axC, textFontName, parameterTickLabelFontSize, xTickValues, xTickLabels)
    if plotExternalFieldLogY
        set(axC, 'YScale', 'log')
        ylim(axC, panelCYLimits(datasetIdx, :))
        yticks(axC, panelCYTicks{datasetIdx})
        yticklabels(axC, panelCYTickLabels{datasetIdx})
    end
    xlim(axC, parameterXLimits)

    if datasetIdx == length(files)
        xlabel(axC, [labelFont 'Number of neurons {\itN}'], ...
            'Interpreter', 'tex', ...
            'FontSize', parameterAxisLabelFontSize)
    else
        xticklabels(axC, repmat({''}, size(get(axC, 'XTick'))))
    end
end

%% d. Exact inferred interaction strength and double-well approximation

for datasetIdx = 1:length(files)
    fprintf('Figure 5: plotting panel d for %s.\n', datasetNames(datasetIdx))
    axD = axes('Parent', parameterFigure, 'Position', panelDPositions(datasetIdx, :));
    hold(axD, 'on')

    datasetCurveData = cacheData.parameterCurves(datasetIdx);

    for popIdx = 1:numel(datasetCurveData.populationCurves)
        curve = datasetCurveData.populationCurves(popIdx);
        popColor = curve.popColor;

        if isfield(curve, 'sourcePopulation') && ~isempty(curve.sourcePopulation)
            fprintf('  panel d curve %d uses source population %d.\n', ...
                popIdx, curve.sourcePopulation)
        end

        % Circles are exact inversions. Solid and dashed lines show approximations.
        scatter(axD, curve.Nplot, curve.lambdaFit, parameterMarkerSize, ...
            'Marker', 'o', ...
            'MarkerEdgeColor', markerEdgeColor, ...
            'MarkerFaceColor', popColor, ...
            'MarkerFaceAlpha', markerFaceAlpha, ...
            'LineWidth', markerLineWidth)
        plot(axD, curve.Nplot, curve.lambdaTheory, ...
            'Color', popColor, ...
            'LineWidth', theoryLineWidth)
        plot(axD, [min(curve.Nplot), max(curve.Nplot)], curve.lambdaInfinity * [1, 1], ...
            '--', ...
            'Color', popColor, ...
            'LineWidth', infinityLineWidth)
    end

    if datasetIdx == 1
        exactLegendHandle = scatter(axD, nan, nan, parameterMarkerSize, ...
            'Marker', 'o', ...
            'MarkerEdgeColor', markerEdgeColor, ...
            'MarkerFaceColor', [0.65 0.65 0.65], ...
            'LineWidth', markerLineWidth);
        doubleWellLegendHandle = plot(axD, nan, nan, ...
            'k-', ...
            'LineWidth', theoryLineWidth);
        infinityLegendHandle = plot(axD, nan, nan, ...
            'k--', ...
            'LineWidth', infinityLineWidth);

        legend(axD, [exactLegendHandle, doubleWellLegendHandle, infinityLegendHandle], ...
            {'Exact', 'Double-well', 'N \rightarrow \infty'}, ...
            'Interpreter', 'tex', ...
            'FontName', textFontName, ...
            'FontSize', legendFontSize, ...
            'Location', 'best')
    end

    ylabel(axD, [labelFont 'Interaction strength \lambda'], ...
        'Interpreter', 'tex', ...
        'FontSize', parameterAxisLabelFontSize)

    formatParameterAxis(axD, textFontName, parameterTickLabelFontSize, xTickValues, xTickLabels)
    xlim(axD, parameterXLimits)
    ylim(axD, panelDYLimits)
    yticks(axD, panelDYTicks)
    yticklabels(axD, arrayfun(@num2str, panelDYTicks, 'UniformOutput', false))

    if datasetIdx == length(files)
        xlabel(axD, [labelFont 'Number of neurons {\itN}'], ...
            'Interpreter', 'tex', ...
            'FontSize', parameterAxisLabelFontSize)
    else
        xticklabels(axD, repmat({''}, size(get(axD, 'XTick'))))
    end
end

set(parameterFigure, 'Renderer', 'painters')

%% Separate combined c/d comparison figure

combinedParameterFigure = figure('Name', 'Figure 5 combined panels c and d', ...
    'Color', 'w');

combinedCPosition = [0.13 0.57 0.78 0.34];
combinedDPosition = [0.13 0.12 0.78 0.34];
combinedCYLimits = [-10, -1e-4];
combinedCYTicks = [-10, -1, -1e-1, -1e-2, -1e-3, -1e-4];
combinedCYTickLabels = {'-10^{1}', '-10^{0}', '-10^{-1}', ...
    '-10^{-2}', '-10^{-3}', '-10^{-4}'};
combinedNValues = [];

% Combined c panel: all external-field curves on one axis.
axCCombined = axes('Parent', combinedParameterFigure, ...
    'Position', combinedCPosition);
hold(axCCombined, 'on')
fprintf('Figure 5: plotting combined panel c.\n')

for datasetIdx = 1:length(files)
    datasetCurveData = cacheData.parameterCurves(datasetIdx);

    for popIdx = 1:numel(datasetCurveData.populationCurves)
        curve = datasetCurveData.populationCurves(popIdx);
        popColor = curve.popColor;
        combinedNValues = [combinedNValues; curve.Nplot(:)]; %#ok<AGROW>

        scatter(axCCombined, curve.Nplot, curve.hFit, combinedParameterMarkerSize, ...
            'Marker', 'o', ...
            'MarkerEdgeColor', markerEdgeColor, ...
            'MarkerFaceColor', popColor, ...
            'MarkerFaceAlpha', markerFaceAlpha, ...
            'LineWidth', markerLineWidth)
        plot(axCCombined, curve.Nplot, curve.hTheory, ...
            'Color', popColor, ...
            'LineWidth', theoryLineWidth)

        if ~plotExternalFieldLogY
            plot(axCCombined, [min(curve.Nplot), max(curve.Nplot)], ...
                curve.hInfinity * [1, 1], ...
                '--', ...
                'Color', popColor, ...
                'LineWidth', infinityLineWidth)
        end
    end
end

ylabel(axCCombined, [labelFont 'External field {\ith}'], ...
    'Interpreter', 'tex', ...
    'FontSize', combinedParameterAxisLabelFontSize)
formatParameterAxis(axCCombined, textFontName, combinedParameterTickLabelFontSize, ...
    xTickValues, xTickLabels)

if plotExternalFieldLogY
    set(axCCombined, 'YScale', 'log')
    ylim(axCCombined, combinedCYLimits)
    yticks(axCCombined, combinedCYTicks)
    yticklabels(axCCombined, combinedCYTickLabels)
end

% Combined d panel: all interaction-strength curves on one axis.
axDCombined = axes('Parent', combinedParameterFigure, ...
    'Position', combinedDPosition);
hold(axDCombined, 'on')
fprintf('Figure 5: plotting combined panel d.\n')

for datasetIdx = 1:length(files)
    datasetCurveData = cacheData.parameterCurves(datasetIdx);

    for popIdx = 1:numel(datasetCurveData.populationCurves)
        curve = datasetCurveData.populationCurves(popIdx);
        popColor = curve.popColor;

        scatter(axDCombined, curve.Nplot, curve.lambdaFit, combinedParameterMarkerSize, ...
            'Marker', 'o', ...
            'MarkerEdgeColor', markerEdgeColor, ...
            'MarkerFaceColor', popColor, ...
            'MarkerFaceAlpha', markerFaceAlpha, ...
            'LineWidth', markerLineWidth)
        plot(axDCombined, curve.Nplot, curve.lambdaTheory, ...
            'Color', popColor, ...
            'LineWidth', theoryLineWidth)
        plot(axDCombined, [min(curve.Nplot), max(curve.Nplot)], ...
            curve.lambdaInfinity * [1, 1], ...
            '--', ...
            'Color', popColor, ...
            'LineWidth', infinityLineWidth)
    end
end

exactLegendHandle = scatter(axDCombined, nan, nan, combinedParameterMarkerSize, ...
    'Marker', 'o', ...
    'MarkerEdgeColor', markerEdgeColor, ...
    'MarkerFaceColor', [0.65 0.65 0.65], ...
    'LineWidth', markerLineWidth);
doubleWellLegendHandle = plot(axDCombined, nan, nan, ...
    'k-', ...
    'LineWidth', theoryLineWidth);
infinityLegendHandle = plot(axDCombined, nan, nan, ...
    'k--', ...
    'LineWidth', infinityLineWidth);

legend(axDCombined, [exactLegendHandle, doubleWellLegendHandle, infinityLegendHandle], ...
    {'Exact', 'Double-well', 'N \rightarrow \infty'}, ...
    'Interpreter', 'tex', ...
    'FontName', textFontName, ...
    'FontSize', legendFontSize, ...
    'Location', 'best')

ylabel(axDCombined, [labelFont 'Interaction strength \lambda'], ...
    'Interpreter', 'tex', ...
    'FontSize', combinedParameterAxisLabelFontSize)
xlabel(axDCombined, [labelFont 'Number of neurons {\itN}'], ...
    'Interpreter', 'tex', ...
    'FontSize', combinedParameterAxisLabelFontSize)
formatParameterAxis(axDCombined, textFontName, combinedParameterTickLabelFontSize, ...
    xTickValues, xTickLabels)
ylim(axDCombined, panelDYLimits)
yticks(axDCombined, panelDYTicks)
yticklabels(axDCombined, arrayfun(@num2str, panelDYTicks, 'UniformOutput', false))

combinedNValues = combinedNValues(isfinite(combinedNValues) & combinedNValues > 0);
if ~isempty(combinedNValues)
    xlim(axCCombined, [min(combinedNValues), max(combinedNValues)])
    xlim(axDCombined, [min(combinedNValues), max(combinedNValues)])
end

set(combinedParameterFigure, 'Renderer', 'painters')

%% e. Full-statistic double-well trajectories in h-lambda space

axE = axes('Parent', mainFigure, 'Position', panelEPosition);
hold(axE, 'on')
fprintf('Figure 5: plotting panel e.\n')

plot(axE, [trajectoryHZeroPlot, trajectoryHZeroPlot], ...
    [1, trajectoryLambdaYLim(2)], ...
    'r', ...
    'LineWidth', trajectoryCriticalLineWidth)
plot(axE, linspace(trajectoryHXLim(1), trajectoryHZeroPlot, 500), ...
    zeros(1, 500), ...
    'Color', trajectoryReferenceLineColor, ...
    'LineWidth', 2)

for datasetIdx = 1:length(files)
    [Data, ~, ~] = loadDataFile(files(datasetIdx));
    [~, mus, chis] = summarizeData(Data);

    fullMus = lastofarray(mus);
    fullChis = lastofarray(chis);

    fullMus = fullMus(:).';
    fullChis = fullChis(:).';

    validTargets = isfinite(fullMus) ...
        & isfinite(fullChis) ...
        & fullMus < 0 ...
        & abs(fullMus) < 1 ...
        & fullChis > 0;

    baseColor = hex2rgb(datasetColors(datasetIdx));
    lowNColor = baseColor + (1 - baseColor) * trajectoryLowNLightenAmount;
    highNColor = trajectoryHighNColor;

    fprintf('  panel e: plotting %d full %s recordings.\n', ...
        nnz(validTargets), char(datasetNames(datasetIdx)))

    for popIdx = find(validTargets)
        mTarget = fullMus(popIdx);
        chiTarget = fullChis(popIdx);

        u0 = sqrt(chiTarget ./ trajectoryNValues + mTarget.^2);
        hTheory = nan(size(trajectoryNValues));
        lambdaTheory = nan(size(trajectoryNValues));

        goodTheory = isfinite(u0) ...
            & u0 > 0 ...
            & u0 < 1 ...
            & abs(mTarget ./ u0) < 1;

        lambdaTheory(goodTheory) = atanh(u0(goodTheory)) ./ u0(goodTheory);
        hTheory(goodTheory) = atanh(mTarget ./ u0(goodTheory)) ...
            ./ (trajectoryNValues(goodTheory) .* u0(goodTheory));

        goodPlot = isfinite(hTheory) ...
            & isfinite(lambdaTheory) ...
            & hTheory < 0 ...
            & lambdaTheory > 0;

        for nIdx = 1:length(trajectoryNValues) - 1
            if goodPlot(nIdx) && goodPlot(nIdx + 1)
                colorWeight = (log10(mean(trajectoryNValues(nIdx:nIdx + 1))) ...
                    - log10(trajectoryNMin)) ...
                    / (log10(trajectoryNMax) - log10(trajectoryNMin));
                colorWeight = max(0, min(1, colorWeight));

                if colorWeight <= trajectoryColorMidpoint
                    colorWeightLocal = colorWeight / trajectoryColorMidpoint;
                    currentColor = (1 - colorWeightLocal) * lowNColor ...
                        + colorWeightLocal * baseColor;
                else
                    colorWeightLocal = (colorWeight - trajectoryColorMidpoint) ...
                        / (1 - trajectoryColorMidpoint);
                    currentColor = (1 - colorWeightLocal) * baseColor ...
                        + colorWeightLocal * highNColor;
                end

                plot(axE, hTheory(nIdx:nIdx + 1), ...
                    lambdaTheory(nIdx:nIdx + 1), ...
                    'Color', currentColor, ...
                    'LineWidth', trajectoryLineWidth)
            end
        end
    end

end

scatter(axE, trajectoryHZeroPlot, 1, trajectoryCriticalPointSize, ...
    'Marker', 'o', ...
    'MarkerFaceColor', 'r', ...
    'MarkerEdgeColor', 'r')

axis(axE, 'square')
set(axE, 'XScale', 'log')
xlim(axE, trajectoryHXLim)
ylim(axE, trajectoryLambdaYLim)
xticks(axE, trajectoryHTickValues)
xticklabels(axE, trajectoryHTickLabels)
xtickangle(axE, 0)
yticks(axE, trajectoryLambdaTicks)
yticklabels(axE, arrayfun(@num2str, trajectoryLambdaTicks, 'UniformOutput', false))

xlabel(axE, [labelFont 'External field {\ith}'], ...
    'Interpreter', 'tex', ...
    'FontSize', axisLabelFontSize)
ylabel(axE, [labelFont 'Interaction strength \lambda'], ...
    'Interpreter', 'tex', ...
    'FontSize', axisLabelFontSize)

for datasetIdx = 1:length(files)
    baseColor = hex2rgb(datasetColors(datasetIdx));
    lowNColor = baseColor + (1 - baseColor) * trajectoryLowNLightenAmount;
    highNColor = trajectoryHighNColor;

    for segmentIdx = 1:trajectoryLegendSegmentCount
        colorWeight = (segmentIdx - 0.5) / trajectoryLegendSegmentCount;

        if colorWeight <= trajectoryColorMidpoint
            colorWeightLocal = colorWeight / trajectoryColorMidpoint;
            legendColor = (1 - colorWeightLocal) * lowNColor ...
                + colorWeightLocal * baseColor;
        else
            colorWeightLocal = (colorWeight - trajectoryColorMidpoint) ...
                / (1 - trajectoryColorMidpoint);
            legendColor = (1 - colorWeightLocal) * baseColor ...
                + colorWeightLocal * highNColor;
        end

        plot(axE, trajectoryLegendX(segmentIdx:segmentIdx + 1), ...
            trajectoryLegendY(datasetIdx) * [1, 1], ...
            'Color', legendColor, ...
            'LineWidth', trajectoryNLegendLineWidth)
    end
end

set(axE, 'TickDir', 'both')
axE.FontName = textFontName;
axE.FontSize = tickLabelFontSize;
axE.TickLabelInterpreter = 'tex';
box(axE, 'on')

set(mainFigure, 'Renderer', 'painters')

%% Separate linear h-axis version of panel e

figure('Name', 'Figure 5 panel e linear inset', ...
    'Color', 'w')
hold on
fprintf('Figure 5: plotting separate linear panel e inset.\n')

plot([0, 0], [1, trajectoryLambdaYLim(2)], ...
    'r', ...
    'LineWidth', trajectoryCriticalLineWidth)
plot(linspace(linearTrajectoryHXLim(1), linearTrajectoryHXLim(2), 500), ...
    zeros(1, 500), ...
    'Color', trajectoryReferenceLineColor, ...
    'LineWidth', 2)

for datasetIdx = 1:length(files)
    [Data, ~, ~] = loadDataFile(files(datasetIdx));
    [~, mus, chis] = summarizeData(Data);

    fullMus = lastofarray(mus);
    fullChis = lastofarray(chis);

    fullMus = fullMus(:).';
    fullChis = fullChis(:).';

    validTargets = isfinite(fullMus) ...
        & isfinite(fullChis) ...
        & fullMus < 0 ...
        & abs(fullMus) < 1 ...
        & fullChis > 0;

    baseColor = hex2rgb(datasetColors(datasetIdx));
    lowNColor = baseColor + (1 - baseColor) * trajectoryLowNLightenAmount;
    highNColor = trajectoryHighNColor;

    fprintf('  linear panel e: plotting %d full %s recordings.\n', ...
        nnz(validTargets), char(datasetNames(datasetIdx)))

    for popIdx = find(validTargets)
        mTarget = fullMus(popIdx);
        chiTarget = fullChis(popIdx);

        u0 = sqrt(chiTarget ./ trajectoryNValues + mTarget.^2);
        hTheory = nan(size(trajectoryNValues));
        lambdaTheory = nan(size(trajectoryNValues));

        goodTheory = isfinite(u0) ...
            & u0 > 0 ...
            & u0 < 1 ...
            & abs(mTarget ./ u0) < 1;

        lambdaTheory(goodTheory) = atanh(u0(goodTheory)) ./ u0(goodTheory);
        hTheory(goodTheory) = atanh(mTarget ./ u0(goodTheory)) ...
            ./ (trajectoryNValues(goodTheory) .* u0(goodTheory));

        goodPlot = isfinite(hTheory) ...
            & isfinite(lambdaTheory) ...
            & hTheory < 0 ...
            & lambdaTheory > 0;

        for nIdx = 1:length(trajectoryNValues) - 1
            if goodPlot(nIdx) && goodPlot(nIdx + 1)
                colorWeight = (log10(mean(trajectoryNValues(nIdx:nIdx + 1))) ...
                    - log10(trajectoryNMin)) ...
                    / (log10(trajectoryNMax) - log10(trajectoryNMin));
                colorWeight = max(0, min(1, colorWeight));

                if colorWeight <= trajectoryColorMidpoint
                    colorWeightLocal = colorWeight / trajectoryColorMidpoint;
                    currentColor = (1 - colorWeightLocal) * lowNColor ...
                        + colorWeightLocal * baseColor;
                else
                    colorWeightLocal = (colorWeight - trajectoryColorMidpoint) ...
                        / (1 - trajectoryColorMidpoint);
                    currentColor = (1 - colorWeightLocal) * baseColor ...
                        + colorWeightLocal * highNColor;
                end

                plot(hTheory(nIdx:nIdx + 1), ...
                    lambdaTheory(nIdx:nIdx + 1), ...
                    'Color', currentColor, ...
                    'LineWidth', trajectoryLineWidth)
            end
        end
    end
end

scatter(0, 1, trajectoryCriticalPointSize, ...
    'Marker', 'o', ...
    'MarkerFaceColor', 'r', ...
    'MarkerEdgeColor', 'r')

axis square
xlim(linearTrajectoryHXLim)
ylim(linearTrajectoryLambdaYLim)
xticks(linearTrajectoryHTicks)
yticks(trajectoryLambdaTicks)
yticklabels(arrayfun(@num2str, trajectoryLambdaTicks, 'UniformOutput', false))

xlabel([labelFont 'External field {\ith}'], ...
    'Interpreter', 'tex', ...
    'FontSize', axisLabelFontSize)
ylabel([labelFont 'Interaction strength \lambda'], ...
    'Interpreter', 'tex', ...
    'FontSize', axisLabelFontSize)

set(gca, 'TickDir', 'both')
ax = gca;
ax.FontName = textFontName;
ax.FontSize = tickLabelFontSize;
ax.TickLabelInterpreter = 'tex';
box on
set(gcf, 'Renderer', 'painters')

%% Local functions

function saveFigure5Cache(cacheFile, cacheData, cacheSettings, cacheVersion, isComplete)
    cacheData.isComplete = isComplete;
    save(cacheFile, 'cacheData', 'cacheSettings', 'cacheVersion', '-v7.3')
end

function tf = isFigure5CacheCompatible(cachedSettings, currentSettings)
    tf = true;
    fields = fieldnames(currentSettings);

    for idx = 1:numel(fields)
        fieldName = fields{idx};

        if ~isfield(cachedSettings, fieldName)
            if strcmp(fieldName, 'selectedPopulationByDataset') ...
                    && all(isnan(currentSettings.(fieldName)))
                continue
            end

            tf = false;
            return
        end

        if ~isequaln(cachedSettings.(fieldName), currentSettings.(fieldName))
            tf = false;
            return
        end
    end
end

function panelB = computePanelBData(files, datasetColors, fullMarkerSizeRange)
    fullNRange = fullDataNRange(files);
    panelB = struct('mus', {}, 'chis', {}, 'fullNs', {}, 'markerSizes', {}, 'color', {});

    for i = 1:length(files)
        fprintf('  full-data panel: loading %s.\n', files(i))
        [Data, ~, ~] = loadDataFile(files(i));
        [Nss, mus, chis] = summarizeData(Data);

        fullNs = lastofarray(Nss);
        panelB(i).mus = lastofarray(mus);
        panelB(i).chis = lastofarray(chis);
        panelB(i).fullNs = fullNs;
        panelB(i).markerSizes = markerSizeFromN(fullNs, fullNRange, fullMarkerSizeRange);
        panelB(i).color = hex2rgb(datasetColors(i));
    end
end

function parameterCurves = computeAllParameterCurveData(files, datasetNames, datasetColors, parameterPopulationMode, selectedPopulationByDataset, cacheFile, cacheSettings, cacheVersion, cacheData)
    parameterCurves = cacheData.parameterCurves;

    for datasetIdx = 1:length(files)
        if numel(parameterCurves) >= datasetIdx ...
                && isfield(parameterCurves(datasetIdx), 'populationCurves') ...
                && ~isempty(parameterCurves(datasetIdx).populationCurves)
            fprintf('  panels c/d: using cached %s data.\n', datasetNames(datasetIdx))
            continue
        end

        fprintf('  panels c/d: computing %s (%d of %d).\n', ...
            datasetNames(datasetIdx), datasetIdx, length(files))

        [populationCurves, selectedPopulation, availablePopulations] = computeDatasetParameterCurveData( ...
            files(datasetIdx), datasetColors(datasetIdx), datasetNames(datasetIdx), ...
            parameterPopulationMode, selectedPopulationByDataset(datasetIdx));

        parameterCurves(datasetIdx).datasetName = datasetNames(datasetIdx);
        parameterCurves(datasetIdx).populationCurves = populationCurves;
        parameterCurves(datasetIdx).selectedPopulation = selectedPopulation;
        parameterCurves(datasetIdx).availablePopulations = availablePopulations;

        cacheData.parameterCurves = parameterCurves;
        saveFigure5Cache(cacheFile, cacheData, cacheSettings, cacheVersion, false)
        fprintf('  panels c/d: saved partial cache after %s.\n', datasetNames(datasetIdx))
    end
end

function [populationCurves, selectedPopulation, availablePopulations] = computeDatasetParameterCurveData(fileName, datasetColor, datasetName, parameterPopulationMode, requestedPopulation)
    [Data, ~, ~] = loadDataFile(fileName);
    [Nss, mus, chis] = summarizeData(Data);

    Nss = asPopulationMatrix(Nss);
    mus = asPopulationMatrix(mus);
    chis = asPopulationMatrix(chis);

    nPopulations = size(Nss, 2);
    availablePopulations = nPopulations;
    selectedPopulation = nan;
    sourcePopulationIndices = 1:nPopulations;
    baseColor = hex2rgb(datasetColor);

    if parameterPopulationMode == "single"
        if isfinite(requestedPopulation)
            selectedPopulation = round(requestedPopulation);

            if selectedPopulation < 1 || selectedPopulation > nPopulations
                error('Requested population %d is outside the %s population range 1-%d.', ...
                    selectedPopulation, datasetName, nPopulations)
            end
        else
            counts = sum(isfinite(Nss), 1);
            maxCount = max(counts);
            candidatePopulations = find(counts == maxCount);
            selectedPopulation = candidatePopulations(randi(numel(candidatePopulations)));
            fprintf('    panels c/d: %s has %d populations tied with %d valid N values.\n', ...
                datasetName, numel(candidatePopulations), maxCount)
        end

        fprintf('    panels c/d: selected %s population %d of %d.\n', ...
            datasetName, selectedPopulation, nPopulations)

        Nss = Nss(:, selectedPopulation);
        mus = mus(:, min(selectedPopulation, size(mus, 2)));
        chis = chis(:, min(selectedPopulation, size(chis, 2)));
        nPopulations = 1;
        sourcePopulationIndices = selectedPopulation;
    end

    populationCurves = struct( ...
        'Nplot', {}, ...
        'hFit', {}, ...
        'lambdaFit', {}, ...
        'hTheory', {}, ...
        'lambdaTheory', {}, ...
        'hInfinity', {}, ...
        'lambdaInfinity', {}, ...
        'sourcePopulation', {}, ...
        'popColor', {});

    for popIdx = 1:nPopulations
        fprintf('    panels c/d: %s population %d of %d.\n', ...
            datasetName, sourcePopulationIndices(popIdx), availablePopulations)

        Nuse = Nss(:, popIdx);
        muse = mus(:, min(popIdx, size(mus, 2)));
        chiuse = chis(:, min(popIdx, size(chis, 2)));

        valid = isfinite(Nuse) & isfinite(muse) & isfinite(chiuse) ...
            & abs(muse) < 1 & chiuse > 0;

        Nuse = Nuse(valid);
        muse = muse(valid);
        chiuse = chiuse(valid);

        if numel(Nuse) < 2
            fprintf('      skipping: fewer than two valid points.\n')
            continue
        end

        [Nuse, sortIdx] = sort(Nuse);
        muse = muse(sortIdx);
        chiuse = chiuse(sortIdx);

        hFit = nan(size(Nuse));
        lambdaFit = nan(size(Nuse));

        for pointIdx = 1:numel(Nuse)
            fprintf('      hlambda point %d of %d, N = %.0f.\n', ...
                pointIdx, numel(Nuse), Nuse(pointIdx))

            try
                [hFit(pointIdx), lambdaFit(pointIdx)] = hlambda( ...
                    muse(pointIdx), chiuse(pointIdx), round(Nuse(pointIdx)));
            catch
                hFit(pointIdx) = nan;
                lambdaFit(pointIdx) = nan;
            end
        end

        u0 = sqrt(chiuse ./ Nuse + muse.^2);
        hTheory = nan(size(muse));
        lambdaTheory = nan(size(muse));

        goodTheory = isfinite(u0) ...
            & u0 > 0 ...
            & abs(u0) < 1 ...
            & abs(muse ./ u0) < 1;

        lambdaTheory(goodTheory) = atanh(u0(goodTheory)) ./ u0(goodTheory);
        hTheory(goodTheory) = atanh(muse(goodTheory) ./ u0(goodTheory)) ...
            ./ (Nuse(goodTheory) .* u0(goodTheory));

        good = isfinite(hFit) & isfinite(lambdaFit) ...
            & isfinite(hTheory) & isfinite(lambdaTheory);

        if ~any(good)
            fprintf('      skipping: no valid exact/theory pairs.\n')
            continue
        end

        curveIdx = numel(populationCurves) + 1;
        populationCurves(curveIdx).Nplot = Nuse(good);
        populationCurves(curveIdx).hFit = hFit(good);
        populationCurves(curveIdx).lambdaFit = lambdaFit(good);
        populationCurves(curveIdx).hTheory = hTheory(good);
        populationCurves(curveIdx).lambdaTheory = lambdaTheory(good);
        populationCurves(curveIdx).hInfinity = 0;

        finalMu = muse(find(isfinite(muse), 1, 'last'));
        populationCurves(curveIdx).lambdaInfinity = atanh(abs(finalMu)) / max(abs(finalMu), eps);
        populationCurves(curveIdx).sourcePopulation = sourcePopulationIndices(popIdx);
        populationCurves(curveIdx).popColor = baseColor;
    end
end

function A = asPopulationMatrix(A)
    if isvector(A)
        A = A(:);
    end
end

function formatSquareAxis(ax, textFontName, tickLabelFontSize)
    axis(ax, 'square')
    box(ax, 'on')
    set(ax, 'TickDir', 'both')
    ax.FontName = textFontName;
    ax.FontSize = tickLabelFontSize;
    ax.LineWidth = 1.2;
    ax.TickLabelInterpreter = 'tex';
end

function formatSubAxis(ax, textFontName, tickLabelFontSize)
    box(ax, 'on')
    set(ax, 'TickDir', 'both')
    ax.FontName = textFontName;
    ax.FontSize = tickLabelFontSize;
    ax.LineWidth = 1.2;
    ax.TickLabelInterpreter = 'tex';
    xlim(ax, [-1, 1])
end

function formatParameterAxis(ax, textFontName, tickLabelFontSize, xTickValues, xTickLabels)
    box(ax, 'on')
    set(ax, 'TickDir', 'both')
    set(ax, 'XScale', 'log')
    ax.FontName = textFontName;
    ax.FontSize = tickLabelFontSize;
    ax.LineWidth = 1.2;
    ax.XTick = xTickValues;
    ax.XTickLabel = xTickLabels;
    ax.TickLabelInterpreter = 'tex';
end

function setParameterAxisXLim(ax, datasetCurveData)
    if isempty(datasetCurveData.populationCurves)
        return
    end

    allN = [];

    for idx = 1:numel(datasetCurveData.populationCurves)
        allN = [allN; datasetCurveData.populationCurves(idx).Nplot(:)]; %#ok<AGROW>
    end

    allN = allN(isfinite(allN) & allN > 0);

    if isempty(allN)
        return
    end

    xlim(ax, [min(allN), max(allN)])
end

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

function y = logGammaLanczos(x)
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
    y = zeros(size(x));

    maskReflect = x < 0.5;
    xr = x(maskReflect);

    if any(maskReflect)
        y(maskReflect) = log(pi) ...
            - log(abs(sin(pi * xr))) ...
            - logGammaLanczos(1 - xr);
    end

    maskMain = ~maskReflect;
    xm = x(maskMain);

    if any(maskMain)
        xm = xm - 1;
        a = p(1) * ones(size(xm));

        for i = 2:length(p)
            a = a + p(i) ./ (xm + i - 1);
        end

        t = xm + g + 0.5;

        y(maskMain) = 0.5 * log(2 * pi) ...
            + (xm + 0.5) .* log(t) ...
            - t ...
            + log(a);
    end
end
