%% KEEP PLOTS

% Standalone script for the six-panel Figure 5 layout.
% The top row contains the mean-field statistics region, neural data in
% statistics space, and two free-energy panels. The bottom row compares
% exact inferred parameters with the double-well approximation for all
% available populations in each dataset.

clear
clc

plotsDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(plotsDir);
addpath(repoRoot)
addpath(fullfile(repoRoot, 'Allen'))
addpath(fullfile(repoRoot, 'Hippo'))
addpath(fullfile(repoRoot, 'Stringer'))

%% Figure settings

files = ["Allenhldata.mat", "stringerhldata.mat", "hippomuchidata.mat"];
datasetNames = ["Allen", "Stringer", "Hippocampus"];
datasetColors = ["#2676ad", "#2f682c", "#4f4cc4"];

textFontName = 'Helvetica';
labelFont = ['\fontname{' textFontName '}'];
axisLabelFontSize = 22;
subAxisLabelFontSize = 18;
tickLabelFontSize = 14;
legendFontSize = 12;

referenceLineColor = "#939598";
markerSize = 55;
markerEdgeColor = [0 0 0];
markerLineWidth = 0.5;
markerFaceAlpha = 0.75;
fullMarkerSizeRange = [85 240];
statisticsXLim = [-1, -0.8];
statisticsYLim = [0.05, 70];
upperBoundShadeColor = [1, 0, 0];
upperBoundShadeAlpha = 0.15;
upperBoundLineWidth = 2;

freeEnergyNs = [20, 50, 100, 200, 1000];
freeEnergyColors = ["#f9cd6b", "#edb34f", "#c2591c", "#ba340d", "#b41a02"];
freeEnergyInfinityColor = [0 0 0];
freeEnergyLineWidth = 1.4;
freeEnergyMu = -0.962;
freeEnergyChi = 0.24;
fixedParameterN = 1000;
muGrid = linspace(-0.999, 0.999, 1000);

theoryLineWidth = 1.5;
infinityLineWidth = 1.0;
xTickValues = 10.^(0:4);
xTickLabels = arrayfun(@(x) sprintf('10^{%d}', x), 0:4, 'UniformOutput', false);

%% Figure 2 color scheme for parameter-space tiling

tileSideLength = 0.05;
numPtsPerSide = 20;
lambdaValues = -0.001:tileSideLength:3;
hTileValues = -1 + tileSideLength/2:tileSideLength:-tileSideLength/2;
jColorValues = [
    0.5765    0.5843    0.5961
    0.5020    0.1765    0.1765
    0.7020         0         0
    0.8941    0.2039         0
    0.9725    0.5373         0
    1.0000    0.6745    0.0667
    1.0000    0.7745    0.1667
];
hSaturationFloor = 0.15;

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
cacheSettings.tileSideLength = tileSideLength;
cacheSettings.numPtsPerSide = numPtsPerSide;
cacheSettings.lambdaValues = lambdaValues;
cacheSettings.hTileValues = hTileValues;
cacheSettings.jColorValues = jColorValues;
cacheSettings.hSaturationFloor = hSaturationFloor;

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

    if ~isfield(cacheData, 'panelA')
        fprintf('Figure 5: computing panel a mean-field tiles.\n')
        cacheData.panelA = computePanelATileData( ...
            hTileValues, lambdaValues, tileSideLength, numPtsPerSide, ...
            jColorValues, hSaturationFloor);
        saveFigure5Cache(cacheFile, cacheData, cacheSettings, cacheVersion, false)
    else
        fprintf('Figure 5: using cached panel a data.\n')
    end

    if ~isfield(cacheData, 'panelB')
        fprintf('Figure 5: computing panel b full-data points.\n')
        cacheData.panelB = computePanelBData(files, datasetColors, fullMarkerSizeRange);
        saveFigure5Cache(cacheFile, cacheData, cacheSettings, cacheVersion, false)
    else
        fprintf('Figure 5: using cached panel b data.\n')
    end

    fprintf('Figure 5: computing panels e and f parameter curves.\n')
    if ~isfield(cacheData, 'parameterCurves')
        cacheData.parameterCurves = struct('datasetName', {}, 'populationCurves', {});
    end

    cacheData.parameterCurves = computeAllParameterCurveData( ...
        files, datasetNames, datasetColors, ...
        cacheFile, cacheSettings, cacheVersion, cacheData);

    saveFigure5Cache(cacheFile, cacheData, cacheSettings, cacheVersion, true)
    fprintf('Figure 5: saved complete computed data to %s.\n', cacheFile)
end

%% Layout

figure('Name', 'Figure 5 panels', ...
    'Color', 'w', ...
    'Position', [100 100 1500 1000])

mainLayout = tiledlayout(4, 3, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

%% a. Mean-field single-minimum statistics region

axA = nexttile(mainLayout, 1, [2 1]);
hold(axA, 'on')
fprintf('Figure 5: plotting panel a.\n')

for k = 1:numel(cacheData.panelA.mTiles)
    fill(axA, cacheData.panelA.mTiles{k}, cacheData.panelA.chiTiles{k}, ...
        cacheData.panelA.tileColors(k, :), ...
        'EdgeColor', 'none')
end

xlim(axA, statisticsXLim)
ylim(axA, statisticsYLim)
set(axA, 'YScale', 'log')
xlabel(axA, [labelFont 'Average activity {\itm}'], ...
    'Interpreter', 'tex', ...
    'FontSize', axisLabelFontSize)
ylabel(axA, [labelFont 'Susceptibility \chi'], ...
    'Interpreter', 'tex', ...
    'FontSize', axisLabelFontSize)
formatSquareAxis(axA, textFontName, tickLabelFontSize)

%% b. Full neural data relative to the single-minimum region

axB = nexttile(mainLayout, 2, [2 1]);
hold(axB, 'on')
fprintf('Figure 5: plotting panel b.\n')

x = linspace(statisticsXLim(1), statisticsXLim(2), 1000);
independentBound = 1 - x.^2;
upperBound = x .* (1 - x.^2) ./ (x - atanh(x) .* (1 - x.^2));
validShade = isfinite(upperBound) & upperBound > 0 & upperBound < statisticsYLim(2);

fill(axB, [x(validShade), fliplr(x(validShade))], ...
    [upperBound(validShade), statisticsYLim(2) * ones(1, nnz(validShade))], ...
    upperBoundShadeColor, ...
    'FaceAlpha', upperBoundShadeAlpha, ...
    'EdgeColor', 'none')
plot(axB, x, independentBound, ...
    'LineWidth', 2, ...
    'Color', referenceLineColor)
plot(axB, x(validShade), upperBound(validShade), ...
    'LineWidth', upperBoundLineWidth, ...
    'Color', upperBoundShadeColor)

for i = 1:numel(cacheData.panelB)
    scatter(axB, cacheData.panelB(i).mus, cacheData.panelB(i).chis, ...
        cacheData.panelB(i).markerSizes, ...
        'Marker', 'o', ...
        'MarkerEdgeColor', markerEdgeColor, ...
        'MarkerFaceColor', cacheData.panelB(i).color, ...
        'MarkerFaceAlpha', markerFaceAlpha, ...
        'LineWidth', markerLineWidth)
end

xlim(axB, statisticsXLim)
ylim(axB, statisticsYLim)
set(axB, 'YScale', 'log')
xlabel(axB, [labelFont 'Average activity {\itm}'], ...
    'Interpreter', 'tex', ...
    'FontSize', axisLabelFontSize)
ylabel(axB, [labelFont 'Susceptibility \chi'], ...
    'Interpreter', 'tex', ...
    'FontSize', axisLabelFontSize)
formatSquareAxis(axB, textFontName, tickLabelFontSize)

%% c-d. Free-energy panels in the top-right tile

axC = nexttile(mainLayout, 3);
hold(axC, 'on')
fprintf('Figure 5: computing and plotting panel c.\n')

[fixedH, fixedLambda] = hlambda(freeEnergyMu, freeEnergyChi, fixedParameterN);
freeEnergyLabels = strings(1, numel(freeEnergyNs) + 1);

for k = 1:numel(freeEnergyNs)
    N = freeEnergyNs(k);

    score = (N * fixedH * muGrid + N * fixedLambda * muGrid.^2 / 2 ...
        + logGammaLanczos(N + 1) ...
        - logGammaLanczos(N * (1 + muGrid) / 2 + 1) ...
        - logGammaLanczos(N * (1 - muGrid) / 2 + 1)) / N;

    y = normalizedFreeEnergy(score);

    plot(axC, muGrid, y, ...
        'Color', freeEnergyColors(k), ...
        'LineWidth', freeEnergyLineWidth)

    freeEnergyLabels(k) = sprintf('N = %d', N);
end

pPlus = (1 + muGrid) / 2;
pMinus = (1 - muGrid) / 2;
entropyTerm = pPlus .* log(pPlus) + pMinus .* log(pMinus);
scoreInf = fixedLambda * muGrid.^2 / 2 + fixedH * muGrid - entropyTerm;
yInf = normalizedFreeEnergy(scoreInf);

plot(axC, muGrid, yInf, ...
    'Color', freeEnergyInfinityColor, ...
    'LineWidth', freeEnergyLineWidth)
freeEnergyLabels(end) = 'N = \infty';

legend(axC, freeEnergyLabels, ...
    'Interpreter', 'tex', ...
    'FontName', textFontName, ...
    'FontSize', 9, ...
    'Location', 'best')

ylabel(axC, [labelFont 'Free energy {\itf}(\mu)'], ...
    'Interpreter', 'tex', ...
    'FontSize', subAxisLabelFontSize)
formatSubAxis(axC, textFontName, tickLabelFontSize)

axD = nexttile(mainLayout, 6);
hold(axD, 'on')
fprintf('Figure 5: computing and plotting panel d.\n')

freeEnergyLabels = strings(1, numel(freeEnergyNs) + 1);

for k = 1:numel(freeEnergyNs)
    N = freeEnergyNs(k);
    [h, lambda] = hlambda(freeEnergyMu, freeEnergyChi, N);

    score = (N * h * muGrid + N * lambda * muGrid.^2 / 2 ...
        + logGammaLanczos(N + 1) ...
        - logGammaLanczos(N * (1 + muGrid) / 2 + 1) ...
        - logGammaLanczos(N * (1 - muGrid) / 2 + 1)) / N;

    y = normalizedFreeEnergy(score);

    plot(axD, muGrid, y, ...
        'Color', freeEnergyColors(k), ...
        'LineWidth', freeEnergyLineWidth)

    freeEnergyLabels(k) = sprintf('N = %d', N);
end

lambdaInf = atanh(abs(freeEnergyMu)) / max(abs(freeEnergyMu), eps);
pPlus = (1 + muGrid) / 2;
pMinus = (1 - muGrid) / 2;
entropyTerm = pPlus .* log(pPlus) + pMinus .* log(pMinus);
scoreInf = lambdaInf * muGrid.^2 / 2 - entropyTerm;
yInf = normalizedFreeEnergy(scoreInf);

plot(axD, muGrid, yInf, ...
    'Color', freeEnergyInfinityColor, ...
    'LineWidth', freeEnergyLineWidth)
freeEnergyLabels(end) = 'N = \infty';

legend(axD, freeEnergyLabels, ...
    'Interpreter', 'tex', ...
    'FontName', textFontName, ...
    'FontSize', 9, ...
    'Location', 'best')

xlabel(axD, [labelFont 'Activity \mu'], ...
    'Interpreter', 'tex', ...
    'FontSize', subAxisLabelFontSize)
ylabel(axD, [labelFont 'Free energy {\itf}(\mu)'], ...
    'Interpreter', 'tex', ...
    'FontSize', subAxisLabelFontSize)
formatSubAxis(axD, textFontName, tickLabelFontSize)

%% e-f. Exact inference and double-well approximation for all populations

for datasetIdx = 1:length(files)
    fprintf('Figure 5: plotting panels e/f for %s.\n', datasetNames(datasetIdx))
    axH = nexttile(mainLayout, 6 + datasetIdx);
    hold(axH, 'on')

    axLambda = nexttile(mainLayout, 9 + datasetIdx);
    hold(axLambda, 'on')

    plotAllDatasetParameterCurves(axH, axLambda, cacheData.parameterCurves(datasetIdx), ...
        markerSize, markerEdgeColor, markerLineWidth, markerFaceAlpha, ...
        theoryLineWidth, infinityLineWidth)

    title(axH, datasetNames(datasetIdx), ...
        'FontName', textFontName, ...
        'FontSize', tickLabelFontSize)

    ylabel(axH, [labelFont 'External field {\ith}'], ...
        'Interpreter', 'tex', ...
        'FontSize', subAxisLabelFontSize)
    ylabel(axLambda, [labelFont 'Interaction strength \lambda'], ...
        'Interpreter', 'tex', ...
        'FontSize', subAxisLabelFontSize)
    xlabel(axLambda, [labelFont 'Number of neurons {\itN}'], ...
        'Interpreter', 'tex', ...
        'FontSize', subAxisLabelFontSize)

    formatParameterAxis(axH, textFontName, tickLabelFontSize, xTickValues, xTickLabels)
    formatParameterAxis(axLambda, textFontName, tickLabelFontSize, xTickValues, xTickLabels)
    ylim(axLambda, [0, 2.5])
end

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

        if ~isfield(cachedSettings, fieldName) ...
                || ~isequaln(cachedSettings.(fieldName), currentSettings.(fieldName))
            tf = false;
            return
        end
    end
end

function panelA = computePanelATileData(hTileValues, lambdaValues, tileSideLength, numPtsPerSide, jColorValues, hSaturationFloor)
    [centers, edgeVals, tileHalfLength, tileColors] = makeParameterTiles( ...
        hTileValues, lambdaValues, tileSideLength, numPtsPerSide, ...
        jColorValues, hSaturationFloor);

    panelA.mTiles = {};
    panelA.chiTiles = {};
    panelA.tileColors = [];
    nTiles = size(centers, 1);

    for k = 1:nTiles
        if k == 1 || mod(k, 100) == 0 || k == nTiles
            fprintf('  panel a: mean-field tile %d of %d.\n', k, nTiles)
        end

        h0 = centers(k, 1);
        lambda0 = centers(k, 2);

        [hSquare, lambdaSquare] = squareBoundary(h0, lambda0, edgeVals, tileHalfLength, numPtsPerSide);
        [mVals, chiVals] = meanFieldStatistics(hSquare, lambdaSquare);
        valid = isfinite(mVals) & isfinite(chiVals) & chiVals > 0;

        if all(valid)
            panelA.mTiles{end + 1} = mVals; %#ok<AGROW>
            panelA.chiTiles{end + 1} = chiVals; %#ok<AGROW>
            panelA.tileColors(end + 1, :) = tileColors(k, :); %#ok<AGROW>
        end
    end
end

function panelB = computePanelBData(files, datasetColors, fullMarkerSizeRange)
    fullNRange = fullDataNRange(files);
    panelB = struct('mus', {}, 'chis', {}, 'markerSizes', {}, 'color', {});

    for i = 1:length(files)
        fprintf('  panel b: loading full data for %s.\n', files(i))
        [Data, ~, ~] = loadDataFile(files(i));
        [Nss, mus, chis] = summarizeData(Data);

        fullNs = lastofarray(Nss);
        panelB(i).mus = lastofarray(mus);
        panelB(i).chis = lastofarray(chis);
        panelB(i).markerSizes = markerSizeFromN(fullNs, fullNRange, fullMarkerSizeRange);
        panelB(i).color = hex2rgb(datasetColors(i));
    end
end

function parameterCurves = computeAllParameterCurveData(files, datasetNames, datasetColors, cacheFile, cacheSettings, cacheVersion, cacheData)
    parameterCurves = cacheData.parameterCurves;

    for datasetIdx = 1:length(files)
        if numel(parameterCurves) >= datasetIdx ...
                && isfield(parameterCurves(datasetIdx), 'populationCurves') ...
                && ~isempty(parameterCurves(datasetIdx).populationCurves)
            fprintf('  panels e/f: using cached %s data.\n', datasetNames(datasetIdx))
            continue
        end

        fprintf('  panels e/f: computing %s (%d of %d).\n', ...
            datasetNames(datasetIdx), datasetIdx, length(files))

        parameterCurves(datasetIdx).datasetName = datasetNames(datasetIdx);
        parameterCurves(datasetIdx).populationCurves = computeDatasetParameterCurveData( ...
            files(datasetIdx), datasetColors(datasetIdx), datasetNames(datasetIdx));

        cacheData.parameterCurves = parameterCurves;
        saveFigure5Cache(cacheFile, cacheData, cacheSettings, cacheVersion, false)
        fprintf('  panels e/f: saved partial cache after %s.\n', datasetNames(datasetIdx))
    end
end

function populationCurves = computeDatasetParameterCurveData(fileName, datasetColor, datasetName)
    [Data, ~, ~] = loadDataFile(fileName);
    [Nss, mus, chis] = summarizeData(Data);

    Nss = asPopulationMatrix(Nss);
    mus = asPopulationMatrix(mus);
    chis = asPopulationMatrix(chis);

    nPopulations = size(Nss, 2);
    baseColor = hex2rgb(datasetColor);
    populationCurves = struct( ...
        'Nplot', {}, ...
        'hFit', {}, ...
        'lambdaFit', {}, ...
        'hTheory', {}, ...
        'lambdaTheory', {}, ...
        'hInfinity', {}, ...
        'lambdaInfinity', {}, ...
        'popColor', {});

    for popIdx = 1:nPopulations
        fprintf('    panels e/f: %s population %d of %d.\n', ...
            datasetName, popIdx, nPopulations)

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
        populationCurves(curveIdx).popColor = baseColor;
    end
end

function [centers, edgeVals, tileHalfLength, tileColors] = makeParameterTiles(hTileValues, lambdaValues, tileSideLength, numPtsPerSide, jColorValues, hSaturationFloor)
    centers = [];

    for lambda0 = lambdaValues
        for h0 = hTileValues
            centers = [centers; [h0, lambda0]];
        end
    end

    tileHalfLength = tileSideLength / 2;
    edgeVals = linspace(-tileHalfLength, tileHalfLength, numPtsPerSide);

    centerHs = centers(:, 1);
    centerLambdas = centers(:, 2);
    hNorm = (centerHs - min(centerHs)) / (max(centerHs) - min(centerHs));
    lambdaNorm = (centerLambdas - min(centerLambdas)) / (max(centerLambdas) - min(centerLambdas));
    colorPositions = linspace(0, 1, size(jColorValues, 1));

    tileColors = zeros(size(centers, 1), 3);

    for k = 1:size(centers, 1)
        lambdaColor = interp1(colorPositions, jColorValues, lambdaNorm(k), 'linear');
        hsvColor = rgb2hsv(lambdaColor);
        saturationScale = hSaturationFloor + (1 - hSaturationFloor) * hNorm(k);
        hsvColor(2) = hsvColor(2) * saturationScale;
        tileColors(k, :) = hsv2rgb(hsvColor);
    end
end

function [hSquare, lambdaSquare] = squareBoundary(h0, lambda0, edgeVals, tileHalfLength, numPtsPerSide)
    hSquare = [h0 + edgeVals, ...
        h0 + tileHalfLength * ones(1, numPtsPerSide), ...
        h0 + fliplr(edgeVals), ...
        h0 - tileHalfLength * ones(1, numPtsPerSide)];

    lambdaSquare = [lambda0 - tileHalfLength * ones(1, numPtsPerSide), ...
        lambda0 + edgeVals, ...
        lambda0 + tileHalfLength * ones(1, numPtsPerSide), ...
        lambda0 + fliplr(edgeVals)];
end

function [mVals, chiVals] = meanFieldStatistics(hVals, lambdaVals)
    mVals = nan(size(hVals));
    chiVals = nan(size(hVals));

    for idx = 1:numel(hVals)
        h = hVals(idx);
        lambda = lambdaVals(idx);
        m = solveMeanFieldBranch(h, lambda);
        denominator = 1 - lambda * (1 - m^2);
        chi = (1 - m^2) / denominator;

        if isfinite(m) && isfinite(chi) && denominator > 0
            mVals(idx) = m;
            chiVals(idx) = chi;
        end
    end
end

function m = solveMeanFieldBranch(h, lambda)
    if abs(h) < 1e-12
        m = 0;
        return
    end

    if h > 0
        interval = [0, 0.999999];
    else
        interval = [-0.999999, 0];
    end

    f = @(m) m - tanh(h + lambda * m);

    try
        m = fzero(f, interval);
    catch
        m = fzero(f, sign(h) * 0.5);
    end
end

function y = normalizedFreeEnergy(score)
    peaks = findpeaks(score);

    if numel(peaks) >= 2
        topPeaks = maxk(peaks, 2);
        referenceValue = min(topPeaks);
    else
        referenceValue = max(score);
    end

    y = -(score - referenceValue);
end

function plotAllDatasetParameterCurves(axH, axLambda, datasetCurveData, markerSize, markerEdgeColor, markerLineWidth, markerFaceAlpha, theoryLineWidth, infinityLineWidth)
    for popIdx = 1:numel(datasetCurveData.populationCurves)
        curve = datasetCurveData.populationCurves(popIdx);
        popColor = curve.popColor;

        scatter(axH, curve.Nplot, curve.hFit, markerSize, ...
            'Marker', 'o', ...
            'MarkerEdgeColor', markerEdgeColor, ...
            'MarkerFaceColor', popColor, ...
            'MarkerFaceAlpha', markerFaceAlpha, ...
            'LineWidth', markerLineWidth)
        plot(axH, curve.Nplot, curve.hTheory, ...
            'Color', popColor, ...
            'LineWidth', theoryLineWidth)

        scatter(axLambda, curve.Nplot, curve.lambdaFit, markerSize, ...
            'Marker', 'o', ...
            'MarkerEdgeColor', markerEdgeColor, ...
            'MarkerFaceColor', popColor, ...
            'MarkerFaceAlpha', markerFaceAlpha, ...
            'LineWidth', markerLineWidth)
        plot(axLambda, curve.Nplot, curve.lambdaTheory, ...
            'Color', popColor, ...
            'LineWidth', theoryLineWidth)

        plot(axH, [min(curve.Nplot), max(curve.Nplot)], curve.hInfinity * [1, 1], ...
            '--', ...
            'Color', popColor, ...
            'LineWidth', infinityLineWidth)
        plot(axLambda, [min(curve.Nplot), max(curve.Nplot)], curve.lambdaInfinity * [1, 1], ...
            '--', ...
            'Color', popColor, ...
            'LineWidth', infinityLineWidth)
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
