%% KEEP PLOTS

% Recompute full-data m and chi after independently shuffling each neuron's
% time series. This removes zero-lag covariance while preserving each
% neuron's mean activity.

clear
clc

plotsDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(plotsDir);
addpath(repoRoot)
addpath(fullfile(repoRoot, 'Allen'))
addpath(fullfile(repoRoot, 'Hippo'))
addpath(fullfile(repoRoot, 'Stringer'))

%% Figure settings

datasetNames = ["Allen", "Hippocampus", "Stringer"];
datasetColors = ["#2676ad", "#4f4cc4", "#2f682c"];

textFontName = 'Helvetica';
labelFont = ['\fontname{' textFontName '}'];
axisLabelFontSize = 32;
tickLabelFontSize = 18;

dataMarker = 'o';
markerEdgeColor = [0 0 0];
markerLineWidth = 0.5;
markerFaceAlpha = 0.75;
fullMarkerSizeRange = [85 240];
referenceLineColor = "#939598";
upperBoundShadeColor = [1, 0, 0];
upperBoundShadeAlpha = 0.15;
upperBoundLineWidth = 2;

statisticsXLim = [-1, -0.8];
statisticsYLim = [0.05, 70];
modelYLim = [-1, 2.5];
hZeroPlot = -0.00005;
criticalPointSize = 85;

forceRecompute = false;
useParallelProcessing = true;
numParallelWorkers = 10;  % Use [] to keep an existing/default pool size.
cacheVersion = 2;
cacheFile = fullfile(plotsDir, 'ShuffledFullDataFigure1_cache.mat');

%% Compute or load shuffled full-data statistics

if ~forceRecompute && exist(cacheFile, 'file') == 2
    cached = load(cacheFile, 'shuffledStats', 'cacheVersion');

    if isfield(cached, 'cacheVersion') && cached.cacheVersion == cacheVersion
        shuffledStats = cached.shuffledStats;
        fprintf('Loaded shuffled full-data cache from %s\n', cacheFile)
    else
        fprintf('Cache version changed. Recomputing shuffled full-data statistics.\n')
        shuffledStats = computeShuffledFullDataStats( ...
            repoRoot, datasetNames, datasetColors, ...
            useParallelProcessing, numParallelWorkers);
        save(cacheFile, 'shuffledStats', 'cacheVersion', '-v7.3')
    end
else
    fprintf('Computing shuffled full-data statistics.\n')
    shuffledStats = computeShuffledFullDataStats( ...
        repoRoot, datasetNames, datasetColors, ...
        useParallelProcessing, numParallelWorkers);
    save(cacheFile, 'shuffledStats', 'cacheVersion', '-v7.3')
    fprintf('Saved shuffled full-data cache to %s\n', cacheFile)
end

shuffledStats = applyDatasetColors(shuffledStats, datasetNames, datasetColors);
[shuffledStats, addedIndependentFits] = addIndependentModelFits(shuffledStats);

if addedIndependentFits
    save(cacheFile, 'shuffledStats', 'cacheVersion', '-v7.3')
    fprintf('Saved independent h/lambda fits to %s\n', cacheFile)
end

allNs = [shuffledStats.N];
fullNRange = [min(allNs), max(allNs)];

%% Plot shuffled full-data m/chi and h/lambda panels

figure('Name', 'Shuffled full-data Figure 1 panels', ...
    'Color', 'w')
panelLayout = tiledlayout(1, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

%% Shuffled full-data m and chi

axStats = nexttile(panelLayout, 1);
hold(axStats, 'on')

x = linspace(statisticsXLim(1), statisticsXLim(2), 1000);
independentBound = 1 - x.^2;
upperBound = x.*(1 - x.^2)./(x - atanh(x).*(1 - x.^2));
validShade = isfinite(upperBound) & upperBound > 0 & upperBound < statisticsYLim(2);

fill(axStats, [x(validShade), fliplr(x(validShade))], ...
    [upperBound(validShade), statisticsYLim(2) * ones(1, nnz(validShade))], ...
    upperBoundShadeColor, ...
    'FaceAlpha', upperBoundShadeAlpha, ...
    'EdgeColor', 'none')
plot(axStats, x, independentBound, ...
    'LineWidth', 2, ...
    'Color', referenceLineColor)
plot(axStats, x(validShade), upperBound(validShade), ...
    'LineWidth', upperBoundLineWidth, ...
    'Color', upperBoundShadeColor)

for i = 1:numel(shuffledStats)
    markerSize = markerSizeFromN(shuffledStats(i).N, fullNRange, fullMarkerSizeRange);
    scatter(axStats, shuffledStats(i).m, shuffledStats(i).chi, markerSize, ...
        'Marker', dataMarker, ...
        'MarkerEdgeColor', markerEdgeColor, ...
        'MarkerFaceColor', shuffledStats(i).color, ...
        'MarkerFaceAlpha', markerFaceAlpha, ...
        'LineWidth', markerLineWidth)
end

xlim(axStats, statisticsXLim)
ylim(axStats, statisticsYLim)
axis(axStats, 'square')
set(axStats, 'YScale', 'log')
xlabel(axStats, [labelFont 'Average activity {\itm}'], ...
    'Interpreter', 'tex', ...
    'FontSize', axisLabelFontSize)
ylabel(axStats, [labelFont 'Correlation \chi'], ...
    'Interpreter', 'tex', ...
    'FontSize', axisLabelFontSize)
formatFigureAxis(axStats, textFontName, tickLabelFontSize)

%% Shuffled full-data h and lambda

axModel = nexttile(panelLayout, 2);
hold(axModel, 'on')

plot(axModel, hZeroPlot + 0 * linspace(1, modelYLim(2)), ...
    linspace(1, modelYLim(2)), ...
    'r', ...
    'LineWidth', 2)
plot(axModel, linspace(-2.5, hZeroPlot), ...
    0 * linspace(1, modelYLim(2)), ...
    'Color', referenceLineColor, ...
    'LineWidth', 2)

for i = 1:numel(shuffledStats)
    markerSize = markerSizeFromN(shuffledStats(i).N, fullNRange, fullMarkerSizeRange);
    scatter(axModel, shuffledStats(i).h, shuffledStats(i).lambda, markerSize, ...
        'Marker', dataMarker, ...
        'MarkerEdgeColor', markerEdgeColor, ...
        'MarkerFaceColor', shuffledStats(i).color, ...
        'MarkerFaceAlpha', markerFaceAlpha, ...
        'LineWidth', markerLineWidth)
end

scatter(axModel, hZeroPlot, 1, criticalPointSize, ...
    'Marker', dataMarker, ...
    'MarkerFaceColor', 'r', ...
    'MarkerEdgeColor', 'r')

axis(axModel, 'square')
xlim(axModel, [-2.5, 0])
ylim(axModel, modelYLim)
% xscale(axModel, 'log')
% xticks(axModel, [-1, -0.1, -0.01, -0.001, -0.0001, hZeroPlot])
% xticklabels(axModel, {'-10^{0}', '-10^{-1}', '-10^{-2}', '-10^{-3}', ...
%     '-10^{-4}', '    0'})
% xtickangle(axModel, 0)
xlabel(axModel, [labelFont 'External field {\ith}'], ...
    'Interpreter', 'tex', ...
    'FontSize', axisLabelFontSize)
ylabel(axModel, [labelFont 'Interaction strength \lambda'], ...
    'Interpreter', 'tex', ...
    'FontSize', axisLabelFontSize)
formatFigureAxis(axModel, textFontName, tickLabelFontSize)

set(gcf, 'Renderer', 'painters')

%% Plot independent full-data m/chi and h/lambda panels

figure('Name', 'Independent full-data Figure 1 panels', ...
    'Color', 'w')
independentLayout = tiledlayout(1, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

%% Independent full-data m and chi

axIndependentStats = nexttile(independentLayout, 1);
hold(axIndependentStats, 'on')

x = linspace(statisticsXLim(1), statisticsXLim(2), 1000);
independentBound = 1 - x.^2;
upperBound = x.*(1 - x.^2)./(x - atanh(x).*(1 - x.^2));
validShade = isfinite(upperBound) & upperBound > 0 & upperBound < statisticsYLim(2);

fill(axIndependentStats, [x(validShade), fliplr(x(validShade))], ...
    [upperBound(validShade), statisticsYLim(2) * ones(1, nnz(validShade))], ...
    upperBoundShadeColor, ...
    'FaceAlpha', upperBoundShadeAlpha, ...
    'EdgeColor', 'none')
plot(axIndependentStats, x, independentBound, ...
    'LineWidth', 2, ...
    'Color', referenceLineColor)
plot(axIndependentStats, x(validShade), upperBound(validShade), ...
    'LineWidth', upperBoundLineWidth, ...
    'Color', upperBoundShadeColor)

for i = 1:numel(shuffledStats)
    markerSize = markerSizeFromN(shuffledStats(i).N, fullNRange, fullMarkerSizeRange);
    scatter(axIndependentStats, shuffledStats(i).m, shuffledStats(i).chiIndependent, markerSize, ...
        'Marker', dataMarker, ...
        'MarkerEdgeColor', markerEdgeColor, ...
        'MarkerFaceColor', shuffledStats(i).color, ...
        'MarkerFaceAlpha', markerFaceAlpha, ...
        'LineWidth', markerLineWidth)
end

xlim(axIndependentStats, statisticsXLim)
ylim(axIndependentStats, statisticsYLim)
axis(axIndependentStats, 'square')
set(axIndependentStats, 'YScale', 'log')
xlabel(axIndependentStats, [labelFont 'Average activity {\itm}'], ...
    'Interpreter', 'tex', ...
    'FontSize', axisLabelFontSize)
ylabel(axIndependentStats, [labelFont 'Correlation \chi'], ...
    'Interpreter', 'tex', ...
    'FontSize', axisLabelFontSize)
formatFigureAxis(axIndependentStats, textFontName, tickLabelFontSize)

%% Independent full-data h and lambda

axIndependentModel = nexttile(independentLayout, 2);
hold(axIndependentModel, 'on')

plot(axIndependentModel, hZeroPlot + 0 * linspace(1, modelYLim(2)), ...
    linspace(1, modelYLim(2)), ...
    'r', ...
    'LineWidth', 2)
plot(axIndependentModel, linspace(-2.5, hZeroPlot), ...
    0 * linspace(1, modelYLim(2)), ...
    'Color', referenceLineColor, ...
    'LineWidth', 2)

for i = 1:numel(shuffledStats)
    markerSize = markerSizeFromN(shuffledStats(i).N, fullNRange, fullMarkerSizeRange);
    scatter(axIndependentModel, shuffledStats(i).hIndependent, shuffledStats(i).lambdaIndependent, markerSize, ...
        'Marker', dataMarker, ...
        'MarkerEdgeColor', markerEdgeColor, ...
        'MarkerFaceColor', shuffledStats(i).color, ...
        'MarkerFaceAlpha', markerFaceAlpha, ...
        'LineWidth', markerLineWidth)
end

scatter(axIndependentModel, hZeroPlot, 1, criticalPointSize, ...
    'Marker', dataMarker, ...
    'MarkerFaceColor', 'r', ...
    'MarkerEdgeColor', 'r')

axis(axIndependentModel, 'square')
xlim(axIndependentModel, [-2.5, 0])
ylim(axIndependentModel, modelYLim)
xlabel(axIndependentModel, [labelFont 'External field {\ith}'], ...
    'Interpreter', 'tex', ...
    'FontSize', axisLabelFontSize)
ylabel(axIndependentModel, [labelFont 'Interaction strength \lambda'], ...
    'Interpreter', 'tex', ...
    'FontSize', axisLabelFontSize)
formatFigureAxis(axIndependentModel, textFontName, tickLabelFontSize)

set(gcf, 'Renderer', 'painters')

%% Local functions

function shuffledStats = computeShuffledFullDataStats(repoRoot, datasetNames, datasetColors, useParallelProcessing, numParallelWorkers)
    shuffledStats = struct( ...
        'datasetName', {}, ...
        'sourceName', {}, ...
        'm', {}, ...
        'chi', {}, ...
        'chiIndependent', {}, ...
        'N', {}, ...
        'h', {}, ...
        'lambda', {}, ...
        'hIndependent', {}, ...
        'lambdaIndependent', {}, ...
        'color', {});

    useParallelProcessing = prepareParallelPool(useParallelProcessing, numParallelWorkers);

    shuffledStats = appendAllenStats( ...
        shuffledStats, repoRoot, datasetNames(1), datasetColors(1), ...
        useParallelProcessing);
    shuffledStats = appendHippoStats(shuffledStats, repoRoot, datasetNames(2), datasetColors(2));
    shuffledStats = appendStringerStats( ...
        shuffledStats, repoRoot, datasetNames(3), datasetColors(3), ...
        useParallelProcessing);
end

function shuffledStats = appendAllenStats(shuffledStats, repoRoot, datasetName, datasetColor, useParallelProcessing)
    dataFolder = fullfile(repoRoot, 'Allen', 'Allen_data_share');
    listing = dir(fullfile(dataFolder, '*.mat'));
    fileNames = string({listing.name});
    numFiles = numel(fileNames);
    baseColor = hex2rgb(datasetColor);
    entries = cell(numFiles, 1);

    if useParallelProcessing
        parfor fileIdx = 1:numFiles
            entries{fileIdx} = processAllenShuffledFile( ...
                dataFolder, fileNames(fileIdx), datasetName, baseColor, ...
                fileIdx, numFiles);
        end
    else
        for fileIdx = 1:numFiles
            entries{fileIdx} = processAllenShuffledFile( ...
                dataFolder, fileNames(fileIdx), datasetName, baseColor, ...
                fileIdx, numFiles);
        end
    end

    for fileIdx = 1:numel(entries)
        if ~isempty(entries{fileIdx})
            shuffledStats(end + 1) = entries{fileIdx}; %#ok<AGROW>
        end
    end
end

function entry = processAllenShuffledFile(dataFolder, fileName, datasetName, baseColor, fileIdx, numFiles)
    entry = [];
    fprintf('Allen shuffled full data %d of %d: %s\n', fileIdx, numFiles, fileName)

    try
        raw = load(fullfile(dataFolder, fileName));
    catch ME
        warning('Skipping %s because it could not be loaded: %s', ...
            fileName, ME.message)
        return
    end

    if ~isfield(raw, 'X')
        warning('Skipping %s because it does not contain X.', fileName)
        return
    end

    spiking_patterns = 2.*raw.X - 1;

    if isfield(raw, 'num_cells') && size(spiking_patterns, 1) ~= raw.num_cells
        warning('Skipping %s because spike matrix size does not match num_cells.', ...
            fileName)
        return
    end

    [m, chi, chiIndependent, N] = shuffledFullStatistics(spiking_patterns);
    [h, lambda] = fitShuffledModelPoint(m, chi, N);

    entry = makeStatsEntry(datasetName, string(fileName), ...
        m, chi, chiIndependent, N, h, lambda, baseColor);
end

function shuffledStats = appendHippoStats(shuffledStats, repoRoot, datasetName, datasetColor)
    dataFolder = fullfile(repoRoot, 'Hippo');
    baseColor = hex2rgb(datasetColor);

    fprintf('Hippocampus shuffled full data.\n')

    spiketime = load(fullfile(dataFolder, "Binary_1485neurons_30Hz_share.csv"));
    numnuerons = max(spiketime(:,1));
    totaltime = max(spiketime(:,2));

    spikes = zeros(totaltime, numnuerons);

    for i = 1:size(spiketime, 1)
        neuronIdx = spiketime(i, 1);
        timeIdx = spiketime(i, 2);
        spikes(timeIdx, neuronIdx) = 1;
    end

    spiking_patterns = 2.*spikes.' - 1;
    [m, chi, chiIndependent, N] = shuffledFullStatistics(spiking_patterns);
    [h, lambda] = fitShuffledModelPoint(m, chi, N);

    shuffledStats(end + 1) = makeStatsEntry( ... %#ok<AGROW>
        datasetName, "Binary_1485neurons_30Hz_share.csv", ...
        m, chi, chiIndependent, N, h, lambda, baseColor);
end

function shuffledStats = appendStringerStats(shuffledStats, repoRoot, datasetName, datasetColor, useParallelProcessing)
    dataFolder = fullfile(repoRoot, 'Stringer', 'data');
    infoFile = fullfile(repoRoot, 'Stringer', 'data_info.csv');
    dataInfo = readtable(infoFile, 'TextType', 'string');
    baseColor = hex2rgb(datasetColor);
    fileNames = dataInfo.file_name;
    stimulusTypes = dataInfo.stimulus;
    numRows = height(dataInfo);
    entries = cell(numRows, 1);

    if useParallelProcessing
        parfor rowIdx = 1:numRows
            entries{rowIdx} = processStringerShuffledRow( ...
                dataFolder, fileNames(rowIdx), stimulusTypes(rowIdx), ...
                datasetName, baseColor, rowIdx, numRows);
        end
    else
        for rowIdx = 1:numRows
            entries{rowIdx} = processStringerShuffledRow( ...
                dataFolder, fileNames(rowIdx), stimulusTypes(rowIdx), ...
                datasetName, baseColor, rowIdx, numRows);
        end
    end

    for rowIdx = 1:numel(entries)
        if ~isempty(entries{rowIdx})
            shuffledStats(end + 1) = entries{rowIdx}; %#ok<AGROW>
        end
    end
end

function entry = processStringerShuffledRow(dataFolder, filename, stimulusType, datasetName, baseColor, rowIdx, numRows)
    entry = [];
    fprintf('Stringer shuffled full data %d of %d: %s (%s)\n', ...
        rowIdx, numRows, filename, stimulusType)

    try
        ExpData = matfile(fullfile(dataFolder, filename + ".mat"));
        Stim = ExpData.stim;
    catch ME
        warning('Skipping %s because it could not be loaded: %s', filename, ME.message)
        return
    end

    if stimulusType == "resp"
        Response = Stim.resp;
    elseif stimulusType == "spont"
        Response = Stim.spont;
    else
        warning('Skipping %s because stimulus type "%s" is unknown.', filename, stimulusType)
        return
    end

    offneuron = find(mean(Response)==0);
    Response(:,offneuron) = [];

    [normResponse, mu, sigma ]= zscore(Response,1,1); %#ok<ASGLU>
    binary = (Response - mu)> 2*sigma;
    spiking_patterns = 2.*binary.' - 1;

    [m, chi, chiIndependent, N] = shuffledFullStatistics(spiking_patterns);
    [h, lambda] = fitShuffledModelPoint(m, chi, N);

    entry = makeStatsEntry( ...
        datasetName, filename + "_" + stimulusType, ...
        m, chi, chiIndependent, N, h, lambda, baseColor);
end

function useParallelProcessing = prepareParallelPool(useParallelProcessing, numParallelWorkers)
    if ~useParallelProcessing
        return
    end

    try
        pool = gcp('nocreate');

        if isempty(pool)
            if isempty(numParallelWorkers)
                parpool;
            else
                parpool(numParallelWorkers);
            end
        end
    catch ME
        warning('Parallel pool could not be started: %s. Running serially.', ME.message)
        useParallelProcessing = false;
    end
end

function [m, chi, chiIndependent, N] = shuffledFullStatistics(spiking_patterns)
    spiking_patterns = double(spiking_patterns);
    N = size(spiking_patterns, 1);
    num_bins = size(spiking_patterns, 2);

    neuronMeans = mean(spiking_patterns, 2);
    m = mean(neuronMeans);
    chiIndependent = mean(1 - neuronMeans.^2);

    shuffled_patterns = spiking_patterns;

    for neuronIdx = 1:N
        shuffled_patterns(neuronIdx, :) = shuffled_patterns(neuronIdx, randperm(num_bins));
    end

    populationActivity = sum(shuffled_patterns, 1);
    sumPairMeans = mean(populationActivity.^2);
    chi = (sumPairMeans - sum(neuronMeans).^2)/N;
end

function [h, lambda] = fitShuffledModelPoint(m, chi, N)
    h = nan;
    lambda = nan;

    if ~isfinite(m) || ~isfinite(chi) || ~isfinite(N) ...
            || abs(m) >= 1 || chi <= 0
        return
    end

    try
        [h, lambda] = hlambda(m, chi, N);
    catch ME
        warning('Skipping shuffled h/lambda fit for N = %.0f: %s', N, ME.message)
    end
end

function entry = makeStatsEntry(datasetName, sourceName, m, chi, chiIndependent, N, h, lambda, color)
    entry.datasetName = datasetName;
    entry.sourceName = sourceName;
    entry.m = m;
    entry.chi = chi;
    entry.chiIndependent = chiIndependent;
    entry.N = N;
    entry.h = h;
    entry.lambda = lambda;
    entry.hIndependent = nan;
    entry.lambdaIndependent = nan;
    entry.color = color;
end

function [shuffledStats, addedIndependentFits] = addIndependentModelFits(shuffledStats)
    addedIndependentFits = false;

    for idx = 1:numel(shuffledStats)
        needsFit = ~isfield(shuffledStats, 'hIndependent') ...
            || ~isfield(shuffledStats, 'lambdaIndependent') ...
            || isempty(shuffledStats(idx).hIndependent) ...
            || isempty(shuffledStats(idx).lambdaIndependent) ...
            || ~isfinite(shuffledStats(idx).hIndependent) ...
            || ~isfinite(shuffledStats(idx).lambdaIndependent);

        if ~needsFit
            continue
        end

        [hIndependent, lambdaIndependent] = fitShuffledModelPoint( ...
            shuffledStats(idx).m, shuffledStats(idx).chiIndependent, shuffledStats(idx).N);
        shuffledStats(idx).hIndependent = hIndependent;
        shuffledStats(idx).lambdaIndependent = lambdaIndependent;
        addedIndependentFits = true;
    end
end

function shuffledStats = applyDatasetColors(shuffledStats, datasetNames, datasetColors)
    for idx = 1:numel(shuffledStats)
        datasetIdx = find(datasetNames == string(shuffledStats(idx).datasetName), 1);

        if ~isempty(datasetIdx)
            shuffledStats(idx).color = hex2rgb(datasetColors(datasetIdx));
        end
    end
end

function markerSizes = markerSizeFromN(Ns, nRange, sizeRange)
    if nRange(1) == nRange(2)
        markerSizes = mean(sizeRange) * ones(size(Ns));
        return
    end

    scaledN = (Ns - nRange(1)) ./ diff(nRange);
    markerSizes = sizeRange(1) + scaledN .* diff(sizeRange);
end

function formatFigureAxis(ax, textFontName, tickLabelFontSize)
    box(ax, 'on')
    set(ax, 'TickDir', 'both')
    ax.FontName = textFontName;
    ax.FontSize = tickLabelFontSize;
    ax.TickLabelInterpreter = 'tex';
end

function rgb = hex2rgb(hex)
    hex = char(erase(string(hex), "#"));
    rgb = sscanf(hex, '%2x%2x%2x', [1 3]) / 255;
end
