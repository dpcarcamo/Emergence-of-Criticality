%% KEEP PLOTS

% Plot double-well parameter trajectories for fixed full-data statistics.
% Each recording uses its final/full m and chi values. The double-well
% approximation is evaluated while N is swept from NMin to NMax.

clear
clc

plotsDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(plotsDir);
addpath(repoRoot)
addpath(fullfile(repoRoot, 'Allen'))
addpath(fullfile(repoRoot, 'Hippo'))
addpath(fullfile(repoRoot, 'Stringer'))

%% Plot settings

files = [
    fullfile(repoRoot, "Allen", "Allenhldata.mat")
    fullfile(repoRoot, "Stringer", "stringerhldata.mat")
    fullfile(repoRoot, "Hippo", "hippomuchidata.mat")
];
datasetNames = ["Allen", "Stringer", "Hippocampus"];
datasetColors = ["#2676ad", "#2f682c", "#4f4cc4"];

NMin = 1e2;
NMax = 1e6;
nNValues = 50;
NValues = unique(round(logspace(log10(NMin), log10(NMax), nNValues)));

trajectoryLineWidth = 1.6;
lowNLightenAmount = 0.84;
highNDarkenAmount = 0.08;

referenceLineColor = "#939598";
hZeroPlot = -10^-6.5;
hXLim = [-0.1, hZeroPlot];
hTickExponents = -1:-1:-6;
hTickValues = [-10.^hTickExponents, hZeroPlot];
hTickLabels = arrayfun(@(x) sprintf('-10^{%d}', x), hTickExponents, ...
    'UniformOutput', false);
hTickLabels{end + 1} = '    0';
linearHXLim = [-0.1, 0];
linearHTicks = -0.1:0.02:0;

lambdaYLim = [0, 3];
lambdaTicks = 0:1:3;
criticalLineWidth = 2;
criticalPointSize = 85;

axisLabelFontSize = 32;
tickLabelFontSize = 18;
textFontName = 'Helvetica';
labelFont = ['\fontname{' textFontName '}'];

%% Plot fixed-statistic double-well trajectories

figure('Name', 'Full-statistic double-well trajectories', ...
    'Color', 'w')
hold on

plot([hZeroPlot, hZeroPlot], [1, lambdaYLim(2)], ...
    'r', ...
    'LineWidth', criticalLineWidth)
plot(linspace(hXLim(1), hZeroPlot, 500), ...
    zeros(1, 500), ...
    'Color', referenceLineColor, ...
    'LineWidth', 2)

for datasetIdx = 1:length(files)
    [Data] = loadDataFile(files(datasetIdx));
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
    lowNColor = lightenColor(baseColor, lowNLightenAmount);
    highNColor = darkenColor(baseColor, highNDarkenAmount);

    fprintf('Double-well trajectories: plotting %d full %s recordings.\n', ...
        nnz(validTargets), char(datasetNames(datasetIdx)))

    for popIdx = find(validTargets)
        mTarget = fullMus(popIdx);
        chiTarget = fullChis(popIdx);

        u0 = sqrt(chiTarget ./ NValues + mTarget.^2);
        hTheory = nan(size(NValues));
        lambdaTheory = nan(size(NValues));

        goodTheory = isfinite(u0) ...
            & u0 > 0 ...
            & u0 < 1 ...
            & abs(mTarget ./ u0) < 1;

        lambdaTheory(goodTheory) = atanh(u0(goodTheory)) ./ u0(goodTheory);
        hTheory(goodTheory) = atanh(mTarget ./ u0(goodTheory)) ...
            ./ (NValues(goodTheory) .* u0(goodTheory));

        goodPlot = isfinite(hTheory) ...
            & isfinite(lambdaTheory) ...
            & hTheory < 0 ...
            & lambdaTheory > 0;

        for nIdx = 1:length(NValues) - 1
            if goodPlot(nIdx) && goodPlot(nIdx + 1)
                colorWeight = (log10(mean(NValues(nIdx:nIdx + 1))) - log10(NMin)) ...
                    / (log10(NMax) - log10(NMin));
                colorWeight = max(0, min(1, colorWeight));
                currentColor = (1 - colorWeight) * lowNColor ...
                    + colorWeight * highNColor;

                plot(hTheory(nIdx:nIdx + 1), ...
                    lambdaTheory(nIdx:nIdx + 1), ...
                    'Color', currentColor, ...
                    'LineWidth', trajectoryLineWidth)
            end
        end
    end
end

scatter(hZeroPlot, 1, criticalPointSize, ...
    'Marker', 'o', ...
    'MarkerFaceColor', 'r', ...
    'MarkerEdgeColor', 'r')

axis square
xscale log
xlim(hXLim)
ylim(lambdaYLim)
xticks(hTickValues)
xticklabels(hTickLabels)
xtickangle(0)
yticks(lambdaTicks)
yticklabels(arrayfun(@num2str, lambdaTicks, 'UniformOutput', false))

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

%% Plot fixed-statistic double-well trajectories with a linear h axis

figure('Name', 'Full-statistic double-well trajectories linear h axis', ...
    'Color', 'w')
hold on

plot([0, 0], [1, lambdaYLim(2)], ...
    'r', ...
    'LineWidth', criticalLineWidth)
plot(linspace(linearHXLim(1), linearHXLim(2), 500), ...
    zeros(1, 500), ...
    'Color', referenceLineColor, ...
    'LineWidth', 2)

for datasetIdx = 1:length(files)
    [Data] = loadDataFile(files(datasetIdx));
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
    lowNColor = lightenColor(baseColor, lowNLightenAmount);
    highNColor = darkenColor(baseColor, highNDarkenAmount);

    fprintf('Double-well linear trajectories: plotting %d full %s recordings.\n', ...
        nnz(validTargets), char(datasetNames(datasetIdx)))

    for popIdx = find(validTargets)
        mTarget = fullMus(popIdx);
        chiTarget = fullChis(popIdx);

        u0 = sqrt(chiTarget ./ NValues + mTarget.^2);
        hTheory = nan(size(NValues));
        lambdaTheory = nan(size(NValues));

        goodTheory = isfinite(u0) ...
            & u0 > 0 ...
            & u0 < 1 ...
            & abs(mTarget ./ u0) < 1;

        lambdaTheory(goodTheory) = atanh(u0(goodTheory)) ./ u0(goodTheory);
        hTheory(goodTheory) = atanh(mTarget ./ u0(goodTheory)) ...
            ./ (NValues(goodTheory) .* u0(goodTheory));

        goodPlot = isfinite(hTheory) ...
            & isfinite(lambdaTheory) ...
            & hTheory < 0 ...
            & lambdaTheory > 0;

        for nIdx = 1:length(NValues) - 1
            if goodPlot(nIdx) && goodPlot(nIdx + 1)
                colorWeight = (log10(mean(NValues(nIdx:nIdx + 1))) - log10(NMin)) ...
                    / (log10(NMax) - log10(NMin));
                colorWeight = max(0, min(1, colorWeight));
                currentColor = (1 - colorWeight) * lowNColor ...
                    + colorWeight * highNColor;

                plot(hTheory(nIdx:nIdx + 1), ...
                    lambdaTheory(nIdx:nIdx + 1), ...
                    'Color', currentColor, ...
                    'LineWidth', trajectoryLineWidth)
            end
        end
    end
end

scatter(0, 1, criticalPointSize, ...
    'Marker', 'o', ...
    'MarkerFaceColor', 'r', ...
    'MarkerEdgeColor', 'r')

axis square
xlim(linearHXLim)
ylim(lambdaYLim)
xticks(linearHTicks)
yticks(lambdaTicks)
yticklabels(arrayfun(@num2str, lambdaTicks, 'UniformOutput', false))

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

%% Helper functions

function [Data] = loadDataFile(fileName)
    obj = load(char(fileName));
    Data = obj.Data;
    Data(Data == 0) = nan;
end

function [Nss, mus, chis] = summarizeData(Data)
    Nss = squeeze(mean(squeeze(Data(3, :, :, :)), 1, 'omitnan'));
    mus = squeeze(mean(squeeze(Data(1, :, :, :)), 1, 'omitnan'));
    chis = squeeze(mean(squeeze(Data(2, :, :, :)), 1, 'omitnan'));
end

function rgb = hex2rgb(hex)
    hex = char(erase(string(hex), "#"));
    rgb = sscanf(hex, '%2x%2x%2x', [1 3]) / 255;
end

function cOut = lightenColor(cIn, amt)
    cOut = cIn + (1 - cIn) * amt;
end

function cOut = darkenColor(cIn, amt)
    cOut = cIn * (1 - amt);
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
