%% KEEP PLOTS

% Plots data points in (mu, chi) and stored data points in (h, lambda).

clear
clc

plotsDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(plotsDir);
addpath(repoRoot)
addpath(fullfile(repoRoot, 'Allen'))
addpath(fullfile(repoRoot, 'Hippo'))
addpath(fullfile(repoRoot, 'Stringer'))

files = ["Allenhldata.mat", "hippomuchidata.mat", "stringerhldata.mat"];
markers = {'square', 'o', '^'};
colors = ["#1f4fd6", "#c65100", "#167a43"];
markerSize = 70;

% -1 denotes the final/full data point.
%figureSets = {[-1]};
figureSets = {[20 100 1000 -1]};

%% Plot mu and chi from data

for f = 1:length(figureSets)
    figure
    hold on

    currentSet = figureSets{f};

    for i = 1:length(files)
        [Data, ~, ~] = loadDataFile(files(i));
        [Nss, mus, chis] = summarizeData(Data);

        baseColor = hex2rgb(colors(i));
        [c20, c100, c1000, cFull] = dataAmountColors(baseColor);
        dim = datasetSearchDim(i);

        for val = currentSet
            if val == -1
                scatter(lastofarray(mus), lastofarray(chis), markerSize, ...
                    'Marker', markers{i}, ...
                    'MarkerEdgeColor', cFull, ...
                    'LineWidth', 1.5)
            elseif val == 0
                continue
            else
                [~, I] = min(abs(Nss - val), [], dim, "linear");
                scatter(mus(I), chis(I), markerSize, ...
                    'Marker', markers{i}, ...
                    'MarkerEdgeColor', colorForN(val, c20, c100, c1000, cFull), ...
                    'LineWidth', 1.5)
            end
        end
    end

    x = linspace(-1, 0, 1000);
    plot(x, 1 - x.^2, 'LineWidth', 2, 'Color', 'blue')

    ylim([0.05, 70])
    xlim([-1, -0.8])
    axis square
    xlabel('\mu', 'FontSize', 18)
    ylabel('\chi', 'FontSize', 18)
    set(gca, 'YScale', 'log')
    set(gca, 'TickDir', 'both')
    box on
end

%% Plot h and lambda from data

for f = 1:length(figureSets)
    figure
    hold on

    plot(-0.00005 + 0 * linspace(1, 2), linspace(1, 2.5), 'r', 'LineWidth', 2)
    plot(linspace(-1.4, -0.00001), 0 * linspace(1, 2.5), 'b', 'LineWidth', 2)

    currentSet = figureSets{f};

    for i = 1:length(files)
        [Data, hs, ls] = loadDataFile(files(i));
        [Nss, ~, ~] = summarizeData(Data);

        baseColor = hex2rgb(colors(i));
        [c20, c100, c1000, cFull] = dataAmountColors(baseColor);
        dim = datasetSearchDim(i);

        for val = currentSet
            if val == -1
                scatter(lastofarray(hs), lastofarray(ls), markerSize, ...
                    'Marker', markers{i}, ...
                    'MarkerEdgeColor', cFull, ...
                    'LineWidth', 1.5)
            elseif val == 0
                continue
            else
                [~, I] = min(abs(Nss - val), [], dim, "linear");
                scatter(hs(I), ls(I), markerSize, ...
                    'Marker', markers{i}, ...
                    'MarkerEdgeColor', colorForN(val, c20, c100, c1000, cFull), ...
                    'LineWidth', 1.5)
            end
        end
    end

    axis square
    xlabel('h', 'FontSize', 20)
    ylabel('\lambda', 'FontSize', 20)
    xscale log
    box on
    set(gca, 'TickDir', 'both')
    xlim([-1, -0.00005])
    ylim([0, 2.5])
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
