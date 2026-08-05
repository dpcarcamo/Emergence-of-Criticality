%% KEEP PLOTS

% Used to make the cartoon plot of the jacobian. Plots tiles in h-lambda
% space and sees where they land and distort in mu chi. Also makes plots of
% the Jacobian for changing h and N.
% Figure mapping: this script generates the h-lambda cartoon, Jacobian
% diagnostics, and mapped mu-chi panels at multiple N values.

clear
clc

currentPath = pwd;
folders = split(currentPath, '\');
newPath = join(folders(1:length(folders)-1),"\");
addpath(newPath{1});
addpath(strcat(newPath{1}, '\Allen'))
addpath(strcat(newPath{1}, '\Hippo'))
addpath(strcat(newPath{1}, '\Stringer'))



%% Plot solid tiles in (h,lambda) space and their images in (mu,chi) space


% --- Plot settings ---
n = 20;          % system size for the Jacobian versus h panel
tileSideLength = 0.05;    % side length for touching squares in (h,lambda) space
numPtsPerSide = 20;      % number of points along each square edge
mappedNValues = [20 100 1000];
jacobianNValues = unique(round(logspace(log10(10), log10(1000), 120)));
lambdaValues = -0.001:tileSideLength:3;
hTileValues = -1 + tileSideLength/2:tileSideLength:-tileSideLength/2;
jacobianLambdaValues = [0.2:0.5:3];
axisLabelFontSize = 34;
tickLabelFontSize = 20;
textFontName = 'Helvetica';
labelFont = ['\fontname{' textFontName '}'];
jColorValues = [
    0.5020/4    0.1765/4    0.1765/4
    0.5020    0.1765    0.1765
    0.7020         0         0
    0.8941    0.2039         0
    0.9725    0.5373         0
    1.0000    0.6745    0.0667
    1.0000    0.7745    0.1667
];
hGreyColor = [0.65 0.65 0.65];
hGreyStrength = 0.95;

centers = [];


for j = lambdaValues
    for i = hTileValues
        centers = [centers; [i,j]];
    end
end

tileHalfLength = tileSideLength/2;
edgeVals = linspace(-tileHalfLength, tileHalfLength, numPtsPerSide);

% Extract center coordinates
centerHs = centers(:,1);
centerLambdas = centers(:,2);
jacobianReferenceH = -0.001;

% Normalize h and lambda independently to [0,1]
hNorm = (centerHs - min(centerHs)) / (max(centerHs) - min(centerHs));
lNorm = (centerLambdas - min(centerLambdas)) / (max(centerLambdas) - min(centerLambdas));
jColorPositions = linspace(0, 1, size(jColorValues, 1));

cols = zeros(size(centers,1),3);

for k = 1:size(centers,1)

    jColor = interp1(jColorPositions, jColorValues, lNorm(k), 'linear');
    greyWeight = hGreyStrength * (1 - hNorm(k));
    cols(k,:) = (1 - greyWeight) * jColor + greyWeight * hGreyColor;
end

figure('Name', 'Jacobian cartoon panels', 'Color', 'w');
tiledlayout(2, 3, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

% --- Plot in (h, lambda) space ---
nexttile(1)
hold on
%set(gca, 'TickDir', 'both')
box on

for k = 1:size(centers,1)
    h0 = centers(k,1);
    lambda0 = centers(k,2);

    % Square boundary in (h, lambda)
    hSquare = [h0 + edgeVals, ...
        h0 + tileHalfLength * ones(1, numPtsPerSide), ...
        h0 + fliplr(edgeVals), ...
        h0 - tileHalfLength * ones(1, numPtsPerSide)];
    lambdaSquare = [lambda0 - tileHalfLength * ones(1, numPtsPerSide), ...
        lambda0 + edgeVals, ...
        lambda0 + tileHalfLength * ones(1, numPtsPerSide), ...
        lambda0 + fliplr(edgeVals)];

    fill(hSquare, lambdaSquare, cols(k,:), ...
        'EdgeColor', 'none');
end

xlabel([labelFont 'External field {\ith}'], ...
    'Interpreter', 'tex', ...
    'FontSize', axisLabelFontSize)
ylabel([labelFont 'Interaction strength \lambda'], ...
    'Interpreter', 'tex', ...
    'FontSize', axisLabelFontSize)

axis equal
axis square
xlim([-1, 0])
ylim([0,3])

ax = gca;

ax.FontName = textFontName;
ax.FontSize = tickLabelFontSize;      % tick labels
ax.LineWidth = 1.5;

% --- Plot Jacobian versus h ---
nexttile(2)
hold on
hs = linspace(-2,-0.001, 1000); 

for l = jacobianLambdaValues

    Jacs = [];
    for h = hs

    Jac = Jacobian(h,l,n);
    Jacs = [Jacs,Jac];

    end

    y = (l - min(centerLambdas)) / (max(centerLambdas) - min(centerLambdas));
    color = interp1(jColorPositions, jColorValues, y, 'linear');

    plot(hs, Jacs, color=color, LineWidth=2)

end


xlabel([labelFont 'External field {\ith}'], ...
    'Interpreter', 'tex', ...
    'FontSize', axisLabelFontSize)
ylabel([labelFont 'Jacobian |{\bfJ}|'], ...
    'Interpreter', 'tex', ...
    'FontSize', axisLabelFontSize)
%yscale log

axis square
ax = gca;
box on

ax.FontName = textFontName;
ax.FontSize = tickLabelFontSize;      % tick labels
ax.LineWidth = 1.5;    % thicker axis lines
ax.TickLabelInterpreter = 'tex';

set(gcf, 'Renderer', 'painters');
xscale log
xlim([-1,-0.001])
xticks([-1, -0.1, -0.01, -0.001])
xticklabels({'-10^{0}', '-10^{-1}', '-10^{-2}', '-10^{-3}'})
xtickangle(0)
%

% --- Plot Jacobian versus N at the h value closest to zero ---
nexttile(3)
hold on

for lambda0 = jacobianLambdaValues
    Jacs = zeros(size(jacobianNValues));

    for idx = 1:length(jacobianNValues)
        Jacs(idx) = Jacobian(jacobianReferenceH, lambda0, jacobianNValues(idx));
    end

    y = (lambda0 - min(centerLambdas)) / (max(centerLambdas) - min(centerLambdas));
    color = interp1(jColorPositions, jColorValues, y, 'linear');

    plot(jacobianNValues, Jacs, color=color, LineWidth=2)
end

xlabel([labelFont 'System size {\itN}'], ...
    'Interpreter', 'tex', ...
    'FontSize', axisLabelFontSize)
ylabel([labelFont 'Jacobian |{\bfJ}|'], ...
    'Interpreter', 'tex', ...
    'FontSize', axisLabelFontSize)
axis square
ax = gca;
box on
ax.FontName = textFontName;
ax.FontSize = tickLabelFontSize;      % tick labels
ax.LineWidth = 1.5;    % thicker axis lines
ax.TickLabelInterpreter = 'tex';
xscale log
yscale log
xlim([10, 1000])
xticks([10 100 1000])
xticklabels({'10^{1}', '10^{2}', '10^{3}'})
xtickangle(0)

% --- Plot mapped filled curves in (mu, chi) space ---
for nIdx = 1:length(mappedNValues)
    currentN = mappedNValues(nIdx);

    nexttile(3 + nIdx)
    hold on
    %set(gca, 'TickDir', 'both')
    box on
    mappedMuMin = -1;
    mappedMuMax = 0;
    mappedChiMin = 0.05;
    mappedChiMax = 5;

    for k = 1:size(centers,1)
        h0 = centers(k,1);
        lambda0 = centers(k,2);

        % Square boundary in (h, lambda)
        hSquare = [h0 + edgeVals, ...
            h0 + tileHalfLength * ones(1, numPtsPerSide), ...
            h0 + fliplr(edgeVals), ...
            h0 - tileHalfLength * ones(1, numPtsPerSide)];
        lambdaSquare = [lambda0 - tileHalfLength * ones(1, numPtsPerSide), ...
            lambda0 + edgeVals, ...
            lambda0 + tileHalfLength * ones(1, numPtsPerSide), ...
            lambda0 + fliplr(edgeVals)];

        % Map boundary into (mu, chi)
        muVals = zeros(size(hSquare));
        chiVals = zeros(size(hSquare));

        for i = 1:length(hSquare)
            [muVals(i), chiVals(i)] = muChiExact2Spin(hSquare(i), lambdaSquare(i), currentN);
        end

        finiteMu = muVals(isfinite(muVals));
        finiteChi = chiVals(isfinite(chiVals) & chiVals > 0);

        if ~isempty(finiteMu)
            mappedMuMin = min(mappedMuMin, min(finiteMu));
            mappedMuMax = max(mappedMuMax, max(finiteMu));
        end

        if ~isempty(finiteChi)
            mappedChiMin = min(mappedChiMin, min(finiteChi));
            mappedChiMax = max(mappedChiMax, max(finiteChi));
        end

        fill(muVals, chiVals, cols(k,:), ...
            'EdgeColor', 'none');
    end

    % Independent line
    x = linspace(-1,-0.0001,500);
    plot(x, 1-x.^2, 'Color', [0.5765    0.5843    0.5961], 'LineWidth', 3)

    xlabel([labelFont 'Average activity {\itm}'], ...
        'Interpreter', 'tex', ...
        'FontSize', axisLabelFontSize)

    if nIdx == 1
        ylabel([labelFont 'Susceptibility \chi'], ...
            'Interpreter', 'tex', ...
            'FontSize', axisLabelFontSize)
    else
        ylabel('')
    end

    axis square
    xlim([mappedMuMin, mappedMuMax])
    % ylim([mappedChiMin, 1.05 * mappedChiMax])
    ylim([0,1000])
    ax = gca;

    yscale log

    ax.FontName = textFontName;
    ax.FontSize = tickLabelFontSize;      % tick labels
    ax.LineWidth = 1.5;    % thicker axis lines
    ax.TickLabelInterpreter = 'tex';
end

set(gcf, 'Renderer', 'painters');
