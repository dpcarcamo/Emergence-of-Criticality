%% KEEP PLOTS

% Used to make the cartoon plot of the jacobian. Plots circles in h-lambda
% space and sees where they land and distort in mu chi. Also makes plots of
% the susceptibility and Jacobian for changing h at fixed lambda. 
% Figure mapping: Fig. 2a-b are the circle mapping panels, Fig. 2c is the
% susceptibility versus h, and Fig. 2d is the Jacobian versus h.

clear
clc

currentPath = pwd;
folders = split(currentPath, '\');
newPath = join(folders(1:length(folders)-1),"\");
addpath(newPath{1});
addpath(strcat(newPath{1}, '\Allen'))
addpath(strcat(newPath{1}, '\Hippo'))
addpath(strcat(newPath{1}, '\Stringer'))



%% Plot solid circles in (h,lambda) space and their images in (mu,chi) space


% --- Plot settings ---
n = 20;          % system size for muChiExact2Spin
r = 0.03;        % radius of each circle in (h,lambda) space
numPts = 400;    % number of points along each circle boundary
axisLabelFontSize = 34;
tickLabelFontSize = 20;
textFontName = 'Helvetica';
labelFont = ['\fontname{' textFontName '}'];
jColorValues = [
    0.5020    0.1765    0.1765
    0.7020         0         0
    0.8941    0.2039         0
    0.9725    0.5373         0
    1.0000    0.6745    0.0667
];
hGreyColor = [0.65 0.65 0.65];
hGreyStrength = 0.85;

centers = [];


for j = 0.2:0.25:1.95
    for i = -0.1:-0.1:-1.95
        centers = [centers; [i,j]];
    end
end

theta = linspace(0, 2*pi, numPts);

% Extract center coordinates
hs = centers(:,1);
ls = centers(:,2);

% Normalize h and lambda independently to [0,1]
hNorm = (hs - min(hs)) / (max(hs) - min(hs));
lNorm = (ls - min(ls)) / (max(ls) - min(ls));
jColorPositions = linspace(0, 1, size(jColorValues, 1));

cols = zeros(size(centers,1),3);

for k = 1:size(centers,1)

    jColor = interp1(jColorPositions, jColorValues, lNorm(k), 'linear');
    greyWeight = hGreyStrength * (1 - hNorm(k));
    cols(k,:) = (1 - greyWeight) * jColor + greyWeight * hGreyColor;
end

figure('Name', 'Jacobian cartoon panels', 'Color', 'w');
tiledlayout(2, 2, ...
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

    % Circle boundary in (h, lambda)
    hCircle = h0 + r*cos(theta);
    lambdaCircle = lambda0 + r*sin(theta);

    fill(hCircle, lambdaCircle, cols(k,:), ...
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
xlim([-2, 0])
ylim([0,2])

ax = gca;

ax.FontName = textFontName;
ax.FontSize = tickLabelFontSize;      % tick labels
ax.LineWidth = 1.5;

% --- Plot mapped filled curves in (mu, chi) space ---
nexttile(2)
hold on
%set(gca, 'TickDir', 'both')
box on

for k = 1:size(centers,1)
    h0 = centers(k,1);
    lambda0 = centers(k,2);

    % Circle boundary in (h, lambda)
    hCircle = h0 + r*cos(theta);
    lambdaCircle = lambda0 + r*sin(theta);

    % Map boundary into (mu, chi)
    muVals = zeros(size(theta));
    chiVals = zeros(size(theta));

    for i = 1:length(theta)
        [muVals(i), chiVals(i)] = muChiExact2Spin(hCircle(i), lambdaCircle(i), n);
    end

    fill(muVals, chiVals, cols(k,:), ...
        'EdgeColor', 'none');
end

x = linspace(-1,-0.0001,500);
plot(x, 1-x.^2, 'Color', [0.5765    0.5843    0.5961], 'LineWidth', 1.5)

xlabel([labelFont 'Average activity \langle\mu\rangle'], ...
    'Interpreter', 'tex', ...
    'FontSize', axisLabelFontSize)
ylabel([labelFont 'Susceptibility \chi'], ...
    'Interpreter', 'tex', ...
    'FontSize', axisLabelFontSize)

axis square
xlim([-1,-0])
ylim([0.05, 5])
ax = gca;

ax.FontName = textFontName;
ax.FontSize = tickLabelFontSize;      % tick labels
ax.LineWidth = 1.5;    % thicker axis lines
ax.TickLabelInterpreter = 'tex';

set(gcf, 'Renderer', 'painters');
%yscale log
%%
showLabels = false;
hs = linspace(-2,-0.001, 1000); 

nexttile(3)
hold on
for l = 0.2:0.25:1.95

    chis = [];
    for h = hs

    [mu, chi] = muChiExact2Spin(h,l,n);
    chis = [chis,chi];

    end

    y = (l - min(ls)) / (max(ls) - min(ls));
    color = interp1(jColorPositions, jColorValues, y, 'linear');

    plot(hs, chis, color=color, LineWidth=2)

    if showLabels
        text(hs(end), chis(end), sprintf('  J = %.2f', l), ...
            'FontSize', 14, ...
            'Color', color, ...
            'VerticalAlignment', 'middle')
    end

end


xlabel([labelFont 'External field {\ith}'], ...
    'Interpreter', 'tex', ...
    'FontSize', axisLabelFontSize)
ylabel([labelFont 'Susceptibility \chi'], ...
    'Interpreter', 'tex', ...
    'FontSize', axisLabelFontSize)

axis square
ax = gca;
%yscale log
xscale log

ax.FontName = textFontName;
ax.FontSize = tickLabelFontSize;      % tick labels
ax.LineWidth = 1.5;    % thicker axis lines

set(gcf, 'Renderer', 'painters');
box on
xlim([-1,-0.001])
xticks([-1, -0.1, -0.01, -0.001])
xticklabels({'-10^{0}', '-10^{-1}', '-10^{-2}', '-10^{-3}'})
xtickangle(0)
ax.TickLabelInterpreter = 'tex';
%

hs = linspace(-2,-0.001, 1000); 

nexttile(4)
hold on
for l = 0.2:0.25:1.95

    Jacs = [];
    for h = hs

    Jac = Jacobian(h,l,n);
    Jacs = [Jacs,Jac];

    end

    y = (l - min(ls)) / (max(ls) - min(ls));
    color = interp1(jColorPositions, jColorValues, y, 'linear');

    plot(hs, Jacs, color=color, LineWidth=2)

    if showLabels
        text(hs(end), Jacs(end), sprintf('  J = %.2f', l), ...
            'FontSize', 14, ...
            'Color', color, ...
            'VerticalAlignment', 'middle')
    end

end


xlabel([labelFont 'External field {\ith}'], ...
    'Interpreter', 'tex', ...
    'FontSize', axisLabelFontSize)
ylabel([labelFont 'Jacobian determinant'], ...
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
