%% KEEP PLOTS

% Used to make the cartoon plot of the jacobian. Plots circles in h-lambda
% space and sees where they land and distort in mu chi. Also makes plots of
% the suscpetibility and Jacobian for changing h at fixed lambda. 

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

grayLight = [0.35 0.35 0.35];
grayDark = [0.75 0.75 0.75];

colorLight = [0.3098    0.4667    0.8392]*1.15;  % nice solid blue
colorDark = [0.7176    0.4902    0.7804];  % lighter blue

cols = zeros(size(centers,1),3);

for k = 1:size(centers,1)

    x = hNorm(k);   % left to right
    y = lNorm(k);   % bottom to top

    % bottom edge: gray ramp
    bottomColor = (1-x)*grayDark + x*grayLight;

    % top edge: colored ramp
    topColor = (1-x)*colorDark + x*colorLight;

    % interpolate vertically
    cols(k,:) = (1-y)*bottomColor + y*topColor;
end

figure;

% --- Plot in (h, lambda) space ---
subplot(1,2,1)
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

xlabel('h', 'FontSize', 22)
ylabel('\lambda', 'FontSize', 22)

axis equal
axis square
xlim([-2, 0])
ylim([0,2])

ax = gca;

ax.FontSize = 16;      % tick labels
ax.LineWidth = 1.5;

% --- Plot mapped filled curves in (mu, chi) space ---
subplot(1,2,2)
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
plot(x, 1-x.^2, 'b', 'LineWidth', 1.5)

xlabel('\mu', 'FontSize', 28)
ylabel('\chi', 'FontSize', 28)

axis square
xlim([-1,-0])
ylim([0.05, 5])
ax = gca;

ax.FontSize = 20;      % tick labels
ax.LineWidth = 1.5;    % thicker axis lines

set(gcf, 'Renderer', 'painters');
%yscale log
%%
showLabels = true;
hs = linspace(-2,-0.001, 1000); 

figure 
hold on
for l = 0.2:0.25:1.95

    chis = [];
    for h = hs

    [mu, chi] = muChiExact2Spin(h,l,n);
    chis = [chis,chi];

    end

    x = 1;   % left to right
    y = (l - min(ls)) / (max(ls) - min(ls));   % bottom to top

    % bottom edge: gray ramp
    bottomColor = (1-x)*grayDark + x*grayLight;

    % top edge: colored ramp
    topColor = (1-x)*colorDark + x*colorLight;

    % interpolate vertically
    color = (1-y)*bottomColor + y*topColor;

    plot(hs, chis, color=color, LineWidth=2)

    if showLabels
        text(hs(end), chis(end), sprintf('  \\lambda = %.2f', l), ...
            'FontSize', 14, ...
            'Color', color, ...
            'VerticalAlignment', 'middle')
    end

end


xlabel('h', 'FontSize', 28)
ylabel('\chi', 'FontSize', 28)

axis square
ax = gca;
%yscale log
xscale log

ax.FontSize = 20;      % tick labels
ax.LineWidth = 1.5;    % thicker axis lines

set(gcf, 'Renderer', 'painters');
box on
xlim([-1,-0.001])
%

hs = linspace(-2,-0.001, 1000); 

figure 
hold on
for l = 0.2:0.25:1.95

    Jacs = [];
    for h = hs

    Jac = Jacobian(h,l,n);
    Jacs = [Jacs,Jac];

    end

    x = 1;   % left to right
    y = (l - min(ls)) / (max(ls) - min(ls));   % bottom to top

    % bottom edge: gray ramp
    bottomColor = (1-x)*grayDark + x*grayLight;

    % top edge: colored ramp
    topColor = (1-x)*colorDark + x*colorLight;

    % interpolate vertically
    color = (1-y)*bottomColor + y*topColor;

    plot(hs, Jacs, color=color, LineWidth=2)

    if showLabels
        text(hs(end), Jacs(end), sprintf('  \\lambda = %.2f', l), ...
            'FontSize', 14, ...
            'Color', color, ...
            'VerticalAlignment', 'middle')
    end

end


xlabel('h', 'FontSize', 22)
ylabel('Jacobian', 'FontSize', 22)
%yscale log

axis square
ax = gca;
box on

ax.FontSize = 20;      % tick labels
ax.LineWidth = 1.5;    % thicker axis lines

set(gcf, 'Renderer', 'painters');
xlim([-1,-0.001])

xscale log
