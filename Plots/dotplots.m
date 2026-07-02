%% KEEP PLOTS

%% Box in (mu,chi) and image in (h,lambda)
% Figure mapping: this script generates the Fig. 3c grid and mapped box
% boundaries in statistics space and model space.
%% Precompute inverse map data
clear; clc; 

Ns = [20, 100, 1000];

mus  = linspace(-0.9,-0.02,24);
chis = linspace(1,3.5,24); 

edgePts = 50;

muMin = min(mus);
muMax = max(mus);
chiMin = min(chis);
chiMax = max(chis);

% Interior grid
hs = zeros(length(mus), length(chis), length(Ns));
ls = zeros(length(mus), length(chis), length(Ns));

% Box edges
muBottom  = linspace(muMin, muMax, edgePts);
chiBottom = chiMin * ones(size(muBottom));

chiRight = linspace(chiMin, chiMax, edgePts);
muRight  = muMax * ones(size(chiRight));

muTop  = linspace(muMax, muMin, edgePts);
chiTop = chiMax * ones(size(muTop));

chiLeft = linspace(chiMax, chiMin, edgePts);
muLeft  = muMin * ones(size(chiLeft));

muEdges  = {muBottom,  muRight,  muTop,  muLeft};
chiEdges = {chiBottom, chiRight, chiTop, chiLeft};

hEdges = cell(length(Ns),4);
lEdges = cell(length(Ns),4);

for a = 1:length(Ns)

    N = Ns(a);

    fprintf('Computing N = %d\n', N)

    % Map grid points
    for j = 1:length(mus)
        for k = 1:length(chis)

            fprintf('N index %d, mu index %d, chi index %d\n', a, j, k)

            mu = mus(j);
            chi = chis(k);

            [h, lambda] = hlambda(mu, chi, N);

            hs(j,k,a) = h;
            ls(j,k,a) = lambda;
        end
    end

    % Map box edges
    for e = 1:4
        hEdge = zeros(size(muEdges{e}));
        lEdge = zeros(size(muEdges{e}));

        for q = 1:length(muEdges{e})
            [e,q]
            [hTemp, lTemp] = hlambda(muEdges{e}(q), chiEdges{e}(q), N);
            hEdge(q) = hTemp;
            lEdge(q) = lTemp;
        end

        hEdges{a,e} = hEdge;
        lEdges{a,e} = lEdge;
    end
end

%%
edgeColors = {
    [0 0 1 0.5], ...
    [0 0.6 0 0.5], ...
    [1 0 0 0.5], ...
    [0 0.8 0.8 0.5]
};

dotSize = 20;
edgeWidth = 3;

muMin = min(mus);
muMax = max(mus);
chiMin = min(chis);
chiMax = max(chis);


for a = 1:length(Ns)

    N = Ns(a);

    % Original box in (mu,chi)
    if a == 1
    figure
    hold on
    set(gca, 'TickDir', 'both')
    box on

    for e = 1:4
        plot(muEdges{e}, chiEdges{e}, ...
            'Color', edgeColors{e}, ...
            'LineWidth', edgeWidth);
    end

    % Interior points only
    mus_in  = mus(2:end-1);
    chis_in = chis(2:end-1);

    [MU, CHI] = meshgrid(mus_in, chis_in);
    scatter(MU(:), CHI(:), dotSize, 'k', 'filled');

    xlabel('\mu', 'FontSize', 22)
    ylabel('\chi', 'FontSize', 22)
    title(['(\mu,\chi), N = ', num2str(N)], 'FontSize', 14)

    ax = gca;
    ax.FontSize = 22;

    axis square
    xlim([muMin-0.05, muMax+0.05])
    ylim([chiMin-0.1, chiMax+0.1])
    end

    % Transformed box in (h,lambda)
    figure
    hold on
    

    hs_in = hs(2:end-1, 2:end-1, a);
    ls_in = ls(2:end-1, 2:end-1, a);

    scatter(hs_in(:), ls_in(:), dotSize, 'k', 'filled', MarkerFaceAlpha=0.6);

    for e = 1:4
        plot(hEdges{a,e}, lEdges{a,e}, ...
            'Color', edgeColors{e}, ...
            'LineWidth', edgeWidth);
    end

    

    xlabel('h', 'FontSize', 30)
    ylabel('\lambda', 'FontSize', 30)
    title(['(h,\lambda), N = ', num2str(N)], 'FontSize', 14)

    axis square

    % xscale log
    % xlim([-1,-0.0005])
    % xticks([-1, -0.1, -0.01, -0.001])
    % xticklabels({'-10^0','-10^{-1}','-10^{-2}','-10^{-3}' })

    ylim([0,3])
    set(gca, 'TickDir', 'both')
    ax = gca;
    ax.FontSize = 26;
    box on
end
