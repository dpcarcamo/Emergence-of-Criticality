%% KEEP PLOTS
% File to plot the double well theory


%% Plot one population from each dataset: refit h/lambda vs theory

files = ["Allenhldata.mat", "hippomuchidata.mat", "stringerhldata.mat"];

markers = {'square', 'o', '^'};
colors = ["#2676ad","#1b9671","#2f682c"];

markerSize = 20;

figure(1); clf; hold on
figure(2); clf; hold on
figure(3); clf; hold on
figure(4); clf; hold on

for i = 1:length(files)

    obj = load(files(i));

    Data = obj.Data;
    Data(Data == 0) = nan;

    Nss  = squeeze(mean(squeeze(Data(3,:,:,:)), 1, 'omitnan'));
    mus  = squeeze(mean(squeeze(Data(1,:,:,:)), 1, 'omitnan'));
    chis = squeeze(mean(squeeze(Data(2,:,:,:)), 1, 'omitnan'));

    % Pick one population
    if i == 2
        % Hippocampus only has one population
        Nuse   = Nss(:);
        muse   = mus(:);
        chiuse = chis(:);
    else
        % Allen/Stringer: columns are populations, rows are changing N
        counts = sum(isfinite(Nss), 1);
        [~, popIdx] = max(counts);

        Nuse   = Nss(:, popIdx);
        muse   = mus(:, popIdx);
        chiuse = chis(:, popIdx);
    end

    % Clean selected data
    valid = isfinite(Nuse) & isfinite(muse) & isfinite(chiuse) ...
          & abs(muse) < 1 & chiuse > 0;

    Nuse   = Nuse(valid);
    muse   = muse(valid);
    chiuse = chiuse(valid);

    % Refit h and lambda using hlambda
    hFit = nan(size(Nuse));
    lFit = nan(size(Nuse));

    for j = 1:length(Nuse)
        N   = round(Nuse(j));
        mu  = muse(j);
        chi = chiuse(j);

        try
            [hFit(j), lFit(j)] = hlambda(mu, chi, N);
        catch
            hFit(j) = nan;
            lFit(j) = nan;
        end
    end

    % Theory from exactly the same selected data
    u0 = sqrt(chiuse ./ Nuse + muse.^2);

    hTheory = nan(size(muse));
    lTheory = nan(size(muse));

    goodTheory = isfinite(u0) ...
              & u0 > 0 ...
              & abs(u0) < 1 ...
              & abs(muse ./ u0) < 1;

    lTheory(goodTheory) = atanh(u0(goodTheory)) ./ u0(goodTheory);

    hTheory(goodTheory) = atanh(muse(goodTheory) ./ u0(goodTheory)) ...
                        ./ (Nuse(goodTheory) .* u0(goodTheory));

    % Keep only points where both fit and theory are valid
    good = isfinite(hFit) & isfinite(lFit) ...
         & isfinite(hTheory) & isfinite(lTheory);

    Nuse = Nuse(good);
    muse = muse(good);
    chiuse = chiuse(good);
    hFit = hFit(good);
    lFit = lFit(good);
    hTheory = hTheory(good);
    lTheory = lTheory(good);

    % Sort by N
    [Nuse, idx] = sort(Nuse);
    muse = muse(idx);
    chiuse = chiuse(idx);
    hFit = hFit(idx);
    lFit = lFit(idx);
    hTheory = hTheory(idx);
    lTheory = lTheory(idx);

    % h versus N
    figure(1)

    plot(Nuse, hFit, ...
        'Marker', markers{i}, ...
        'Color', colors(i), ...
        'LineWidth', 1.5, ...
        'MarkerSize', markerSize)

    plot(Nuse, hTheory, ...
        'Marker', 'none', ...
        'Color', colors(i), ...
        'LineWidth', 1.5, ...
        'MarkerSize', markerSize*0.7)

    figure(2)

    plot(muse, lFit, ...
        'Marker', markers{i}, ...
        'Color', colors(i), ...
        'LineWidth', 1.5, ...
        'MarkerSize', markerSize)

    plot(muse, lTheory, ...
        'Marker', '*', ...
        'Color', colors(i), ...
        'LineWidth', 1.5, ...
        'MarkerSize', markerSize*0.7)

    figure(3)

    plot(Nuse, lFit, ...
        'Marker', markers{i}, ...
        'Color', colors(i), ...
        'LineWidth', 1.5, ...
        'MarkerSize', markerSize)

    plot(Nuse, lTheory, ...
        'Marker', 'none', ...
        'Color', colors(i), ...
        'LineWidth', 1.5, ...
        'MarkerSize', markerSize*0.7)

    figure(4)

    plot(chiuse, lFit, ...
        'Marker', markers{i}, ...
        'Color', colors(i), ...
        'LineWidth', 1.5, ...
        'MarkerSize', markerSize)

    plot(chiuse, lTheory, ...
        'Marker', '*', ...
        'Color', colors(i), ...
        'LineWidth', 1.5, ...
        'MarkerSize', markerSize*0.7)

end

% Format h versus N
figure(1)
set(gca, 'XScale', 'log')
xlabel('N', 'FontSize', 18)
ylabel('h', 'FontSize', 18)
axis square
box on
ax = gca;

ax.FontSize = 20;      % tick labels
set(gca, 'TickDir', 'both')

yscale log

legend({'Mouse Brain', 'Theory', ...
        'Mouse Hippocampus', 'Theory', ...
        'Mouse Visual Cortex', 'Theory'}, ...
        'FontSize', 12, 'Location', 'best')

figure(2)
xlabel('\mu(N)', 'FontSize', 18)
ylabel('\lambda', 'FontSize', 18)
axis square
box on

ax = gca;

ax.FontSize = 20;      % tick labels
set(gca, 'TickDir', 'both')

yscale log

figure(3)
xlabel('N', 'FontSize', 18)
ylabel('\lambda', 'FontSize', 18)
axis square
box on
ax = gca;

ax.FontSize = 20;      % tick labels
set(gca, 'TickDir', 'both')

%yscale log
xscale log
ylim([0,2.5])

legend({'Mouse Brain', 'Theory', ...
        'Mouse Hippocampus', 'Theory', ...
        'Mouse Visual Cortex', 'Theory'}, ...
        'FontSize', 12, 'Location', 'best')

figure(4)
xlabel('\chi(N)', 'FontSize', 18)
ylabel('\lambda', 'FontSize', 18)
ax = gca;
ylim([0.1,2.5])

ax.FontSize = 20;      % tick labels
set(gca, 'TickDir', 'both')

axis square
box on

%yscale log
xscale log

