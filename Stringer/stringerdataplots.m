% Plot Stringer h/lambda trajectories and mu/chi summaries.

load("hslstringer.mat")

%% Inspect one stored h/lambda trajectory

figure
hold on
for i = 30
    
    
    hs = hsls(1,:,i);
    ls = hsls(2,:,i);
    Ns = hsls(3,:,i);

    if sum(hs(:) < -10) > 0
        disp(i)

    end

    if hs(1) == 0
        disp(i)
    end

    plot(hs, ls, '-o', 'MarkerSize', 10)
    axis square
    xlabel('h')
    ylabel('\lambda')
end

%% Plot averaged mu and chi versus population size

Data(Data == 0) = nan;

figure
hold on
yyaxis left
plot(mean(squeeze(Data(3,:,:)),1, 'omitnan'), mean(squeeze(Data(1,:,:)),1, 'omitnan')) 

set(gca, 'XScale', 'log')
title('Mean Field as function of population size')
xlabel('Population Size')
ylabel('Mean Field')
axis square

yyaxis right
plot(mean(squeeze(Data(3,:,:)),1, 'omitnan'), mean(squeeze(Data(2,:,:)),1, 'omitnan'))


set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
title('Mouse Visual')
xlabel('Population Size')
ylabel('Susceptibility')
axis square
