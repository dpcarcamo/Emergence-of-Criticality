% Parse hippocampus spike-time events into mu, chi, and N summaries.
% Spatial nearest-neighbor clusters are built from the cell position table.

clear
clc

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(scriptDir);
addpath(repoRoot)

%% Load sparse spike events and build a time-by-neuron binary matrix

spiketime = load(fullfile(scriptDir, "Binary_1485neurons_30Hz_share.csv"));

numnuerons = max(spiketime(:,1));
totaltime = max(spiketime(:,2));

spikes = zeros(totaltime, numnuerons);

for i = 1:size(spiketime, 1)
    neuronIdx = spiketime(i, 1);
    timeIdx = spiketime(i, 2);
    spikes(timeIdx, neuronIdx) = 1;
end

%% Convert binary activity to plus/minus spin variables

spiking_patterns = 2.*spikes.' - 1;

num_bins = size(spiking_patterns,2);
num_nuerons = size(spiking_patterns,1);

datacorr = spiking_patterns*spiking_patterns.'/num_bins;
datamean = mean(spiking_patterns.');

%% Load cell positions and compute spatial distances

pos = load(fullfile(scriptDir, "cell_positions.csv"));
distanceSMetric =  squareform(pdist(pos));

%% Build mu/chi/N data array from nearest-neighbor clusters

count = 1;
Data = zeros(3,numnuerons,38);
keep = 1:num_nuerons;
for j = unique([round(logspace(log10(2),log10(numnuerons), 40)), [20, 100, 1000] ])
    i = count;

    centerspins = [];
    spins = keep;    
    keep = [];
    for k = 1:round(numnuerons/j)
        
        
        spinidx = randi(length(spins));
        randspin = spins(spinidx);
        spins(spinidx) = [];
        centerspins = [centerspins, randspin];
        keep = [keep, randspin];

        distanceMetric = distanceSMetric(randspin,:);
        [sortDistances, I] = sort(distanceMetric);
        region   = (distanceMetric <=  sortDistances(j));
        
        reduceddatamean = datamean(region);
        reduceddatacorr = datacorr(region, region);

            
        nuspi = length(reduceddatamean);

        Data(1,k, i) = sum(reduceddatamean)/nuspi;
        Data(2, k, i) = (sum(reduceddatacorr,'all') - sum(reduceddatamean).^2)/nuspi;
        Data(3,k,i) = j;


    end
    count = count + 1;
    fprintf('Processed hippocampus cluster size %d\n', j);
end
%% Clean placeholder zeros

Data(Data == 0) = nan;

%% Diagnostic mu/chi plot versus cluster size

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
title('Mouse Hippocampus')
xlabel('Population Size')
ylabel('Susceptibility')
axis square

%% Fit Curie-Weiss h/lambda from averaged mu/chi curves

figure
hold on


mus = mean(squeeze(Data(1,:,:)),1, 'omitnan');
chis = mean(squeeze(Data(2,:,:)),1, 'omitnan');
Ns = mean(squeeze(Data(3,:,:)),1, 'omitnan');
hs = [];
ls = [];

for j = 1:38
    
    N = Ns(j);
    mu = mus(j);
    chi = chis(j);
    [h, lambda] = hlambda(mu,chi,N);
    hs = [hs, h];
    ls = [ls, lambda];
    

end

%% Plot h/lambda trajectory

plot(hs, ls, '-o', 'MarkerSize', 10)
axis square
xlabel('h')
ylabel('lambda')
