% Parse Allen spike recordings into population-level mu, chi, and N arrays.
% Also fits equivalent Curie-Weiss h/lambda values from the averaged
% mu/chi curves.

clear
clc

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(scriptDir);
dataFolder = fullfile(scriptDir, 'Allen_data_share');

addpath(repoRoot);
addpath(dataFolder)

listing = dir(fullfile(dataFolder, '*.mat'));

%% Build mu/chi/N data array from nearest-neighbor spatial clusters

Data = zeros(3,5000, 40, 45); 
for namenum = 1:length(listing)

    fprintf('Processing Allen recording %d of %d\n', namenum, length(listing));

    filename = listing(namenum).name;
    load(fullfile(dataFolder, filename));
   
    xs = cell_data.anterior_posterior_ccf_coordinate;
    ys = cell_data.dorsal_ventral_ccf_coordinate;
    zs = cell_data.left_right_ccf_coordinate;
    
    spiking_patterns = 2.*X-1;
    
    num_bins = size(spiking_patterns,2);
    num_nuerons = size(spiking_patterns,1);

    if num_nuerons ~= num_cells
        warning('Skipping %s because spike matrix size does not match num_cells.', filename);
        continue
    end
    
    datacorr = spiking_patterns*spiking_patterns.'/num_bins;
    datamean = mean(spiking_patterns.');
        
    
    Locations = [xs, ys, zs];
    
    datacorr_pseudo = datacorr;    
    
    % Spatial distance matrix for nearest-neighbor cluster selection.
    distanceSMetric = squareform(pdist(Locations));
    
    
    count = 1;
    keep = 1:num_nuerons;
    for j = unique([round(logspace(log10(2),log10(num_nuerons), 40)), [20, 100, 1000]])
        if j > num_nuerons
            continue
        end
        i = count;
    
        centerspins = [];
        spins = keep;    
        keep = [];
        for k = 1:round(num_nuerons/j)
            
            
            spinidx = randi(length(spins));
            randspin = spins(spinidx);
            spins(spinidx) = [];
            centerspins = [centerspins, randspin];
            keep = [keep, randspin];
    
            distanceMetric = distanceSMetric(randspin,:);
            [sortDistances, I] = sort(distanceMetric);
            region   = (distanceMetric <=  sortDistances(j));
            
            reduceddatamean = datamean(region);
            reduceddatacorr = datacorr_pseudo(region, region);
            nuspi = length(reduceddatamean);
    
            Data(1,k,i, namenum) = sum(reduceddatamean)/nuspi;
            Data(2, k, i, namenum) = (sum(reduceddatacorr,'all') - sum(reduceddatamean).^2)/nuspi;
            Data(3,k,i, namenum) = j;
    
    
        end
        count = count + 1;
        j;
    end

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
plot(mean(squeeze(Data(3,:,:)),1, 'omitnan'), mean(squeeze(Data(2,:,:)),1, 'omitnan'), '-o')


set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
title('Mouse Visual')
xlabel('Population Size')
ylabel('Susceptibility')
axis square

%% Fit Curie-Weiss h/lambda for each averaged Allen curve

hs = zeros(size(Data,3), size(Data,4));
ls = zeros(size(Data,3), size(Data,4));

figure
hold on

for k = 1:size(Data,4)
    mus = mean(squeeze(Data(1,:,:,k)),1, 'omitnan');
    chis = mean(squeeze(Data(2,:,:,k)),1, 'omitnan');
    Ns = mean(squeeze(Data(3,:,:,k)),1, 'omitnan');
    for j = 1:size(Data,3)
        
        N = Ns(j);
        mu = mus(j);
        chi = chis(j);
        if isnan(N)
            continue
        end
        [h, lambda] = hlambda(mu,chi,N);
        hs(j,k) = h;
        ls(j,k) = lambda;
        
    end
end

%% Diagnostic h/lambda plots

Nss = squeeze(mean(squeeze(Data(3,:,:,:)), 1, 'omitnan'));

Nss(Nss==0)=nan;
hs(hs==0) = nan;

valid = Nss < 1200;

figure
plot(hs(valid), ls(valid), '-b', 'MarkerSize', 10)
axis square
xlabel('h')
ylabel('\lambda')


figure
plot(Nss(valid), hs(valid), '-b', 'MarkerSize', 10)
axis square
xlabel('N')
ylabel('h')


figure

plot(Nss(valid), ls(valid), '-b', 'MarkerSize', 10)
axis square
xlabel('N')
ylabel('\lambda')

