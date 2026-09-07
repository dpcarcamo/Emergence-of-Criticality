% Script to analyze all different iterations of networks.
clear
clc

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(scriptDir);
dataFolder = fullfile(scriptDir, 'data');
infoFile = fullfile(scriptDir, 'data_info.csv');
outputFile = fullfile(scriptDir, 'stringerhldata2.mat');

addpath(repoRoot)
addpath(dataFolder)

dataInfo = readtable(infoFile, 'TextType', 'string');

randtrial = 50;
maxSubpopulations = max(dataInfo.num_cells);
maxSizes = 45;
Data = nan(3, maxSubpopulations, maxSizes, height(dataInfo));

%% Parse recordings into mu, chi, and N arrays

for namenum = 1:height(dataInfo)

    filename = dataInfo.file_name(namenum);
    stimulusType = dataInfo.stimulus(namenum);
    fprintf('Processing Stringer row %d of %d: %s (%s)\n', ...
        namenum, height(dataInfo), filename, stimulusType)
    
    ExpData = matfile(fullfile(dataFolder, filename + ".mat"));
    
    Stim = ExpData.stim;
    
    if stimulusType == "resp"
        Response = Stim.resp;
    elseif stimulusType == "spont"
        Response = Stim.spont;
    else
        error('Unknown stimulus type "%s" for %s.', stimulusType, filename)
    end
    
    Locations = ExpData.med;
    dbinfo = ExpData.db;
    mouse_name = dbinfo.mouse_name; %#ok<NASGU>
    
    offneuron = find(mean(Response)==0);
    Response(:,offneuron) = [];
    Locations(offneuron,:) = [];
    
    [normResponse, mu, sigma ]= zscore(Response,1,1); %#ok<ASGLU>
    
    
    % Binarizing Data
    binary = (Response - mu)> 2*sigma;
    %imshow(binary)
    
    
    % Filter Data
    spiking_patterns = binary.';
    spiking_patterns = 2.*spiking_patterns - 1;
    num_bins = size(spiking_patterns,2);
    num_nuerons = size(spiking_patterns,1);
    subpopulationSizes = unique([ ...
        round(logspace(log10(2), log10(num_nuerons), 45))]);
    
    datacorr = spiking_patterns*spiking_patterns.'/num_bins;
    datamean = mean(spiking_patterns.');
    
    datacorr_pseudo = datacorr;
    
    
    
    
    % Distance
    
    distanceSMetric = squareform(pdist(Locations));
    
    %
    
    
    count = 1;
    for j = subpopulationSizes
        i = count;

        centerspins = [];
        for k = 1:round(num_nuerons/j)
            
            spins = 1:num_nuerons;
            spins(centerspins) = [];
            randspin = spins(randi(length(spins)));
            centerspins = [centerspins, randspin];
            

            distanceMetric = distanceSMetric(randspin,:);
            [sortDistances, I] = sort(distanceMetric); %#ok<ASGLU>
            region   = (distanceMetric <=  sortDistances(j));
            
            reduceddatamean = datamean(region);
            reduceddatacorr = datacorr_pseudo(region, region);
            reducedLocation = Locations(region,:); %#ok<NASGU>

            % Calculate mu and chi for each spatial subpopulation.
            nuspi = length(reduceddatamean);
            Data(1,k,i,namenum) = sum(reduceddatamean)/nuspi;
            Data(2,k,i,namenum) = (sum(reduceddatacorr,'all') - sum(reduceddatamean).^2)/nuspi;
            Data(3,k,i,namenum) = j;
           
        end
        count = count + 1
        j
    end
    
end



%% Fit Curie-Weiss h/lambda values from averaged mu/chi curves

Nss = squeeze(mean(squeeze(Data(3,:,:,:)), 1, 'omitnan'));
hs = nan(size(Nss));
ls = nan(size(Nss));

for k = 1:size(Data,4)
    fprintf('Fitting h/lambda for Stringer row %d of %d\n', k, size(Data,4))

    mus = mean(squeeze(Data(1,:,:,k)),1, 'omitnan');
    chis = mean(squeeze(Data(2,:,:,k)),1, 'omitnan');
    Ns = mean(squeeze(Data(3,:,:,k)),1, 'omitnan');
    
    for j = 1:size(Data,3)
        
        N = Ns(j);
        mu = mus(j);
        chi = chis(j);
        
        if isnan(N) || isnan(mu) || isnan(chi)
            continue
        end
        
        try
            [h, lambda] = hlambda(mu,chi,N);
            hs(j,k) = h;
            ls(j,k) = lambda;
        catch ME
            warning('Skipping h/lambda fit for row %d, size index %d: %s', ...
                k, j, ME.message)
        end
        
    end
end

%% Save processed Stringer data

save(outputFile, 'Data', 'hs', 'ls', 'Nss', 'dataInfo', '-v7.3')
fprintf('Saved processed Stringer data to %s\n', outputFile)
