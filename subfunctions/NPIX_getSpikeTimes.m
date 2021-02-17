function tSp = NPIX_getSpikeTimes (dire, clust, SR)
%
% extracts spike times from phy data
%
% ---------------------------------------------------------------------------------------------------------------------
% SYNTAX
% [tSp] = NPIX_getSpikeTimes(dire, clust, SR)
%
% ---------------------------------------------------------------------------------------------------------------------
% DESCRIPTION
% NPIX_getSpikeTimes.m extracts spike times from data that was clustered with 
% kilosort and curated with phy
%
% inputs:
% dire       -  folder of where the files are located
% clust      -  double precicion number that defines the cluster for which spike times should be extracted
% SR         -  sampling rate of original recordings (i.e. can be determined using 'read_sglxMeta.m')
%
% outputs:
% tSp       -  vector with spiketimes
%
% ---------------------------------------------------------------------------------------------------------------------
% author:  V. Hofmann
% last changes: 2019-10-07
% ---------------------------------------------------------------------------------------------------------------------


%    dire = 'D:\DATA_NPIX\2019-07-29_NPIX_ELL_RapheStim\Steps\CAR';
%    clust = 505;
%   SR = 30000;


%% load data 
    %spike_clusters = readNPY([dire, '\spike_clusters.npy']);        % [DESKTOP] from PHY
    spike_clusters = readNPY([dire, '/spike_clusters.npy']);        % [LAPTOP] from PHY
    spike_clusters=double(spike_clusters);
    
    %spike_ind = readNPY([dire, '\spike_times.npy']);              % [DESKTOP] loads spiketimes
    spike_ind = readNPY([dire, '/spike_times.npy']);              % [LAPTOP] loads spiketimes
    spike_ind = double(spike_ind);
    
    %cluster_info = tdfread ([dire, '\cluster_info.tsv']);           % [DESKTOP] loads cluster infos
    cluster_info = tdfread ([dire, '/cluster_info.tsv']);           % [LAPTOP] loads cluster infos
        
%% sort for spike times of interst

    % for i=1:length(cluster_info.group)
    %    a(i) = strcmp(cluster_info.group(i,:), 'good ');
    % end
    % a=cluster_info.sorted>=1; % for step raster & psth
    % a=cluster_info.sorted==1; % for step raster & psth
    
    %indices = cluster_info.id(clust); % clustered data 
    %indices=sort(unique(spike_clusters),'ascend'); % any detected thing

    spike_ind = spike_ind(spike_clusters == clust);
       
%% convert indices to spike times
    tSp = spike_ind*SR;

end