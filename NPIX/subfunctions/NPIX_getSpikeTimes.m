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


%% load data 
% set correct delimiter
if ispc
    pLoc = '\';
elseif ismac
    pLoc = '/';
end

spike_clusters = readNPY([dire, pLoc, 'spike_clusters.npy']);               % from PHY
spike_clusters=double(spike_clusters);

spike_ind = readNPY([dire, pLoc, 'spike_times.npy']);                       % loads spiketimes
spike_ind = double(spike_ind);

cluster_info = tdfread([dire, pLoc, 'cluster_info.tsv']);                  % loads cluster infos

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