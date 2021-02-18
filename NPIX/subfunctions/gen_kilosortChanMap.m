function gen_kilosortChanMap (channels)
%
% generates a channel map based on specified channels for the use in kilosort 
%
% ---------------------------------------------------------------------------------------------------------------------
% SYNTAX
% gen_kilosortChanMap (chans)
%
% ---------------------------------------------------------------------------------------------------------------------
% DESCRIPTION
% NPIX_kilosortChanMap.m generates a channel map-file that is stored and can be used to spike sort data that was recorded with SpikeGLX on a
% subset of channels.
%
% inputs:
% chans      -  0-based vector specifying the channel numbers that were saved during recording
%
% outputs:
% the function will request a target folder for saving the results file
% with the filename reflecting minimum and maximum channel number that were
% recorded. 
%
% ---------------------------------------------------------------------------------------------------------------------
% author:  V. Hofmann
% last changes: 2019-11-13
% ---------------------------------------------------------------------------------------------------------------------



%%
channels = channels(:)+1;                               % bring channels into column vector
chanMap = channels-min(channels);                     % assign channel numbers
chanMap0ind = channels-min(channels)-1;               % assign channel numbers with 0 based indexing
connected = true (length(channels),1);                  % indicate which channels were connected
shankInd = ones(size(channels));

%%
% answer = inputdlg('Enter name for Channel map');    % create a name for the probe layout
% name = answer{1};
name = ['neuropixPhase3B2_kilosortChanMap_chan', num2str(min(channels)), 'to', num2str(max(channels))];


%%
% set up coordinates for full electrodes
nSitesTot = 384;
x = [43; 11; 59; 27];
x = repmat(x, nSitesTot/length(x),1);

y = [20; 0];
y = cumsum(repmat(y, nSitesTot/length(y),1));

% select relevant coordinates 
xcoords = x(channels);
ycoords = y(channels);

%%
dire = uigetdir('D:\DATA_NPIX');        % select a folder to save chan map to 
savename = [dire, '\neuropixPhase3B2_kilosortChanMap_chan', num2str(min(channels)-1), 'to', num2str(max(channels)-1), '.mat'];
save(savename, 'channels', 'chanMap', 'chanMap0ind', 'connected', 'shankInd', 'name', 'xcoords', 'ycoords')


end