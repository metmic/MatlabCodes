function NPIX_RawData
%
% Converts Phy and spikeGLX outputs to .mat format
%
% -------------------------------------------------------------------------
% DESCRIPTION
% NPIX_RawData extracts the NPIX data from spikeGLX and Phy files and saves
% it for analysis.
%
% -------------------------------------------------------------------------
% REQUIRED SUBFUNCTIONS
% read_sglxMeta.m
% NPIX_getClusterInfo.m
% NPIX_getSpikeTimes.m
%
% -------------------------------------------------------------------------
% author:  V. Hofmann, M. Haggard
% update: M. Metzen
% last changes: 2021-02-15
% -------------------------------------------------------------------------
%%
clc
close all
clear
% for analysis, it might be important to use the correct delimiter based on the OS
%----- chose back- or forward slash based on OS ----------------------------------------
if ispc
    pLoc = '\';
elseif ismac
    pLoc = '/';
end

% set save directory
disp('^^^^^^ set save dir ^^^^^^')
sdire = uigetdir();
% get the directory with the SpikeGLX bin files (here it is the folder with
% the kilosorted data
disp('^^^^^^ get kilosort dir ^^^^^^')
bdir = uigetdir();
FileName = dir([bdir, pLoc, '*.bin']);
% current dirctory path
dire = pwd; 
% get sampling rate from SGLX recording
[meta] = read_sglxMeta([bdir, pLoc, FileName(1).name]);                     % get meta data from SpikeGLX recording
SR = str2double(meta.imSampRate);                                           % get imec sample rate

% get cluster info (connect ID and channel info from phy to STAs/neuron IDs)
info = NPIX_getClusterInfo(bdir);                                           % read in info about sorted clusters ('cluster view')

% get spike times
disp('^^^^^^ get spike times ^^^^^^')
tSp = cell(1,length(info.id));                                              % pre-allocate variable
for n = 1:length(info.id)
    tSp{n} = NPIX_getSpikeTimes (bdir, info.id(n), 1/SR);                   % extract spike times and store in tSP
end

% if a label (i.e., 1 for perfectly sorted; 2 for well sorted; 3 for crappy sorted, etc)
% was used during curation with phy, this can now be used to extract neurons of certain quality; adjust as needed 

% take neurons with only a certain label
% tSp = tSp(info.sorted==1);                                                % select for high quality sorted clusters based on marker mannually assigned in phy2
% Ch = info.channel(info.sorted==1);                                        % main channel # for each neuron for STA section
% take neurons with more than one lables
tSp = tSp(info.sorted>=1 & info.sorted<=3);
Ch = info.channel(info.sorted>=1 & info.sorted<=3);                         % main channel # for each neuron for STA section

% re-order ascending channel: deep to superficial / caudal to rostral
% depending on NPIX placement
[Ch, idxCh] = sort(Ch);
tSp = tSp(idxCh);
info.id = id(idxCh);
info.sorted = info.sorted(idxCh);

N = length(tSp);                                                            % number of high quality units
disp(['you have ' num2str(N) ' neurons in this dataset'])
% Impose a forced refractory period 
refract = 0.002;                                                            % duration of imposed refractory period in [s]
disp('imposing forced refractory period...')
for n = 1:N
    tSp{n} = tSp{n}(diff(tSp{n})>=refract);
end
clear refract n

%% get recording durations to evaluate individual stimuli
disp('^^^^^^ get individual stimulus onsets ^^^^^^')
% if t-files have been concatenated, we need to know the durations of these
% single t-files in order to have the exact start times for each stimulus block
% (stored in tEdges); make sure to only take the bin files that have been
% concatenated!
disp('^^^^^^ get original SpikeGLX bin dir ^^^^^^')
dirBin = uigetdir();
% files contained in directory with extension *.bin
dB = dir([dirBin, pLoc, '*.bin']);
% tool to select those bin files that have been concatenated
fnB = {dB.name};
[indx,tf] = listdlg('PromptString','select BIN files that have been concatenated','SelectionMode','multiple','ListString',fnB);
files = dB(indx,tf);
clear dB fnB indx tf

nChansTotal = str2double(meta.nSavedChans);                                 % total recording channels stored in file 

nSampsTotal = NaN(1,length(files));                                         % pre-allocate
for I = 1:length(files)
    d = dir([dirBin,pLoc, files(I).name]);                                  % get file size
    nSampsTotal(I) = d.bytes/nChansTotal/2;                                 % calculate number of samples per channel        
end
dur = nSampsTotal./SR;                                                      % calculate duration of recordings
tEdges = [[0,cumsum(dur(1:end-1))]', cumsum(dur)'];                         % calculate start and end times for each file in concatenated bin (s)

clear files nChansTotal I d dur SR

%% Save for analysis
disp('^^^^^^ save extracted data ^^^^^^')
name = input('choose save name: ','s');
save([sdire, pLoc, name '_extracted.mat'], 'Ch', 'info', 'meta', 'N', 'tEdges', 'tSp', 'nSampsTotal')
