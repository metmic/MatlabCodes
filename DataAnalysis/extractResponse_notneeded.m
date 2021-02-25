function [neuronResponse, behaviorRespons, stimulus] = extractResponse(input1,input2)

% this code extracts neuronal and behavioral responses from NPIX recordings
% it uses the sorted neurnal data as well as the corresponding spike2 file(s)
% and stores everything for later in depth analysis

% Maso Feb 2021
%%
clear
clc
%----- chose delimiter based on OS ----------------------------------------
if ispc
    pLoc = '\';
elseif ismac
    pLoc = '/';
end
%----- set save directory -------------------------------------------------
disp('^^^^^^ set save directory ^^^^^^')
sdir = uigetdir;
% current directory
dire = pwd;
%----- load the raw data file (i.e., sorted and curated) ------------------
%----- load the data file into "data" (easier to clear later if multiple files are analyzed)
disp('^^^^^^ get directory with spike sorted data file ^^^^^^')
dirD = uigetdir;
disp('^^^^^^ get spike sorted data file ^^^^^^')
data = load([dirD pLoc uigetfile([dirD pLoc '*.mat'])]);
% sampling interval of final files
dt = 1/2000;
% number of neurons
N = data.N;
%----- selector for analysis ----------------------------------------------
% here we set a dialogue that allows to "jump" to a specific stimulus section
fn = {'01: baseline'...
    '02: noise'...
    '03: AMs'...
    '04: envelopes'...
    '05: adaptation'...
    '06: nothing'};
DataAnalysis = listdlg('PromptString','results to plot:','SelectionMode','multiple','ListString',fn,'ListSize',[300,200]);
disp(['***** working on ' fn{DataAnalysis} ' *****'])
%%
if DataAnalysis == 1
    %----------------------------------------------------------------------
    % BASELINE DATA:
    %----------------------------------------------------------------------
    %------ Binary for Baseline Firing Rate -------------------------------
    disp('^^^^^^ extract baseline data ^^^^^^')
    % usually baseline is at the beginning of a file
    durBL = input('duration of baseline (s): '); % in sec
    % load spike2 datafile
    disp('^^^^^^ get spike2 dir ^^^^^^')
    dirS2 = uigetdir();
    disp('****** pick spike2 file with baseline ******')
    dN = dir([dirS2, pLoc, '*.mat']);
    fnN = {dN.name};
    [indx,tf] = listdlg('PromptString','Spike2:','SelectionMode','multiple','ListString',fnN);
    files = dN(indx,tf);
    clear d fn indx tf
    % loads baseline spike2 file into dataS2 structure
    dataS2 = load([dirS2, pLoc' files(1).name]);
    marker = dataS2.Ch9.times; 
    
    % on- and offset of baseline recording in SGLX file
    nStim = input('enter order of baseline recording in file: ');
    if nStim == 1
        % if baseline is the 1st recording, then we take the 1st batch
        tEdgesHere = [1,durBL+1]; % includes the extra second from spikeGLX
    end
    % initiate a binary matrix to store baseline activity
    bin_BL = zeros (N, durBL/dt);
    % initiate vector for baseline firing rates
    FR_BL = NaN(N,1);
    disp(['^^^ you have ' num2str(N) ' neurons in this dataset... ^^^'])
    % loop over neurons to extract individual binaries and firing rates 
    for n = 1:N
       idx = round(data.tSp{n}(data.tSp{n}>tEdgesHere(1,1) & data.tSp{n}<=tEdgesHere(1,2))/dt);
       idx = idx(idx>0 & idx<=length(bin_BL));
       bin_BL(n,idx) = 1;
       % get a time dependent firing rate in Hz
       [~, ~, FR_filt] = ConvolveBin (bin_BL(n,:), dt, .1, .1, .01, 2);
       FR_BL(n,1) = mean(FR_filt);
       clear idx FR_filt
    end
    clear n
    neuronResponse.binary = bin_BL;
    neuronResponse.FR = FR_BL;

    
end