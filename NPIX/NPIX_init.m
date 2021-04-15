% guideline and codes for NPIX processing
% introduction to clustering in kilosort and manual curation with phy
% prior to analysis

% Maso Feb 2021
%% GET SAMPLES FOR EACH SpikeGLX t-File FROM THW META FILE
% extracts the number of samples and the samplingrate for each NPIX t-file
% this is needed if concatenated NPIX files are analyzed

% varargin = search pattern for subfolder names (i.e. 'T*')
[nSamps, SR_meta] = getNSamps_batch('ELL*');
%% GENERATE NEW CHANNEL MAP IF SAVED # OF NPIX CHANNELS IS DIFFERENT TO 385
% generates a new channel map if only a subset of channels are saved from NPIX
% chans: number of channels in subset; has to be continuous
nCh = input('enter # of NPIX channels: ');
chans = 0:nCh;
gen_kilosortChanMap(chans);
%% PRE-PROCESS RAW NPIX DATA [filtering, setting bad channels to NaN]
clc
% syntax: NPIXpreproc(channels, desired SR, original SR,flag);
% channels: number of channels in recording
% desired SR: desired SR for re-sampling
% original SR: sampling rate used in spike GLX
% flag: (0) for custom high and band pass filtering; uses EOD frequency for
%           cutoff; lets set channels without neurons to NaN and loc filter
%           settings [use for recordings in TS]
%       (1) for default car filtering 

% Nchan = numel(chans);
Nchan = 385;
% SR = 30000.554317548747; very first NPIX recording with un-calibrated probe [Maso]
% if you need to resample the raw data (i.e., for spike sorting in Spike2),
% enter new sampling rate (e.g., desSR = 10000;); if sorted in kilosort, leave
% blank)
desSR = []; % 10000;
flag = 0;
NPIXpreproc(Nchan,[],[],flag);

%% ----- GUIDELINES -------------------------------------------------------
% after pre-pocessing:
% --> make sure everything is setup properly: anaconda with python, CUDA, Visual Studio Community 2013 
% (look at their websites for instructions:
% kilosort: https://github.com/MouseLand/Kilosort
% phy: https://github.com/cortex-lab/phy and https://phy.readthedocs.io/en/latest/
% better use an earlier version of kilosort, as the latest version might
% have not been fully functional

% FIRSTLY: concatenate files (if needed or if you want to keep track of neuron IDs)
% (1) copy filtered files (output from NPIXpreproc) that should be concatenated into a new folder (name it i.e., "conc")
% (2) open cmd (terminal) and navigate to the folder "conc"
% (3) type: copy /b *.ap.bin conc.bin
%     [this concatenates all files ending with *.ap.bin into a new file
%     conc.bin; to make sure all files are being concatenated chronologically,
%     add/change numbers if neccesarry]
% (4) rename one of the meta files such that it has the same name as the concatenated
%     file and open in editor
% (5) type: "dir" in the terminal (it should still be in the path with the concatenated file)
% (6) note the new filesize of the concatenated file and use this in the
%     open metafile as new number in the line starting with "fileSizeBytes"
% (7) if data was re-sampled, you have to enter the new sampling rate in
%     line starting with "imSampRate" as well
%% ----- KILOSORT ---------------------------------------------------------
% SECONDLY: use kilosort to get isolated neurons
% (1) open matlab, navigate to the file to be sorted and type "kilosort" in the command window
%     a new window opens (kilosort GUI)
kilosort
% (2) select the file
% (3) select the correct NPIX probe (we have 3B1 staggered)
% (4) press "run" [this takes a while]
%% ------------------------------------------------------------------------
% after kilosort is done (no errors!), it saves a bunch of new files in the
% same folder
% because curation with phy changes only some of these files, it would be good
% to create a backup folder and back up these: ALL *.tsv files and spike_clusters.npy
% do the backup after EACH phy session in case merging and splitting needs
% to be undone!
%% ------------------------------------------------------------------------
% curation with phy
% (1) open a terminal (cmd) with administartor rights and navigate to the
%     folder with kilosorted data
% (2) type: "conda activate phy2" and then: "phy template-gui params.py"
%     phy opens and you can start looking at the sorted data
%% ----- SPIKE EXTRACTION AND ANALYSIS PREPARATION ------------------------
% After you are done with the spikesorting, you can start analysing the data.
% Therefore, you need to extract the spiketimes and align all the responses
% with the specific stimuli, extract the stimulus dipole and the behavioral
% responses from the associated Spike2 files.
% First thing to do is extract the spiketimes of each sorted neuron:
% check the code for some more details on what to extract
NPIX_RawData

% once the spiketimes of each cluster / neuron have been extracted, you can
% proceed with the analysis

NPIX_PreAnalysis
% this code creates structures that contain the spiketimes, behavioral responses
% as well as the corresponding stimuli
