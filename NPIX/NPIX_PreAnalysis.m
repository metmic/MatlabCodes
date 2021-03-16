% NPIX_PreAnalysis
%
% compiled script for analysis of some NPIX data
%
% ---------------------------------------------------------------------------------------------------------------------
% SYNTAX
% NPIX_Analysis
%
% ---------------------------------------------------------------------------------------------------------------------
% DESCRIPTION
% NPIX_Analysis.m calculates firing properties of neuron clusters sorted in phy2 under baseline conditions, and in response
% to certain stimuli
% ---------------------------------------------------------------------------------------------------------------------
% REQUIRED SUBFUNCTIONS
% readShankMap.m
% SpikeCountCorr.m
% distFunc.m
% STA.m
% PhaseHist.m
% tVarCorr.m
% FRvsFR.m
% reSample.m
%
% ---------------------------------------------------------------------------------------------------------------------
% author:  V. Hofmann, M. Haggard
% update: M. Metzen
% last changes: 2021-02-15
% ---------------------------------------------------------------------------------------------------------------------
%%
clear
clc
% for analysis, it might be important to use the correct delimiter based on the OS
%----- chose back- or forward slash based on OS ----------------------------------------
if ispc
    pLoc = '\';
elseif ismac
    pLoc = '/';
end
%----- set save directory -------------------------------------------------
disp('^^^^^^ set save directory ^^^^^^')
sdir = uigetdir;
% gets the current directory
dire = pwd;
%----- load the raw data file (i.e., sorted and curated) ------------------
%----- load the data file into "data" (easier to clear later if multiple files are analyzed)
disp('^^^^^^ get directory with spike extracted data file ^^^^^^')
dirD = uigetdir;
disp('^^^^^^ get extracted spikes file ^^^^^^')
data = load([dirD pLoc uigetfile([dirD pLoc '*_extracted.mat'])]);
% sampling interval of final files (usually set to 2k)
dt = 1/2000;
% number of extracted neurons in current datafile
N = data.N;
%----- selector for analysis ----------------------------------------------
% here we set a dialogue that allows to "jump" to a specific stimulus section
fn = {'01: baseline'...
    '02: noise'...
    '03: AMs'...
    '04: envelopes'...
    '05: adaptation'...
    '06: nothing'};
DataAnalysis = listdlg('PromptString','stimulus to process:','SelectionMode','multiple','ListString',fn,'ListSize',[150,200]);
%%
while DataAnalysis<numel(fn)
    
    if DataAnalysis == 1
        disp(['***** working on ' fn{DataAnalysis} ' *****'])
        %----------------------------------------------------------------------
        % BASELINE DATA (no stimulus!):
        %----------------------------------------------------------------------
        %------ Binary for Baseline Firing Rate -------------------------------
        disp('^^^^^^ get baseline data ^^^^^^')
        % usually baseline is at the beginning of a file
        durBL = input('duration of baseline (s): ');                            % in sec
        sglxT = 1;                                                              % enter the extra time from SpikeGLX in sec
        % on- and offset of baseline recording
        StimInd = input('enter index of baseline recording in current dataset (1 if it is the first): ');
        if StimInd == 1
            startP = sglxT;
        end
        
        tEdgesHere = [startP,durBL+startP];                                               % includes the extra time from spikeGLX
        % initiate a binary matrix to store baseline activity
        bin_BL = zeros (N, round(tEdgesHere(1,2)/dt));
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
        
        % euclidian distances between relevant channels
        combSp = nchoosek(1:N, 2);                                              % all possible pairwise combinations of neurons
        disp('^^^^^^ distance between neurons ^^^^^^')
        [~, Pos] = readShankMap (data.meta);                                    % read the imec shank map from meta file
        distBL = distFunc(combSp, data.Ch, Pos);
        % get neuron IDs and info
        Ch = data.Ch;
        info = data.info;
        
        name = input('set save name: ','s');
        disp('saving dataset baseline')
        save([sdir pLoc name '_baseline.mat'],'bin_BL','combSp','FR_BL','Pos','distBL','durBL','dt','dirD','info','Ch')
        clearvars -except pLoc DataAnalysis fn data sdir dire dirD dt N
        DataAnalysis = listdlg('PromptString','results to plot:','SelectionMode','multiple','ListString',fn,'ListSize',[150,200]);
        %%
    elseif DataAnalysis == 2
        disp(['***** working on ' fn{DataAnalysis} ' *****'])
        %----------------------------------------------------------------------
        % NOISE STIMULATION:
        %----------------------------------------------------------------------
        noiseDur = input('duration of noise (s): ');                            % in sec
        sglxT = 1;                                                              % enter the extra time from SpikeGLX in sec
        % this part extracts spiketimes and the corresponding noise stimulus
        % and does some initial analysis (STA, ON- OFF-type)
        disp('^^^^^^ process white noise data ^^^^^^')
        % load spike2 datafile
        disp('^^^^^^ get spike2 dir ^^^^^^')
        dirS2 = uigetdir();
        disp('****** pick noise spike2 file ******')
        dN = dir([dirS2, pLoc, '*.mat']);
        fnN = {dN.name};
        [indx,tf] = listdlg('PromptString','Spike2:','SelectionMode','multiple','ListString',fnN);
        files = dN(indx,tf);
        clear d fnN indx tf
        dataS2 = load([dirS2, pLoc' files(1).name]);
        
        %     figure(007); plot(dataS2.Ch3.times, dataS2.Ch3.values)
        %     pause
        %     close 007                                                          % duration of noise stimulus; in sec
        tStim = dataS2.Ch9.times;                                                          % get stim start times in spike2 using the stim triggers
        nStim = numel(tStim);
        idxStim = round(tStim(nStim)/dt:(tStim(nStim)+noiseDur)/dt);
        idxStim(end) = [];
        vStim = dataS2.Ch3.values(idxStim+1);                                              % extract stimulus waveform
        dipole = dataS2.Ch1.values(round(tStim(nStim)/dataS2.Ch3.interval):round((tStim(nStim)+noiseDur)/dataS2.Ch3.interval)-1);
        figure(007); plot(((1:length(vStim)).*dataS2.Ch3.interval)-dataS2.Ch3.interval, vStim)    % visualizes the noise stimulation as function of time
        pause
        close 007
        
        tStimNew = tStim-tStim(1);
        clear tStim
        tEdges_here = [(tStimNew(nStim))+sglxT , (tStimNew(nStim))+noiseDur+sglxT];         % the extra time (sglxT) is because SpikeGLX adds a second to the beginning/end of each file
        
        preTrg = 0.1;                                                               % pre-trigger time for STA
        postTrg = 0.1;                                                              % post-trigger time for STA
        
        Type = NaN (1,N);                                                           % preallocate memory
        vSTA = NaN (N, (preTrg+postTrg)/dataS2.Ch3.interval+1);                            % preallocate memory
        bin_noise = zeros(N,noiseDur/dt);
        disp('***** STA analysis *****')
        for n = 1:N
            tSp_here = data.tSp{n}(data.tSp{n}>tEdges_here(1) & data.tSp{n}<tEdges_here(2))-tEdges_here(1);            % extract relevant spikes during stimulus period
            idx = round(tSp_here/dt);
            idx = idx(idx>0 & idx<=length(bin_noise));
            bin_noise(n,idx) = 1;
            % [vSTA, tSTA, Type, SL] = STA (spT, vStim, dtStim, varargin)
            [vSTA(n,:), tSTA, Type(n)] = STA(tSp_here, vStim, dataS2.Ch3.interval, 1, preTrg, postTrg, 0.01, 0.007);   % calculate STA & evaluate response type (On- vs Off-type)
            sgtitle(['displaying STA of neuron # ' num2str(n)])
            pause
            close
            clear tSp_here idx
        end
        disp(['***** there are ' num2str(sum(Type == 0)) ' OFF type and ' num2str(sum(Type == 1)) ' ON type neurons in the dataset *****'])
        clear NPIX gblstim nFile nStim tStim idxStim tEdges_here preTrg postTrg n
        
        %----- Determine pair type, same type pairs: 1, opposite type pairs: 0 ----
        combSp = nchoosek(1:N, 2);
        [~, Pos] = readShankMap (data.meta);
        pairType = NaN(length(combSp),1);
        for c = 1:length(combSp)
            if (Type(combSp(c,1))==0 && Type(combSp(c,2))==0) || (Type(combSp(c,1))==1 && Type(combSp(c,2))==1)
                pairType(c,1) = 1;
            else
                pairType(c,1) = 0;
            end
        end
        % Indices in combSp and r_n etc. for same and opposite type pairs
        sameType = find(pairType==1);
        oppType = find(pairType==0);
        clear c
        % %     %----- Visualize all STAs -------------------------------------------------
        % %     figure; hold on; box on;
        % %     plot(tSTA, vSTA); plot([0 0],[-1.5 1.5], '--k'); plot([-0.05 0.05],[0 0], '--k');
        % %     leg = strsplit(num2str(1:N));
        % %     title('STA', 'FontSize', 18);
        
        % get neuron IDs and info
        Ch = data.Ch;
        info = data.info;
        name = input('chose save name noise: ','s');
        disp('saving the noise dataset')
        save([sdir pLoc name '_noise.mat'],'bin_noise','tSTA','combSp','noiseDur','dipole','oppType','pairType','sameType','Type','vSTA','vStim','dt','Pos','dirD','info','Ch')
        
        clearvars -except pLoc DataAnalysis fn data sdir dire dirD dt N Type
        DataAnalysis = listdlg('PromptString','results to plot:','SelectionMode','multiple','ListString',fn,'ListSize',[150,200]);
        %%
    elseif DataAnalysis == 3
        %----------------------------------------------------------------------
        % AMPLITUDE MODULATIONS
        %----------------------------------------------------------------------
        % this code extracts spiketimes, behavior and associated AM stimuli and
        % stores them in the structure "AMtim"
        % it also does some initial analysis: VS, phase per neuron per AMf
        sglxT = 1;                                                              % enter the extra time from SpikeGLX in sec
        disp('^^^^^^ analysing AMs ^^^^^^')
        disp('^^^^^ get spike2 dir ^^^^^')
        dirS2 = uigetdir();
        disp('pick AM spike2 file')
        d = dir([dirS2, pLoc, '*.mat']);
        fnAM = {d.name};
        [indx,tf] = listdlg('PromptString','spike2 files with AMs:','SelectionMode','multiple','ListString',fnAM);
        files = d(indx,tf);
        clear d fnAM indx tf
        dataS2 = load([dirS2, pLoc' files(1).name]);
        
        % generate vector with AM frequencies used during recording
        AMf = logspace(log10(1),log10(256),9);
        % if the AM series are between baseline and the noise and contain 9
        % different frequencies, then the indices within the spike2 files are
        % 2 - 11;
        nStim = 2:11;
        tStim = dataS2.Ch9.times;                                               % get stim start times in spike2 using the stim triggers
        tStimNpix = tStim-tStim(1);
        
        % loop over AM stimuli to extract dipole, global stim, behavior, spikes
        for I = 1:numel(AMf)
            % initiate binary matrix for current AM frequency
            tEdges_here = [(tStimNpix(nStim(I)))+sglxT , (tStimNpix(nStim(I+1))-2)+sglxT];  % the extra sglxT time is because SpikeGLX adds a second to the beginning/end of each file
            bin_AM = zeros(N,ceil(abs(diff(tEdges_here))/dt)-1);
            
            % extract dipole signal and time
            Tdipole = dataS2.Ch1.values(round(tStim(nStim(I))/dataS2.Ch1.interval):round((tStim(nStim(I+1))-2)/dataS2.Ch1.interval)-1);                           % get dipole recording
            Tdipoletime = dataS2.Ch1.times(round(tStim(nStim(I))/dataS2.Ch1.interval):round((tStim(nStim(I+1))-2)/dataS2.Ch1.interval)-1);                       % get time vector of dipole
            Tdipoletime = Tdipoletime-Tdipoletime(1);
            % extract "generated" stimulus and time
            Tstim = dataS2.Ch3.values(round(tStim(nStim(I))/dataS2.Ch3.interval):round((tStim(nStim(I+1))-2)/dataS2.Ch3.interval)-1);                           % get dipole recording
            Tstimtime = dataS2.Ch3.times(round(tStim(nStim(I))/dataS2.Ch3.interval):round((tStim(nStim(I+1))-2)/dataS2.Ch3.interval)-1);                       % get time vector of dipole
            Tstimtime = Tstimtime-Tstimtime(1);
            % extract EODf (aka, behavior) and time
            TEOD = dataS2.Ch8.values(round(tStim(nStim(I))/dataS2.Ch8.interval):round((tStim(nStim(I+1))-2)/dataS2.Ch8.interval)-1);                              % get EOD trace
            TEODtime = dataS2.Ch8.times(round(tStim(nStim(I))/dataS2.Ch8.interval):round((tStim(nStim(I+1))-2)/dataS2.Ch8.interval)-1);                           % get time vector of EOD
            TEODtime = TEODtime-TEODtime(1);
            
            % get 0xings of AM stimulus
            mpd = round((1/dt/AMf(I))-(1/dt/AMf(I)*0.2)); % determine the cycleduration and subtract 10% as mi dist between 2 succ. peaks
            [peaks,t0Xings] = findpeaks(Tstim,'minpeakdistance',mpd,'minpeakheight',0.5);
            t0Xings = t0Xings(peaks>mean(Tstim));
            t0Xings = t0Xings(t0Xings<length(Tstim));
            t0Xings = (t0Xings*dt);
            
            % loop over neurons
            for n = 1:N
                tSp_here = data.tSp{n}(data.tSp{n}>tEdges_here(1) & data.tSp{n}<tEdges_here(2))-tEdges_here(1);            % extract relevant spikes during stimulus period
                idx = round(tSp_here/dt);
                idx = idx(idx>0 & idx<=length(bin_AM));
                bin_AM(n,idx) = 1;
                clear tSp_here idx
            end
            AMstim{I}.tEdges = tEdges_here;
            AMstim{I}.binary = bin_AM;
            AMstim{I}.dipole = Tdipole;
            AMstim{I}.t_dipole = Tdipoletime;
            AMstim{I}.stim = Tstim;
            AMstim{I}.t_stim = Tstimtime;
            AMstim{I}.EOD = TEOD;
            AMstim{I}.t_EOD = TEODtime;
            AMstim{I}.t0Xings = t0Xings;
            
            clear tEdges_here bin_AM Tdipole Tdipoletime TEOD TEODtime t0Xings mpd
        end
        
        disp('****** compute phase, vector strength, and stats for each neuron ******')
        binwidth = 30;
        warning off
        for I = 1:length(AMf)
            phasesN = NaN(500, N);
            VSsN = NaN(1, N);
            FR_modN = NaN(1,N);
            for n = 1:N
                % Phase Histogram and Vector Strength
                tSp_here = data.tSp{n}(data.tSp{n}>AMstim{1,I}.tEdges(1)...
                    & data.tSp{n}<AMstim{1,I}.tEdges(2))-AMstim{1,I}.tEdges(1); % extract relevant spikes during stimulus period
                
                [vS, vPhi, ~, ~, yHist] = phase_histogram(tSp_here, AMstim{1,I}.t0Xings, binwidth);
                % FR Modulation
                ySp = yHist./(binwidth* (numel(AMstim{1,I}.t0Xings)-1));
                %             % Bootstrapping to get SEM of phase and VS at each location
                %             nboot = 100; % Use test_bootstrap.m to test for how many boot straps is sufficient.
                %             boot_phases = NaN(nboot,1);
                %             boot_VSs = NaN(nboot,1);
                %             boot_FR_mod = NaN(nboot,1);
                
                %             for i = 1:nboot
                %                 resamp_tSp = reSample(tSp_here,AMstim{1,I}.t0Xings);
                %                 [vPhi, vS, ~,yHist, ~, ~] = PhaseHist (resamp_tSp, AMstim{1,I}.t0Xings);
                %                 binwidth = 0.25/length(yHist);                                  % binwidth in PhaseHist.m function - duration of one trial (s) / # bins
                %                 boot_phases(i,1) = vPhi;
                %
                %                 boot_VSs(i,1) = vS;
                %                 ySp = yHist./(binwidth*n_trials);
                %                 boot_FR_mod(i,1) = max(ySp)-min(ySp);                           % FR modulation (spikes/second)
                %                 clear resamp_tSp vPhi vS yHist ySp
                %             end
                
                %             phases(n,I) = mean(boot_phases);   phases_SEM(n,I) = std(boot_phases);
                %             VSs(n,I) = mean(boot_VSs);         VS_SEM(n,I) = std(boot_VSs);
                %             FR_mod(n, I) = mean(boot_FR_mod);  FR_mod_SEM(n,I)= std(boot_FR_mod);
                
                phasesN(1:length(vPhi),n) = vPhi;
                VSsN(:,n) = vS;
                FR_modN(:,n) = max(ySp)-min(ySp);
                clear vPhi vS ySp tSp_here yHist
            end
            AMstim{I}.phases = phasesN;
            AMstim{I}.VS = VSsN;
            AMstim{I}.FR_mod = FR_modN;
            clear phasesN VSsN FR_modN
        end
        [~, Pos] = readShankMap (data.meta);
        % get neuron IDs and info
        Ch = data.Ch;
        info = data.info;
        name = input('set savename for AMs: ','s');
        disp('saving AM dataset')
        save([sdir pLoc name '_AMs.mat'],'AMstim','Type','AMf','Pos','dirD','dt','info','Ch')
        
        clearvars -except pLoc DataAnalysis fn data sdir dire dirD dt N Type
        DataAnalysis = listdlg('PromptString','results to plot:','SelectionMode','multiple','ListString',fn,'ListSize',[300,200]);
        %%
    elseif DataAnalysis == 4
        disp(['***** working on ' fn{DataAnalysis} ' *****'])
        %----------------------------------------------------------------------
        % ENVELOPES
        %----------------------------------------------------------------------
        sglxT = 1;                                                              % enter the extra time from SpikeGLX in sec
        indEnvStim = input('enter position of envelope stimuli (i.e.: [2 3]): ');
        quest = input('sorted with KS (1 for YES)? ');
        questEnv = input('how many times were all envelopes repeated? ');   % if two or more envelope freuquency runs are analyzed, enter the number here
        % load envelope stimuli
        disp('***** load envelope stimulus file *****')
        dirStim = uigetdir();
        load([dirStim, pLoc, 'envelope_stims.mat']);
        
        disp('***** get spike2 dir *****')
        dirS2 = uigetdir();
        disp('***** pick spike2 file(s) of concatenated analysis file *****')
        d = dir([dirS2, pLoc, '*.mat']);
        fnE = {d.name};
        [indx,tf] = listdlg('PromptString','spike2 files:','SelectionMode','multiple','ListString',fnE);
        files = d(indx,tf);
        clear d fnE indx tf
        for I = 1:size(files,1)
            dataS2{I} = load([dirS2, pLoc' files(I).name]);
        end
        Rec = input('position of 1st envelope NPIX recording: '); % enter position in concat file of first envelope run
        
        % envelope frequencies (in Hz)
        envF = [0.05 0.1 0.25 0.5 0.75 1];
        % number of cycles per envelope frequencies
        cyclVect = repmat(20,numel(envF),1)';
        envDurations = cyclVect.*(1./envF);
        
        % get single envelope stimuli in separate cells and for each env frequency run
        envB = [];
        for I = 1:numel(envF)
            if I == 1
                temp = env(1:cyclVect(I)*1/envF(I)/dt);
            else
                start = sum(envDurations(1:I-1))+(2*(I-1));
                temp = env((start/dt)+1:(start/dt)+1 + ((cyclVect(I)*1/envF(I)/dt)-1));
            end
            envStim{I} = temp;
            envB = [envB; temp; zeros(2*2000,1)];
            clear temp
        end
        clear start
        
        nStim = 1:numel(envF)+1;
        tStim = [];
        for I = indEnvStim
            tStimT = dataS2{I}.Ch9.times;                                                % get stim start times in spike2 using the stim triggers
            tStim = [tStim tStimT];
            clear tStimT
        end
        
        % check if all envelope onsets are correct and delete false ones if neccessary
        q = 1;
        figure; tiledlayout(questEnv,1)
        for I = indEnvStim
            nexttile
            plot(dataS2{I}.Ch3.times,dataS2{I}.Ch3.values);hold on;plot(tStim(:,q),zeros(size(tStim,1),1),'*')
            q = q+1;
        end
        delet = [];
        questB = input('Delete timepoints (YES = 1)? ');
        close all
        if questB == 1
            delet = input('enter timepoint indices to delete; format [number number]: ');
            tStim(delet) = [];
            q = 1;
            figure; tiledlayout(questEnv,1)
            for I = indEnvStim
                nexttile
                plot(dataS2{I}.Ch3.times,dataS2{I}.Ch3.values);hold on;plot(tStim(:,q),zeros(size(tStim,1),1),'*')
                q = q+1;
            end
            pause
            close
        end
        
        % add end-time
        q = 1;
        for I = indEnvStim
            tStimT(:,q) = [tStim(:,q); dataS2{I}.Ch3.times(end)];
            tStimEnv(:,q) = tStimT(:,q)-tStimT(1,q);
            q = q+1;
        end
        tStim = tStimT;clear tStimT
         
        % add time of prev NPIX recordings due to concatenation
        RecOri = Rec;
        for I = 1:questEnv
            if quest == 1
                if Rec == 1
                    pretime = 0;
                else
                    pretime = data.nSampsTotal(Rec-1)./str2double(data.meta.imSampRate);
                end
                tStimNeurons(:,I) = tStimEnv(:,I)+data.tEdges(Rec,1);
            else
                tStimNeurons(:,I) = tStimEnv(:,I); % extra second for single t-files (for data that is NOT concatenated)
            end
            Rec = Rec+1;
        end
        Rec = RecOri;
        
        % correct for start point of next envelope if one envelope was faulty and delete the faulty envelope
        tStimNeuronsCorrected = tStimNeurons;
        q = 1;
        for J = indEnvStim
            % behavior
            disp('^^^^^^ conpute behavioral gain, phase, offset and get scaling factor and contrast ^^^^^^')
            [~, behavior.EODf, behavior.EODf_low, behavior.gain, behavior.phase, behavior.offset, scalingFactor{q}, StimContrast{q}] = ...
                TuningEnvelopeBN(dataS2{J}, envF, envStim, delet);
            
% %             if questB == 1
% %                 envFD = envF(delet);
% %                 tStimEnvcorrected(delet+1:end) = tStimEnvcorrected(delet+1:end)-(10*(1/envFD));
% %                 stimend = abs(length(env)/2000-tStimEnvcorrected(7));
% %                 tStimEnvcorrected(7) = tStimEnvcorrected(7)-stimend;
% %                 
% %                 tStimcorrected(delet) = [];
% %                 tStimEnvcorrected(delet) = [];
% %                 tStimNeuronsCorrected(delet) = [];
% %                 envF(delet) = [];
% %                 cycldur(delet) = [];
% %                 neurons{delet} = [];
% %                 dataS2.Ch9.times(delet) = [];
% %             end
            
            % extract binaries for each envelope frequency and compute gain, phase, VS and zscore
            disp('^^^^^^ conpute neuronal gain, phase, VS and zscore ^^^^^^')
            for I = 1:numel(envF)
                cycl = cyclVect(I);
                disp([' working on envelope frequency ' num2str(envF(I)) 'Hz'])
                tEdges_here = [];
                tEdges_here = [(tStimNeuronsCorrected(nStim(I),q))+sglxT , (tStimNeuronsCorrected(nStim(I),q)+(cycl*(1/envF(I))))+sglxT]; % the extra 1 is because SpikeGLX adds a second to the beginning/end of each file
                bin_envelope = zeros(N,ceil(abs(diff(tEdges_here))/dt));
                % loop over neurons
                for n = 1:N
                    tSp_here = data.tSp{n}(data.tSp{n}>tEdges_here(1) & data.tSp{n}<tEdges_here(2))-tEdges_here(1);            % extract relevant spikes during stimulus period
                    idx = round(tSp_here/dt);
                    idx = idx(idx>0 & idx<=length(bin_envelope));
                    bin_envelope(n,idx) = 1;
                    clear tSp_here idx
                end
                if size(bin_envelope,2) > ceil((1/envF(I))*cycl/dt)
                    bin_envelope(:,ceil( (1/envF(I))*cycl/dt )+1:end) = [];
                end
                
                disp('^^^^^^ tuning ^^^^^^')
                for II = 1:N
                    output = EnvTuning(envStim{I}, bin_envelope(II,:), envF(I), 1/dt, cycl, scalingFactor{1,q}(I), 0);
                    neurons{I}.zscore(:,II) = output(2:end,14);
                    neurons{I}.gain(:,II) = output(2:end,5);
                    %----- if bursts and isolated spikes are needed -----------------------
                    %         neurons{I}.gainB(:,II) = outputM(2:end,6);
                    %         neurons{I}.gainI(I,II) = outputM(2:end,7);
                    
                    neurons{I}.phase(:,II) = output(2:end,8);
                    %         neurons{I}.phaseB(:,II) = outputM(2:end,9);
                    %         neurons{I}.phaseI(:,II) = outputM(2:end,10);
                    
                    neurons{I}.VS(:,II) = output(2:end,11);
                    %         neurons{I}.VSB(:,II) = outputM(2:end,12);
                    %         neurons{I}.VSI(:,II) = outputM(2:end,13);
                    clear output
                end
                neurons{I}.binary = bin_envelope;
                clear tEdges_heret0Xings mpd Tenvelope bin_envelope
            end
            
            % get positions of channels
            [~, EnvRun{q}.Pos] = readShankMap (data.meta);
            % get neuron IDs and info
            EnvRun{q}.Ch = data.Ch;
            EnvRun{q}.info = data.info;
            
            EnvRun{q}.neurons = neurons;
            EnvRun{q}.behavior = behavior;
            clear neurons behavior
            q = q+1;
        end
        
        disp('***** saving envelope dataset *****')
        name = input('enter save name: ','s');
        save([sdir pLoc name '_envelopes.mat'],'EnvRun','envF','envStim','scalingFactor','StimContrast','dirD','dt')
        clearvars -except pLoc DataAnalysis fn data sdir dire dirD dt N Type
        DataAnalysis = listdlg('PromptString','results to plot:','SelectionMode','multiple','ListString',fn,'ListSize',[150,200]);
        %%
    end
    
end
disp('^^^^^^ nothing more to do... ^^^^^^')