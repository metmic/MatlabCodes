function [envStim, EODf, EODf_low, gain, phase, offset, scalingFactor, StimContrast] = TuningEnvelopeBN(input1, input2, input3, input4)
%----------------
% inputs:
% input1: spike2 data containing EOD trace, EOD0Xings, dipole, etc
% input2: envelope freuqencies
% input3: envelope stimuli
% input4: if stimuli where deleted
%---------------

envF = input2;clear input2
% number of cycles per envelope frequencies
cyclVect = repmat(20,numel(envF),1)';
envDurations = cyclVect.*(1./envF);

% sampling rate behavior; make sure, the channel number (Ch) corresponds
% with the EOD trace in spike2!
SR = round(1/input1.Ch8.interval);
% sampling rate glbl stim; check channel number
SRglbl = 1/input1.Ch3.interval;
% sampling rate dipole; check channel number
SRdipole = 1/input1.Ch1.interval;
% on- offsets of stimulus; this gives the total duration for all envelope
% frequencies
marker = input1.Ch31.times;
% get onsets for each envelope frequency; there is a 2sec pause between
% each envelope stimulus
disp('^^^^^^^^^^^^^ get stimulus onsets ^^^^^^^^^^^^^^^^')
envON(1) = marker(1);
for I = 2:numel(envF)
    envON(I) = envON(I-1)+2+envDurations(I-1);
end

envON(input4) = [];
envF(input4) = [];

figure;
plot(input1.Ch3.times,input1.Ch3.values,'k');hold on;plot(envON,zeros(numel(envON),1),'rd','MarkerFaceColor','r')
title('stimulus onsets')
disp(' press enter ')
pause
close all

%% compute envelope gain, phase and offset for behavior
clc
% get baseline EOD amplitude; 5sec without stimulus
disp('^^^^^^^^^^^^^ get baseline EOD amplitude ^^^^^^^^^^^^^^^^')
blEODamp = baselineEOD(input1.Ch1.values(1:round(marker(1)*SRdipole)-1)*10, 1/SRdipole);
disp('******************************************')
disp('*** compute gain/phase/offset; BEHAVIOR ***')
disp('******************************************')
disp(' ')
% compute for each single envelope frequency
for I = 1:numel(envF)
    close all
    clear EODfT EODf_lowT StimEnv temp*
    % current envelope frequency and cycle duration
    CurrEnvFreq = envF(I);
    cycl = cyclVect(I);
    disp('**************************')
    disp(['****** ' num2str(CurrEnvFreq) 'Hz ******'])
    disp('**************************')
    clear temp*
    % get EOD, global stimulus and dipole epoch for each frequency
    % EOD trace
    tempEOD = input1.Ch8.values(round(envON(I)*SR):round((envON(I)+envDurations(I))*SR)-1);
    % get mean EODf just before stimulus
    disp('^^^^^^^^^^^^^ get mean EODf before stimulus ^^^^^^^^^^^^^^^^')
    tempEODbefore = input1.Ch8.values(round((envON(I)-1.5)*SR):round(envON(I)*SR)-1);
    meanEODf = EODftimedependent(tempEODbefore, SR, SRglbl);
    [~, meanEODf] = chirpremove(meanEODf, 10, SR, 0);
    meanEODf = nanmean(meanEODf(500:2500));
    clear tempEODbefore
    
%     % global stimulus
%     tempGLBL = input1.Ch3.values(round(envON(I)*SRglbl):round((envON(I)+envDurations(I))*SRglbl)-1);
    % dipole
    tempDipole = input1.Ch1.values(round(envON(I)*SRdipole):round((envON(I)+envDurations(I))*SRdipole)-1);
    tempDipole = tempDipole.*10; % convert into mV/cm (factor 10 if dipole is 1mm apart)
    
    % get stimulius contrast, meanEODf and amplitude of dipole
    disp('^^^^^^^^^^^^^ get contrast and dipole amplitude ^^^^^^^^^^^^^^^^')
    StimPer = [1 numel(tempDipole)];
    [StimContrast(I), ~, Dipamp(I), envStimTemp] = getContrast (tempDipole, input3{I}, 1/SRdipole, blEODamp, CurrEnvFreq, StimPer);
    envStim{I} = envStimTemp;
    clear StimPer envStimTemp
    
    % downdample dipole to 2k
    tempDipole = decimate(tempDipole,SRdipole/SRglbl);
    
    % plot EOD trace and threshold to get the EOD_0Xings
    disp('^^^^^^^^^^^^^ extract EODf ^^^^^^^^^^^^^^^^')
    tempfEOD = EODftimedependent(tempEOD, SR, SRglbl);

    % chrip removal
    disp('^^^^^^^^^^^^^ remove chirps ^^^^^^^^^^^^^^^^')
    tempEODnochirps = chirpremove(tempfEOD, CurrEnvFreq, SRglbl, 1);
    tempEODnochirps = tempEODnochirps+(meanEODf-mean(tempEODnochirps));
    
    % filter fEOD
    cutoffA = CurrEnvFreq/50;
    cutoffB = CurrEnvFreq/10000;
    [B,A] = butter(2, cutoffA);
    EODfT = filtfilt(B,A,[tempEODnochirps; tempEODnochirps; tempEODnochirps]);
    EODf{I} = EODfT(numel(tempEODnochirps)+1:end-numel(tempEODnochirps));

    [B1,A1] = butter(1,cutoffB);
    EODf_lowT = filtfilt(B1,A1,[EODf{I}; EODf{I}; EODf{I}]);
    EODf_low{I} = EODf_lowT(numel(EODf{I})+1:end-numel(EODf{I}));
   
    %----------------------------------------------------------------------
    scalingFactor(I) = (Dipamp(I)); % divide gain by this
    disp('^^^^^^^^^^^^^ compute gain / phase / gain ^^^^^^^^^^^^^^^^')
    [gain(I), phase(I), offset(I)] = gainphaseoffset(input3{I},  EODf{I}, meanEODf, scalingFactor(I), CurrEnvFreq, SRglbl);
end

