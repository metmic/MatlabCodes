% TUNING BEHAVIOR
% code to compute behavioral gain to different envlope frequencies

% Maso 02/2021
%% initial settings and load matlab converted spike2 data
clear 
close all
clc
% use delimiter based on OS
if ispc
    pLoc = '\';
elseif ismac
    pLoc = '/';
end
% envelope frequencies (in Hz)
envF = [0.05 0.1 0.25 0.5 0.75 1];
% number of cycles per envelope frequencies
cyclVect = repmat(20,numel(envF),1)';
envDurations = cyclVect.*(1./envF);
% load the *.mat converted spike2 file
[file, curDir] = uigetfile('*.mat');
data = load([curDir, file]);
% savename
% name = file(end-8:end-4);
name = '_pre';
% sampling rate behavior; make sure, the channel number (Ch) corresponds
% with the EOD trace in spike2!
SR = round(1/data.Ch8.interval);
% sampling rate glbl stim; check channel number
SRglbl = 1/data.Ch3.interval;
% sampling rate dipole; check channel number
SRdipole = 1/data.Ch1.interval;
% on- offsets of stimulus; this gives the total duration for all envelope
% frequencies
marker = data.Ch31.times;
% get onsets for each envelope frequency; there is a 2sec pause between
% each envelope stimulus
disp('^^^^^^^^^^^^^ get stimulus onsets ^^^^^^^^^^^^^^^^')
envON(1) = marker(1);
for I = 2:numel(envF)
    envON(I) = envON(I-1)+2+envDurations(I-1);
end
figure;
plot(data.Ch3.times,data.Ch3.values,'k');hold on;plot(envON,zeros(numel(envON),1),'rd','MarkerFaceColor','r')
title('stimuli onsets')
disp(' press enter ')
pause
close all

%% compute envelope gain and phase for behavior
clc
% get baseline EOD amplitude; 5sec without stimulus
disp('^^^^^^^^^^^^^ get baseline EOD amplitude ^^^^^^^^^^^^^^^^')
blEODamp = baselineEOD(data.Ch1.values(1:round(marker(1)*SRdipole)-1)*10, 1/SRdipole);
disp('******************************************')
disp('*** compute gain/phase/offset; BEHAVIOR ***')
disp('******************************************')
disp(' ')
% compute for each single envelope frequency
for I = 1:numel(envF)
    close all
    % current envelope frequency and cycle duration
    CurrEnvFreq = envF(I);
    cycl = cyclVect(I);
    disp('**************************')
    disp(['****** ' num2str(CurrEnvFreq) 'Hz ******'])
    disp('**************************')
    clear temp*
    % get EOD, global stimulus and dipole epoch for each frequency
    % EOD trace
    tempEOD = data.Ch8.values(round(envON(I)*SR):round((envON(I)+envDurations(I))*SR)-1);
    % get mean EODf just before stimulus
    disp('^^^^^^^^^^^^^ get mean EODf before stimulus ^^^^^^^^^^^^^^^^')
    tempEODbefore = data.Ch8.values(round((envON(I)-1.5)*SR):round(envON(I)*SR)-1);
    meanEODf = EODftimedependent(tempEODbefore, SR, SRglbl);
    [~, meanEODf] = chirpremove(meanEODf, CurrEnvFreq, SR, 0);
    meanEODf = nanmean(meanEODf(500:2500));
    clear tempEODbefore
    
    % global stimulus
    tempGLBL = data.Ch3.values(round(envON(I)*SRglbl):round((envON(I)+envDurations(I))*SRglbl)-1);
    % dipole
    tempDipole = data.Ch1.values(round(envON(I)*SRdipole):round((envON(I)+envDurations(I))*SRdipole)-1);
    tempDipole = tempDipole.*10; % convert into mV/cm (factor 10 if dipole is 1mm apart)
    
    % get stimulius contrast, meanEODf and amplitude of dipole
    disp('^^^^^^^^^^^^^ get contrast and dipole amplitude ^^^^^^^^^^^^^^^^')
    StimPer = [1 numel(tempDipole)];
    [StimContrast(I), ~, Dipamp(I), ~] = getContrast (tempDipole, 1/SRdipole, blEODamp, CurrEnvFreq, StimPer);
    clear StimPer
    
    % downdample dipole to 2k
    tempDipole = decimate(tempDipole,SRdipole/SRglbl);
    
    % plot EOD trace and threshold to get the EOD_0Xings
    disp('^^^^^^^^^^^^^ extract EODf ^^^^^^^^^^^^^^^^')
    tempfEOD = EODftimedependent(tempEOD, SR, SRglbl);

    % chrip removal
    disp('^^^^^^^^^^^^^ remove chirps ^^^^^^^^^^^^^^^^')
    tempEODnochirps = chirpremove(tempfEOD, CurrEnvFreq, SRglbl, 1);
    
    %----------------------------------------------------------------------
    scalingFactor(I)=1./(Dipamp(I)); % multiply gain by this
    disp('^^^^^^^^^^^^^ compute gain / phase / gain ^^^^^^^^^^^^^^^^')
    [gain(I), phase(I), offset(I)] = gainphaseoffset(tempDipole, tempEODnochirps, meanEODf, scalingFactor(I), CurrEnvFreq, SRglbl);
end
% save results
save([curDir, pLoc,'results' name '.mat'],'gain','phase','offset','envF','scalingFactor')
%% plot data
clear
close all
clc

pre = load('results_pre.mat');
post = load('results_pre.mat');
% post = load('results_post.mat');

k = .3;
gainLim(2) = max([pre.gain post.gain]);
gainLim(1) = min([pre.gain post.gain]);
gainLim(2) = gainLim(2)+(k*gainLim(2));
gainLim(1) = gainLim(1)-(k*gainLim(1));

offsetLim(2) = max([pre.offset post.offset]);
offsetLim(1) = min([pre.offset post.offset]);
offsetLim(2) = offsetLim(2)+(k*offsetLim(2));
offsetLim(1) = offsetLim(1)-(k*offsetLim(1));


figure;clf
tiledlayout(3,1)

nexttile
plot(pre.envF,pre.gain,'kd-','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','k')
hold on
plot(post.envF,post.gain,'rd-','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','r')
hold off
ylim(gainLim);xlim([0.04 1.1])
logx; logy
ylabel('gain (Hz/mV/cm)')
legend('pre','post')

nexttile
plot(pre.envF,pre.phase,'kd-','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','k')
hold on
plot(post.envF,post.phase,'rd-','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','r')
hold off
ylim([-180 0]);xlim([0.04 1.1])
logx
ylabel('phase (rad)')

nexttile
plot(pre.envF,pre.offset,'kd-','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','k')
hold on
plot(post.envF,post.offset,'rd-','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','r')
hold off
ylim(offsetLim);xlim([0.04 1.1])
logx
logx;
ylabel('offset (Hz)')
xlabel('envelope frequency (Hz)')