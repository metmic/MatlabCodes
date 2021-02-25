function [output_matrix,stimulus] = EnvTuningNPIXDUO(stimulus, binary, Envf, SR, numCylc, scalingFactor, plotFig)
% envelope tuning

%----- INPUT VARIABLES ----------------------------------------------------
% stimulus = scaled envelope stimulus
% binary = binary matrix neurons
% Envf = vector containing the envelope frequencies in Hz
% SR = sampling frequency binary
% numCylc = number of full cycles per envelope frequency
% scalingFactor = number to scale stimulusO up to mV/cm
% plotFig = set 1 if tuning should be plotted

%----- OUTPUT VARIABLES ---------------------------------------------------
% output_matrix = matrix containg the results
%       col 1: upper frequency carrier
%       col 2: firing rate
%       col 3: burst fraction
%       col 4: envelope frequency
%       col 5 - 7: gain (all spikes - bursts - isolated spikes)
%       col 8 - 10: phase (all spikes - bursts - isolated spikes)
%       col 11 - 13: vector strength (all spikes - bursts - isolated spikes)
%       col 14: z-score VS 

% Maso 2018
%==========================================================================

%----- burst threshold ----------------------------------------------------
ttr = 0.01;
%----- get data -----------------------------------------------------------
envDur = (numCylc.*(1./Envf))*SR;
time = 1/SR:1/SR:size(stimulus,1)/SR;
%--------------------------------------------------------------------------
%----- initialize output matrix -------------------------------------------
output_matrix = 1:14;
%----- analyze ------------------------------------------------------------
% disp('analyze responses')
%----- convert into spiketimes---------------------
spiketimes1 = (find(binary == 1))./SR;
%-----  get the zero=crossings of the envelope stimulus ---------------
mpd = round((50/(1/Envf)))-1;
[peaks,OXings_env] = findpeaks(stimulus,'minpeakdistance',mpd);
OXings_env = OXings_env(peaks>mean(stimulus));
%     figure;plot(stimulus);hold on; plot(OXings_env,mean(stimulus)*ones(length(OXings_env)),'*r');hold off
%     pause(1)
OXings_env=OXings_env(OXings_env<length(stimulus));
OXings_env=OXings_env./SR;    % OXings in sec
%----- separate bursts and iso spikes ---------------------------------
spiketimes = spiketimes1;             % takes all spiketimes since we have selected the start and end of stim
[spiketimes_bursts,spiketimes_isolated] = SepBurstIso(spiketimes,ttr);
%----- get the measures gain and phase --------------------------------
% all spikes
%     [vector_strength,z_stat,gain,gainraw,phase]
[VSA,z_stat,gainA,gainrawA,phaseA] = GainPhaseMI(spiketimes,OXings_env,Envf,scalingFactor,plotFig);
%     ind_int=find(z_stat>=3.5); % only plot vector strength estimates that are statistically significant (corresponds to p<0.05);
% % % % %     % bursts
% % % % %     if isempty(spiketimes_bursts)==1
% % % % %         VSB=NaN;
% % % % %         gainB=NaN;
% % % % %         gainrawB=NaN;
% % % % %         phaseB=NaN;
% % % % %     else
% % % % %         [VSB,~,gainB,gainrawB,phaseB]=GainPhaseMI(spiketimes_bursts,OXings_env,env_freq,scalingFactor(ii),ii);
% % % % %     end
% % % % %     % isolated spikes
% % % % %     [VSI,~,gainI,gainrawI,phaseI]=GainPhaseMI(spiketimes_isolated,OXings_env,env_freq,scalingFactor(ii),ii);
%----- Assign Values into output matrix -------------------------------
K = 1;
output(K,1)=15;
output(K,2)=NaN;
output(K,3)=length(spiketimes_bursts)/length(spiketimes);
output(K,4)=Envf;
output(K,5)=gainA; % all spikes
output(K,6)=NaN;%gainB; % bursts
output(K,7)=NaN;%gainI; % isolated spikes

output(K,8)=phaseA;
output(K,9)=NaN;%phaseB;
output(K,10)=NaN;%phaseI;

output(K,11)=VSA; % vector strength all spikes
output(K,12)=NaN;%VSB; % vector strength bursts
output(K,13)=NaN;%VSI; % vector strength isolated spikes
output(K,14)=z_stat;

K=K+1;
output_matrix=[output_matrix; output];
clear binaryTemp spiketimes1 sTemp
close all
