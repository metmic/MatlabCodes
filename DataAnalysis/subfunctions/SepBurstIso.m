function [spiketimes_bursts,spiketimes_isolated]=SepBurstIso(spiketimes_use,ttr)
% this code separates spikes in a burst and isolated spikes using a
% burst-threshold "ttr"
% code assumes that the spiketimes are in sec. 

isis=diff(spiketimes_use); % compute the ISI sequence of the spiketrain
isis(isis<0.001)=[];
threshold=ttr; % set the threshold to 10 msec, change as necessary

% we'll loop over the isis
spiketimes_bursts=[];
spiketimes_events=[];  % set up the appropriate indices
spiketimes_isolated=spiketimes_use; % initially all spikes are isolated
I=2; % loop variable for the isis
while I<=length(isis)-1 % this is a "while" loop 
    if isis(I-1)>threshold && isis(I)<=threshold % we are at the beginning of a burst
    spiketimes_events=[spiketimes_events spiketimes_use(I)]; % add the current spiketime to "spiketimes_events" as it belongs 
    % to a burst
    spiketimes_bursts=[spiketimes_bursts spiketimes_use(I) spiketimes_use(I+1)]; % also add it to the "indices_bursts"
    dummy=I; % dummy variable which we set to the current value of I 
    while isis(I)<=threshold && I<=length(isis)-1 % check if the next ISI is also less than threshold
    spiketimes_bursts=[spiketimes_bursts spiketimes_use(I) spiketimes_use(I+1)];
    I=I+1;
    end
    spiketimes_isolated=setdiff(spiketimes_isolated,spiketimes_use(dummy:I)); % remove the indices belonging to bursts from isolated
    if I>dummy
        I=I-1; % go back one step
    end
    end
I=I+1; % this is where we increment the variable
end

spiketimes_bursts=unique(spiketimes_bursts); % eliminate any duplicates
spiketimes_bursts=spiketimes_bursts(:);
spiketimes_isolated=unique(spiketimes_isolated); % eliminate any duplicates
spiketimes_events=unique(spiketimes_events); % eliminate any duplicates