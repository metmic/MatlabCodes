function [EODf, EOD_chirp_filt] = chirpremove(fEOD, Freq, SR, flag)
% code to remove chirps from EOD trace and detrend

% filter EODtrace
[BB,AA] = butter(2,(100*Freq)/(SR*.5));
EOD_chirp_filt = filtfilt(BB,AA,[fEOD; fEOD; fEOD]);
EOD_chirp_filt = EOD_chirp_filt(numel(fEOD)+1:2*numel(fEOD));

% % % set threshhold
figure(001);plot(EOD_chirp_filt(200:end-200));title('threshold chirps; if no chirps then select below EOD trace ouside window')
axis tight
ylim([min(EOD_chirp_filt(200:end-200))-10 max(EOD_chirp_filt(200:end-200))+10])
[~,thresh] = ginput(1);close 001
% use fixed threshhold
% thresh = 0.1;

EOD_chirp_filtO = EOD_chirp_filt;
% replace chirps with extrapolated data between each start and stop point
if thresh < min(EOD_chirp_filt(200:end-200))
else
    [~,loc]=findpeaks(EOD_chirp_filt,'MinPeakHeight',thresh);
    % chirp duration
    durSC = 250;
    loc(loc<durSC) = [];
    loc(loc>numel(EOD_chirp_filt)-durSC) = [];
    for II = 1:length(loc)
        if loc(II) < durSC
            % get EOD frequency at beginning and end of chirp
            start = 1;
            endc = EOD_chirp_filt(round(1+durSC));
            % duration of chirp event
            dur = length((EOD_chirp_filt(1:1+durSC)));
            % create a dummy vector containing the mean EOD frequency for
            % the duration of the chirp event
            extrapoldata = mean([start endc])*ones(dur,1);
            % replace chirp with dummy vector
            EOD_chirp_filt(((1:1+durSC))) = extrapoldata;
        else
            % get EOD frequency at beginning and end of chirp
            start = EOD_chirp_filt(round(loc(II)-durSC));
            endc = EOD_chirp_filt(round(loc(II)+durSC));
            % duration of chirp event
            dur = length((EOD_chirp_filt(loc(II)-durSC:loc(II)+durSC)));
            % create a dummy vector containing the mean EOD frequency for
            % the duration of the chirp event
            extrapoldata = mean([start endc])*ones(dur,1);
            % replace chirp with dummy vector
            EOD_chirp_filt(((loc(II)-durSC:loc(II)+durSC))) = extrapoldata;
        end
        clear extrapoldata start endc dur
    end
    clear loc
end
% % figure(666);
% % plot(EOD_chirp_filtO);hold on;plot(EOD_chirp_filt,'LineWidth',1.3);legend('w/ chirps','w/o chirps')
% % % title('press key')
% % pause(1)
% % close 666
EODf = [];
if flag == 1
    % filter fEOD to remove non-stationarity
    [B1,A1] = butter(1,0.01/1000);
    EODfnochirpsDetrended = filtfilt(B1,A1,[EOD_chirp_filt; EOD_chirp_filt; EOD_chirp_filt]);
    EODfnochirpsDetrended = EODfnochirpsDetrended(numel(EOD_chirp_filt)+1:2*numel(EOD_chirp_filt));
    EODf = EOD_chirp_filt-EODfnochirpsDetrended+(nanmean(EOD_chirp_filt(2000:2000+(1/Freq*12*SR))));
end


