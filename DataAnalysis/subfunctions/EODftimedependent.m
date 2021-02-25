function EODf = EODftimedependent(Response, SR, SRnew)

% figure;plot(Response(1:1000))
% [~,EODthresh] = ginput(1);
EODthresh = 0.1;
close all
% extract the zero-crossings of the EOD trace
EODOXings = find([Response' -100]>=EODthresh & [-100 Response']<EODthresh)/(SR);
ind = find(diff(EODOXings)<1/1500);           % excludes frequencies higher than 1500Hz
EODOXings = EODOXings(setdiff(1:1:length(EODOXings),ind));
% EOD frequency series
fEOD = 1./diff(EODOXings); %
oldaxis = EODOXings(1:end-1);
newaxis = 1/SR:1/SR:numel(Response)/SR;
EODf = interp1(oldaxis,fEOD,newaxis,'spline');   % interpolate to 50k
sinds = find(isnan(EODf) == 1);
sinds1 = find(isnan(EODf) == 0);
EODf(sinds) = mean(EODf(sinds1));
EODf = EODf(:);
% downsample to 2k
EODf = decimate(EODf,SR/SRnew);