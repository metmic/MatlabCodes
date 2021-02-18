function medianTrace = applyCAR(filename, hpf, outputDir)
%
% performs common average reference (CAR) subtraction
%
% ---------------------------------------------------------------------------------------------------------------------
% SYNTAX
% medianTrace = applyCAR(filename, hpf, outputDir)
%
% ---------------------------------------------------------------------------------------------------------------------
% DESCRIPTION
% applyCAR.m performs a CAR subtraction on input data. An optional HPF can be added to further reduce
% contamination of spike data by the EOD artifact. 
%
% inputs:
% filename    -  filename of *.bin file to be processed. Filename should include path information
% channels    -  vector with channel numbers to be included in processing (indexed as 1-based indeces)
% hpf         -  swich command for hpf. if hpf == 1 a 3rd order butterworth hpf will be applied to the data
% outputDir   -  optional output directory for stored files. If empty files will be stored in current directory
%
% outputs:
% medianTrace -  medianTrace used for CAR subtraction
%
% comments:
% a few variations of this codes are still in place but deactivated: 
% 1. L. 110-117 a scaled CAR (median trace is scaled according to variance in target channel before subtraction)
% 2. L. 49 & 124-158 a spatial (local) subtractive filter
% 3. L. 159-168 original CAR (without selection of target channels)
% 
%
% ---------------------------------------------------------------------------------------------------------------------
% author:  V. Hofmann (based on code by N. Steinmetz)
% last changes: 2019-08-23
% ---------------------------------------------------------------------------------------------------------------------


%% read shank data
copy_sglxMeta(filename);
[meta] = read_sglxMeta(filename);
[sMap, Pos] = readShankMap (meta); % required for local spatial filters
sr = round(str2double(meta.imSampRate));


%%
chunkSize = round(sr)*60;

d = dir(filename);
nChansTotal = str2double(meta.nSavedChans); %nChansTotal = 385; 
nSampsTotal = d.bytes/nChansTotal/2;
nChunksTotal = ceil(nSampsTotal/chunkSize);

fid = []; fidOut = [];
try
    
    [pathstr, name, ext] = fileparts(filename);
    fid = fopen(filename, 'r');
    idx = strfind(name, '.imec');
    if nargin < 3
        outputFilename  = [pathstr filesep name(1:idx-1) '_CAR' name(idx:end) ext];
        mdTraceFilename = [pathstr filesep name(1:idx-1) '_CAR' name(idx:end) '_medianTrace.mat'];
    else
        outputFilename  = [outputDir filesep name(1:idx-1) '_CAR' name(idx:end) ext];
        mdTraceFilename = [outputDir filesep name(1:idx-1) '_CAR' name(idx:end) '_medianTrace.mat'];
    end
    fidOut = fopen(outputFilename, 'w');
    
    % theseInds = 0;
    chunkInd = 1;
    medianTrace = zeros(1, nSampsTotal);
    while 1
        fprintf(1, 'chunk %d/%d\n', chunkInd, nChunksTotal);
        dat = fread(fid, [nChansTotal chunkSize], '*int16');
        
        if ~isempty(dat)
            
            %% subtract median of each channel
            dat = bsxfun(@minus, dat, median(dat,2));
%[pxx, freq] = pwelch(double(dat(77,:)),bartlett(2048),1024,2048, sr);
%figure; loglog(freq, pxx)
            
            if hpf == 1
                %% LP/HP filter
                locutoff = 100;
                hicutoff = 2500;
                [B,A] = butter (3, [locutoff/(0.5*sr), hicutoff/(0.5*sr)], 'bandpass');
                
                    %if size(dat,2)>2*sr+size(dat,2)
                    tmp = NaN(size(dat,1), size(dat,2)+2*sr);
                    for i = 1:nChansTotal-1
                        tmp(i,:) = filtfilt(B,A, double([dat(i,1:sr),dat(i,:),dat(i,1:sr)]'));
                    end
                    tmp(nChansTotal,:)=[dat(nChansTotal,1:sr), dat(nChansTotal,:), dat(nChansTotal,1:sr)];
                    clear dat
                    dat = int16(tmp(:, sr+1:size(tmp,2)-sr));
            end
            
            
            %% CAR across selected channel range only
            %             tm = median(dat(channels,:),1);                     % CAR across active channels
            %             
            %             %subtract median (unscaled)
            %             datCAR = NaN(size(dat));
            %             datCAR(channels,:) = bsxfun(@minus, dat(channels,:), tm);
            %             
            %             % subtract median (scaled)
            %             %             tmAll = NaN(size(dat));
            %             %             for i = 1:size(channels,2)
            %             %                 scaleVar = (std(double(dat(channels(i),:))))./std(double(tm));  % get relative scaling of var in target an average channel
            %             %                 tmAll (channels(i),:) = tm.*scaleVar;
            %             %             end
            %             %             datCAR = dat-int16(tmAll);
            %             %             tm = nanmedian(tmAll);
            % 
            %             % write data
            %             fwrite(fidOut, datCAR, 'int16');
            %             medianTrace((chunkInd-1)*chunkSize+1:(chunkInd-1)*chunkSize+numel(tm)) = tm;

            
            %% local spatial averaging
                        FSi = 20; % inner diameter of filter smallest electrode distance is 25.6 mum
                        FSo = 60; % outer diameter of filter
                                               
                        dat_filt = NaN (size(dat)); % create data matrix
            
                        for i = 1:nChansTotal-1 % i = channelRange % loop through channels of interst
            
                            for ii = 1:size(sMap,1) % calculate distance between current channel and all other channels
                                dist(ii) = sqrt((Pos(i,1)-Pos(ii,1)).^2  + (Pos(i,2)-Pos(ii,2)).^2); % calculate euclidian distance of electrodes
                            end
            
                            % calculate "ranked distance" (i.e. closest set of electrodes = 1, second closest set of electrodes = 2, ...)
                            %             [sdist,~] = sort(unique(dist));
                            %             d_rank = NaN(size(dist));
                            %             for iii = 1:length(sdist)
                            %                 d_rank(dist == sdist(iii))= find(sdist==sdist(iii))-1;
                            %             end
            
            
                            % subtract spatial average
                            idx = dist>=FSi & dist<=FSo; % select channels within perimeter
                            % idx = d_rank>=FSi & d_rank<=FSo;
                            dat_filt(i,:) = bsxfun(@minus, dat(i,:), median(dat(idx,:)));     % subtract median of each channel
                        end
            
            
                        tm = nanmedian(dat_filt,1);
                        dat_filt(isnan(dat_filt))=0;
            
                        % write data
                        fwrite(fidOut, dat_filt, 'int16');
                        medianTrace((chunkInd-1)*chunkSize+1:(chunkInd-1)*chunkSize+numel(tm)) = tm;
            
            
            %% original CAR filter
            %       %         theseInds = theseInds(end):theseInds(end)+chunkSize-1;
            %       dat = bsxfun(@minus, dat, median(dat,2));     % subtract median of each channel
            %       tm = median(dat,1);
            %       dat = bsxfun(@minus, dat, tm);                % subtract median of each time point
            %
            %       % write data
            %       fwrite(fidOut, dat, 'int16');
            %       medianTrace((chunkInd-1)*chunkSize+1:(chunkInd-1)*chunkSize+numel(tm)) = tm;
            
        else
            break
        end
        
        chunkInd = chunkInd+1;
    end
    
    save(mdTraceFilename, 'medianTrace', '-v7.3');
    fclose(fid);
    fclose(fidOut);
    
catch me
    
    if ~isempty(fid)
        fclose(fid);
    end
    
    if ~isempty(fidOut)
        fclose(fidOut);
    end
    
    rethrow(me)
    
end


end