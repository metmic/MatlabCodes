function medianTrace = applyCARtoDatORI(filename, nChansTotal, goodChans, EODf, outputDir, SRo, desSR, sMap, Pos, FSi, FSo, flag)
%                        applyCARtoDatNEWori([path file], Tchan, Bchan, EODf, pwd, SRo, desSR, sMap, Pos, FSi, FSo, flag);
%                           applyCARtoDatNEW(filename, nChansTotalN, goodChans, EODf, outputDir, SRo, desSR, sMap, Pos, FSi, FSo, flag)
% Subtracts median of each channel, then subtracts median of each time
% point.
%
% filename should include the extension
% outputDir is optional, by default will write to the directory of the input file
%
% should make chunk size as big as possible so that the medians of the
% channels differ little from chunk to chunk.

%%
copy_sglxMeta(filename);
if isempty(desSR)
    % keep sGLX sampling
    SRXk = SRo;
    chunkSize = round(SRXk)*60;
    % TR is for channel distance calculation used in the loc filtering
    TR = 0;
else
    % re-sample
    SRXk = desSR;
    chunkSize = 1000000;
    TR = 1;
    [P,Q] = rat(SRXk/SRo);
end

fid = []; fidOut = [];
d = dir(filename);
nChansTotalO = 385;
nSampsTotal = ceil(d.bytes/nChansTotal/2);%volker
nChunksTotal = ceil(nSampsTotal/chunkSize);
allChans = 1:nChansTotalO;

if flag == 1
    disp('high pass filter')
    [B, A] = butter(2,(300)/(SRo*0.5),'high');
else
    disp('band pass filter')
    fcutlow = 300;   %low cut frequency in Hz
    fcuthigh = EODf*3;   %high cut frequency in Hz
    [B, A] = butter(2,[fcutlow,fcuthigh]./(SRo/2),'bandpass');
end

try
    [pathstr, name, ext] = fileparts(filename);
    fid = fopen(filename, 'r');
    idx = strfind(name, '.imec');
%     if nargin < 7
%         outputFilenameXk  = [pathstr filesep name(1:idx-1) '_FiXk' name(idx:end) ext];
%         mdTraceFilename = [pathstr filesep name(1:idx-1) '_Fi' name(idx:end) '_medianTrace.mat'];
%     else
    outputFilenameXk  = [outputDir filesep name(1:idx-1) '_FiXk' name(idx:end) ext];
    mdTraceFilename = [outputDir filesep name(1:idx-1) '_Fi' name(idx:end) '_medianTrace.mat'];
%     end
    fidOutXk = fopen(outputFilenameXk, 'w');
    
    % theseInds = 0;
    chunkInd = 1;
    medianTrace = zeros(1, ceil(nSampsTotal));
    while 1
        fprintf(1, 'chunk %d/%d\n', chunkInd, nChunksTotal);
        dat = fread(fid, [nChansTotal chunkSize], '*int16');%old
        if ~isempty(dat)
            % copy good channels
            dat = double(dat);
            datF = NaN(size(dat));
            %% filter each raw trace using high or band oass filter
            for I = 1:nChansTotal
                if I<=goodChans
                    temp = dat(I,:);
                    tempF = filtfilt(B,A,[temp temp temp]);
                    tempT = tempF(numel(temp)+1:2*numel(temp));
                    datF(I,:) = tempT;
                    clear temp tempF tempT
                else
                    datF(I,:) = dat(I,:);
                end
            end
            %% loc filter
            dat_filt = datF; % create data matrix
            if (goodChans) == 1
                if SRXk == SRo
                    datFXk(1,:) = dat_filt;
                else
                    % resample
                    tempTXk = resample(dat_filt,P,Q);
                    datFXk(1,:) = tempTXk;
                    clear tempTXk
                end
            else
                for i = 1:goodChans % loop through channels of interest
                    dist = [];
                    for ii = 1:size(sMap,1)-TR % calculate distance between current channel and all other channels
                        dist(ii) = sqrt((Pos(i,1)-Pos(ii,1)).^2  + (Pos(i,2)-Pos(ii,2)).^2); % calculate euclidian distance of electrodes
                    end
                    % subtract spatial average
                    idx = (dist)>=FSi & (dist)<=FSo; % select channels within perimeter
                    dat_filt(i,:) = bsxfun(@minus, datF(i,:), nanmedian(datF(idx,:),1) );     % subtract median of each channel
                    if SRXk == SRo
                        datFXk(i,:) = dat_filt(i,:);
                    else
                        %resample
                        tempTXk = resample(dat_filt(i,:),P,Q);
                        datFXk(i,:) = tempTXk;
                        clear tempTXk
                    end
                end
                % fix for KS2.5
                datFXk(goodChans+1:nChansTotalO,:) = 0;
            end
            tm = nanmedian(datF,1);
            datFXk = int16(datFXk);
            %% write data
            fwrite(fidOutXk, datFXk, 'int16');  
            medianTrace((chunkInd-1)*chunkSize+1:(chunkInd-1)*chunkSize+numel(tm)) = tm;
%             if chunkInd == nChunksTotal
%                 medianTrace(:,chunkInd*chunkSize+size(datFXk,2)+1:end) = [];
%             end
            clear datTemp datFXk dat_filt datF dat tm
        else
            break
        end
        chunkInd = chunkInd+1;
    end
    %   save(mdTraceFilename, 'medianTrace', '-v7.3');
    fclose(fid);
    fclose(fidOutXk);
    
catch me
    if ~isempty(fid)
        fclose(fid);
    end
    
    if ~isempty(fidOut)
        fclose(fidOut);
    end
    
    if ~isempty(fidOutXk)
        fclose(fidOutXk);
    end
    rethrow(me) 
end

end