function [gain, phase, offset] = gainphaseoffset(stimulus, response, meanResponse, scalingFactor, Freq, SR)

% gain/phase/offset measure

    figure(201);plot(response)
    title('select most stationary part of response')
    ylim([meanResponse-10 meanResponse+10])
    % choose whether to take entire trace or subsection;
    % this is needed if the response is "inconsistent" or if tracking is only
    % good for one portion
    choice = menu('Take entire trace?','YES','NO');
    switch choice
        case 1
            TempStimulus = stimulus;
            TempResponse = response;
        case 2
            % define the part you want to analyze using a crosshair
            XPart = ginput(2);
            XPart = round(XPart);
            if XPart(1)<0
                XPart(1)=1;
            end
            if XPart(2) > numel(response)
                XPart(2) = numel(response);
            end
            TempStimulus = stimulus(XPart(1):XPart(2));
            TempResponse = response(XPart(1):XPart(2));
    end
    close 201
    
    % get 0Xings of the envelope stimulus
    minpkdist = floor((numel(TempStimulus)/(numel(TempStimulus)/SR*Freq))- ...
        (numel(TempStimulus)/(numel(TempStimulus)/SR*Freq)*0.1));
    if (minpkdist<numel(TempStimulus))
        [~,peaks] = findpeaks(TempStimulus,'minpeakheight',0, ...
            'minpeakdistance',minpkdist);
    else
        minpkdist = numel(TempStimulus)-1;
        [~,peaks] = findpeaks(TempStimulus,'minpeakheight',0, ...
            'minpeakdistance',minpkdist);
    end
    peaks = peaks(peaks<numel(TempStimulus));
    if numel(peaks)>1
        if length(peaks)==2
            maxdiff = max(diff(peaks(1:end)));
        else
            maxdiff = ceil((numel(TempStimulus)/(numel(TempStimulus)/SR*Freq)));
        end
        peaks(1) = peaks(2)-maxdiff;
        if peaks(1)<0
            peaks(1) = 1;
        end
        indices3 = find(diff(peaks)<maxdiff/2);
        indices2 = setdiff(1:1:numel(peaks),indices3);
        peaks = peaks(indices2);
    end
    newpeaks = zeros(numel(peaks),1);
    for J = 1:numel(newpeaks)
        newpeaks(J) = peaks(1)+(J-1)*maxdiff;
    end
    peaks = newpeaks;
    clear newpeaks
    %---------------------------------------------------------------------
    psth = zeros(maxdiff,1);
    
    if Freq>=1
        MM = 10;
    else
        MM = 1;
    end
    if numel(peaks)>1
        for J = 1:numel(peaks)-MM
            temp = TempResponse(peaks(J):peaks(J)+maxdiff-1);
            temp = temp';
            psth = psth+temp;
        end
        psth = psth./((numel(peaks)-MM));
    end
    chq1 = diff(psth);
    maxchq1 = max(chq1);
    minchq1 = min(chq1);
    [~,inds_chq1_pos] = min(abs(chq1-maxchq1));
    [~,inds_chq1_neg] = min(abs(chq1-minchq1));
    if inds_chq1_pos>inds_chq1_neg
        psth(inds_chq1_pos+1:end,:) = (psth(inds_chq1_pos+1:end))- ...
            (abs(chq1(inds_chq1_pos)));
    elseif inds_chq1_pos<inds_chq1_neg
        psth(inds_chq1_neg+1:end) = ...
            psth(inds_chq1_neg+1:end)+(abs(chq1(inds_chq1_neg)));
    else
    end
    clear chq1 inds_chq1_pos chq_inds_neg ff maxchq1 minchq1
    %---------------------------------------------------------------------
    % Nonlinear fitting
    % the function funENVELOPE generates a sine/cosine wave to fit your PSTH
    edges = (1:1:numel(psth))./SR;
    edges = edges-edges(1);
    funENVELOPE = @(beta0,x) beta0(1)+beta0(2)*cos(2*pi.*Freq.*x+beta0(3));
    [beta,r,~,Sigma] = nlinfit((1:1:numel(psth))./SR,psth',funENVELOPE,...
        [nanmean(psth) (max(psth)-min(psth)/2) 0]);
    yfit = nlpredci(funENVELOPE,(1:1:numel(psth))./SR,beta,r,'Covar',Sigma);
    SSE = sum((psth-yfit).^2);
    SST = sum((psth-mean(psth)).^2);
    R = sqrt(1-SSE/SST);
    figure(999);plot((1:1:length(psth))./SR,psth);
    beta(3) = mod(beta(3),2*pi);
    if beta(3)<0
        beta(2) = -beta(2);
        beta(3) = beta(3)+pi;
    end
    if beta(3)<-pi
        beta(3) = beta(3)+2*pi;
    end
    if beta(3)>pi
        beta(3) = beta(3)-2*pi;
    end
    if beta(3)>0
        beta(3) = beta(3)-2*pi;
    end
    % by using the command "hold on" the last generated figure will be updated
    % with a second plot
    hold on;plot(edges,funENVELOPE(beta,edges),'r');
    pause(1)
    close 999
    % now calculate gain, phase and offset of the response to envelope stimuli
    gain = abs((beta(2))-(beta(1)))*scalingFactor;
    offset = abs(nanmean(TempResponse)-meanResponse);
    time_estim = (1:1:numel(psth))*SR';
    if time_estim(end)>=1/Freq-(1/Freq/10)
        pha = (beta(3)/2/pi*360);
        if pha<-300
            pha = pha+300;
        end
        phase = pha;                              % phase
        clear pha time_estim
    else
        phase = NaN;
        clear time_estim
    end