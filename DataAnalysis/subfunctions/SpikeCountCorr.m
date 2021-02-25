function [r_RAW, r_SIG, r_NOI, T, num_RAW, den_RAW, num_SIG, den_SIG, num_NOI, den_NOI] = SpikeCountCorr (bin1, bin2, varargin)
%
% spike count correlation of input spike trains
%
% ---------------------------------------------------------------------------------------------------------------------
% SYNTAX
% [C_RAW, C_SIG, C_NOI, num_RAW, T, den_RAW, num_SIG, den_SIG, num_NOI, den_NOI] = SpikeCountCorr (bin1, bin2, T, N, dt)
%
% ---------------------------------------------------------------------------------------------------------------------
% DESCRIPTION
% SpikeCountCorr.m is a calculates spike count correlations (raw, signal
% and noise) for two input binary spike trains.
%
% inputs:
% bin1       -  input binary of neuron 1; NxM matrix with N stimulus repetitions and M representing time according to resolution defined by dt
% bin2       -  input binary of neuron 2
% T          -  timescales for which pUnit activity should be calculated (default = logspace(log10(.001), log10(1.5), 15))
% dt         -  time increment of input-binaries (default = 1/2000)
% N          -  number of shuffles for the calculation of signal correlations (default = 20)
%
% outputs:
% C_*        -  correlation coefficient (RAW, SIG, NOI). Signal
%               correlations are based on N shuffles
% num_*      -  numerator of correlation coefficient
% den_*      -  denominator of correlation coefficient
%
% ---------------------------------------------------------------------------------------------------------------------
% author:  V. Hofmann
% last changes: 2018-12-10
% ---------------------------------------------------------------------------------------------------------------------


%% check inputs
if isempty(varargin)
    T = logspace(log10(.001), log10(1.5), 15);
    dt = 1/2000;
    N = 20;
elseif length(varargin) == 1
    T = varargin{1};
    N = 20;
    dt = 1/2000;
elseif length(varargin) == 2
    T = varargin{1};
    N =  varargin{2};
    dt = 1/2000;
elseif length(varargin) == 3
    T = varargin{1};
    N =  varargin{2};
    dt = varargin{3};
elseif length(varargin) >=3
    error('SpikeCountCorr.m : wrong number of input arguments')
end


%% pre-allocate variables
r_RAW = NaN(1,length(T)); r_SIG = NaN(1,length(T)); r_NOI = NaN(1,length(T));
num_RAW = NaN(1,length(T)); den_RAW = NaN(1,length(T)); num_SIG = NaN(1,length(T));
den_SIG = NaN(1,length(T)); num_NOI = NaN(1,length(T)); den_NOI = NaN(1,length(T));


for kk = 1:length(T)                                        % loop through time windows
    %% GENERATE SPIKE COUNTS
    iT = round(1:round(T(kk)/dt):length(bin1)+1);           % indices of spikecount window in N1/N2
    
    cSN1 = NaN(size(bin1,1) ,length(iT)-1);                 % preallocate cSN1
    cSN2 = cSN1;                                            % preallocate cSN2
    
    for ii = 1:size(bin1,1)                                 % loop through binaries
        for i = 1:length(iT)-1                              % loop through indicies
            cSN1(ii,i) = nansum(bin1(ii,iT(i):iT(i+1)-1));  % spike counts for bin1
            cSN2(ii,i) = nansum(bin2(ii,iT(i):iT(i+1)-1));  % spike counts for bin1
        end
    end
    
    
    %% RAW CORRELATIONS
    covM = cov(cSN1,cSN2);                                  % cov and var of spike counts
    r_RAW(kk)=covM(2,1)/sqrt(covM(1,1)*covM(2,2));          % pearson correlation coefficient
    num_RAW(kk)=covM(2,1);                                  % numerator
    den_RAW(kk)=sqrt(covM(1,1)*covM(2,2));                  % denominator
    clear covM
    
    
    %% SIGNAL CORRELAITONS
    num = NaN(1,N);
    denom = NaN(1,N);
    
    for f=1:N                                               % shuffle N times
        order1 = randperm(size(cSN2,1), size(cSN2,1));      % generate a new order for episodes
        order2 = randperm(size(cSN2,1), size(cSN2,1));
        cSN1_reshuf = cSN1(order1, :);                      % assign order to bianry
        cSN2_reshuf = cSN2(order2, :);
        
        covM = cov(cSN1_reshuf,cSN2_reshuf);                % cov and var of shuffled spike counts
        num(f) = covM(2,1);                                 % numerator
        denom(f) = sqrt(covM(1,1)*covM(2,2));               % denominator
    end
    
    r_SIG(kk) = mean(num)/mean(denom);                      % pearson correlation coefficient for shuffled spike counts form average cov and var
    num_SIG(kk) = mean(num);                                % average numerator
    den_SIG(kk) = mean(denom);                              % average denominator
    clear covM
    
    
    %% NOISE CORRELATIONS
    mcSN1 = repmat(mean(cSN1),size(cSN1,1),1);              % calculate a mean across all episodes, and reproduce it by the number of episode for subtraction
    mcSN2 = repmat(mean(cSN2),size(cSN2,1),1);
    
    dcSN1 = cSN1 - mcSN1;                                   % calculate residuals of spike counts
    dcSN2 = cSN2 - mcSN2;
    
    covM = cov(dcSN1, dcSN2);                               % cov and var of spike count residuals
    r_NOI(kk)=covM(2,1)/sqrt(covM(1,1)*covM(2,2));          % pearson correlation coefficient for spike count residuals
    num_NOI(kk) = covM(2,1);                                % numerator
    den_NOI(kk) = sqrt(covM(1,1)*covM(2,2));                % denominator
    clear covM
    
    
    %%
    clear iT cSN1 cSN2 mcSN1 mcSN1 dcSN1 dcSN2
end

