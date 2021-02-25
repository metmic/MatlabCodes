function [Contr, fEOD, amp, envStim] = getContrast (Dipole, Stim, dt, blEODamp, envF, varargin)

if size(varargin,2) == 1
    idxStim = varargin{1};                 % indices of stimulus periods to be analyzed
    idxBL = [round(1./dt), round(5./dt)];  % use second 1-5 of recording for baseline measurement
elseif size(varargin,2) == 2
    idxStim = varargin{1};                 % indices of stimulus periods to be analyzed
    idxBL = [round(1./dt), round(5./dt)];  % use second 1-5 of recording for baseline measurement
elseif size(varargin,2) >= 2
    error('getContrast.m: wrong number of input arguments')
end

%% get a rough estimate of the EOD frequency
% do FFT and find the peak
nfft = 5000;
[p, f] = pwelch(Dipole-mean(Dipole), bartlett(nfft), .5*nfft, nfft, 1/dt);
fEOD = f(p==max(p));

%% extract envelope signal of entire trace
HighCut = fEOD;
[b,a] = butter(3, HighCut./(.5/dt), 'low');                                 % LowPass filter at EOD frequency in order to get rid of HF artifacts
Dipole_filt = filtfilt(b,a, Dipole);                                        % filtered dipole signal
env = envelope(Dipole_filt, envF*5000, 'peaks');                                  % extract envelope signal 

%% visualization whole trace 
% figure; hold on;
% plot(Dipole_filt, 'k')
% plot(env, 'r')
% title (['low pass filtered at ', num2str(fEOD), ' Hz'])

%% calculate baseline EOD amplitude
% bin_bl = ones (length(env),1);                                            % generate a binary of non-stimulated times
% for i=1:size(idxStim,1)
%     bin_bl(idxStim(i,1):idxStim(i,2))=0;
% end
% tdead = (1/dt)*2;                                                         % exclude 1st 2 seconds from analysis
% bin_bl(1:tdead) = 0;
% 
% aEODbl = mean(env(logical(bin_bl)));                                      % baseline EOD amplitude

% aEODbl = mean(env(idxBL(1):idxBL(2)));                                    % baseline EOD amplitude
aEODbl = blEODamp;                                                          % baseline EOD amplitude

%% get std of EOD during stimulation 
for k = 1:size(idxStim,1)
    stdEODstim(1,k) = (std(env(idxStim(k,1):idxStim(k,2))));                % caluclate contrast for dipole 1
end

%% calculate contrasts for all stimulus episodes
Contr = (stdEODstim/aEODbl)*100;
%% calculate dipole amplitude for all stimulus episodes
figure;plot(Dipole_filt);
title('select peak and trough of envelope modulation')
[~, Y] = ginput(2);
close all

mult = abs(diff(Y)) / abs(diff([max(Stim) min(Stim)]));
offset = (sum(Y)/2);

figure;
yyaxis left; plot(1/10000:1/10000:numel(Dipole)/10000,Dipole_filt);ylim([0 5]);
yyaxis right; plot(1/2000:1/2000:numel(Stim)/2000,(Stim*mult)+offset,'LineWidth',1.5);ylim([0 5])
title('amplitude agrees most of the time?')
quest = input('match [1: yes; 0: no]?: ');
while quest == 0
    close all
    disp(['### manually chose new values; mult: ' num2str(mult) '; offset: ' num2str(offset) ' ###'])
    mult = input('enter number to multiply: ');
    offset = input('enter offset: ');
    figure;
    yyaxis left;plot(Dipole_filt);ylim([0 5]);
    yyaxis right;plot((Stim*mult)+offset,'LineWidth',1.5);ylim([0 5]);
    quest = input('match [1: yes; 0: no]?: ');
end
envStim = Stim*mult+offset;
amp = abs(max(envStim)-min(envStim));
end