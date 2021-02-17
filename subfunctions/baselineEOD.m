function aEODbl = baselineEOD(Dipole, dt)


%% get a rough estimate of the EOD frequency
% do FFT and find the peak
nfft = 5000;
[p, f] = pwelch(Dipole-mean(Dipole), bartlett(nfft), .5*nfft, nfft, 1/dt);
fEOD = f(p==max(p));
%% extract envelope signal of entire trace
HighCut = fEOD;
[b,a] = butter(3, HighCut./(.5/dt), 'low');                  % LowPass filter at EOD frequency in order to get rid of HF artifacts
Dipole_filt = filtfilt(b,a, Dipole);                         % filtered dipole signal
env = envelope(Dipole_filt, 250, 'peaks');                   % extract envelope signal 
aEODbl = nanmean(env);                        % baseline EOD amplitude