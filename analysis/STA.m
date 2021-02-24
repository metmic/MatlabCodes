function [vSTA, tSTA, Type, SL] = STA (spT, vStim, dtStim, varargin)
%function STA_Volker (spT, vStim, dtStim, varargin)
%
% spike triggered average 
%
% ---------------------------------------------------------------------------------------------------------------------
% SYNTAX
% [vSTA, tSTA, Type] = STA (tSp, vStim, dt, vis, lead, lag, window, offset)
%
% ---------------------------------------------------------------------------------------------------------------------
% DESCRIPTION
% STA.m calculates a spike triggered average of the inputs. Furthermore it evaluates the type of neuorn (ON vs OFF) 
% based on the slope of the STA in an evaluation window 
% and noise) for two input binary spike trains.
%
% inputs:
% tSp        -  spiketimes relative to the stimulus waveform
% vStim      -  stimulus waveform
% dt         -  1/sampling rate of vStim
%
% optional inputs: 
% vis        -  swich on/off the visualization of the STA ('1' for on, '0' for off)
%               (default = on)
% lead       -  lead time in s of STA before 0 (default = 0.1 s)
% lag        -  lag time in s of STA after time 0 (default = 0.1 s)
% window     -  window width in s for STA evaluation (default = 0.01 s)
% offset     -  negative offset of evaluation window in s relative to time 0
%               implemented to account of synaptic delay (default = 0.007 s for ELL)
%
% outputs:
% vSTA       -  waveform of STA
% tSTA       -  time vector corresponding to vSTA
% Type       -  '1' if ON-type cell; '0' if OFF-type cell
% SL         -  sum of slope in evaluation window
%
% ---------------------------------------------------------------------------------------------------------------------
% author:  V. Hofmann
% last changes: 2019-01-15
% ---------------------------------------------------------------------------------------------------------------------
% tSp=spT{2}; %which unit
tSp = spT;
dt = dtStim;

%% check inputs
if isempty(varargin)
    vis = 1;
    lead = 0.1;
    lag = 0.1;
    window = 0.01;
    offset = 0.007;
elseif length(varargin) == 1
    vis = varargin{1};
    lead = 0.1;
    lag = 0.1;
    window = 0.01;
    offset = 0.007;
elseif length(varargin) == 2
    vis = varargin{1};
    lead = varargin{2};
    lag = 0.1;
    window = 0.01;
    offset = 0.007;
elseif length(varargin) == 3
    vis = varargin{1};
    lead = varargin{2};
    lag = varargin{3};
    window = 0.01;
    offset = 0.007;
elseif length(varargin) == 4
    vis = varargin{1};
    lead = varargin{2};
    lag = varargin{3};
    window = varargin{4};
    offset = 0.007;
elseif length(varargin) == 5
    vis = varargin{1};
    lead = varargin{2};
    lag = varargin{3};
    window = varargin{4};
    offset = varargin{5};
elseif length(varargin) >= 6
    error('wrong number of input arguments')
end


%% calculate indices
leadIDX = lead/dt;                                                % number of indices for lead time
lagIDX = lag/dt;                                                  % number of indices of lag time
spIDX = round(tSp/dt);                                            % indicies of spike times

%% eliminate spike with overlapping ends 
spIDX = spIDX(spIDX>leadIDX & spIDX< length(vStim)-lagIDX);       % exclude spikes that are too close to the start and end


%% gather waveforms
trg_vStim = NaN(length(spIDX), leadIDX+lagIDX+1);
for i=1:length(spIDX)
    trg_vStim (i,:) = vStim(spIDX(i)-lead/dt:spIDX(i)+lag/dt);
end

%% calculate spike triggered average
vSTA = mean(trg_vStim);

delay = 0.00;                                                    % add delay if need to shift the STA on way or the other
tSTA = -1*lead:dt:lag;
tSTA = tSTA+delay;


%% determine type
slope = [NaN, diff(vSTA)];
wind = [offset*-1-window/2 offset*-1+window/2];
windIDX = leadIDX+wind/dt;

windSlope = slope(windIDX(1):windIDX(2));

if sum(windSlope)>0
    Type = 1;
elseif sum(windSlope)<0
    Type = 0;
else
    Type = NaN;
end

SL = sum(windSlope);


%% visualize STA
if vis == 1
    h = figure; 
    subplot 211; hold on; box on;
%     plot([0 0], [min(vSTA)*1.1 max(vSTA)*1.1], '--k')
    plot([0 0], [min(vSTA) max(vSTA)], '--k')
    plot([-lead lag], [0 0], '--k')
    plot(tSTA,vSTA, 'k', 'LineWidth', 2)
    plot(tSTA(windIDX(1):windIDX(2)),vSTA(windIDX(1):windIDX(2)), 'r', 'LineWidth', 2)
%     set(gca, 'YLim', [min(vSTA)*1.1 max(vSTA)*1.1], 'XLim', [-lead+delay lag+delay])
    set(gca, 'YLim', [min(vSTA) max(vSTA)], 'XLim', [-lead+delay lag+delay])
    xlabel('time (s)')
    ylabel('amplitude')


    subplot 212; hold on; box on;
    patch([offset*-1-window/2 offset*-1+window/2  offset*-1+window/2  offset*-1-window/2], 1.1*[min(slope) min(slope) max(slope) max(slope)], '--k', 'FaceColor', 'k', 'FaceAlpha', 0.1)
    plot(tSTA,slope, 'k', 'LineWidth', 2)
    plot(tSTA(windIDX(1):windIDX(2)), slope(windIDX(1):windIDX(2)), 'r', 'LineWidth', 2)
%     set(gca, 'YLim', [min(slope)*1.1 max(slope)*1.1], 'XLim', [-lead+delay lag+delay])
    set(gca, 'YLim', [min(slope) max(slope)], 'XLim', [-lead+delay lag+delay])
    xlabel('time (s)')
    ylabel('slope of STA')

    text(0, max(slope)/3*2, ['\leftarrow sum slope = ', num2str(sum(windSlope)), '\newline     cell Type = ', num2str(Type)])
%    pause
%    close (h)
end
end