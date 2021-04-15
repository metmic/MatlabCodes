function [fr_box, fr_gauss, fr_filt] = ConvolveBin (bin, dt, varargin)
%
% ConvolveBin calculates firing rate of binary
%
% ---------------------------------------------------------------------------------------------------------------------
% SYNTAX
% [fr_box, fr_filt] = ConvolveBin (bin, dt, w, sigma, cutoff, order)
%
% ---------------------------------------------------------------------------------------------------------------------
% DESCRIPTION
% ConvolveBin calculates firing rate of binary by 
% 1. boxcar convolution 
% 2. gaussian convolution
% 3. low pass filtering the binary
%
% inputs:
% bin      -    binary with spike data
% dt       -    1/SR of bin
% w        -    [optional] width of boxcar in s (default = 100 ms)
% sigma    -    [optional] sigma of gaussian kernel (default = 100 ms)
% cutoff   -    [opitonal] cutoff for low frequency filter (default = 10 Hz)
% order    -    [optional] order of the butterworth filter
%
% outputs:
% FR_box  -     box car convolved firing rate
% FR_gauss -    gaussian convolved firing rate
% FR_filt  -    filtered firing rate
%
% ---------------------------------------------------------------------------------------------------------------------
% author:  V. Hofmann
% last changes: 2020-02-24
% ---------------------------------------------------------------------------------------------------------------------


%% check inputs
if isempty (varargin)
    w=.1;
    sigma=0.1;
    cutoff = 10;
    order = 2;
elseif length(varargin) == 1
    w = varargin{1};
    sigma = 0.1;
    cutoff = 10;
    order = 2;
elseif length(varargin) == 2
    w = varargin{1};
    sigma = varargin{2};
    cutoff = 10;
    order = 2;
elseif length(varargin) == 3
    w = varargin{1};
    sigma = varargin{2};
    cutoff = varargin{3};
    order = 2;
elseif length(varargin) == 4
    w = varargin{1};
    sigma = varargin{2};
    cutoff = varargin{3};
    order = varargin{4};
else
    error('ConvolveBin.m : too many input arguments')
end


%% Boxcar convolution
boxc = ones(1, w/dt)/w;                                              % generate box and divide by width to normalize area to 1
fr_box = conv(bin, boxc, 'same');                                    % convolve with cutting the edges


%% Gausian convolution
kernel_x = - 1.5*sigma+dt : dt :  1.5*sigma; sigma = 0.5*sigma;         % create a x-vector for the gauss kernel and fold sigma
kernel = 1/(sqrt(2*pi)* sigma ) * exp(-(kernel_x).^2 /(2*sigma^2));     % create gaussian kernel

fr_gauss = conv(bin, kernel, 'same');                                   % convolve



%% filtered firing rate
cutoff = (2*dt).*cutoff;                                             % nyquist frequency of signal * cutoff

[a,b]= butter(order, cutoff);                                            % create filter
fr_filt = filtfilt(a,b, [bin,bin,bin])./dt;                          % apply filter bin is 3x concatenated to avoid edge effects
                                                                     % divide by dt to get real FR
                                                                    
fr_filt(fr_filt<0)=0;                                                % set negative values to 0
fr_filt = fr_filt(1+(size(fr_filt,2)/3):2*(size(fr_filt,2)/3));      % select only center piece of concatenated data


%% visualize
% figure; box on; hold on
% plot(0:dt:(length(bin)*dt)-dt, fr_filt, 'k', 'LineWidth', 2)
% plot(0:dt:(length(bin)*dt)-dt, fr_box, 'r', 'LineWidth', 2)
% plot(0:dt:(length(bin)*dt)-dt, fr_gauss, 'b', 'LineWidth', 2)
% 
% xlabel('time (s)'); ylabel('FR (Hz)')
% [mean(fr_filt), mean(fr_box), mean(fr_gauss)]


end