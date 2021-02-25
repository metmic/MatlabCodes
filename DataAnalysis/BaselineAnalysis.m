function [r_BL, rawC, pValue, r_RAW_corrected, T] = BaselineAnalysis

% do some computations
close all
clear
clc
% for analysis, it might be important to use the correct delimiter based on the OS
%----- chose back- or forward slash based on OS ----------------------------------------
if ispc
    pLoc = '\';
elseif ismac
    pLoc = '/';
end

disp('^^^^^^ set save directory ^^^^^^')
sdir = uigetdir;
%------ load data -------------------------------
disp('^^^^^^ get directory with baseline data ^^^^^^')
dirD = uigetdir;
data = load([dirD pLoc uigetfile([dirD pLoc '*_baseline.mat'])]);
N = numel(data.Ch);
disp(['^^^ you have ' num2str(N) ' neurons in this dataset... ^^^'])

%------ do some analysis ----------------------------------------------
disp('^^^ compute baseline correlations ^^^')
combSp = nchoosek(1:N, 2);                                              % all possible pairwise combinations of neurons
% time window vector for baseline correlations
WinMin = 0.002;                                                         % min time window; in sec
WinMax = 2;                                                             % max time window; in sec
incr = 10;                                                              % increment
% generate correlation time window vector
T = logspace(log10(WinMin),log10(WinMax),incr);
%     T = [.002, .005, .01, .02, .025, .05, .0625, .075, .125, .2 , .3, .4, .5, .75, 1, 1.5, 2];
    
r_BL = NaN(numel(combSp), numel(T));                                        % pre-allocate
for I = 1:size(combSp,1)
    % BL spike count correlation
    r_BL(I,:) = SpikeCountCorr (data.bin_BL(combSp(I,1),:), data.bin_BL(combSp(I,2),:), T);
end

% baseline correlations as a function of physical distance
for I = 1:incr
    [r,p] = corrcoef(data.distBL, abs(r_BL(:,I)));
    rawC(I) = r(2);
    pValue(I) = p(2);
    clear r p
end
    
% plot data with highest distance correlation
clear mdl xs ys
T_plot = find(abs(rawC) == max(abs(rawC)));
mdl = fitlm(data.distBL, abs(r_BL(:,T_plot)));
% fit datapoints
p = polyfit(data.distBL, abs(r_BL(:,T_plot)),1);
xs = 0:1:ceil(max(data.distBL)); ys = polyval(p, xs) ;

%------ correlation as a function of timescale ----------------------------
r_RAW_corrected = r_BL(:,1:end-1);
r_RAW_corrected(mean(r_BL,2)<0,1:incr-1) = r_RAW_corrected(mean(r_BL,2)<0,1:incr-1)*-1;

m = nanmean(r_RAW_corrected);
E1 = nanstd(r_RAW_corrected)./sqrt(size(r_RAW_corrected,1));
E2 = nanstd(r_RAW_corrected);

% plot results
figure; clf;
tiledlayout(2,1)
nexttile
hold all; box on
plot(data.distBL, abs(r_BL(:,T_plot)), '.k', 'MarkerSize', 10) % actual data
xlabel('distance (\mum)', 'FontSize', 14)
ylabel('|r_{baseline}|', 'FontSize', 14)
hold on
plot(xs, ys, 'r', 'LineWidth', 2) % fit
hold off
set(gca, 'YLim', [0 1]);
title({'BL correlations as function of distance'
    ['r = ', num2str(rawC(T_plot)), ', r^2 = ', num2str(mdl.Rsquared.Adjusted)] 
    ['pairs = ', num2str(numel(data.distBL)), ', p-value = ', num2str(pValue(T_plot))]}, 'FontSize', 12)

nexttile
hold all; box on; 
plot([min(T), max(T(1:end-1))], [0 0], '--k')
patch([T(1:end-1), fliplr(T(1:end-1))], [m+E2, fliplr(m-E2)], 'k', 'FaceAlpha', '0.3', 'EdgeColor', 'none')
patch([T(1:end-1), fliplr(T(1:end-1))], [m+E1, fliplr(m-E1)], 'b', 'EdgeColor', 'none')
plot(T(1:end-1), m, 'w', 'LineWidth', 2)
xlabel('timescale (ms)', 'FontSize', 14); ylabel('|r_{baseline}|', 'FontSize', 14);
set(gca, 'XScale', 'log', 'XLim', [min(T(1:end-1)) max(T(1:end-1))])
title('baseline spike count corr.', 'FontSize', 12);
