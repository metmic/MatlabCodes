function [vPhi, vS, xHist, yHist, xSine, ySine] = PhaseHist (tSp, OXings, varargin)
%
% Phase histogram
%
% ---------------------------------------------------------------------------------------------------------------------
% SYNTAX
%  [vPhi, vS, xHist, yHist, xSine, ySine] = PhaseHist(tSp, OXings, vis)
%
% ---------------------------------------------------------------------------------------------------------------------
% DESCRIPTION
% PhaseHistogram.m calculates and plots a phase histogram, PSTH and raster plot for the inputs
%
% inputs:
% tSp           -   vector with spike times
% OXings        -   vector with zero-crossings of stimulus
% vis           -   if varargin is not empty, data will be visualized. Default is no visualization
%
% outputs:
% vPhi          -   phase of mean vector
% vS            -   strength of mean vector. If vector strength is not significant vS will return 0
% xHist         -   x-values of PSTH
% yHist         -   y-values of PSTH
% xSine         -   x-values of fitted sine wave
% ySine         -   y-values of fitted sine wave
%
%
% ---------------------------------------------------------------------------------------------------------------------
% author:  V. Hofmann (2018-12-17)
% last changes: 2020-03-04
% ---------------------------------------------------------------------------------------------------------------------



%% phase histogram
phiSp = NaN(size(tSp));

for jj = 1:length(tSp)                                                                          % determine phase for each spike
    IDXstim = find(OXings <= tSp(jj), 1, 'last');
    phiSp(jj) = 2 * pi * (tSp(jj) - OXings(IDXstim))/(OXings(IDXstim + 1) - OXings(IDXstim));
end

nBins = 36;                                                                             % define number of bins
xHist = linspace(0,2*pi, nBins+1)';                                                     % define phase bins
yHist = histcounts(phiSp, xHist);                                                       % extract value inside bin
xHist=xHist(1:end-1)+2*pi/nBins/2;                                                      % transform bin edges to bin centers


vS = sqrt(mean(cos(phiSp)).^2 + mean(sin(phiSp)).^2);                                           % calculate vector strength from phase histogram

if vS^2*length(phiSp)<4.5                                                                       % set VS to 0 if below significance
    vS=NaN;
end

[b,~,~,~,~]=nlinfit(xHist,yHist',@(b,x)(b(1)+(b(2)-b(1)).*sin(x+b(3))),[min(yHist) max(yHist) max(xHist(yHist==max(yHist)))]);   % fit sinewave to determine phase

xSine = deg2rad(1:1:360);
ySine = (b(1)+(b(2)-b(1)).*sin(xSine+b(3)));
%ySine = (ySine./range(ySine))-min((ySine./range(ySine)));
vPhi = rad2deg(xSine(ySine==max(ySine)));               % get average phase
phiSp = phiSp(:);


%% visualize

if ~isempty(varargin) & varargin{1} ~=0
    figure;
    subplot(5,1,1)
    plot(0:360, sin(deg2rad(0:360)), 'k', 'LineWidth', 2)
    set(gca, 'XLim', [0 360])
    
    subplot(5,1,2:3); hold on; box on
    for i=1:length(OXings)-1
        phiSpEp=phiSp(tSp>OXings(i) & tSp<=OXings(i+1));
        plot([phiSpEp, phiSpEp]', [ones(length(phiSpEp),1)*i, ones(length(phiSpEp),1)*(i+1)]', 'k-')
    end
    set(gca, 'XLim', deg2rad([0, 360]), 'XTick', '', 'YTick','');
    
    subplot(5,1,4:5)
    hold on, box on,
    bar(rad2deg(xHist), yHist, 'FaceColor', [.7 .7 .7])
    plot(rad2deg(xSine), ySine, '--', 'Color', [1 0 0], 'LineWidth', 2)
   
    plot(rad2deg(xSine), ((ySine./range(ySine))-min((ySine./range(ySine))))*max(ySine), '-', 'Color', [.5 .5 .5], 'LineWidth', 2)
    plot([vPhi vPhi], [0 max(yHist)*1.15], '-k', 'LineWidth', 2)
    text(vPhi+5, max(yHist)*0.95, {['PH: ',num2str(vPhi)],['VS: ',num2str(vS)]})
    set(gca, 'XLim', [0 360])
end




end
