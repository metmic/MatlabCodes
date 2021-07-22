function [vector_strength,z_stat,gain,gainraw,phase]=GainPhaseMI(spiketimes,OXings_env,env_freq,scaling_factor,plotFig)
% code to compute gain, phase, offset and MI from neural recording to env
% sines
spiketimes=spiketimes(spiketimes>=OXings_env(1) & spiketimes<=OXings_env(end));
%     firingrate=length(spiketimes)/spiketimes(end);
phases=zeros(length(spiketimes),1);
for GH=1:length(spiketimes)
    maxi=find(OXings_env<=spiketimes(GH), 1, 'last' );
    if ~isempty(maxi)
        if maxi==length(OXings_env)
            phases(GH)=(spiketimes(GH)-OXings_env(maxi))*2*pi/(OXings_env(maxi)-OXings_env(maxi-1));
        else
            phases(GH)=(spiketimes(GH)-OXings_env(maxi))*2*pi/(OXings_env(maxi+1)-OXings_env(maxi));
        end
    else
        phases(GH)=(spiketimes(GH)-(OXings_env(1)-(OXings_env(2)-OXings_env(1))))*2*pi/(OXings_env(2)-OXings_env(1));
    end
end

edges=linspace(0,2*pi,21);
for I=1:length(phases)
    if phases(I)>(2*pi)
        phases(I)=phases(I)-(2*pi);
    end
end

values=histc(phases,edges);
values=values(1:end-1);
edges=edges(1:end-1);
vector_strength=sqrt(mean(cos(phases)).^2+mean(sin(phases)).^2);
z_stat=length(phases)*vector_strength.^2;
%     ind_int=find(z_stat>=3.5) % only plot vector strength estimates that are statistically significant (corresponds to p<0.05);
% % %     if z_stat >= 3.5 % only plot vector strength estimates that are statistically significant (corresponds to p<0.05);
% % %     else
% % %         vector_strength = NaN;
% % %     end
% %     values=values/(length(OXings_env)-1)/((edges(2)-edges(1))/2/pi/env_freq(ii));
values=values/(length(OXings_env)-1)/((edges(2)-edges(1))/2/pi/env_freq);

ind_PH=(1:1:length(edges));
edges_new=edges(ind_PH);
values_new=values(ind_PH);
indZero = (values_new>0);

if sum(indZero) >= 2 
    % output is: mean - max - min - phase of fitted function
     [BETA1,~,~,~,~]=nlinfit(edges_new,values_new',@sinewave,[max(values_new) min(values_new) env_freq max(edges(values_new==max(values_new)))]);    
    if plotFig == 1
        figure(777);clf;bar(edges,values);hold on; plot(edges,sinewave(BETA1,edges),'r');hold off;
        pause;close 777
    end
    %     dbstop if error
    
    
%     % get gain and phase from sine fit
%     yFit = sinewave(BETA1,edges);
%     yFit(yFit<0) = 0;
%     gain = abs(max(yFit)-min(yFit))/scaling_factor;
%     gainraw = abs(max(yFit)-min(yFit));
%     indmaxFit=find(yFit==max(yFit));
%     BETA1(4)=edges(indmaxFit);
    
    BETA1(4) = mod(BETA1(4),2*pi);
    while BETA1(4)<0
        BETA1(4) = BETA1(4)+2*pi;
    end
    while BETA1(4) > 2*pi
       BETA1(4) = BETA1(4)-2*pi;
    end
    
    gain = abs(BETA1(2)-BETA1(1))*2/scaling_factor;
    gainraw = abs(BETA1(2)-BETA1(1))*2;

    phase = rad2deg(BETA1(4));
else
    gain = NaN;
    gainraw = NaN;
    phase = NaN;
end
K=1;
    
