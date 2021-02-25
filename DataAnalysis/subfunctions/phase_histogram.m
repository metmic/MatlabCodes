function [vector_strength, phase, bins, edges, yfit] = phase_histogram(spiketimes, trigtimes, number_bins)
spiketimes = spiketimes(spiketimes > trigtimes(1) & spiketimes < trigtimes(end));
phase = zeros(length(spiketimes), 1);

if length(spiketimes)<length(trigtimes)*0.25
    phase = NaN(length(spiketimes), 1);
    yfit = NaN;
    vector_strength = NaN;
    bins = NaN;
    edges = NaN;
else
    % loop over the spiketimes
    for I = 1:length(spiketimes)
        ind = find(trigtimes <= spiketimes(I), 1, 'last' );
        phase(I) = 2 * pi * (spiketimes(I) - trigtimes(ind))...
            / (trigtimes(ind+1) - trigtimes(ind));
    end
    edges = linspace(0, 2*pi, number_bins);
    bins = histc(phase, edges);
    % % figure; hold on; box on;
    % % bar(edges, bins);
    % fit
    % edges=edges-edges(1);
    [beta,r,J,Sigma]=nlinfit(edges,bins',@sinewave,[min(bins) max(bins) max(edges)]);
    yfit = sinewave(beta,edges);
    % % hold on;
    % % plot(edges,yfit)
    
    % vector strength
    vector_strength = sqrt(mean(cos(phase)).^2 + mean(sin(phase)).^2);
end
% disp(['VS: ' num2str(vector_strength)])
