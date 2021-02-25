function [dist] = distFunc(comb, channel, Pos)
% Calculates the euclidiean distance between pairs of channels along
% the combination input
combCh = NaN(size(comb));                                                   % pre-allocate variable (channel# in each pair)
dist = NaN(size(comb,1),1);                                                 % pre-allocate variable (euclidian distance of main channels in each pair)

for c=1:length(comb)
    combCh(c,1:2) = [channel(comb(c,1)), channel(comb(c,2))]; % get channel# of each neuron in each pair
end

for c = 1:length(combCh)
    x1=Pos(combCh(c,1),1);                                                  % xpos of Ch1 in pair i
    y1=Pos(combCh(c,1),2);                                                  % ypos of Ch1 in pair i
    x2=Pos(combCh(c,2),1);                                                  % xpos of Ch2 in pair i
    y2=Pos(combCh(c,2),2);                                                  % ypos of Ch2 in pair i

    dist(c) = sqrt((x1-x2)^2 + (y1-y2)^2);                                  % euclidian distance between Ch1 and Ch2 in pair i in [mum]
end
end
