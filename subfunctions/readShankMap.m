function [sMap, Pos] = readShankMap (meta)
%
% converts meta.snsShankMap to double matrix, calulates indexed and physical positions of electrodes
%
% ---------------------------------------------------------------------------------------------------------------------
% SYNTAX
% [meta] = read_sglxMeta(binName)
%
% ---------------------------------------------------------------------------------------------------------------------
% DESCRIPTION
% readShankMap.m reads meta data from struct array and calculates electrode positions
% to caclulate physical distance between electrode i and electrode ii use: 
% d = sqrt((Pos(i,1)-Pos(ii,1)).^2  + (Pos(i,2)-Pos(ii,2)).^2)
%
% inputs:
% meta       -  struct array generated with read_sglxMeta.m
%
% outputs:
% sMap       -  matrix of the shank map containing 4 cloumns per channel:
%               col 1 - shank#; col 2 - col#; col 3 - row# ; col 4 - on/off flag (spikeGLX visualization);        
% Pos        -  physical position in mum of electrodes with integer numbers
%               representing positions (according to manual)
%
% ---------------------------------------------------------------------------------------------------------------------
% author:  V. Hofmann
% last changes: 2019-08-23
% ---------------------------------------------------------------------------------------------------------------------
%
%% read out entries
k = strfind(meta.snsShankMap, '(');
kdiff = [diff(k), length(meta.snsShankMap)+1-k(end)];

sMap_cell = cell(1,length(k));
for i=1:length(k)
    sMap_cell{i} = meta.snsShankMap(k(i)+1:k(i)+kdiff(i)-2);
end
clear k kdiff

%% convert entries to double matrix
sMap = NaN (length(sMap_cell)-1,4);

for i = 2:length(sMap_cell)
    k = [0,strfind(sMap_cell{i},':'), length(sMap_cell{i})+1];
    for ii = 1:length(k)-1
        sMap(i-1,ii) = str2double(sMap_cell{i}(k(ii)+1:k(ii+1)-1));
    end
end

%% convert from 0-based to common MATLAB indexing
sMap(:,1:3) = sMap(:,1:3)+1;



%% create Coordinates 
Pos (mod(sMap(:,3), 2) & sMap(:,2)==1,1) = 2;      % even y and left x
Pos (mod(sMap(:,3), 2) & sMap(:,2)==2,1) = 4;     % even y and right x
Pos (~mod(sMap(:,3), 2) & sMap(:,2)==1,1) = 1;    % odd y and left x
Pos (~mod(sMap(:,3), 2) & sMap(:,2)==2,1) = 3;    % odd y and right x

Pos = [Pos, sMap(:,3)];

Pos(:,1) = Pos(:,1).*16-16;
Pos(:,2) = Pos(:,2).*20-20;

% figure; plot(Pos(:,1), Pos(:,2), '*'); axis equal; % plot probe layout in indexed position

end