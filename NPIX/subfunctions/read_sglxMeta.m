function [meta] = read_sglxMeta(binName)
%
% reads *.meta file into struct array
%
% ---------------------------------------------------------------------------------------------------------------------
% SYNTAX
% [meta] = read_sglxMeta(binName)
%
% ---------------------------------------------------------------------------------------------------------------------
% DESCRIPTION
% read_sglMeta reads meta file of defined *.bin file
%
% inputs:
% binName    -  filename of *.bin file in question. Filename should include pathinformation
%
% outputs:
% meta       -  struct array containing data of meta file
%
%
% ---------------------------------------------------------------------------------------------------------------------
% author:  V. Hofmann
% last changes: 2019-08-23
% ---------------------------------------------------------------------------------------------------------------------

%% Create the matching metafile name
[pathstr,name,~] = fileparts(binName);
metaName = strcat(name, '.meta');

%% Parse ini file into cell entries C{1}{i} = C{2}{i}
fid = fopen(fullfile(pathstr, metaName), 'r');
C = textscan(fid, '%[^=] = %[^\r\n]');
% -------------------------------------------------------------
fclose(fid);

%% New empty struct
meta = struct();

% Convert each cell entry into a struct entry
for i = 1:length(C{1})
    tag = C{1}{i};
    if tag(1) == '~'
        % remake tag excluding first character
        tag = sprintf('%s', tag(2:end));
    end
    meta = setfield(meta, tag, C{2}{i});
end
end