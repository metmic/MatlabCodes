function copy_sglxMeta(oldname, varargin)
%
% copies *.meta file corresponding to *.bin file under new name
%
% ---------------------------------------------------------------------------------------------------------------------
% SYNTAX
% copy_sglxMeta(oldname)
% copy_sglxMeta(oldname, newname)
%
% ---------------------------------------------------------------------------------------------------------------------
% DESCRIPTION
% copy_sglxMeta.m copies the *.meta file corresponding to the defined *.bin file under a new filename
%
% inputs:
% oldname    -  filename of *.bin file in question. Filename should include pathinformation
% newname    -  filename that should be used for target *.meta file
%               if newname is not defined copy_sglxMeta.m changes input file name according to
%               ELLRec_02_g0_t0.imec0.ap.bin' => ELLRec_02_g0_t0_CAR.imec0.ap.bin'
%
% ---------------------------------------------------------------------------------------------------------------------
% author:  V. Hofmann
% last changes: 2019-08-23
% ---------------------------------------------------------------------------------------------------------------------


%% Input evaluation
if nargin>2
    error('copy_sglxMeta: wrong number of input arguments')
elseif nargin == 2
    newname = varargin{1};
elseif nargin == 1
    [pathstr, name, ~] = fileparts(oldname);
    idx = strfind(name, '.imec');
    newname =  [pathstr filesep name(1:idx-1) '_CAR' name(idx:end) '.meta'];
end

%% Create the matching metafile name
    [pathstr,name,~] = fileparts(oldname);
    metaName = strcat(name, '.meta');

    
%% Read original meta file
    fid = fopen(fullfile(pathstr, metaName), 'r');
    C = textscan(fid, '%[^=] = %[^\r\n]');
    fclose(fid);
    
    
%% write copy of meta file with new file name      
    fid = fopen(newname, 'a+');
    for i=1:size(C{1},1)
       fprintf(fid, [C{1}{i,:}, '=', C{2}{i,:}, '\n']);
    end
    fclose (fid);
    
    
end 