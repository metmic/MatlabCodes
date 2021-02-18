function [nSamps, SR] = getNSamps_batch(varargin)
%
% calculates number of samples from IMEC recorded data (with SpikeGLX) for
% files within a user defined folder (i.e. before files are concatinated
% for clustering / spike sorting). 
%
% ---------------------------------------------------------------------------------------------------------------------
% SYNTAX
% [nSamp, SR] = getNSamps_batch(varargin)
%
% ---------------------------------------------------------------------------------------------------------------------
% DESCRIPTION
% getNSamps_batch.m reads "*ap.bin" and "ap.meta" files (i.e. SpikeGLX
% recorded IMEC data) that are found within a userdefined folder and
% calculates the sampling rate and number of samples of each file. 
%
% inputs:
% GUI prompted folder has to be selected
%
% UPDATE: files don't have to be all in one folder; select parent folder
% that contains subfolders with files of interest and use varagin as the
% pattern to look for in subfolder names
%
% varargin = search pattern for subfolder names (i.e. 'T*')
%
% outputs:
% nSamp      -  vector with number of samples in all files in the file list
% SR         -  sampling rates of the RAW files
%
% a file "*folder-name_Batch.mat" will be saved to the target folder
%
% ---------------------------------------------------------------------------------------------------------------------
% author:  V. Hofmann; update Maso
% last changes: 2019-08-06


%% get folder to analyze
% Get a list of all files and folders in current folder.
folders = dir(pwd);
% Get a logical vector that tells which is a directory.
dirFlags = [folders.isdir];
% Extract only those that are directories.
subFolders = folders(dirFlags);
subFolders(1:2) = [];
% Take only subFolders that contain data
for I = 1:size(subFolders,1)
    sNames{I} = subFolders(I).name;
end
B = varargin;
TF = strncmpi(sNames,B,1);

dire = pwd;
% dire = uigetdir('D:\Data_NPIX');

%% loop over subFolders
nSamps = [];
SR = [];
k = 1;
for F = 1:length(TF)
    if TF(F) == 1
        cd(subFolders(F).name)
% %         % extract number of samples in bin file
% %         files = dir('*ap.bin');
% %         d = dir(files.name);
% %         nSamps(k,1) = d.bytes/385/2;
% %         clear files
% %         % extract real sampling rate of file
% %         files = dir('*ap.meta');
% %         meta = ReadMeta (files.name,pwd);
% %         SR(k) = str2double(meta.imSampRate);
% %         clear meta files
    %% extract number of samples in files
    files = dir('*ap.bin');
    
    for I = 1:length(files)
        d = dir(files(I).name);
        nSamps(I,1) = d.bytes/385/2;
        clear d
    end
    clear files
%     clearvars -except dire nSamps
    
    
    %% extract real sampling rates of files
    files = dir('*ap.meta');
    
    for i = 1:length(files)
        meta = ReadMeta (files(i).name, pwd);
        SR(i) = str2double(meta.imSampRate);
        clear meta
    end
% %         k = k+1;
        cd(dire)
    else
    end
end
%%
savename = [dire, '\', dire(max(strfind(dire, '\'))+1:end), '_Batch.mat'];
save(savename, 'nSamps', 'SR')


end