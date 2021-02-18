function applyCAR_batch
%
% applies common average reference (CAR) subtraction to a set of *.bin files
%
% subfunction: applyCAR.m
%
% ---------------------------------------------------------------------------------------------------------------------
% SYNTAX
% applyCAR_batch
%
% ---------------------------------------------------------------------------------------------------------------------
% DESCRIPTION
% applyCAR_batch.m performs a CAR subtraction on a batch of input data files. 
%
% inputs:
% script will ask for a target folder with files to be processed via user interface
%
% ---------------------------------------------------------------------------------------------------------------------
% author:  V. Hofmann
% last changes: 2019-08-26
% ---------------------------------------------------------------------------------------------------------------------


%%

dire = uigetdir('D:\RFProject\');

%%
%dire = uigetdir;
files = dir([dire, '\*.ap.bin']);

disp (['found ', num2str(length(files)), ' files in folder']);
for i = 1:length(files)
    disp (['processing ', files(i).name])
    filename = [dire, '\', files(i).name];
    applyCAR(filename, 1);
end
disp (['Successfully applied CAR to ', num2str(length(files)), ' files']);