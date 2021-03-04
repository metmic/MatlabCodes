function NPIXpreproc(Tchan,desSR,SR,Vflag)

% NPIX recodings preprocess

%----- INPUT VARIABLES
% Tchan = number of total channels
% desSR = desired sampling rate
% SR = original SR
% Vflag = use of ELL data with spiky channels saved
% Maso Oct 2019 || update Mar 2020
%%
if ispc
    pLoc = '\';
elseif ismac
    pLoc = '/';
end
if Vflag == 0
    % set parameters for filtering
    defaultValue={'', '', '20','80','2'};
    answer = inputdlg({'start of bad channels','EOD frequency (Hz)','inner diameter (um) of Loc filter; Loc 1 is 20','outer diameter (um) of Loc filter; Loc 2 is 55','what filter (1: high pass; 2: band pass)'},...
        'Enter filtering parameters', [1 50], defaultValue);
    Bchan = str2num(answer{1});
    EODf = str2num(answer{2});
    FSi = str2num(answer{3});
    FSo = str2num(answer{4});
    flag = str2num(answer{5});
   
    files = dir([pwd, pLoc, '*.ap.bin']);

    disp (['found ', num2str(length(files)), ' files in folder']);
    for i = 1:length(files)
        disp (['processing ', files(i).name])
        filename = [pwd, pLoc, files(i).name];
        %      % select file to analyze
        %     [file,path] = uigetfile('*.bin');
        % load meta file
        [meta] = read_sglxMeta(filename);
        [sMap, Pos] = readShankMap (meta); % required for local spatial filters
        if nargin == 4
            SRo = str2num(getfield(meta,'imSampRate')); % actual sampling rate of the recording
        else
            SRo = SR;
        end
        %--- set bad channels to NaN and apply HP & Loc filter --------------------
        applyCARtoDatORI(filename, Tchan, Bchan, EODf, pwd, SRo, desSR, sMap, Pos, FSi, FSo, flag);
    end
    disp (['Successfully applied CAR to ', num2str(length(files)), ' files']);

    elseif Vflag == 1 % for car processing of data
    disp('ELL data')
    chans = 0:1:Tchan;
    applyCAR_batch
end