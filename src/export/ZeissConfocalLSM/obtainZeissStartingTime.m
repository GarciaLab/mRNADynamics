function StartingTime = obtainZeissStartingTime(Folder, LSMIndex, LSMMeta2, NDigits, FileCreationDate)
  % Get the starting time of this acquisition
  % This is different if I have an LSM or CZI file
  LSMDir = dir([Folder, filesep, '*.lsm']);     %Zeiss confocal, old LSM format
  CZIDir = dir([Folder, filesep, '*.czi']);     %Zeiss confocal, new CZI format
  
  if ~isempty(LSMDir)
    StartingTime = LSMMeta2.get(['TimeStamp #', iIndex(1, NDigits)]); %SEANCHANGE
  elseif ~isempty(CZIDir)
    TimeStampString = LSMMeta2.get('Global Information|Image|T|StartTime #1');

    %Get the number of days since 1/1/0000
    if isempty(TimeStampString)
        disp('StartTime was not present in CZI metadata. Using file creation date as start date');
        % JP we "fake" an experiment start date since CZI does not appear
        % to have this info, just relative timestamps.
        % FileCreationDate is a datenum, we need to reformat it
        % as the code downstream expects it ('2020-11-26T00:00:00')
        TimeStamp = datetime(FileCreationDate, 'ConvertFrom', 'datenum');
    else
        TimeStamp = datetime(TimeStampString(1:19),'InputFormat','yyyy-MM-dd''T''HH:mm:ss');
    end
    %Convert the time to seconds. We're ignoring the seconds.
    StartingTime = datenum(TimeStamp)*24*60*60; %SEANCHANGE from StartingTime(LSMIndex)
  else
    error('Wrong file format. The code should not have gotten this far.')
  end
end
