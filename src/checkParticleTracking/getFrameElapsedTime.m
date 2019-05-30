function ElapsedTime = getFrameElapsedTime(FrameInfo, numFrames)
  %% JP 03/17/2019: why is this inside a try/end block that swallows all errors?
  try
    if isfield(FrameInfo, 'FileMode')
      if strcmp(FrameInfo(end).FileMode, 'TIF')

        for j = 1:numFrames
          ElapsedTime(j) = etime(datevec(FrameInfo(j).TimeString), datevec(FrameInfo(1).TimeString));
        end

      elseif strcmp(FrameInfo(end).FileMode, 'LSM') || strcmp(FrameInfo(end).FileMode, 'LSMExport') || ...
          strcmp(FrameInfo(end).FileMode, 'OMETIFF') || strcmp(FrameInfo(end).FileMode, 'LIFExport') || strcmp(FrameInfo(end).FileMode, 'LAT')

        for j = 1:numFrames
          ElapsedTime(j) = FrameInfo(j).Time - FrameInfo(1).Time;
        end

      else
        error('File mode not supported. Cannot extract time information. Include format in ExportDataForLivemRNA.m')
      end

    else
      warning('No FileMode information found. Assuming that this is TIF from the 2-photon.')

      for j = 1:numFrames
        ElapsedTime(j) = etime(datevec(FrameInfo(j).TimeString), datevec(FrameInfo(1).TimeString));
      end

    end
  end
  
  ElapsedTime = ElapsedTime / 60; %Time is in minutes
end