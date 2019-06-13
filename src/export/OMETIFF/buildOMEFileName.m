function FileName = buildOMEFileName(Prefix, FramesFormat, CurrentNFrame, ChannelsFormat, CurrentNChannel, ZSlicesFormat, CurrentZSlice)
  FileName = [Prefix, '_', sprintf(FramesFormat, str2double(CurrentNFrame)), '_z',...
    sprintf(ZSlicesFormat, str2double(CurrentZSlice)), '_ch', sprintf(ChannelsFormat, str2double(CurrentNChannel)), '.tif'];
end