function hisImage = openHistoneImage(Prefix, PreProcPath, CurrentFrame, NDigits)

HistoneImageFileNamePrefix = [PreProcPath, filesep, Prefix, filesep, Prefix];
HistoneImageFileNameSuffix = [iIndex(CurrentFrame, NDigits), '.tif'];

try
    hisImage = imread([HistoneImageFileNamePrefix, '-His_', HistoneImageFileNameSuffix]);
catch
    
    try
        hisImage = imread([HistoneImageFileNamePrefix, '_His_', HistoneImageFileNameSuffix]);
    catch
        hisImage = 0;
    end
    
end

end