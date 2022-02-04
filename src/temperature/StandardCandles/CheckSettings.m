function SeriesMatched = CheckSettings(Prefix)


liveExperiment = LiveExperiment(Prefix);
prefix_regex = '(?<datefolder>20[0-9][0-9]-[0-9][0-9]-[0-9][0-9])-(?<subfolder>[A-Za-z0-9_-]+)';
foldersplit = regexp(liveExperiment.Prefix, prefix_regex, 'names');

rawDataFolder = [fileparts(fileparts(fileparts( liveExperiment.preFolder))),...
    filesep, 'RawDynamicsData', filesep, foldersplit.datefolder, filesep, foldersplit.subfolder, filesep];


[LIFImages, LIFMeta] = loadLIFFile(rawDataFolder);
[NSeries, NFrames, NSlices,...
        NPlanes, NChannels, Frame_Times] = getFrames(LIFMeta, {liveExperiment.experimentType});
SubFolderLineAverageRegex = '_LA(?<LANumber>[0-9]+)';
LARegexNames = regexp(foldersplit.subfolder, SubFolderLineAverageRegex, 'names');
try
    LA_refvalue = str2num(LARegexNames.LANumber);
catch
   LA_refvalue = NaN;
end
SubFolderZoomRegex = '_Zoom(?<ZoomNumber>[0-9]+)';
ZoomRegexNames = regexp(foldersplit.subfolder, SubFolderZoomRegex, 'names');
Zoom_refvalue = str2num(ZoomRegexNames.ZoomNumber);
SubfolderDimensionsRegex = '_(?<InDim>[0-9]+)x(?<OutDim>[0-9]+)_';
DimensionsRegexNames = regexp(foldersplit.subfolder, SubfolderDimensionsRegex, 'names');
try
    InDimensions_refvalue = str2num(DimensionsRegexNames.InDim);
    OutDimensions_refvalue = str2num(DimensionsRegexNames.OutDim);
catch
    InDimensions_refvalue = NaN;
    OutDimensions_refvalue = NaN;
end

SubfolderzStepsRegex =  '_z[sS]tep(?<zSteps>[0-9]+)_';
zStepsRegexNames = regexp(foldersplit.subfolder, SubfolderzStepsRegex, 'names');
zSteps_refvalue = str2num(zStepsRegexNames.zSteps);
zStepsRegex = '<DimensionDescription BitInc="[0-9]+" BytesInc="[0-9]+" Unit="(?<UnitVar>[a-zA-Z]+)" Length="(?<LengthVar>[0-9\.e-]+)" Origin="(?<OriginVar>[0-9\.e-]+)" NumberofElements="(?<SlicesVar>[0-9]+)" DimID="3"';
zStepsXMLRegex = 'ValidBeginStack="1" ValidEndStack="1" Begin="(?<BeginVar>[0-9eE\-\.]+)" End="(?<EndVar>[0-9eE\-\.]+)" Sections="(?<SectionsVar>[0-9eE\-\.]+)"';


[XMLFolder, seriesPropertiesXML, seriesXML] = getSeriesFiles(rawDataFolder);
NSeries = length(seriesXML);



for series = 1:NSeries
    xmltext = fileread([XMLFolder,filesep,seriesXML(series).name]);
    LineAccumExpressionobj = '(?<=" Line_Accumulation=").*?(?=")';
    LineAccumPossiblePositions = regexp(xmltext, LineAccumExpressionobj, 'match');
    try
        LineAccumulations(series) = str2double(LineAccumPossiblePositions{1}); 
    catch
        LineAccumulations(series) = NaN;
    end
    
    LineAvgExpressionobj = '(?<=" LineAverage=").*?(?=")';
    LineAveragePossiblePositions = regexp(xmltext, LineAvgExpressionobj, 'match');
    try
        LineAverages(series) = str2double(LineAveragePossiblePositions{1}); 
    catch
        LineAverages(series) = NaN;
    end
    
    FrameAvgExpressionobj = '(?<=" FrameAverage=").*?(?=")';
    FrameAveragePossiblePositions = regexp(xmltext, FrameAvgExpressionobj, 'match');
    try
        FrameAverages(series) = str2double(FrameAveragePossiblePositions{1}); 
    catch
        FrameAverages(series) = NaN;
    end
    
    ZoomExpressionobj = '(?<=" Zoom=").*?(?=")';
    ZoomPossiblePositions = regexp(xmltext, ZoomExpressionobj, 'match');
    try
       Zoom(series) = str2double(ZoomPossiblePositions{1}); 
    catch
       Zoom(series) = NaN;
    end
    
    OutDimensionExpressionobj = '(?<=" OutDimension=").*?(?=")';
    OutDimensionPossiblePositions = regexp(xmltext, OutDimensionExpressionobj, 'match');
    try
       OutDimensions(series) = str2double(OutDimensionPossiblePositions{1}); 
    catch
        OutDimensions(series) = NaN;
    end
    
    InDimensionExpressionobj = '(?<=" InDimension=").*?(?=")';
    InDimensionPossiblePositions = regexp(xmltext, InDimensionExpressionobj, 'match');
    try
       InDimensions(series) = str2double(InDimensionPossiblePositions{1}); 
    catch
        InDimensions(series) = NaN;
    end
    
    ScanSpeedExpressionobj = '(?<=" ScanSpeed=").*?(?=")';
    ScanSpeedPossiblePositions = regexp(xmltext, ScanSpeedExpressionobj, 'match');
    try
       ScanSpeeds(series) = str2double(ScanSpeedPossiblePositions{1}); 
    catch
        ScanSpeeds(series) = NaN;
    end
    
    zStepsXMLRegexNames= regexp(xmltext, zStepsXMLRegex,'names', 'once');
    try
        zStepVar(series) = round(10^8*abs(str2num(zStepsXMLRegexNames.BeginVar)-str2num(zStepsXMLRegexNames.EndVar))/(str2num(zStepsXMLRegexNames.SectionsVar)-1));
    catch
        zStepVar(series) = NaN;
    end
end

SeriesMatched = true(1, NSeries);

if ~isnan(zSteps_refvalue)
    if ~all(round(zStepVar,2) == round(zSteps_refvalue,2))
        SeriesMatched(round(zStepVar,2) ~= round(zSteps_refvalue,2)) = false;
    end
elseif ~all(round(zStepVar,2) == round(mode(zStepVar),2))
    SeriesMatched(round(zStepVar,2) ~= round(mode(zStepVar),2)) =false;
end
  
if ~all(round(ScanSpeeds) == round(mode(ScanSpeeds)))
    SeriesMatched(round(ScanSpeeds) ~= round(mode(ScanSpeeds))) = false;
end

if ~isnan(LA_refvalue)
    if ~all(round(LineAverages) == round(LA_refvalue))
        SeriesMatched(round(LineAverages) ~= round(LA_refvalue)) = false;
    end
elseif ~all(round(LineAverages) == round(mode(LineAverages)))
    SeriesMatched(round(LineAverages) ~= round(mode(LineAverages))) = false;
end

if ~isnan(Zoom_refvalue)
    if ~all(round(Zoom,2) == round(Zoom_refvalue,2)) 
        SeriesMatched(round(Zoom,2) ~= round(Zoom_refvalue,2))  = false;
    end
elseif ~all(round(Zoom,2) == round(mode(Zoom),2))
    SeriesMatched(round(Zoom,2) ~= round(mode(Zoom),2)) = false;
end

if ~isnan(InDimensions_refvalue)
    if ~all(round(InDimensions,2)== round(InDimensions_refvalue,2))
        SeriesMatched(round(InDimensions,2)~= round(InDimensions_refvalue,2)) = false;
    end
elseif ~all(round(InDimensions,2) == round(mode(InDimensions),2))
    SeriesMatched(round(InDimensions,2) ~= round(mode(InDimensions),2)) = false;
end

if ~isnan(OutDimensions_refvalue)
    if ~all(round(OutDimensions,2)== round(OutDimensions_refvalue,2))
        SeriesMatched(round(OutDimensions,2) ~= round(OutDimensions_refvalue,2)) = false;
    end
elseif ~all(round(OutDimensions,2) == round(mode(OutDimensions),2))
    SeriesMatched(round(OutDimensions,2) ~= round(mode(OutDimensions),2)) = false;
end

if ~all(round(LineAccumulations) == round(mode(LineAccumulations)))
    SeriesMatched(round(LineAccumulations) ~= round(mode(LineAccumulations))) = false;
end

if ~all(round(FrameAverages) == round(mode(FrameAverages)))
    SeriesMatched(round(FrameAverages) ~= round(mode(FrameAverages))) = false;
end

SeriesMatchedFinal = true(1, sum(NFrames));
CurrentIndex = 0;
for i = 1:NSeries
    SeriesMatchedFinal(CurrentIndex+1:CurrentIndex+NFrames(i)) = SeriesMatched(i);
end

SeriesMatched = SeriesMatchedFinal;
    


