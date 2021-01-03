% Added PreferredFileName so we can automate testing and bypass the user prompt when there are many files available.
function FrameInfo = processLIFExportMode(rawDataFolder, ProjectionType, Channels,...
    Prefix, PreProcFolder, PreferredFileNameForTest,...
    nuclearGUI, skipExtraction, skipNuclearProjection)

disp('Exporting movie file...');

cleanupObj = onCleanup(@myCleanupFun);

liveExperiment = LiveExperiment(Prefix);

resultsFolder = liveExperiment.resultsFolder;

moviePrecision = 'uint16';
hisPrecision = 'uint16';

if skipExtraction
  FrameInfo = [];
end

%Load the reference histogram for the fake histone channel
load('ReferenceHist.mat', 'ReferenceHist');

markandfind = false;
if ~skipExtraction
    %%
    %this section is being deprecated
    try
        %Loads file and metadata
        [XMLFolder, seriesPropertiesXML, seriesXML] = getSeriesFiles(rawDataFolder);
    catch
        XMLFolder = '';
        seriesPropertiesXML = '';
        seriesXML = '';
    end

    try
        if contains(seriesPropertiesXML(1).name, 'Mark_and_Find')
            markandfind = true;
        end
    catch % do nothing
    end
    %%
    [LIFImages, LIFMeta] = loadLIFFile(rawDataFolder);

    %Obtains frames information
    [NSeries, NFrames, NSlices,...
        NPlanes, NChannels, Frame_Times] = getFrames(LIFMeta);

    %use the old method(exported from lasx) if the files are exported
    %already. if they're not, just use bioformats. the lasx method is being
    %deprecated.
    if ~isempty(XMLFolder)
        timeStampRetrievalMethod = 'lasx';
    else
        timeStampRetrievalMethod = 'bioformats';
    end

    if sum(NFrames)~=0

        switch timeStampRetrievalMethod

            case 'manual'

                xml_file = [liveExperiment.rawFolder, filesep, 'lifMeta.xml'];

                if ~exist(xml_file, 'file')
                    generateLIFMetaDataXML(Prefix, xml_file);
                end

                Frame_Times = getTimeStampsFromLifXML(xml_file);

            case 'lasx'

                Frame_Times = obtainFrameTimes(XMLFolder, seriesPropertiesXML,...
                    NSeries, NFrames, NSlices, NChannels);

            case 'bioformats'

                Frame_Times = getFrameTimesFromBioFormats(LIFMeta, NSlices);

            otherwise, error('what?');

        end
    end

    [InitialStackTime, zPosition] = getFirstSliceTimestamp(NSlices,...
        NSeries, NPlanes, NChannels, Frame_Times, XMLFolder, seriesXML);

    FrameInfo = recordFrameInfo(NFrames, NSlices, InitialStackTime, LIFMeta, zPosition);

    if markandfind
        FrameInfo = repmat(FrameInfo, NSeries, 1);
    end

    save([resultsFolder, filesep, 'FrameInfo.mat'], 'FrameInfo', '-v6');

    %Find the flat field (FF) information
    LIFExportMode_flatFieldImage(LIFMeta,...
        rawDataFolder, PreProcFolder, Prefix, PreferredFileNameForTest);

    if sum(NFrames) == 0
        NFrames = ~NFrames;
    end

    % this function exports tif z stacks
    
    % 1 was added as an input for 'zslicesPadding' in exportTifStacks by
    % Clay
    exportTifStacks(LIFImages, 'LIF', NChannels, NFrames, NSlices, Prefix, ...
        moviePrecision, hisPrecision, nuclearGUI,...
        ProjectionType, Channels, ReferenceHist, skipNuclearProjection,1)  

    if nuclearGUI

        chooseAnaphaseFrames(...
            Prefix, 'ProjectionType', ProjectionType,...
            'ReferenceHist', ReferenceHist);

    end

    try close(waitbarFigure); catch; end

end


disp('Movie files exported.');


end