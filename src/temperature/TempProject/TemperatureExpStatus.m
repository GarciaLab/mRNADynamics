classdef TemperatureExpStatus
    %liveProject Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Prefix = '';
        ProjectName = '';
        Zoom = [];
        Temp_sp = [];
        Temp_obs = [];
        Temp_bath = [];
        Temp_scope = [];
        Region = '';
        
        hasFullEmbryoImages = false;
        hasFullMatchStack = false;
        hasFullMatchPost = false;
        hasExportedTIFs = false;
        hasExportedDataForLivemRNA = false;
        hasStitchedFullEmbryoImages = false;
        hasTrackedNuclei = false;
        hasCheckedNucleiSegmentation = false;
        hasFlexIntegrateSchnitzFluo = false;
        hasCheckedNuclearTracking = false;
        hasFoundAPAxis = false;
        hasAddedParticlePosition = false;
        hasCheckedDivisionTimes = false;
        hasCorrectedDivisionTime = false;
        hasCompiledNuclearProtein = false;
        hasChosenCytoplasmFrames = false;
        hasFilteredHistoneWeka = false;
        hasSegmentedCytoplasm = false;
        
        include_set = false;
        
        
        
        
        MinAP = [];
        MaxAP = [];
        
        StitchedEmbryoImagesDate = [];
        TrackedNucleiDate = [];
        CheckedNucleiSegmentationDate = [];
        FlexIntegrateFluoDate = [];
        CheckedNuclearTrackingDate = [];
        FoundAPAxisDate = [];
        AddedParticlePositionDate = [];
        CheckedDivisionTimesDate = [];
        CorrectedDivisionTimesDate = [];
        CompiledNuclearProteinDate = []'
        SetDefaultFluoDate = [];
        ConvertedCompiledNucleiDate = [];
        ExtractedNuclearFluoProfilesAPDate = [];
        PlottedNuclearFluoProfilesAPDate = [];
        ChosenCytoplasmFramesDate = [];
        FilteredHistoneWekaDate = [];
        SegmentedCytoplasmDate = [];
        
        
        
        
        
    end
    
    properties (Hidden = true)
        Comments = '';
    end
    
    
    methods
        % Constructor
        function this = TemperatureExpStatus(dataType,Prefix)

            dataStatusFilename = 'DataStatus.*';    %This naming convention is now enforced inside findDataStatus.m
            
            % Get all Dropbox/Results folders
            [~, ~, ~, ~, ~, ~, ~, ~, ~, allDropboxFolders] =  DetermineLocalFolders;
            
            % Find all DataStatus.xlsx files in all DropboxFolders
            dataStatusFolders = findDataStatus(allDropboxFolders);
            
            %Look in all DataStatus.xlsx files to find the tab specified by the dataType
            %user input
            dataStatusWithDataTypeFolder = findDataStatusTab(dataStatusFolders, dataType);
            
            %Redefine the DropboxFolder according to the DataStatus.xlsx we'll use
            dropboxFolder = dataStatusWithDataTypeFolder;
            
            %Load the contents of the DataStatus.XLSX tab we just found
            dataStatusDir = dir([dropboxFolder,filesep,dataStatusFilename]);
            dataTypeTabContents = readcell([dropboxFolder,filesep,dataStatusDir(1).name], ...
                'Sheet', dataType);
            dataStatusRownames = {};
            for i = 1:size(dataTypeTabContents, 1)
                if ~ismissing(dataTypeTabContents{i, 1})
                    dataStatusRownames{i} = dataTypeTabContents{i, 1};
                else
                    dataStatusRownames{i} = '';
                end
            end
            this.Prefix = Prefix;
            this.ProjectName = dataType;
            
            % Load all Project Prefixes and find
            allPrefixes = getProjectPrefixes(dataType);
            ColIndex = find(ismember(allPrefixes, Prefix))+1;

            %%
            % Add Zoom Info
            ZoomIndex = find(contains(dataStatusRownames, 'Zoom'));
            if isnumeric(dataTypeTabContents{ZoomIndex, ColIndex})
                this.Zoom = double(dataTypeTabContents{ZoomIndex, ColIndex});
            end
            % Add Temp_sp Info
            TempSPIndex = find(contains(dataStatusRownames, 'Temp_sp'));
            if isnumeric(dataTypeTabContents{TempSPIndex, ColIndex})
                this.Temp_sp = double(dataTypeTabContents{TempSPIndex, ColIndex});
            end
            % Add Temp_obs Info
            TempObsIndex = find(contains(dataStatusRownames, 'Temp_obs'));
            if isnumeric(dataTypeTabContents{TempObsIndex, ColIndex})
                this.Temp_obs = double(dataTypeTabContents{TempObsIndex, ColIndex});
            end
            % Add Temp_bath Info
            TempBathIndex = find(contains(dataStatusRownames, 'Temp_bath'));
            if isnumeric(dataTypeTabContents{TempBathIndex, ColIndex})
                this.Temp_bath = double(dataTypeTabContents{TempBathIndex, ColIndex});
            end
            % Add Temp_scope Info
            TempScopeIndex = find(contains(dataStatusRownames, 'Temp_scope'));
            if isnumeric(dataTypeTabContents{TempScopeIndex, ColIndex})
                this.Temp_scope = double(dataTypeTabContents{TempScopeIndex, ColIndex});
            end
            % Add Region Info
            RegionIndex = find(contains(dataStatusRownames, 'Region'));
            this.Region = dataTypeTabContents{RegionIndex, ColIndex};
            
            % Add FullEmbryo & StitchingInfo
            hasFullEmbryoIndex = find(contains(dataStatusRownames, 'hasFullEmbryoImages'));
            hFEstring =  dataTypeTabContents{hasFullEmbryoIndex, ColIndex};
            if contains(lower(hFEstring), 'yes') |  contains(lower(hFEstring), 'done')
                this.hasFullEmbryoImages = 1;
            end
            hasFullMatchStackIndex = find(contains(dataStatusRownames, 'hasFullMatchStack'));
            hFMSstring = dataTypeTabContents{hasFullMatchStackIndex, ColIndex};
            if contains(lower(hFMSstring), 'yes') |  contains(lower(hFMSstring), 'done')
                this.hasFullMatchStack = 1;
            end
            hasFullMatchPostIndex = find(contains(dataStatusRownames, 'hasFullMatchPost'));
            hFMPstring = dataTypeTabContents{hasFullMatchPostIndex, ColIndex};
            if contains(lower(hFMPstring), 'yes') |  contains(lower(hFMPstring), 'done')
                this.hasFullMatchPost = 1;
            end
            hasStitchedEmbryoImagesIndex = find(contains(dataStatusRownames, 'StitchFullEmbryoImages'));
            hSEI = dataTypeTabContents{hasStitchedEmbryoImagesIndex, ColIndex};
            if isdatetime(hSEI)
                hSEIstring = datestr(hSEI);
                hSEIdate = hSEI;
            else
                hSEIstring = hSEI;
                hSEIdate = [];
            end
            if contains(lower(hSEIstring), 'yes') |  contains(lower(hSEIstring), 'done')
                this.hasStitchedFullEmbryoImages = 1;
            elseif contains(lower(hSEIstring), 'no')
                this.hasStitchedFullEmbryoImages = 0;
            elseif isdatetime(hSEI)
                this.hasStitchedFullEmbryoImages = 1;
                this.StitchedEmbryoImagesDate = hSEI;
            elseif ismissing(hSEI)
                this.hasStitchedFullEmbryoImages = 0;
            end
            % Store ExportedTIFs info
            hasExportedTIFsIndex = find(contains(dataStatusRownames, 'ExportedTIFs'));
            if ~isempty(hasExportedTIFsIndex)
                hET = dataTypeTabContents{hasExportedTIFsIndex, ColIndex};
                if isdatetime(hET)
                    hETstring = datestr(hET);
                    hETdate = hET;
                else
                    hETstring = hET;
                    hETdate = [];
                end
                if contains(lower(hETstring), 'yes') |  contains(lower(hETstring), 'done')
                    this.hasExportedTIFs = 1;
                elseif contains(lower(hETstring), 'no')
                    this.hasExportedTIFs = 0;
                elseif isdatetime(hET)
                    this.hasExportedTIFs = 1;
                elseif ismissing(hET)
                    this.hasExportedTIFs = 0;
                end
            end
            % Store ExportedData info
            hasExportedDataForLivemRNAIndex = find(contains(dataStatusRownames, 'ExportDataForLivemRNA'));
            if ~isempty(hasExportedDataForLivemRNAIndex)
                hEDFLm = dataTypeTabContents{hasExportedDataForLivemRNAIndex, ColIndex};
                if isdatetime(hEDFLm)
                    hEDFLmstring = datestr(hEDFLm);
                    hEDFLmdate = hEDFLm;
                else
                    hEDFLmstring = hEDFLm;
                    hEDFLmdate = [];
                end
                if contains(lower(hEDFLmstring), 'yes') |  contains(lower(hEDFLmstring), 'done')
                    this.hasExportedDataForLivemRNA = 1;
                elseif contains(lower(hEDFLmstring), 'no')
                    this.hasExportedDataForLivemRNA = 0;
                elseif isdatetime(hEDFLm)
                    this.hasExportedDataForLivemRNA = 1;
                elseif ismissing(hEDFLm)
                    this.hasExportedDataForLivemRNA = 0;
                end
            end
            
            
            % Store TrackNuclei info
            hasTrackedNucleiIndex = find(contains(dataStatusRownames, 'TrackNuclei'));
            if ~isempty(hasTrackedNucleiIndex)
                hTN = dataTypeTabContents{hasTrackedNucleiIndex, ColIndex};
                if isdatetime(hTN)
                    hTNstring = datestr(hTN);
                    hTNdate = hTN;
                else
                    hTNstring = hTN;
                    hTNdate = [];
                end
                if contains(lower(hTNstring), 'yes') |  contains(lower(hTNstring), 'done')
                    this.hasTrackedNuclei = 1;
                elseif contains(lower(hTNstring), 'no')
                    this.hasTrackedNuclei = 0;
                elseif isdatetime(hTN)
                    this.hasTrackedNuclei = 1;
                    this.TrackedNucleiDate = hTN;
                elseif ismissing(hTN)
                    this.hasTrackedNuclei = 0;
                end
            end
            
            % Store CheckNucleiSegmentation info
            hasCheckedNucleiSegmentationIndex = find(contains(dataStatusRownames, 'CheckNucleiSegmentation'));
            if ~isempty(hasCheckedNucleiSegmentationIndex)
                hCNS = dataTypeTabContents{hasCheckedNucleiSegmentationIndex, ColIndex};
                if isdatetime(hCNS)
                    hCNSstring = datestr(hCNS);
                    hCNSdate = hCNS;
                else
                    hCNSstring = hCNS;
                    hCNSdate = [];
                end
                if contains(lower(hCNSstring), 'yes') |  contains(lower(hCNSstring), 'done')
                    this.hasCheckedNucleiSegmentation = 1;
                elseif contains(lower(hCNSstring), 'no')
                    this.hasCheckedNucleiSegmentation = 0;
                elseif isdatetime(hCNS)
                    this.hasCheckedNucleiSegmentation = 1;
                    this.CheckedNucleiSegmentationDate = hCNS;
                elseif ismissing(hCNS)
                    this.hasCheckedNucleiSegmentation = 0;
                end
            end
            
            % Store CheckNuclearTracking info
            hasCheckedNucleiTrackingIndex = find(contains(dataStatusRownames, 'CheckNuclearTracking'));
            if ~isempty(hasCheckedNucleiTrackingIndex)
                hCNT = dataTypeTabContents{hasCheckedNucleiTrackingIndex, ColIndex};
                if isdatetime(hCNT)
                    this.hasCheckedNuclearTracking = 1;
                    this.CheckedNuclearTrackingDate = hCNT;
                elseif ischar(hCNT) | isstring(hCNT)
                    if contains(lower(hCNT), 'yes') |  contains(lower(hCNT), 'done')
                        this.hasCheckedNuclearTracking = 1;
                    else
                        this.hasCheckedNuclearTracking = 0;
                    end
                elseif ismissing(hCNT)
                    this.hasCheckedNuclearTracking = 0;
                end
            end
            
            % Store FlexIntegrate info
            hasFlexIntegrateSchnitzFluoIndex = find(contains(dataStatusRownames, 'flexIntegrateSchnitzFluo'));
            if ~isempty(hasFlexIntegrateSchnitzFluoIndex)
                hFISF = dataTypeTabContents{hasFlexIntegrateSchnitzFluoIndex, ColIndex};
                if isdatetime(hFISF)
                    this.hasFlexIntegrateSchnitzFluo = 1;
                    this.FlexIntegrateFluoDate = hFISF;
                elseif ischar(hFISF) | isstring(hFISF)
                    if contains(lower(hFISF), 'yes') |  contains(lower(hFISF), 'done')
                        this.hasFlexIntegrateSchnitzFluo = 1;
                    else
                        this.hasFlexIntegrateSchnitzFluo = 0;
                    end
                elseif ismissing(hFISF)
                    this.hasFlexIntegrateSchnitzFluo = 0;
                end
            end
            % Store FindAPAxis info
            hasFoundAPAxisIndex = find(contains(dataStatusRownames, 'FindAPAxisTile'));
            if ~isempty(hasFoundAPAxisIndex)
                hFAPAx = dataTypeTabContents{hasFoundAPAxisIndex, ColIndex};
                if isdatetime(hFAPAx)
                    this.hasFoundAPAxis = 1;
                    this.FoundAPAxisDate = hFAPAx;
                elseif ischar(hFAPAx) | isstring(hFAPAx)
                    if contains(lower(hFAPAx), 'yes') |  contains(lower(hFAPAx), 'done')
                        this.hasFoundAPAxis = 1;
                    else
                        this.hasFoundAPAxis = 0;
                    end
                elseif ismissing(hFAPAx)
                    this.hasFoundAPAxis = 0;
                end
            end
            % Store AddParticlePosition info
            hasAddParticlePositionIndex = find(contains(dataStatusRownames, 'AddParticlePositionTile'));
            if ~isempty(hasAddParticlePositionIndex)
                hAPP = dataTypeTabContents{hasAddParticlePositionIndex, ColIndex};
                if isdatetime(hAPP)
                    this.hasAddedParticlePosition = 1;
                    this.AddedParticlePositionDate = hAPP;
                elseif ischar(hAPP) | isstring(hAPP)
                    if contains(lower(hAPP), 'yes') |  contains(lower(hAPP), 'done')
                        this.hasAddedParticlePosition = 1;
                    else
                        this.hasAddedParticlePosition = 0;
                    end
                elseif ismissing(hAPP)
                    this.hasAddedParticlePosition = 0;
                end
            end
            
            
            % Store CheckDivisionTimes info
            hasCheckedDivisionTimesIndex = find(contains(dataStatusRownames, 'CheckDivisionTimes'));
            if ~isempty(hasCheckedDivisionTimesIndex)
                hChDT = dataTypeTabContents{hasCheckedDivisionTimesIndex, ColIndex};
                if isdatetime(hChDT)
                    this.hasCheckedDivisionTimes = 1;
                    this.CheckedDivisionTimesDate = hChDT;
                elseif ischar(hChDT) | isstring(hChDT)
                    if contains(lower(hChDT), 'yes') |  contains(lower(hChDT), 'done')
                        this.hasCheckedDivisionTimes = 1;
                    else
                        this.hasCheckedDivisionTimes = 0;
                    end
                elseif ismissing(hChDT)
                    this.hasCheckedDivisionTimes = 0;
                end
            end
            
            % Store CorrectDivisionTimes info
            hasCorrectedDivisionTimesIndex = find(contains(dataStatusRownames, 'CorrectDivisionTimes'));
            if ~isempty(hasCorrectedDivisionTimesIndex)
                hCoDT = dataTypeTabContents{hasCorrectedDivisionTimesIndex, ColIndex};
                if isdatetime(hCoDT)
                    this.hasCorrectedDivisionTime = 1;
                    this.CorrectedDivisionTimesDate = hCoDT;
                elseif ischar(hCoDT) | isstring(hCoDT)
                    if contains(lower(hCoDT), 'yes') |  contains(lower(hCoDT), 'done')
                        this.hasCorrectedDivisionTime = 1;
                    else
                        this.hasCorrectedDivisionTime = 0;
                    end
                elseif ismissing(hCoDT)
                    this.hasCorrectedDivisionTime = 0;
                end
            end
            
            % Store CompiledNuclearProtein info
            hasCompiledNuclearProteinIndex = find(contains(dataStatusRownames, 'FilterCompileNuclei'));
            if ~isempty(hasCompiledNuclearProteinIndex)
                hCNP = dataTypeTabContents{hasCompiledNuclearProteinIndex, ColIndex};
                if isdatetime(hCNP)
                    this.hasCompiledNuclearProtein = 1;
                    this.CompiledNuclearProteinDate = hCNP;
                elseif ischar(hCNP) | isstring(hCNP)
                    if contains(lower(hCNP), 'yes') |  contains(lower(hCNP), 'done')
                        this.hasCompiledNuclearProtein = 1;
                    else
                        this.hasCompiledNuclearProtein = 0;
                    end
                elseif ismissing(hCNP)
                    this.hasCompiledNuclearProtein = 0;
                end
            end
       
            % Store ChooseCytoplasmFrames info
            hasChooseCytoplasmFramesIndex = find(contains(dataStatusRownames, 'ChooseCytoplasmFrames'));
            if ~isempty(hasChooseCytoplasmFramesIndex)
                hCCF = dataTypeTabContents{hasChooseCytoplasmFramesIndex, ColIndex};
                if isdatetime(hCCF)
                    this.hasChosenCytoplasmFrames = 1;
                    this.ChosenCytoplasmFramesDate = hCCF;
                elseif ischar(hCCF) | isstring(hCCF)
                    if contains(lower(hCCF), 'yes') |  contains(lower(hCCF), 'done')
                        this.hasChosenCytoplasmFrames = 1;
                    else
                        this.hasChosenCytoplasmFrames = 0;
                    end
                elseif ismissing(hCCF)
                    this.hasChosenCytoplasmFrames = 0;
                end
            end
            
            % Store filterHistoneWeka info
            hasFilterHistoneWekaIndex = find(contains(dataStatusRownames, 'filterHistoneWeka'));
            if ~isempty(hasFilterHistoneWekaIndex)
                hFHW = dataTypeTabContents{hasFilterHistoneWekaIndex, ColIndex};
                if isdatetime(hFHW)
                    this.hasFilteredHistoneWeka = 1;
                    this.FilteredHistoneWekaDate = hFHW;
                elseif ischar(hFHW) | isstring(hFHW)
                    if contains(lower(hFHW), 'yes') |  contains(lower(hFHW), 'done')
                        this.hasFilteredHistoneWeka = 1;
                    else
                        this.hasFilteredHistoneWeka = 0;
                    end
                elseif ismissing(hFHW)
                    this.hasFilteredHistoneWeka = 0;
                end
            end
            
            % Store segmentCytoplasm info
            hasSegmentedCytoplasmIndex = find(contains(dataStatusRownames, 'segmentCytoplasm'));
            if ~isempty(hasSegmentedCytoplasmIndex)
                hSC = dataTypeTabContents{hasSegmentedCytoplasmIndex, ColIndex};
                if isdatetime(hSC)
                    this.hasSegmentedCytoplasm = 1;
                    this.SegmentedCytoplasmDate = hSC;
                elseif ischar(hSC) | isstring(hSC)
                    if contains(lower(hSC), 'yes') |  contains(lower(hSC), 'done')
                        this.hasSegmentedCytoplasm = 1;
                    else
                        this.hasSegmentedCytoplasm = 0;
                    end
                elseif ismissing(hSC)
                    this.hasSegmentedCytoplasm = 0;
                end
            end
            
            % Store whether set should be included in analysis
            IncludeSetIndex = find(contains(dataStatusRownames, 'include_set'));
            if ~isempty(IncludeSetIndex)
                IncludeSetCell = dataTypeTabContents{IncludeSetIndex, ColIndex};
                if isnumeric(IncludeSetCell)
                    if IncludeSetCell == 1
                        this.include_set = 1;
                    else
                        this.include_set = 0;
                    end
                elseif ismissing(IncludeSetCell)
                    this.include_set = 0;
                end
            end
            
            % Store any comments
            CommentsIndex = find(contains(dataStatusRownames, 'Comments'));
            if ~isempty(CommentsIndex)
                CommentsCell = dataTypeTabContents{CommentsIndex, ColIndex};
                if ~ismissing(CommentsCell)
                    this.Comments = CommentsCell;
                end
            end
            [this.MinAP, this.MaxAP] = getAPrange(this);
            
            
        end
        
        
        function [MinAP, MaxAP] = getAPrange(this)
            if this.hasCompiledNuclearProtein
                liveExperiment = LiveExperiment(this.Prefix);
                load([liveExperiment.resultsFolder, filesep, 'HealthSummary.mat']);
                MinAP = floor(HealthSummary.APboundaries(1)/liveExperiment.APResolution)*liveExperiment.APResolution;
                
                MaxAP = floor(HealthSummary.APboundaries(2)/liveExperiment.APResolution)*liveExperiment.APResolution;
            else
                MinAP = [];
                MaxAP = [];
            end
   
        end
        
    end
end