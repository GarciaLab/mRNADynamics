function CompiledNuclei = SetDefaultFluo(Prefix, FluoLabel)
% ConvertCompiledNucleiToTableArray(Prefix, varargin)
%
% DESCRIPTION
% This function converts the struct array output by the fucntion 
% CompileNuclearProtein to a table array format. 
%
% ARGUMENTS
% Prefix: prefix string of the data set to analyze. 
% varargin: A cell in which the first element nis the prefix string of the data set
%           to analyze. Subsequent elements can be the options below.
% 
% Author (contact): Gabriella Martini (martini@berkeley.edu)
% Created: 9/26/20
% Last Updated: 9/26/20
%
% Documented by: Gabriella Martini (martini@berkeley.edu)
%% 

padding = 1;
liveExperiment = LiveExperiment(Prefix);

[~,~,DefaultDropboxFolder,~,~]=...
    DetermineLocalFolders;


%Load the folder information
[SourcePath, FISHPath, DefaultDropboxFolder, DropboxFolder, MS2CodePath, PreProcPath,...
configValues, movieDatabasePath] = DetermineAllLocalFolders(Prefix);

[Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
nc9, nc10, nc11, nc12, nc13, nc14, CF] = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder);


load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat']);
load([DropboxFolder,filesep,Prefix,filesep,'APDivision.mat']);
load([DropboxFolder,filesep,Prefix,filesep,'CompiledNuclei.mat']);


CNfieldnames = fieldnames(CompiledNuclei);
if nargin == 1
    FluoLabel = chooseFluo(CNfieldnames);
    FluoTimeTraceLabel = [FluoLabel, '_TimeTrace'];
    FluoZTraceLabel = [FluoLabel, '_Z'];
    FluoFrameApprovedLabel = [FluoLabel, '_FrameApproved'];
else
    foundmatch = false;
    for i=1:length(CNfieldnames)
        fn = CNfieldnames{i};
        if contains(fn, FluoLabel)
            foundmatch = true;
            break
        end
    end
    if ~foundmatch 
       disp('Could not find matching fluorescence field. Please select a default fluorescence field.')
       FluoLabel = chooseFluo(CNfieldnames);
    end
    FluoTimeTraceLabel = [FluoLabel, '_TimeTrace'];
    FluoZTraceLabel = [FluoLabel, '_Z'];
    FluoFrameApprovedLabel = [FluoLabel, '_FrameApproved']; 
end
%% 
radiusexpr = '^Fluo(?<radstr>[0-9_]*)um_[A-Za-z]';
radiusString = regexp(FluoLabel, radiusexpr, 'names').radstr;
radiusString = strrep(radiusString, '_', '.');
IntegrationRadius = str2num(radiusString); % in microns
PixelSize = FrameInfo(end).PixelSize; % in microns

%% 
NFrames = length(FrameInfo);
ncFrames = [nc9, nc10, nc11, nc12, nc13, nc14, NFrames];
NZslices = FrameInfo(end).NumberSlices;
xSize = FrameInfo(end).PixelsPerLine;
ySize = FrameInfo(end).LinesPerFrame;
NCDeltas = {[], [], [], [], [], []};
for k=1:length(CompiledNuclei)
    cycle = CompiledNuclei(k).nc;
    anaphaseFrame= CompiledNuclei(k).anaphaseFrame;
    Frames = CompiledNuclei(k).Frames;
    FluoTimeTrace = CompiledNuclei(k).(FluoTimeTraceLabel);
    FluoZ = CompiledNuclei(k).(FluoZTraceLabel);
    FrameApproved = logical(CompiledNuclei(k).(FluoFrameApprovedLabel));
    FrameApproved(FluoTimeTrace == 0) = 0;
    xPos = CompiledNuclei(k).xPos;
    yPos = CompiledNuclei(k).yPos;
    MeanAP = CompiledNuclei(k).MeanAP;
    APbin = uint16(ceil(MeanAP/APResolution));
%     ncFrames = [APDivision(9, APbin), APDivision(10, APbin),...
%         APDivision(11, APbin), APDivision(12, APbin),APDivision(13, APbin),...
%         APDivision(14, APbin), NFrames];
    ncStart = ncFrames(cycle-8);            

    ncOffset = FrameInfo(ncStart).Time/60;
    AllFrameTimes = vertcat(FrameInfo.Time)/60;
    NCFrameTimes = AllFrameTimes(Frames)-ncOffset;

    NucleusFluoDeltas = diff(FluoTimeTrace);
    TimeDeltas = diff(NCFrameTimes);
    NucleusDeltas = NucleusFluoDeltas./(TimeDeltas.');
    NCDeltas{cycle-8} = [NCDeltas{cycle-8} NucleusDeltas];
    
    CompiledNuclei(k).FluoTimeTrace = FluoTimeTrace;
    CompiledNuclei(k).FrameApproved = FrameApproved;
    %CompiledNuclei(k).FluoDeltas = [0 NucleusDeltas];
    
    if isempty(anaphaseFrame)
        anaphaseFrame = ncFrames(cycle-8);
    end
    if anaphaseFrame < min(Frames(FrameApproved))
        CompiledNuclei(k).Flag2 = double(min(Frames(FrameApproved)) - anaphaseFrame+1)/(ncFrames(cycle-7)-anaphaseFrame+1);
    end
    if max(Frames(FrameApproved)) < ncFrames(cycle - 7)
        if (xPos(end) > 2*IntegrationRadius/PixelSize) & (xPos(end) < xSize -  2*IntegrationRadius/PixelSize) & ...
                 (yPos(end) > 2*IntegrationRadius/PixelSize) & (yPos(end) < ySize -  2*IntegrationRadius/PixelSize)
            CompiledNuclei(k).Flag3 =....
                double(ncFrames(cycle - 7)- max(Frames(FrameApproved))+1)/double(ncFrames(cycle-7)-anaphaseFrame+1);
        else
            CompiledNuclei(k).Flag3 =...
                double(ncFrames(cycle - 7)- max(Frames(FrameApproved))+1)/double(ncFrames(cycle-7)-anaphaseFrame+1);
            CompiledNuclei(k).Flag6 = 1;
        end

    end



    if length(Frames(FrameApproved)) < (ncFrames(cycle-7)-anaphaseFrame)
        CompiledNuclei(k).Flag4 = ...
            double((ncFrames(cycle-7)-anaphaseFrame)+1 - length(Frames(FrameApproved)))/double(ncFrames(cycle-7)-anaphaseFrame+1);

    end

    if length(Frames((FluoZ <= NZslices - padding) & (FluoZ > padding) & FrameApproved)) < length(Frames(FrameApproved))
        CompiledNuclei(k).Flag5 =...
            1-length(Frames((FluoZ <= NZslices - padding) & (FluoZ > padding) & FrameApproved))/length(Frames(FrameApproved));

    end

end

%% %Do Additional Flagging Here


for nc=9:14
   if (~isempty(find([CompiledNuclei(:).nc] == nc,1))) &...
      (length(NCDeltas{nc-8}) > 100)
        deltas = abs(NCDeltas{nc-8});
        maxdelta = nanmean(deltas) + 4*nanstd(deltas);
        ncNuclei = find([CompiledNuclei(:).nc] == nc);
        for i=1:length(ncNuclei)
            NCID = ncNuclei(i);
            Fluos = CompiledNuclei(NCID).FluoTimeTrace;
            
            FrameApproved = CompiledNuclei(NCID).FrameApproved;
            OriginalUnapprovedFrames = sum(~FrameApproved);
            AllFrameTimes = vertcat(FrameInfo.Time)/60;
            Times = [FrameInfo(CompiledNuclei(NCID).Frames).Time]/60;
            NucleusFluoDeltas = diff(Fluos);
            TimeDeltas = diff(Times);
            for j=1:length(NucleusFluoDeltas)
                if FrameApproved(j+1)
                    MinIdx = max([1, j-9]);
                    for k =0:(j-MinIdx)
                        if FrameApproved(j-k)
                            TestDelta = abs((Fluos(j+1)-Fluos(j-k))/(Times(j+1)-Times(j-k)));
                            if TestDelta >= maxdelta 
                                 FrameApproved(j+1) = 0;
                            end
                            break
                        end  
                    end
                end
            end
            CompiledNuclei(NCID).FrameApproved = FrameApproved;
            CompiledNuclei(NCID).Flag7 =...
                (sum(~FrameApproved)-OriginalUnapprovedFrames)/length(Frames);

        end
   end
end

%% 

try
    save([DropboxFolder,filesep,Prefix,filesep,'CompiledNuclei.mat'],...
        'CompiledNuclei','-v6');
catch
    save([DropboxFolder,filesep,Prefix,filesep,'CompiledNuclei.mat'],...
        'CompiledNuclei','-v7.3', '-nocompression');
end














