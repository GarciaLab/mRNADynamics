function [CompiledNucleiTable] = ConvertCompiledNucleiToTableArray(Prefix, varargin)
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
% Created: 11/29/19
% Last Updated: 12/04/19
%
% Documented by: Gabriella Martini (martini@berkeley.edu)


%% Get Information about about folders
[~,~,DefaultDropboxFolder,~,~]=...
    DetermineLocalFolders;


%Load the folder information
[SourcePath, FISHPath, DefaultDropboxFolder, DropboxFolder, MS2CodePath, PreProcPath,...
configValues, movieDatabasePath] = DetermineAllLocalFolders(Prefix);

[Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
nc9, nc10, nc11, nc12, nc13, nc14, CF] = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder);



%Load all the information
load([DropboxFolder,filesep,Prefix,filesep,'CompiledNuclei.mat']);
load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat']);
% load([DropboxFolder,filesep,Prefix,filesep,'APDivision.mat']);
% 
% APDivisionTimes = zeros(size(APDivision));
% for i=1:size(APDivision, 1)
%     for j = 1:size(APDivision, 2)
%         if APDivision(i,j) ~= 0
%             APDivisionTimes(i,j) = FrameInfo(APDivision(i,j)).Time;
%         end
%     end
% end


% 
%% 
numNuclei = size(CompiledNuclei, 2);
fns = fieldnames(CompiledNuclei);

varnames =  { 'Prefix','NucleusID', 'Frame', 'FrameNC','FrameAnaphase', 'Time', 'TimeNC','timeSinceAnaphase', 'FrameCount', 'xPos',...
    'yPos', 'Radius', 'cellno', 'nc', 'APPos', 'MeanAP', 'MedianAP', 'APbin','APbinIdx', 'APbinIdxV2',...
    'MeanDV', 'MedianDV', 'ncStart', 'Fluo', 'FrameApproved'};

 
vartypes = { 'string','int64', 'int64', 'int64', 'int64','double', 'double', 'double', 'int64', 'int64',...
    'int64', 'double', 'int64', 'int64','double', 'double', 'double','double','int64', 'int64',...
    'double', 'double','int64', 'double', 'int64'};
FluoTimeStrings = {};
FluoZStrings = {};
FluoApprovedStrings = {};
FlagStrings = {};


for i=1:length(fns)
    fn = fns{i};
    if regexp(fn, '^Fluo[A-za-z_0-9]*_TimeTrace')
        FluoTimeStrings{length(FluoTimeStrings)+1} = fn;
        FluoZStrings{length(FluoZStrings)+1} = strrep(fn, '_TimeTrace', '_Z');
        FluoApprovedStrings{length(FluoApprovedStrings)+1} = strrep(fn, '_TimeTrace', '_FrameApproved');
        varnames{length(varnames) + 1} = fn;
        vartypes{length(vartypes) + 1} = 'double';
    end
end

for i=1:length(fns)
    fn = fns{i};
    if regexp(fn, '^Fluo[A-za-z_0-9]*_TimeTrace')
        varnames{length(varnames) + 1} =  strrep(fn, '_TimeTrace', '_Z');
        vartypes{length(vartypes) + 1} = 'double';
    end
end

for i=1:length(fns)
    fn = fns{i};
    if regexp(fn, '^Fluo[A-za-z_0-9]*_TimeTrace')
        varnames{length(varnames) + 1} = strrep(fn, '_TimeTrace', '_FrameApproved');
        vartypes{length(vartypes) + 1} = 'double';
        
    end
end


for i=1:length(fns)
    fn = fns{i};
    if regexp(fn, '^Flag[A-za-z_0-9]*')
        FlagStrings{length(FlagStrings) + 1} = fn;
        varnames{length(varnames) + 1} = fn;
        vartypes{length(vartypes) + 1} = 'double';
        
    end
end

%% 

% STOPPED HERE 
nvars = length(varnames);

compnucfields = fieldnames(CompiledNuclei(1));
n = 1;
totalSamples = 0;
while n <= numNuclei
    totalSamples = totalSamples + size(CompiledNuclei(n).Frames, 1);
    n = n + 1;
end






AllPrefixes = cell(totalSamples,1);
AllNucleusID = zeros(totalSamples,1, 'uint16');
AllFrame = zeros(totalSamples,1, 'uint16');
AllFrameNC = zeros(totalSamples,1, 'uint16');
AllFrameAnaphase = zeros(totalSamples,1, 'uint16');
AllTime = zeros(totalSamples,1, 'double');
AllTimeNC = zeros(totalSamples,1, 'double');
AlltimeSinceAnaphase = zeros(totalSamples,1, 'double');
AllFrameCount = zeros(totalSamples,1, 'uint16');
AllxPos = zeros(totalSamples,1, 'uint16');
AllyPos = zeros( totalSamples,1, 'uint16');
AllRadius = zeros(totalSamples,1, 'double');
Allcellno = zeros(totalSamples,1, 'uint16');
Allnc = zeros(totalSamples,1, 'uint16');
AllAPPos = zeros(totalSamples,1, 'double');
AllMeanAP = zeros(totalSamples,1, 'double');
AllMedianAP = zeros(totalSamples,1, 'double');
AllAPbin = zeros(totalSamples,1, 'double');
AllAPbinIdx = zeros(totalSamples,1, 'uint8');
AllAPbinIdxV2 = zeros(totalSamples,1, 'uint8');
AllMeanDV = zeros(totalSamples,1, 'double');
AllMedianDV = zeros(totalSamples,1, 'double');
AllncStart = zeros(totalSamples,1, 'uint16');
AllFluo = zeros(totalSamples,1, 'double');
AllFrameApproved = zeros(totalSamples,1, 'uint16');
AllFluoTimeTraces = zeros(totalSamples, length(FluoTimeStrings), 'double');
AllFluoZTraces = zeros(totalSamples, length(FluoTimeStrings), 'uint16');
AllFluoFrameApproved= zeros(totalSamples, length(FluoTimeStrings), 'uint16');
AllFlags = zeros(totalSamples, length(FlagStrings), 'double');

sample = 1;
n = 1;
nc_frames = [nc9, nc10, nc11, nc12, nc13, nc14];
NCDeltas = {[], [], [], [], [], []};
%% 

h=waitbar(0,'Converting Compiled Nuclei');
while n <= numNuclei
     FrameCount = size(CompiledNuclei(n).Frames, 1);
     nc = CompiledNuclei(n).nc;
     schnitz = CompiledNuclei(n).schnitz;
     MeanAP = CompiledNuclei(n).MeanAP;
     MedianAP = CompiledNuclei(n).MedianAP;
     APbin = ceil(MeanAP/APResolution)*APResolution;
     APbinIdx = ceil(MeanAP/APResolution);
     MeanDV = CompiledNuclei(n).MeanDV;
     MedianDV = CompiledNuclei(n).MedianDV;
     ncStart = nc_frames(nc-8);
     ncStart = max(ncStart, 1);
     ncOffset = FrameInfo(ncStart).Time/60;
     FrameTimes = vertcat(FrameInfo.Time)/60;
     NCFrameTimes = FrameTimes(CompiledNuclei(n).Frames)-ncOffset;
      
     NucleusFluoDeltas = diff(CompiledNuclei(n).FluoTimeTrace);
     TimeDeltas = diff(NCFrameTimes);
     NucleusDeltas = NucleusFluoDeltas./(TimeDeltas.');
     NCDeltas{nc-8} = [NCDeltas{nc-8} NucleusDeltas];
     if ~isempty(CompiledNuclei(n).anaphaseFrame)
        anaphaseFrame = CompiledNuclei(n).anaphaseFrame;
     else
         anaphaseFrame = nan;
     end
     
     for f=1:FrameCount
        waitbar(sample/totalSamples,h)
        AllPrefixes{sample} = Prefix;
        %CompiledNucleiTable.Prefix(sample) = Prefix;
        %CompiledNucleiTable.NucleusID(sample) = schnitz; 
        AllNucleusID(sample) = schnitz;
        %CompiledNucleiTable.Frame(sample) = CompiledNuclei(n).Frames(f);
        AllFrame(sample) = CompiledNuclei(n).Frames(f);
        AllFrameNC(sample) = CompiledNuclei(n).Frames(f)-nc_frames(nc-8);
        if ~isempty(CompiledNuclei(n).anaphaseFrame)
            AllFrameAnaphase(sample) = CompiledNuclei(n).Frames(f)-anaphaseFrame;
        else
            AllFrameAnaphase(sample) = nan;
        end
        AllTime(sample) = FrameInfo(AllFrame(sample)).Time/60;
        AllTimeNC(sample) = AllTime(sample)-FrameInfo(nc_frames(nc-8)).Time/60;
        if ~isnan(anaphaseFrame)
            AlltimeSinceAnaphase(sample) = AllTime(sample)-FrameInfo(anaphaseFrame).Time/60;
        else
            AlltimeSinceAnaphase = nan;
        end
        AllFrameCount(sample) = FrameCount;
        AllxPos(sample) = CompiledNuclei(n).xPos(f);
        AllyPos(sample) = CompiledNuclei(n).yPos(f);
        AllRadius(sample) = CompiledNuclei(n).Radius(f);
        Allcellno(sample) = CompiledNuclei(n).cellno(f);
        Allnc(sample) = nc;
        AllAPPos(sample) = CompiledNuclei(n).APpos(f);
        AllMeanAP(sample) = MeanAP;
        AllMedianAP(sample) = MedianAP;
        AllAPbin(sample) = APbin;
        AllAPbinIdx(sample) = APbinIdx;
        AllAPbinIdxV2(sample) = ceil(AllAPPos(sample)/APResolution);
        
        AllMeanDV(sample) = MeanDV;
        AllMedianDV(sample) = MedianDV;
        
        AllncStart(sample) = ncStart;
        
        AllFluo(sample) = CompiledNuclei(n).FluoTimeTrace(f);
        AllFrameApproved(sample) = CompiledNuclei(n).FrameApproved(f);
        for j=1:length(FluoTimeStrings)
            %CompiledNucleiTable.(FluoTimeStrings{j})(sample) = CompiledNuclei(n).(FluoTimeStrings{j})(f);
            AllFluoTimeTraces(sample,j) = CompiledNuclei(n).(FluoTimeStrings{j})(f);
            %CompiledNucleiTable.(FluoZStrings{j})(sample) 
            AllFluoZTraces(sample,j)= CompiledNuclei(n).(FluoZStrings{j})(f);
            %CompiledNucleiTable.(FluoApprovedStrings{j})(sample) 
            AllFluoFrameApproved(sample,j)= CompiledNuclei(n).(FluoApprovedStrings{j})(f);
        end
        
        for j=1:length(FlagStrings)
            %CompiledNucleiTable.(FlagStrings{j})(sample) 
            AllFlags(sample,j)= double(CompiledNuclei(n).(FlagStrings{j}));
        end
        
        
        sample = sample + 1;
     end
     n = n + 1;
end
%% 
close(h)

% CompiledNucleiTable = table('Size', [totalSamples nvars],...
%     'VariableTypes', vartypes,'VariableNames', varnames);
try
    CompiledNucleiTable = table(AllPrefixes, AllNucleusID, AllFrame, AllFrameNC, AllFrameAnaphase,...
        AllTime, AllTimeNC, AlltimeSinceAnaphase, AllFrameCount, AllxPos, AllyPos, AllRadius,...
        Allcellno, Allnc, AllAPPos, AllMeanAP, AllMedianAP, AllAPbin, AllAPbinIdx, AllAPbinIdxV2, AllMeanDV, AllMedianDV,...
        AllncStart, AllFluo, AllFrameApproved, AllFluoTimeTraces, AllFluoZTraces, AllFluoFrameApproved, AllFlags);
catch
    AlltimeSinceAnaphase = AlltimeSinceAnaphase.';
    CompiledNucleiTable = table(AllPrefixes, AllNucleusID, AllFrame, AllFrameNC, AllFrameAnaphase,...
        AllTime, AllTimeNC, AlltimeSinceAnaphase, AllFrameCount, AllxPos, AllyPos, AllRadius,...
        Allcellno, Allnc, AllAPPos, AllMeanAP, AllMedianAP, AllAPbin, AllAPbinIdx,  AllAPbinIdxV2,AllMeanDV, AllMedianDV,...
        AllncStart, AllFluo, AllFrameApproved, AllFluoTimeTraces, AllFluoZTraces, AllFluoFrameApproved, AllFlags);%,...
end
%'VariableNames', varnames);
Nvars = size(CompiledNucleiTable, 2);
CompiledNucleiTableSub = table2array(CompiledNucleiTable(:,Nvars-3));
CompiledNucleiTableSub = array2table(CompiledNucleiTableSub);
CompiledNucleiTableSub2= table2array(CompiledNucleiTable(:,Nvars-2:Nvars-1));
CompiledNucleiTableSub2 = array2table(CompiledNucleiTableSub2);
CompiledNucleiTableSub3 = table2array(CompiledNucleiTable(:,Nvars));
CompiledNucleiTableSub3 = array2table(CompiledNucleiTableSub3);
CompiledNucleiTableOld = CompiledNucleiTable;
CompiledNucleiTable = horzcat([CompiledNucleiTableOld(:, 1:Nvars-4) CompiledNucleiTableSub CompiledNucleiTableSub2 CompiledNucleiTableSub3]);
CompiledNucleiTable.Properties.VariableNames = varnames;

%% 


save([DropboxFolder,filesep,Prefix,filesep,'CompiledNucleiTable.mat'],...
        'CompiledNucleiTable');

