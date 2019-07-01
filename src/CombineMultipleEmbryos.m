function saveFolder = CombineMultipleEmbryos(DataType, varargin)
% CombineMultipleEmbryos(DataType,Prefix)
%
% DESCRIPTION
% This function combines multiple embryos into a single embryo that can be
% passed to FitMeanAPSymmetric just like a normal embryo. It primarily
% combines nc13 data, but if the start of nc12 was caught, it can combine
% nc12 data as well.
%
% PARAMETERS
% DataType: This is a string that is identical to the name of the tab in
% dataStatus.xlsx that you wish to analyze.
% Prefix: Prefix of the data set to analyze
%
% OUTPUT
% Minimum required variables for FitMeanAP (nc12, nc13, nc14, NParticlesAP,
%MeanVectorsAP, SDVectorAP, ElapsedTime, APDivision) corresponding to the
%embryos combined
%
% Author (contact): Meghan Turner (meghan_turner@berkeley.edu [Could also be slack])
% Created: 01/01/2016
% Last Updated: 12/31/2016
%
% Commented by: Simon Alamos (simon.alamos@berkeley.edu)

%allData stores the following variables from each of the embryos to be
%averaged:
%APFilter, APbinArea, APbinID, AllTracesVector, CompiledParticles,
%ElapsedTime, EllipsePos, EllipsesOnAP, FilteredParticlesPos, MaxAPIndex,
%MaxCyto, MaxFrame, MeanCyto, MeanOffsetVector, MeanSlopeVectorAP,
%MeanVectorAP, MeanVectorAll, MeanVectorAllAP, MeanVectorAnterior,
%MedianCyto, MinAPIndex, NEllipsesAP, NOffsetParticles, NParticlesAP,
%NParticlesAll, NSlopeAP, NewCyclePos, OnRatioAP, ParticleCountAP,
%ParticleCountProbAP, Prefix, SDCyto, SDOffsetVector, SDSlopeVectorAP,
%SDVectorAP, SDVectorAll, SEVectorAllAP, StemLoopEnd, TotalEllipsesAP,
%nc10, nc11, nc12, nc13, nc14, nc9, ncFilter, ncFilterID, SetName,
%APDivision, schnitzcells, Ellipses, Particles

dim = 2;
for i = 1:length(varargin)
    if strcmpi(varargin{i}, '3D')
        dim = 3;
    end
end

if ischar(DataType)
    [allData, Prefixes, resultsFolder] = LoadMS2Sets(DataType);
else
    allData = DataType;
    DataType = inputname(1);
end

[Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
	Channel1, Channel2,Objective, Power,  DataFolder, DropboxFolderName, Comments,...
    nc9, nc10, nc11, nc12, nc13, nc14, CF, Channel3,prophase,metaphase, anaphase, DVResolution] = getExperimentDataFromMovieDatabase(Prefixes{1}, resultsFolder);

dv = strcmpi(ExperimentAxis, 'DV');
if dv
    allDataParticles = [allData.Particles];
    allDataNuclei = [allData.Nuclei];
    allDataOld = allData;
    allData = allDataParticles;
end

%% Find and save some useful variables

%Find number of embryos being combined
    
allData(cellfun(@isempty,{allData.NParticlesAP})) = []; %don't include embryos that have no particles
numEmbryos = size(allData,2);

%Find the total number of frames for each embryo
numFrames = zeros(1, numEmbryos);
maxAPIndex = zeros(1, numEmbryos);
if dv
    maxDVIndex = zeros(1, numEmbryos);
end

maxTime = zeros(1, numEmbryos);
for embryo = 1:numEmbryos
    numFrames(embryo) = size(allData(embryo).ElapsedTime, 2);
    if ~isempty(allData(embryo).MaxAPIndex)
        maxAPIndex(embryo) = allData(embryo).MaxAPIndex;
    else
        maxAPIndex(embryo) = 41; %assume standard AP resolution if unspecified.
    end
    if ~isempty(allData(embryo).MaxDVIndex)
        maxDVIndex(embryo) = allData(embryo).MaxDVIndex;
    else
        maxDVIndex(embryo) = 20; %assume standard DV resolution if unspecified.
    end
    maxTime(embryo) = allData(embryo).ElapsedTime(numFrames(embryo));
end


allData(isempty(allData(embryo).NParticlesAP)) = [];
numEmbryos = length(allData);

%Store the number of AP bins (this should always be 41).
numAPBins = maxAPIndex(1);
timeStep = median(diff([allData.ElapsedTime]));

if dv
    numDVBins = maxDVIndex(1);
end

%Store all variables to be combined in a single structure.
combinedData = struct('APDivision',{},'ElapsedTime',{},'NParticlesAP',{},...
    'MeanVectorAP',{},'SDVectorAP',{}, ...
    'OnRatioLineageAP',{});

if dv
   combinedData.DVDivision = {};
   combinedData.NParticlesDV = {};
   combinedData.MeanVectorDV = {};
   combinedData.SDVectorDV = {};
   combinedData.OnRatioLineageDV = {};
end

%% ALIGN ElapsedTime AND SHIFT OTHER VARIABLE TO MATCH NEW TIME VECTOR

% Define the new ElapsedTime vector for the combined embryo. This
% ElapsedTime variable has evenly spaces time points (40 seconds, 0.6667
% minutes) and is long enough to accomodate all the data for the longest
% running embryo
% timestep = .6667; %min
maxElapsedTime = max(ceil(maxTime./timeStep));
ElapsedTime = timeStep*(0:maxElapsedTime);

spotChannels = getCoatChannel(Channel1, Channel2, varargin);
if length(spotChannels) == 1
    ch = spotChannels;
else
    error('uh oh too many spot channels.')
end
    

%This inserts the variables to be combined into a structure such that there
%is room to shift the elements around as needed to align by ElapsedTime.
for embryo = 1:numEmbryos
    if ~isempty(allData(embryo).NParticlesAP)
        %ElapsedTime
        combinedData(embryo).ElapsedTime = NaN(1,maxElapsedTime + 1);
        combinedData(embryo).ElapsedTime(1,1:length(allData(embryo).ElapsedTime))...
            = allData(embryo).ElapsedTime;
        %NParticlesAP
        if iscell(allData(embryo).NParticlesAP)
            allData(embryo).NParticlesAP = allData(embryo).NParticlesAP{ch};
            if dim == 3
                allData(embryo).MeanVector3DAP = allData(embryo).MeanVector3DAP{ch};
            else
                allData(embryo).MeanVectorAP = allData(embryo).MeanVectorAP{ch};
                
            end
            allData(embryo).SDVectorAP = allData(embryo).SDVectorAP{ch};
        end
        combinedData(embryo).NParticlesAP = zeros((maxElapsedTime + 1), numAPBins);
        combinedData(embryo).NParticlesAP(1:size(allData(embryo).NParticlesAP,1),:)...
            = allData(embryo).NParticlesAP;
        %MeanVectorAP
        
        combinedData(embryo).MeanVectorAP = NaN((maxElapsedTime + 1), numAPBins);
        if dim == 3
            combinedData(embryo).MeanVectorAP(1:size(allData(embryo).MeanVector3DAP,1),:)...
                = allData(embryo).MeanVector3DAP;
        else
            combinedData(embryo).MeanVectorAP(1:size(allData(embryo).MeanVectorAP,1),:)...
                = allData(embryo).MeanVectorAP;
        end
        %SDVectorAP
        
        combinedData(embryo).SDVectorAP = NaN((maxElapsedTime + 1), numAPBins);
        combinedData(embryo).SDVectorAP(1:size(allData(embryo).SDVectorAP,1),:)...
            = allData(embryo).SDVectorAP;
        %APDivision
        combinedData(embryo).APDivision = allData(embryo).APDivision;
        
    end
end

%Adjust all the individual elapsed time vectors to align with the correct
%bin in the combined ElasedTime vector. A data point is considered to be in
%a time bin if it is within 20 seconds (0.33335 minutes), plus or minus, of
%the correct time (found in ElapsedTime).
for embryo = 1:numEmbryos
    aligned = 0;
    while ~aligned
        difference = combinedData(embryo).ElapsedTime - ElapsedTime;
        if isempty(difference(difference >= timeStep/2))
            aligned = 1;
        else
            index = find((difference >= timeStep/2), 1, 'first');
            %Shift the ElapsedTime vector
            combinedData(embryo).ElapsedTime((index+1):(maxElapsedTime+1)) = ...
                combinedData(embryo).ElapsedTime(index:(maxElapsedTime+1-1));
            combinedData(embryo).ElapsedTime(index) = NaN;
            
            %Shift the NParticlesAP matrix
            combinedData(embryo).NParticlesAP((index+1):(maxElapsedTime+1),:) = ...
                combinedData(embryo).NParticlesAP(index:(maxElapsedTime+1-1),:);
            combinedData(embryo).NParticlesAP(index, :) = 0;
            
            %Shift the MeanVectorAP matrix
            combinedData(embryo).MeanVectorAP((index+1):(maxElapsedTime+1),:) = ...
                combinedData(embryo).MeanVectorAP(index:(maxElapsedTime+1-1),:);
            combinedData(embryo).MeanVectorAP(index, :) = NaN;
            
            %Shift the SDVectorAP matrix
            combinedData(embryo).SDVectorAP((index+1):(maxElapsedTime+1),:) = ...
                combinedData(embryo).SDVectorAP(index:(maxElapsedTime+1-1),:);
            combinedData(embryo).SDVectorAP(index, :) = NaN;
            
            %Adjust APDivision. Every time the other variables are shifted,
            %any AP bin that has a start time frame greater than or equal
            %to the shifted index must be incremented by one to ensure that
            %APDivision still refers to the correct time and data.
            binsAboveIndex = combinedData(embryo).APDivision >= index;
            combinedData(embryo).APDivision = combinedData(embryo).APDivision + ...
                binsAboveIndex;
        end
    end
end

%% COMBINE NParticlesAP, MeanVectorAP, SDVectorAP, AND APDivision

% NParticlesAP, MeanVectorAP, and SDVectorAP are (elapsed time)x(number of
% AP bins) large, where the elapsed time is the combined embryo ElapsedTime
% vector defined just above
NParticlesAP = zeros(maxElapsedTime + 1, max(maxAPIndex));
MeanVectorAP = zeros(maxElapsedTime + 1, max(maxAPIndex));
SDVectorAP = zeros(maxElapsedTime + 1, max(maxAPIndex));
APDivision = zeros(14, numAPBins);

%Keeping track of how many embryos contributed to a certain Mean
%Fluorescence Average
numEmbryosWithMeanFluor = zeros(maxElapsedTime + 1, max(maxAPIndex));


for APBin = 1:max(maxAPIndex)
    %Determine which frame all embryos should be aligned to for this AP bin
    %for nc13
    
    startFrames13 = zeros(1,numEmbryos);
    endFrames13 = zeros(1,numEmbryos);
    startFrames12 = zeros(1,numEmbryos);
    endFrames12 = zeros(1,numEmbryos);
    
    for embryo = 1:numEmbryos
        startFrames12(1,embryo) = combinedData(embryo).APDivision(12, APBin);
        endFrames12(1,embryo) = combinedData(embryo).APDivision(13, APBin) - 1;
        
        startFrames13(1,embryo) = combinedData(embryo).APDivision(13, APBin);
        endFrames13(1,embryo) = combinedData(embryo).APDivision(14, APBin) - 1;
    end
    combinedStart13 = max(startFrames13);
    combinedEnd13 = combinedStart13 + max(endFrames13 - startFrames13);
    
    combinedStart12 = max(startFrames12);
    combinedEnd12 = combinedStart12 + max(endFrames12 - startFrames12);
    
    combinedStart14 = combinedEnd13 + 1;
    
    
    %Only add embryos' data to the combined variables if at least one
    %embryo exists in this AP bin
    if combinedStart13 ~= 0
        %Combine APDivision
        APDivision(12, APBin) = combinedStart12;
        APDivision(13, APBin) = combinedStart13;
        APDivision(14, APBin) = combinedStart14;
        
        %Align and combine the embryos' data for this AP bin
        for embryo = 1:numEmbryos
            begin13 = combinedData(embryo).APDivision(13, APBin);
            end13 = combinedData(embryo).APDivision(14, APBin) - 1;
            length13 = end13 - begin13;
            
            begin12 = combinedData(embryo).APDivision(12, APBin);
            end12 = combinedData(embryo).APDivision(13, APBin) - 1;
            length12 = end12 - begin12;
            
            %Only bother adding the embryo to the combined variables if the
            %embryo exists in this AP bin
            if begin13 ~= 0
                %Combine NParticlesAP for all embryos
                newParticles = combinedData(embryo).NParticlesAP(begin13:end13,APBin);
                NParticlesAP(combinedStart13:(combinedStart13 + length13),APBin) = ...
                    NParticlesAP(combinedStart13:(combinedStart13 + length13),APBin) ...
                    + newParticles;
                
                %Combine MeanVectorAP for all embryos
                newMeans = combinedData(embryo).MeanVectorAP(begin13:end13,APBin);
                newMeans(isnan(newMeans)) = 0;      %Replacing all the NaN's with
                %zeros to be able to add
                %matrices
                newMeans = newMeans .* newParticles;
                newEmbryosWithMeanFluor = newMeans~=0;
                MeanVectorAP(combinedStart13:(combinedStart13 + length13),APBin) = ...
                    MeanVectorAP(combinedStart13:(combinedStart13 + length13),APBin)...
                    + newMeans;
                
                numEmbryosWithMeanFluor(combinedStart13:(combinedStart13 + length13),APBin) = ...
                    numEmbryosWithMeanFluor(combinedStart13:(combinedStart13 + length13),APBin)...
                    + newEmbryosWithMeanFluor;
                
                %Combine SDVectorAP for all embryos
                %This will be done by squaring the SD to get the Variances,
                %averaging all the Variances, and finally taking the square root to
                %regain the SDs.
                newSDs = combinedData(embryo).SDVectorAP(begin13:end13,APBin);
                newSDs(isnan(newSDs)) = 0;          %Replacing all the NaN's with
                %zeros to be able to add
                %matrices
                newVars = newSDs.^2;     %Changing SDs to Variances
                newVars = newVars .* newParticles;
                SDVectorAP(combinedStart13:(combinedStart13 + length13),APBin) = ...
                    SDVectorAP(combinedStart13:(combinedStart13 + length13),APBin)...
                    + newVars;
            end
        end
    end
    
    
    if combinedStart12 ~= 0
        %Combine APDivision
        APDivision(12, APBin) = combinedStart12;
        
        %Align and combine the embryos' data for this AP bin
        for embryo = 1:numEmbryos
            
            
            begin12 = combinedData(embryo).APDivision(12, APBin);
            end12 = combinedData(embryo).APDivision(13, APBin) - 1;
            length12 = end12 - begin12;
            
            %Only bother adding the embryo to the combined variables if the
            %embryo exists in this AP bin
            if begin12 ~= 0
                %Combine NParticlesAP for all embryos
                newParticles = combinedData(embryo).NParticlesAP(begin12:end12,APBin);
                NParticlesAP(combinedStart12:(combinedStart12 + length12),APBin) = ...
                    NParticlesAP(combinedStart12:(combinedStart12 + length12),APBin) ...
                    + newParticles;
                
                %Combine MeanVectorAP for all embryos
                newMeans = combinedData(embryo).MeanVectorAP(begin12:end12,APBin);
                newMeans(isnan(newMeans)) = 0;      %Replacing all the NaN's with
                %zeros to be able to add
                %matrices
                newMeans = newMeans .* newParticles;
                newEmbryosWithMeanFluor = newMeans~=0;
                MeanVectorAP(combinedStart12:(combinedStart12 + length12),APBin) = ...
                    MeanVectorAP(combinedStart12:(combinedStart12 + length12),APBin)...
                    + newMeans;
                
                numEmbryosWithMeanFluor(combinedStart12:(combinedStart12 + length12),APBin) = ...
                    numEmbryosWithMeanFluor(combinedStart12:(combinedStart12 + length12),APBin)...
                    + newEmbryosWithMeanFluor;
                
                %Combine SDVectorAP for all embryos
                %This will be done by squaring the SD to get the Variances,
                %averaging all the Variances, and finally taking the square root to
                %regain the SDs.
                newSDs = combinedData(embryo).SDVectorAP(begin12:end12,APBin);
                newSDs(isnan(newSDs)) = 0;          %Replacing all the NaN's with
                %zeros to be able to add
                %matrices
                newVars = newSDs.^2;     %Changing SDs to Variances
                newVars = newVars .* newParticles;
                SDVectorAP(combinedStart12:(combinedStart12 + length12),APBin) = ...
                    SDVectorAP(combinedStart12:(combinedStart12 + length12),APBin)...
                    + newVars;
            end
        end
    end
    
    
    
end

%Replace all remaining zeros with NaN's (reversing what was done above ^)
MeanVectorAP(MeanVectorAP == 0) = NaN;
SDVectorAP(SDVectorAP == 0) = NaN;

%Take the mean of MeanVectorAP and SDVectorAP
MeanVectorAP = MeanVectorAP./NParticlesAP;
SDVectorAP = SDVectorAP./NParticlesAP;
%Regain the SD (instead of the Variance)
SDVectorAP = sqrt(SDVectorAP);

%Need to check that zeros in SDVectorAP weren't mistakenly changed to NaN
%by the previous lines of code. (This is almost guaranteed to have happened
%based on how the previous part was implemented.)
for m = 1:size(MeanVectorAP, 1)
    for n = 1:size(MeanVectorAP, 2)
        if ~isnan(MeanVectorAP(m,n)) && isnan(SDVectorAP(m,n))
            SDVectorAP(m,n) = 0;
        end
    end
end


%% COMBINE EllipsesOnProbAP
OnRatioLineageAP = zeros(numAPBins,3);
contributingElements = zeros(numAPBins,3);
for embryo = 1:numEmbryos
    if isfield(allData(embryo),'OnRatioLineageAP')
        newRatio = allData(embryo).OnRatioLineageAP;     %Store embryo's ratio data
        newElements = ~isnan(newRatio);             %Store which elements
        %conribute to combined ratio
        newRatio(isnan(newRatio)) = 0;              %Change NaN to zero to be
        %able to add matrices
        
        OnRatioLineageAP = OnRatioLineageAP + newRatio;
        contributingElements = contributingElements + newElements;
    end
end

OnRatioLineageAP = OnRatioLineageAP./contributingElements;
%% DEFINE NEW, COMBINED START FRAMES FOR nc12, nc13, AND n14

%As written right now, nc12 isn't analyzed with this code so I'm just
%saying that the start of nc12 doesn't exist
% nc12 = 0;

nc12 = min(nonzeros(APDivision(12,:)));

%Define nc13 as the earliest nc13 start time found in the combined
%APDivision matrix
nc13 = min(nonzeros(APDivision(13,:)));

%Define nc14 as the earliest nc14 start time found in the combined
%APDivision matrix
nc14 = min(nonzeros(APDivision(14,:)));

%% SAVE DATA THAT DOESN'T NEED TO BE COMBINED OR CHANGED

% APbinID just stores the fraction of the way along the AP axis that
% corresponds to the particular bin. This does not differ at all between
% embryos, as all embryos (should) have the same number of AP bins, 41, and
% AP Bin 1 has an APbinID of 0 and AP Bin 41 has an APbinID of 1
APbinID = allData(1).APbinID;

%% SAVE ALL DATA NEEDED TO RUN 'FitMeanAP*' SCRIPTS

%Minimum required variables for FitMeanAP: nc12, nc13, nc14, NParticlesAP,
%MeanVectorsAP, SDVectorAP, ElapsedTime, APDivision.
saveFolder = [resultsFolder, filesep, DataType, '_combinedEmbryo'];
mkdir(saveFolder);
save([saveFolder,filesep,[DataType,'_Combined_CompiledParticles.mat']],...
    'nc12','nc13', 'nc14','NParticlesAP', 'MeanVectorAP', 'SDVectorAP', ...
    'ElapsedTime','APbinID', 'OnRatioLineageAP');
disp('Combined_CompiledParticles saved');
save([saveFolder,filesep,[DataType,'_Combined_APDivision.mat']],'APDivision');
disp('Combined_APDivision saved');

%Final line of code at which a breakpoint can be placed for debugging
%purposes ('A' is not an actual variable with any important meaning).
%Allows you to view and manipulate the workspace created by this function
%after running the whole thing. Shouldn't do anything, but you can delete
%if you don't need to debug.
A = 0;
