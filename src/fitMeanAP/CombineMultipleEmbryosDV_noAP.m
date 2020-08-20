function saveFolder = CombineMultipleEmbryosDV(DataType, varargin)
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


%accept input as either project string or the dataset itself
if ischar(DataType)
    [allData, Prefixes, resultsFolder] = LoadMS2Sets(DataType);
else
    allData = DataType;
    DataType = inputname(1);
end

thisProject = liveProject(DataType);

firstExperimentPrefix = thisProject.includedExperimentNames{1};
firstExperiment = LiveExperiment(firstExperimentPrefix);

ExperimentAxis = firstExperiment.experimentAxis;
Channel1 = firstExperiment.Channel1;
Channel2 = firstExperiment.Channel2;
Channel3 = firstExperiment.Channel3;

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
        maxDVIndex(embryo) = 19; %assume standard DV resolution if unspecified.
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
if ~dv
    combinedData = struct('APDivision',{},'ElapsedTime',{},'NParticlesAP',{},...
        'MeanVectorAP',{},'SDVectorAP',{}, ...
        'OnRatioLineageAP',{});
else
    combinedData = struct('APDivision',{},'ElapsedTime',{},'NParticlesAP',{},...
        'MeanVectorAP',{},'SDVectorAP',{}, ...
        'OnRatioLineageAP',{}, 'DVDivision',{},'NParticlesDV',{},...
        'MeanVectorDV',{},'SDVectorDV',{}, ...
        'OnRatioLineageDV',{});
end
%% ALIGN ElapsedTime AND SHIFT OTHER VARIABLE TO MATCH NEW TIME VECTOR

% Define the new ElapsedTime vector for the combined embryo. This
% ElapsedTime variable has evenly spaced time points  and is long enough to accomodate all the data for the longest
% running embryo
maxElapsedTime = max(ceil(maxTime./timeStep));
ElapsedTime = timeStep*(0:maxElapsedTime);

spotChannels = getCoatChannel(Channel1, Channel2, Channel3);
if length(spotChannels) == 1
    %     ch = spotChannels;
    ch = 1;
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
            allData(embryo).NParticlesDV = allData(embryo).NParticlesDV{ch};
            
            if dim == 3
                allData(embryo).MeanVector3DAP = allData(embryo).MeanVector3DAP{ch};
                allData(embryo).MeanVector3DDV = allData(embryo).MeanVector3DDV{ch};
                
            else
                allData(embryo).MeanVectorAP = allData(embryo).MeanVectorAP{ch};
                allData(embryo).MeanVectorDV = allData(embryo).MeanVectorDV{ch};
                
            end
            allData(embryo).SDVectorAP = allData(embryo).SDVectorAP{ch};
            allData(embryo).SDVectorDV = allData(embryo).SDVectorDV{ch};
            
        end
        combinedData(embryo).NParticlesAP = zeros((maxElapsedTime + 1), numAPBins);
        combinedData(embryo).NParticlesDV = zeros((maxElapsedTime + 1), numDVBins);
        
        combinedData(embryo).NParticlesAP(1:size(allData(embryo).NParticlesAP,1),:)...
            = allData(embryo).NParticlesAP;
        
        combinedData(embryo).NParticlesDV(1:size(allData(embryo).NParticlesDV,1),:)...
            = allData(embryo).NParticlesDV;
        %MeanVectorAP
        
        combinedData(embryo).MeanVectorAP = NaN((maxElapsedTime + 1), numAPBins);
        combinedData(embryo).MeanVectorDV = NaN((maxElapsedTime + 1), numDVBins);
        
        if dim == 3
            combinedData(embryo).MeanVectorAP(1:size(allData(embryo).MeanVector3DAP,1),:)...
                = allData(embryo).MeanVector3DAP;
            combinedData(embryo).MeanVectorDV(1:size(allData(embryo).MeanVector3DDV,1),:)...
                = allData(embryo).MeanVector3DDV;
        else
            combinedData(embryo).MeanVectorAP(1:size(allData(embryo).MeanVectorAP,1),:)...
                = allData(embryo).MeanVectorAP;
            combinedData(embryo).MeanVectorDV(1:size(allData(embryo).MeanVectorDV,1),:)...
                = allData(embryo).MeanVectorDV;
        end
        %SDVectorAP
        
        combinedData(embryo).SDVectorAP = NaN((maxElapsedTime + 1), numAPBins);
        combinedData(embryo).SDVectorAP(1:size(allData(embryo).SDVectorAP,1),:)...
            = allData(embryo).SDVectorAP;
        
        combinedData(embryo).SDVectorDV = NaN((maxElapsedTime + 1), numDVBins);
        combinedData(embryo).SDVectorDV(1:size(allData(embryo).SDVectorDV,1),:)...
            = allData(embryo).SDVectorDV;
        
        
        
    end
end

%Adjust all the individual elapsed time vectors to align with the correct
%bin in the combined ElasedTime vector. A data point is considered to be in
%a time bin if it is within  plus or minus, of
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
            
            combinedData(embryo).NParticlesDV((index+1):(maxElapsedTime+1),:) = ...
                combinedData(embryo).NParticlesDV(index:(maxElapsedTime+1-1),:);
            combinedData(embryo).NParticlesDV(index, :) = 0;
            
            %Shift the MeanVectorAP matrix
            combinedData(embryo).MeanVectorAP((index+1):(maxElapsedTime+1),:) = ...
                combinedData(embryo).MeanVectorAP(index:(maxElapsedTime+1-1),:);
            combinedData(embryo).MeanVectorAP(index, :) = NaN;
            
            combinedData(embryo).MeanVectorDV((index+1):(maxElapsedTime+1),:) = ...
                combinedData(embryo).MeanVectorDV(index:(maxElapsedTime+1-1),:);
            combinedData(embryo).MeanVectorDV(index, :) = NaN;
            
            %Shift the SDVectorAP matrix
            combinedData(embryo).SDVectorAP((index+1):(maxElapsedTime+1),:) = ...
                combinedData(embryo).SDVectorAP(index:(maxElapsedTime+1-1),:);
            combinedData(embryo).SDVectorAP(index, :) = NaN;
            
            
            combinedData(embryo).SDVectorDV((index+1):(maxElapsedTime+1),:) = ...
                combinedData(embryo).SDVectorDV(index:(maxElapsedTime+1-1),:);
            combinedData(embryo).SDVectorDV(index, :) = NaN;
            
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


NParticlesDV = zeros(maxElapsedTime + 1, max(maxDVIndex));
MeanVectorDV = zeros(maxElapsedTime + 1, max(maxDVIndex));
SDVectorDV = zeros(maxElapsedTime + 1, max(maxDVIndex));
DVDivision = zeros(14, numDVBins);

%Keeping track of how many embryos contributed to a certain Mean
%Fluorescence Average
numEmbryosWithMeanFluor = zeros(maxElapsedTime + 1, max(maxAPIndex));

numEmbryosWithMeanFluorDV = zeros(maxElapsedTime + 1, max(maxDVIndex));
%%
for APBin = 1:max(maxAPIndex)
    %Determine which frame all embryos should be aligned to for this AP bin
    %for nc13
    
    startFrames13 = zeros(1,numEmbryos);
    endFrames13 = zeros(1,numEmbryos);
    startFrames12 = zeros(1,numEmbryos);
    endFrames12 = zeros(1,numEmbryos);
    
    for embryo = 1:numEmbryos
        
        
        [~,~,~,~, ~,~, ~, ~,~,~,~,~, ~, ~, movieDatabase]...
            = readMovieDatabase(Prefixes{embryo});
        
        [~, ~, ~, ~, ~, ~,...
            ~, ~,~, ~,  ~, ~, ~,...
            nc9, nc10, nc11, nc12, nc13, nc14, ~, ~,~,~, ~, ~] = getExperimentDataFromMovieDatabase(Prefixes{embryo}, movieDatabase);
        ncs = [zeros(1,8), nc9, nc10, nc11, nc12, nc13, nc14];
        
        startFrames12(1,embryo) = ncs(12);
        endFrames12(1,embryo) = ncs(13)-1;
        
        startFrames13(1,embryo) = ncs(13);
        endFrames13(1,embryo) = ncs(14)-1;
    end
    
    combinedStart13 = max(startFrames13);
    combinedEnd13 = combinedStart13 + max(endFrames13 - startFrames13);
    
    combinedStart12 = max(startFrames12);
    combinedEnd12 = combinedStart12 + max(endFrames12 - startFrames12);
    
    combinedStart14 = combinedEnd13 + 1;
    
    
    %Only add embryos' data to the combined variables if at least one
    %embryo exists in this AP bin
    if combinedStart13 ~= 0
        
        
        %Align and combine the embryos' data for this AP bin
        for embryo = 1:numEmbryos
            
            [~,~,~,~, ~,~, ~, ~,~,~,~,~, ~, ~, movieDatabase]...
                = readMovieDatabase(Prefixes{embryo});
            [~, ~, ~, ~, ~, ~,...
                ~, ~,~, ~,  ~, ~, ~,...
                nc9, nc10, nc11, nc12, nc13, nc14, ~, ~,~,~, ~, ~] = getExperimentDataFromMovieDatabase(Prefixes{embryo}, movieDatabase);
            ncs = [zeros(1,8), nc9, nc10, nc11, nc12, nc13, nc14];
            
            begin13 = ncs(13);
            end13 = ncs(14)-1;
            length13 = end13 - begin13;
            
            begin12 = ncs(12);
            end12 = ncs(13)-1;
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
        
        
        %Align and combine the embryos' data for this AP bin
        for embryo = 1:numEmbryos
            
            
            [~,~,~,~, ~,~, ~, ~,~,~,~,~, ~, ~, movieDatabase]...
                = readMovieDatabase(Prefixes{embryo});
            [~, ~, ~, ~, ~, ~,...
                ~, ~,~, ~,  ~, ~, ~,...
                nc9, nc10, nc11, nc12, nc13, nc14, ~, ~,~,~, ~, ~] = getExperimentDataFromMovieDatabase(Prefixes{embryo}, movieDatabase);
            ncs = [zeros(1,8), nc9, nc10, nc11, nc12, nc13, nc14];
            
            begin12 = ncs(12);
            end12 =ncs(13)-1;
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



%%
for DVBin = 1:max(maxDVIndex)
    %Determine which frame all embryos should be aligned to for this DV bin
    %for nc13
    
    startFrames13DV = zeros(1,numEmbryos);
    endFrames13DV = zeros(1,numEmbryos);
    startFrames12DV = zeros(1,numEmbryos);
    endFrames12DV = zeros(1,numEmbryos);
    
    for embryo = 1:numEmbryos
        
        [~,~,~,~, ~,~, ~, ~,~,~,~,~, ~, ~, movieDatabase]...
            = readMovieDatabase(Prefixes{embryo});
        [~, ~, ~, ~, ~, ~,...
            ~, ~,~, ~,  ~, ~, ~,...
            nc9, nc10, nc11, nc12, nc13, nc14, ~, ~,~,~, ~, ~] = getExperimentDataFromMovieDatabase(Prefixes{embryo}, movieDatabase);
        ncs = [zeros(1,8), nc9, nc10, nc11, nc12, nc13, nc14];
        startFrames12DV(1,embryo) = ncs(12);
        endFrames12DV(1,embryo) = ncs(13) - 1;
        
        startFrames13DV(1,embryo) = ncs(13);
        endFrames13DV(1,embryo) = ncs(14) - 1;
    end
    combinedStart13DV = max(startFrames13DV);
    combinedEnd13DV = combinedStart13DV + max(endFrames13DV - startFrames13DV);
    
    combinedStart12DV = max(startFrames12DV);
    combinedEnd12DV = combinedStart12DV + max(endFrames12DV - startFrames12DV);
    
    combinedStart14DV = combinedEnd13DV + 1;
    
    
    %Only add embryos' data to the combined variables if at least one
    %embryo exists in this DV bin
    if combinedStart13DV ~= 0
        
        
        %Align and combine the embryos' data for this DV bin
        for embryo = 1:numEmbryos
            
            [~,~,~,~, ~,~, ~, ~,~,~,~,~, ~, ~, movieDatabase]...
                = readMovieDatabase(Prefixes{embryo});
            [~, ~, ~, ~, ~, ~,...
                ~, ~,~, ~,  ~, ~, ~,...
                nc9, nc10, nc11, nc12, nc13, nc14, ~, ~,~,~, ~, ~] = getExperimentDataFromMovieDatabase(Prefixes{embryo}, movieDatabase);
            ncs = [zeros(1,8), nc9, nc10, nc11, nc12, nc13, nc14];
            begin13DV = ncs(13);
            end13DV = ncs(14) - 1;
            length13DV = end13DV - begin13DV;
            
            begin12DV = ncs(12);
            end12DV = ncs(13) - 1;
            length12DV = end12DV - begin12DV;
            
            %Only bother adding the embryo to the combined variables if the
            %embryo exists in this DV bin
            if begin13DV ~= 0
                %Combine NParticlesDV for all embryos
                newParticlesDV = combinedData(embryo).NParticlesDV(begin13DV:end13DV,DVBin);
                NParticlesDV(combinedStart13DV:(combinedStart13DV + length13DV),DVBin) = ...
                    NParticlesDV(combinedStart13DV:(combinedStart13DV + length13DV),DVBin) ...
                    + newParticlesDV;
                
                %Combine MeanVectorDV for all embryos
                newMeansDV = combinedData(embryo).MeanVectorDV(begin13DV:end13DV,DVBin);
                newMeansDV(isnan(newMeansDV)) = 0;      %Replacing all the NaN's with
                %zeros to be able to add
                %matrices
                newMeansDV = newMeansDV .* newParticlesDV;
                newEmbryosWithMeanFluorDV = newMeansDV~=0;
                MeanVectorDV(combinedStart13DV:(combinedStart13DV + length13DV),DVBin) = ...
                    MeanVectorDV(combinedStart13DV:(combinedStart13DV + length13DV),DVBin)...
                    + newMeansDV;
                
                numEmbryosWithMeanFluorDV(combinedStart13DV:(combinedStart13DV + length13DV),DVBin) = ...
                    numEmbryosWithMeanFluorDV(combinedStart13DV:(combinedStart13DV + length13DV),DVBin)...
                    + newEmbryosWithMeanFluorDV;
                
                %Combine SDVectorDV for all embryos
                %This will be done by squaring the SD to get the Variances,
                %averaging all the Variances, and finally taking the square root to
                %regain the SDs.
                newSDsDV = combinedData(embryo).SDVectorDV(begin13DV:end13DV,DVBin);
                newSDsDV(isnan(newSDsDV)) = 0;          %Replacing all the NaN's with
                %zeros to be able to add
                %matrices
                newVarsDV = newSDsDV.^2;     %Changing SDs to Variances
                newVarsDV = newVarsDV .* newParticlesDV;
                SDVectorDV(combinedStart13DV:(combinedStart13DV + length13DV),DVBin) = ...
                    SDVectorDV(combinedStart13DV:(combinedStart13DV + length13DV),DVBin)...
                    + newVarsDV;
            end
        end
    end
    
    
    if combinedStart12DV ~= 0
        
        
        %Align and combine the embryos' data for this DV bin
        for embryo = 1:numEmbryos
            
            [~,~,~,~, ~,~, ~, ~,~,~,~,~, ~, ~, movieDatabase]...
                = readMovieDatabase(Prefixes{embryo});
            [~, ~, ~, ~, ~, ~,...
                ~, ~,~, ~,  ~, ~, ~,...
                nc9, nc10, nc11, nc12, nc13, nc14, ~, ~,~,~, ~, ~] = getExperimentDataFromMovieDatabase(Prefixes{embryo}, movieDatabase);
            ncs = [zeros(1,8), nc9, nc10, nc11, nc12, nc13, nc14];
            
            begin12DV = ncs(12);
            end12DV = ncs(13) - 1;
            length12DV = end12DV - begin12DV;
            
            %Only bother adding the embryo to the combined variables if the
            %embryo exists in this DV bin
            if begin12DV ~= 0
                %Combine NParticlesDV for all embryos
                newParticlesDV = combinedData(embryo).NParticlesDV(begin12DV:end12DV,DVBin);
                NParticlesDV(combinedStart12DV:(combinedStart12DV + length12DV),DVBin) = ...
                    NParticlesDV(combinedStart12DV:(combinedStart12DV + length12DV),DVBin) ...
                    + newParticlesDV;
                
                %Combine MeanVectorDV for all embryos
                newMeansDV = combinedData(embryo).MeanVectorDV(begin12DV:end12DV,DVBin);
                newMeansDV(isnan(newMeansDV)) = 0;      %Replacing all the NaN's with
                %zeros to be able to add
                %matrices
                newMeansDV = newMeansDV .* newParticlesDV;
                newEmbryosWithMeanFluorDV = newMeansDV~=0;
                MeanVectorDV(combinedStart12DV:(combinedStart12DV + length12DV),DVBin) = ...
                    MeanVectorDV(combinedStart12DV:(combinedStart12DV + length12DV),DVBin)...
                    + newMeansDV;
                
                numEmbryosWithMeanFluorDV(combinedStart12DV:(combinedStart12DV + length12DV),DVBin) = ...
                    numEmbryosWithMeanFluorDV(combinedStart12DV:(combinedStart12DV + length12DV),DVBin)...
                    + newEmbryosWithMeanFluorDV;
                
                %Combine SDVectorDV for all embryos
                %This will be done by squaring the SD to get the Variances,
                %averaging all the Variances, and finally taking the square root to
                %regain the SDs.
                newSDsDV = combinedData(embryo).SDVectorDV(begin12DV:end12DV,DVBin);
                newSDsDV(isnan(newSDsDV)) = 0;          %Replacing all the NaN's with
                %zeros to be able to add
                %matrices
                newVarsDV = newSDsDV.^2;     %Changing SDs to Variances
                newVarsDV = newVarsDV .* newParticlesDV;
                SDVectorDV(combinedStart12DV:(combinedStart12DV + length12DV),DVBin) = ...
                    SDVectorDV(combinedStart12DV:(combinedStart12DV + length12DV),DVBin)...
                    + newVarsDV;
            end
        end
    end
    
    
    
end

%%
%Replace all remaining zeros with NaN's (reversing what was done above ^)
MeanVectorAP(MeanVectorAP == 0) = NaN;
SDVectorAP(SDVectorAP == 0) = NaN;

MeanVectorDV(MeanVectorDV == 0) = NaN;
SDVectorDV(SDVectorDV == 0) = NaN;

%Take the mean of MeanVectorAP and SDVectorAP
MeanVectorAP = MeanVectorAP./NParticlesAP;
SDVectorAP = SDVectorAP./NParticlesAP;
%Regain the SD (instead of the Variance)
SDVectorAP = sqrt(SDVectorAP);

%Take the mean of MeanVectorDV and SDVectorDV
MeanVectorDV = MeanVectorDV./NParticlesDV;
SDVectorDV = SDVectorDV./NParticlesDV;
%Regain the SD (instead of the Variance)
SDVectorDV = sqrt(SDVectorDV);


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

for m = 1:size(MeanVectorDV, 1)
    for n = 1:size(MeanVectorDV, 2)
        if ~isnan(MeanVectorDV(m,n)) && isnan(SDVectorDV(m,n))
            SDVectorDV(m,n) = 0;
        end
    end
end

%% COMBINE EllipsesOnProbAP
OnRatioLineageAP = zeros(numAPBins,3);
contributingElements = zeros(numAPBins,3);

OnRatioLineageDV = zeros(numDVBins,3);
contributingElementsDV = zeros(numDVBins,3);

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

for embryo = 1:numEmbryos
    if isfield(allData(embryo),'OnRatioLineageDV')
        newRatioDV = allData(embryo).OnRatioLineageDV;     %Store embryo's ratio data
        newElementsDV = ~isnan(newRatioDV);             %Store which elements
        %conribute to combined ratio
        newRatioDV(isnan(newRatioDV)) = 0;              %Change NaN to zero to be
        %able to add matrices
        
        OnRatioLineageDV = OnRatioLineageDV + newRatioDV;
        contributingElementsDV = contributingElementsDV + newElementsDV;
    end
end

OnRatioLineageDV = OnRatioLineageDV./contributingElementsDV;

%% DEFINE NEW, COMBINED START FRAMES FOR nc12, nc13, AND n14
ncs12 = [];
ncs13 = [];
ncs14 = [];
for embryo = 1:numEmbryos
    [~,~,~,~, ~,~, ~, ~,~,~,~,~, ~, ~, movieDatabase]...
        = readMovieDatabase(Prefixes{embryo});
    [~, ~, ~, ~, ~, ~,...
        ~, ~,~, ~,  ~, ~, ~,...
        nc9, nc10, nc11, nc12, nc13, nc14, ~, ~,~,~, ~, ~] = getExperimentDataFromMovieDatabase(Prefixes{embryo}, movieDatabase);
    ncs = [zeros(1,8), nc9, nc10, nc11, nc12, nc13, nc14];
    ncs12(embryo) = nc12;
    ncs13(embryo) = nc13;
    ncs14(embryo) = nc14;
end

nc12 = min(nonzeros(ncs12));

%Define nc13 as the earliest nc13 start time found in the combined
%APDivision matrix
nc13 = min(nonzeros(ncs13));

%Define nc14 as the earliest nc14 start time found in the combined
%APDivision matrix
nc14 = min(nonzeros(ncs14));

%% DEFINE NEW, COMBINED START FRAMES FOR nc12, nc13, AND n14

nc12DV = min(nonzeros(ncs12));

%Define nc13 as the earliest nc13 start time found in the combined
%APDivision matrix
nc13DV = min(nonzeros(ncs13));

%Define nc14 as the earliest nc14 start time found in the combined
%APDivision matrix
nc14DV = min(nonzeros(ncs14));

%% SAVE DATA THAT DOESN'T NEED TO BE COMBINED OR CHANGED

% APbinID just stores the fraction of the way along the AP axis that
% corresponds to the particular bin. This does not differ at all between
% embryos, as all embryos (should) have the same number of AP bins, 41, and
% AP Bin 1 has an APbinID of 0 and AP Bin 41 has an APbinID of 1
APbinID = allData(1).APbinID;
DVbinID = allData(1).DVbinID;


%% SAVE ALL DATA NEEDED TO RUN 'FitMeanAP*' SCRIPTS

%Minimum required variables for FitMeanAP: nc12, nc13, nc14, NParticlesAP,
%MeanVectorsAP, SDVectorAP, ElapsedTime, APDivision.
saveFolder = [resultsFolder, filesep, DataType, '_combinedEmbryo'];
mkdir(saveFolder);
saveVars = {};
saveVars = [saveVars, 'nc12','nc13', 'nc14','NParticlesAP', 'MeanVectorAP', 'SDVectorAP', ...
    'ElapsedTime','APbinID', 'OnRatioLineageAP'];
if dv
    saveVars = [saveVars, 'nc12DV','nc13DV', 'nc14DV','NParticlesDV', 'MeanVectorDV', 'SDVectorDV', ...
        'ElapsedTime','DVbinID', 'OnRatioLineageDV', 'combinedStart12DV', 'combinedStart13DV', 'combinedStart14DV'];
end
save([saveFolder,filesep,[DataType,'_Combined_CompiledParticles.mat']],...
    saveVars{:});
disp('Combined_CompiledParticles saved');


%Final line of code at which a breakpoint can be placed for debugging
%purposes ('A' is not an actual variable with any important meaning).
%Allows you to view and manipulate the workspace created by this function
%after running the whole thing. Shouldn't do anything, but you can delete
%if you don't need to debug.
A = 0;
