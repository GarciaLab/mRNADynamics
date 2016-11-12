function CombineMultipleEmbryos(DataType, PrefixName)

%This function combines multiple embryos into a single embryo that can be
%passed to FitMeanAPSymmetric just like a normal embryo. It primarily
%combines nc13 data, but if the start of nc12 was caught, it can combine
%nc12 data as well.
%Author: Meghan Turner
%Created: 10/3/2016
%Updated: 10/26/2016

[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;

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
allData = LoadMS2Sets(DataType,PrefixName);


%% Find and save some useful variables

%Find number of embryos being combined
numEmbryos = size(allData,2);

%Find the total number of frames for each embryo
numFrames = zeros(1, numEmbryos);
maxAPIndex = zeros(1, numEmbryos);
maxTime = zeros(1, numEmbryos);
for i = 1:numEmbryos
    numFrames(i) = size(allData(i).ElapsedTime, 2);
    maxAPIndex(i) = allData(i).MaxAPIndex;
    maxTime(i) = allData(i).ElapsedTime(numFrames(i));
end

%Store the number of AP bins (this should always be 41).
numAPBins = maxAPIndex(1);

%Store all variables to be combined in a single structure. 
combinedData = struct('APDivision',{},'ElapsedTime',{},'NParticlesAP',{},...
                        'MeanVectorAP',{},'SDVectorAP',{}, ...
                        'OnRatioLineageAP',{});
                    
%% ALIGN ElapsedTime AND SHIFT OTHER VARIABLE TO MATCH NEW TIME VECTOR

% Define the new ElapsedTime vector for the combined embryo. This
% ElapsedTime variable has evenly spaces time points (40 seconds, 0.6667
% minutes) and is long enough to accomodate all the data for the longest
% running embryo
maxElapsedTime = max(ceil(maxTime./0.6667));
ElapsedTime = 0.6667*(0:maxElapsedTime);
 
%This inserts the variables to be combined into a structure such that there
%is room to shift the elements around as needed to align by ElapsedTime.
for i = 1:numEmbryos
    %ElapsedTime
    combinedData(i).ElapsedTime = NaN(1,maxElapsedTime + 1);
    combinedData(i).ElapsedTime(1,1:length(allData(i).ElapsedTime))... 
        = allData(i).ElapsedTime;
    %NParticlesAP
    combinedData(i).NParticlesAP = zeros((maxElapsedTime + 1), numAPBins);
    combinedData(i).NParticlesAP(1:size(allData(i).NParticlesAP,1),:)... 
        = allData(i).NParticlesAP;
    %MeanVectorAP
    combinedData(i).MeanVectorAP = NaN((maxElapsedTime + 1), numAPBins);
    combinedData(i).MeanVectorAP(1:size(allData(i).MeanVectorAP,1),:)... 
        = allData(i).MeanVectorAP;
    %SDVectorAP
    combinedData(i).SDVectorAP = NaN((maxElapsedTime + 1), numAPBins);
    combinedData(i).SDVectorAP(1:size(allData(i).SDVectorAP,1),:)... 
        = allData(i).SDVectorAP;
    %APDivision
    combinedData(i).APDivision = allData(i).APDivision;
end

%Adjust all the individual elapsed time vectors to align with the correct
%bin in the combined ElasedTime vector. A data point is considered to be in
%a time bin if it is within 20 seconds (0.33335 minutes), plus or minus, of
%the correct time (found in ElapsedTime).
for i = 1:numEmbryos
    aligned = 0;
    while ~aligned
        difference = combinedData(i).ElapsedTime - ElapsedTime;
        if isempty(difference(difference >= 0.33335))
            aligned = 1;
        else
            index = find((difference >= 0.33335), 1, 'first');
            %Shift the ElapsedTime vector
            combinedData(i).ElapsedTime((index+1):(maxElapsedTime+1)) = ...
                combinedData(i).ElapsedTime(index:(maxElapsedTime+1-1));
            combinedData(i).ElapsedTime(index) = NaN;
            
            %Shift the NParticlesAP matrix
            combinedData(i).NParticlesAP((index+1):(maxElapsedTime+1),:) = ...
                combinedData(i).NParticlesAP(index:(maxElapsedTime+1-1),:);
            combinedData(i).NParticlesAP(index, :) = 0;
          
            %Shift the MeanVectorAP matrix
            combinedData(i).MeanVectorAP((index+1):(maxElapsedTime+1),:) = ...
                combinedData(i).MeanVectorAP(index:(maxElapsedTime+1-1),:);
            combinedData(i).MeanVectorAP(index, :) = NaN;
            
            %Shift the SDVectorAP matrix
            combinedData(i).SDVectorAP((index+1):(maxElapsedTime+1),:) = ...
                combinedData(i).SDVectorAP(index:(maxElapsedTime+1-1),:);
            combinedData(i).SDVectorAP(index, :) = NaN;
            
            %Adjust APDivision. Every time the other variables are shifted,
            %any AP bin that has a start time frame greater than or equal
            %to the shifted index must be incremented by one to ensure that
            %APDivision still refers to the correct time and data.
            binsAboveIndex = combinedData(i).APDivision >= index;
            combinedData(i).APDivision = combinedData(i).APDivision + ... 
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
    for i = 1:numEmbryos
        startFrames13(1,i) = combinedData(i).APDivision(13, APBin);
        endFrames13(1,i) = combinedData(i).APDivision(14, APBin) - 1;
    end
    combinedStart13 = max(startFrames13);
    combinedEnd13 = combinedStart13 + max(endFrames13 - startFrames13);
    combinedStart14 = combinedEnd13 + 1;
   
    %Only add embryos' data to the combined variables if at least one 
    %embryo exists in this AP bin
    if combinedStart13 ~= 0
        %Combine APDivision
        APDivision(13, APBin) = combinedStart13;
        APDivision(14, APBin) = combinedStart14;
        
        %Align and combine the embryos' data for this AP bin
        for i = 1:numEmbryos
            begin13 = combinedData(i).APDivision(13, APBin);
            end13 = combinedData(i).APDivision(14, APBin) - 1;
            length13 = end13 - begin13;
            
            %Only bother adding the embryo to the combined variables if the
            %embryo exists in this AP bin
            if begin13 ~= 0
                %Combine NParticlesAP for all embryos
                newParticles = combinedData(i).NParticlesAP(begin13:end13,APBin);
                NParticlesAP(combinedStart13:(combinedStart13 + length13),APBin) = ...
                    NParticlesAP(combinedStart13:(combinedStart13 + length13),APBin) ...
                    + newParticles;

                %Combine MeanVectorAP for all embryos
                newMeans = combinedData(i).MeanVectorAP(begin13:end13,APBin);
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
                newSDs = combinedData(i).SDVectorAP(begin13:end13,APBin);
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
for i = 1:numEmbryos
    newRatio = allData(i).OnRatioLineageAP;     %Store embryo's ratio data
    newElements = ~isnan(newRatio);             %Store which elements 
                                                 %conribute to combined ratio
    newRatio(isnan(newRatio)) = 0;              %Change NaN to zero to be
                                                 %able to add matrices
                                                 
    OnRatioLineageAP = OnRatioLineageAP + newRatio;
    contributingElements = contributingElements + newElements;
end

OnRatioLineageAP = OnRatioLineageAP./contributingElements;
%% DEFINE NEW, COMBINED START FRAMES FOR nc12, nc13, AND n14

%As written right now, nc12 isn't analyzed with this code so I'm just
%saying that the start of nc12 doesn't exist
nc12 = 0;

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
save([DropboxFolder,filesep,[DataType,'_Combined_CompiledParticles.mat']],...
    'nc12','nc13', 'nc14','NParticlesAP', 'MeanVectorAP', 'SDVectorAP', ...
    'ElapsedTime','APbinID', 'OnRatioLineageAP');
display('Combined_CompiledParticles saved');
save([DropboxFolder,filesep,[DataType,'_Combined_APDivision.mat']],'APDivision');
display('Combined_APDivision saved');

%Final line of code at which a breakpoint can be placed for debugging
%purposes ('A' is not an actual variable with any important meaning). 
%Allows you to view and manipulate the workspace created by this function
%after running the whole thing. Shouldn't do anything, but you can delete
%if you don't need to debug.
A = 0;
