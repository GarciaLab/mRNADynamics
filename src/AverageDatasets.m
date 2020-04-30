%% Averaging multiple datasets
function AverageDatasets(DataType,varargin)
% Author : Yang Joon Kim (yjkim90@berkeley.edu)
% This is a script modified from Meghan's CombineMultipleEmbryos.m script
% Last Updated : 7/17/2019

% DESCRIPTION
% This function has input of "DataType" a tab name in the DataStatus.xls, 
% grabs all datasets in that tab,and calculates 
% 1) Averaged MS2 spot fluorescence (weighted sum) ,Standard
% Deviation, and the total number of MS2 spots from multiple embryos in
% nc12, nc13 and nc14.

% Note. This is assuming that you're interested in nc12, nc13 and nc14.
% If you don't have the whole nc12, you might need to comment that dataset
% out in the DataStatus.xlsx
%
% OPTIONS
% 'NC', N : designate the nuclear cycle to begin avearging, this can be 12,13 or
% 14.
% 'NCadjust' : Fine tune the NC registration

% PARAMETERS
% DataType: This is a string that is identical to the name of the tab in
% dataStatus.xlsx that you wish to analyze.
%
% OUTPUT
% Variables for plotting, or more detailed analysis with the Averaged spot
% fluorescence over time. Save as 'Name of the DataType'.mat file
% (nc12, nc13, nc14, NParticlesAP,MeanVectorsAP, SDVectorAP, ElapsedTime) 
% corresponding to the embryos combined

% EXAMPLE
%AverageDatasets(DataType,'NC',13,'savePath',FilePath);

[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;

Data = LoadMS2Sets(DataType,'dontCompare');

%% Define the options
% Save path option
savePath = '';

for i=1:length(varargin)
    if strcmpi(varargin{i}, 'savePath')
        savePath = varargin{i+1};
    end
end

% Define the NC
NC = 13; % Default
for i = 1:length(varargin)
   if strcmpi(varargin{i}, 'NC')
        NC = varargin{i+1};  
   end
end

% instantaneous fraction on
% whether to use the Number of nuclei or the area of APbin
% The idea is that when this option is ON, we will use the area of APbins
% to calculate the fraction on, rather than the number of nuclei. This is
% for the case where there's no nuclear marker.
SpotDensity = false;
for i = 1:length(varargin)
   if strcmpi(varargin{i}, 'SpotDensity')
        SpotDensity = true;
   end
end

%% Initialize the fields to calculate
% Number of datasets
numEmbryos=length(Data);

%Find the total number of frames for each embryo
numFrames = NaN(1, numEmbryos);
maxAPIndex = NaN(1, numEmbryos);
maxTime = NaN(1, numEmbryos);

% Determine the time frames depending on which NCs that we're interested in
for i = 1:numEmbryos
    if NC==13
        numFrames(i) = size(Data(i).ElapsedTime, 2);
        nc13(i) = Data(i).nc13;
        length_total_NC(i) = numFrames(i) - nc13(i)+1; % length of frames from nc13 to the last frame
        maxAPIndex(i) = 41; %Data(i).MaxAPIndex; % Usually 41, in 2.5% binning
        maxTime(i) = Data(i).ElapsedTime(numFrames(i));
    elseif NC==14
        numFrames(i) = size(Data(i).ElapsedTime, 2);
        nc14(i) = Data(i).nc14;
        length_total_NC(i) = numFrames(i) - nc14(i)+1; % length of frames from nc13 to the last frame
        maxAPIndex(i) = 41; %Data(i).MaxAPIndex; % Usually 41, in 2.5% binning
        maxTime(i) = Data(i).ElapsedTime(numFrames(i));
    elseif NC==12
        numFrames(i) = size(Data(i).ElapsedTime, 2);
        nc12(i) = Data(i).nc12;
        length_total_NC(i) = numFrames(i) - nc12(i)+1; % length of frames from nc13 to the last frame
        maxAPIndex(i) = Data(i).MaxAPIndex; % Usually 41, in 2.5% binning
        maxTime(i) = Data(i).ElapsedTime(numFrames(i));
    end
end

%Store the number of AP bins (this should always be 41).
numAPBins = maxAPIndex(1);

%% Synchronize the vectors as the beginning of the (nc 12), 13, and 14
% This could be edited to include even earlier nuclear cycles in the future.
% For now, nc12,nc13 and nc14 might be good enough.

% There are cases that the lengths of nuclear cycles were different (June,2018)
% To average multiple embryos properly, I will separate each nuclear cycle
% separately, since there's no transcription during the mitosis, it's fine
% to synchronize the datasets based on the beginning of the interphase
% (Use the APDivision.mat)

% Define the new ElapsedTime vector for the combined embryo. 
% The new ElapsedTime should start with the beginning of nc12, and also has
% the length of the frames of nc12 ~ nc14 (of the longest dataset for each APbin), all the
% empty values can be plugged with zeross.

% This ElapsedTime variable has evenly spaced time points estimated from diff(ElapsedTime)
% (This is assumption that we took the data with negligible time between serieses, which is pretty fair)

%% For all AP bins, define the new time vector
if NC==12
    Length12 = NaN(numEmbryos,numAPBins);
    Length13 = NaN(numEmbryos,numAPBins);
    Length14 = NaN(numEmbryos,numAPBins);
elseif NC==13
    Length13 = NaN(numEmbryos,numAPBins);
    Length14 = NaN(numEmbryos,numAPBins);
elseif NC==14
    Length14 = NaN(numEmbryos,numAPBins);
end

% For each AP bin, For all embryos, go through to find the longest nuclear cycle 
% (number of frames)
for j=1:numAPBins

    for i=1:numEmbryos
        % Define the nuclear cycle (In case we start from nc 12)
        if NC==12
            nc12(i,j) = Data(i).APDivision(12,j);
        end
            nc13(i,j) = Data(i).APDivision(13,j);
            nc14(i,j) = Data(i).APDivision(14,j);
            
        % Calculate the number of frames (during nc)
        if NC==12 && nc12(i,j)~=0
            Length12 (i,j) = nc13(i,j) - nc12(i,j);
            Length13 (i,j) = nc14(i,j) - nc13(i,j);
            Length14 (i,j) = numFrames(i) - nc14(i,j);
        elseif NC==13 && nc13(i,j)~=0
            Length13 (i,j) = nc14(i,j) - nc13(i,j);
            Length14 (i,j) = numFrames(i) - nc14(i,j);
        elseif NC==14 && nc14(i,j)~=0
            Length14 (i,j) = numFrames(i) - nc14(i,j);
        end
    end
    
    % Find the maximum length for each cycle
    if NC==12
        numFrames12(j) = max(Length12(:,j));
        numFrames13(j) = max(Length13(:,j));
        numFrames14(j) = max(Length14(:,j));
        TotalFrames(j) = numFrames12(j) + numFrames13(j) + numFrames14(j);
    elseif NC==13
        numFrames13(j) = max(Length13(:,j));
        numFrames14(j) = max(Length14(:,j));
        TotalFrames(j) = numFrames13(j) + numFrames14(j);
    elseif NC==14
        numFrames14(j) = max(Length14(:,j));
        TotalFrames(j) = numFrames14(j);
    end
end
 
% NewFrameLength = max(TotalFrames);

%% Define empty matrices (filled with zeros)
% These zeros that are left all the way to the end should be converted to
% NaNs.

if NC==12
    L12 = max(max(Length12));
elseif NC==13
    L12 = 0;
    L13 = max(max(Length13));
    L14 = max(max(Length14)); % Get the minimum for now, we can fix this better later.
elseif NC==14
    L12 = 0;
    L13 = 0;
    L14 = max(max(Length14)); % Get the minimum for now, we can fix this better later.
end


if NC==12
    MeanVectorAP_12 = zeros(L12,numAPBins,numEmbryos);
    SDVectorAP_12 = zeros(L12,numAPBins,numEmbryos);
    NParticlesAP_12 = zeros(L12,numAPBins,numEmbryos);
    FractionON_Instant_12 = zeros(L12,numAPBins,numEmbryos);
    
    MeanVectorAP_13 = zeros(L13,numAPBins,numEmbryos);
    SDVectorAP_13 = zeros(L13,numAPBins,numEmbryos);
    NParticlesAP_13 = zeros(L13,numAPBins,numEmbryos);
    FractionON_Instant_13 = zeros(L13,numAPBins,numEmbryos);

    MeanVectorAP_14 = zeros(L14,numAPBins,numEmbryos);
    SDVectorAP_14 = zeros(L14,numAPBins,numEmbryos);
    NParticlesAP_14 = zeros(L14,numAPBins,numEmbryos);
    FractionON_Instant_14 = zeros(L14,numAPBins,numEmbryos);
elseif NC==13
    MeanVectorAP_13 = zeros(L13,numAPBins,numEmbryos);
    SDVectorAP_13 = zeros(L13,numAPBins,numEmbryos);
    NParticlesAP_13 = zeros(L13,numAPBins,numEmbryos);
    FractionON_Instant_13 = zeros(L13,numAPBins,numEmbryos);

    MeanVectorAP_14 = zeros(L14,numAPBins,numEmbryos);
    SDVectorAP_14 = zeros(L14,numAPBins,numEmbryos);
    NParticlesAP_14 = zeros(L14,numAPBins,numEmbryos);
    FractionON_Instant_14 = zeros(L14,numAPBins,numEmbryos);
elseif NC==14
    MeanVectorAP_14 = zeros(L14,numAPBins,numEmbryos);
    SDVectorAP_14 = zeros(L14,numAPBins,numEmbryos);
    NParticlesAP_14 = zeros(L14,numAPBins,numEmbryos);
    FractionON_Instant_14 = zeros(L14,numAPBins,numEmbryos);
end


% Total matrices
NewFrameLength = L12+L13+L14;

MeanVectorAP = zeros(NewFrameLength,numAPBins,numEmbryos);
SDVectorAP = zeros(NewFrameLength,numAPBins,numEmbryos);
NParticlesAP = zeros(NewFrameLength,numAPBins,numEmbryos);
FractionON_Instant = zeros(NewFrameLength,numAPBins,numEmbryos);

% Synchornize all fields as all of them starts from nc 13 (or 12)
% This synchronization should be done for each AP bin, since they might
% have different anaphase time point.

for j=1:numAPBins       
    for i=1:numEmbryos
        % First, 
        Data(i).NEllipsesAP(Data(i).NEllipsesAP==0)=inf;
        MeanVectorAPTemp = cell2mat(Data(i).MeanVectorAP); %Data(i).MeanVectorAP;
        SDVectorAPTemp = cell2mat(Data(i).SDVectorAP); %Data(i).SDVectorAP;
        NParticlesAPTemp = cell2mat(Data(i).NParticlesAP); %Data(i).NParticlesAP;
        
        % In case there's no nuclear segmentation/tracking (this was the
        % case for the Runt nulls)
        if SpotDensity || isempty(Data(i).NEllipsesAP)
            % In case there's no histone marker.
            %warning('No Fraction ON data. Check whether there is His channel.')
            warning('Calculating the spot density, not using number of nuclei')
            APbinAreaTemp = Data(i).APbinArea;
            NEllipsesAPTemp = ones(size(MeanVectorAPTemp))./APbinAreaTemp;
        else
            NEllipsesAPTemp =Data(i).NEllipsesAP;  %Data(i).NEllipsesAP;
        end
        
        if NC==12 && nc12(i,j)~=0
            % Sync the fields for each nc
            % NC12
            MeanVectorAP_12(1:L12,j,i) = MeanVectorAPTemp(nc12(i,j):nc12(i,j)+L12-1,j);
            SDVectorAP_12(1:L12,j,i) = SDVectorAPTemp(nc12(i,j):nc12(i,j)+L12-1,j);
            NParticlesAP_12(1:L12,j,i) = NParticlesAPTemp(nc12(i,j):nc12(i,j)+L12-1,j);
            FractionON_Instant_12(1:L12,j,i) = NParticlesAPTemp(nc12(i,j):nc12(i,j)+L12-1,j)./...
                                                NEllipsesAPTemp(nc12(i,j):nc12(i,j)+L12-1,j);
            % NC13
            MeanVectorAP_13(1:L13,j,i) = MeanVectorAPTemp(nc13(i,j):nc13(i,j)+L13-1,j);
            SDVectorAP_13(1:L13,j,i) = SDVectorAPTemp(nc13(i,j):nc13(i,j)+L13-1,j);
            NParticlesAP_13(1:L13,j,i) = NParticlesAPTemp(nc13(i,j):nc13(i,j)+L13-1,j);
            FractionON_Instant_13(1:L13,j,i) = NParticlesAPTEmp(nc13(i,j):nc13(i,j)+L13-1,j)./...
                                                NEllipsesAPTemp(nc13(i,j):nc13(i,j)+L13-1,j);
            % NC14                          
            MeanVectorAP_14(1:numFrames(i)-nc14(i,j),j,i) = MeanVectorAPTemp(nc14(i,j):numFrames(i)-1,j);
            SDVectorAP_14(1:numFrames(i)-nc14(i,j),j,i) = SDVectorAPTemp(nc14(i,j):numFrames(i)-1,j);
            NParticlesAP_14(1:numFrames(i)-nc14(i,j),j,i) = NParticlesAPTemp(nc14(i,j):numFrames(i)-1,j);
            FractionON_Instant_14(1:numFrames(i)-nc14(i,j),j,i) = NParticlesAP(nc14(i,j):numFrames(i)-1,j)./...
                                                NEllipsesAP(nc14(i,j):numFrames(i)-1,j);

        elseif NC==13 && nc13(i,j)~=0 && nc14(i,j)~=0
%             MeanVectorAP(1:numFrames(i)-nc13(i)+1,:,i) = Data(i).MeanVectorAP(nc13(i):numFrames(i),:);
%             SDVectorAP(1:numFrames(i)-nc13(i)+1,:,i) = Data(i).SDVectorAP(nc13(i):numFrames(i),:);
%             NParticlesAP(1:numFrames(i)-nc13(i)+1,:,i) = Data(i).NParticlesAP(nc13(i):numFrames(i),:);
%             FractionOn_Instant(1:numFrames(i)-nc13(i)+1,:,i) = Data(i).NParticlesAP(nc13(i):numFrames(i),:)./Data(i).NEllipsesAP(nc13(i):numFrames(i),:);

            % NC13
            MeanVectorAP_13(1:L13,j,i) = MeanVectorAPTemp(nc13(i,j):nc13(i,j)+L13-1,j);
            SDVectorAP_13(1:L13,j,i) = SDVectorAPTemp(nc13(i,j):nc13(i,j)+L13-1,j);
            NParticlesAP_13(1:L13,j,i) = NParticlesAPTemp(nc13(i,j):nc13(i,j)+L13-1,j);
            FractionON_Instant_13(1:L13,j,i) = NParticlesAPTemp(nc13(i,j):nc13(i,j)+L13-1,j)./...
                                               NEllipsesAPTemp(nc13(i,j):nc13(i,j)+L13-1,j);
            % NC14                          
            MeanVectorAP_14(1:numFrames(i)-nc14(i,j),j,i) = MeanVectorAPTemp(nc14(i,j):numFrames(i)-1,j);
            SDVectorAP_14(1:numFrames(i)-nc14(i,j),j,i) = SDVectorAPTemp(nc14(i,j):numFrames(i)-1,j);
            NParticlesAP_14(1:numFrames(i)-nc14(i,j),j,i) = NParticlesAPTemp(nc14(i,j):numFrames(i)-1,j);
            FractionON_Instant_14(1:numFrames(i)-nc14(i,j),j,i) = NParticlesAPTemp(nc14(i,j):numFrames(i)-1,j)./...
                                               NEllipsesAPTemp(nc14(i,j):numFrames(i)-1,j);
        elseif NC==14 && nc14(i,j)~=0
            % NC14                          
            MeanVectorAP_14(1:numFrames(i)-nc14(i,j),j,i) = MeanVectorAPTemp(nc14(i,j):numFrames(i)-1,j);
            SDVectorAP_14(1:numFrames(i)-nc14(i,j),j,i) = SDVectorAPTemp(nc14(i,j):numFrames(i)-1,j);
            NParticlesAP_14(1:numFrames(i)-nc14(i,j),j,i) = NParticlesAPTemp(nc14(i,j):numFrames(i)-1,j);
            FractionON_Instant_14(1:numFrames(i)-nc14(i,j),j,i) = NParticlesAPTemp(nc14(i,j):numFrames(i)-1,j)./...
                                               NEllipsesAPTemp(nc14(i,j):numFrames(i)-1,j);
        end
    end
end

% % Take the most frequent value of dT from the ElapsedTime. It's because the
% % dT can be different in case we stopped and restarted the movie.
deltaT = mode(diff(Data(1).ElapsedTime)); 
ElapsedTime = deltaT*(0:NewFrameLength-1);

%% Concatenate the MeanVectorAP_12,13,and 14 into MeanVectorAP 
% (same for SD, NParticles, and FractionON_Instant
if NC==12
    MeanVectorAP = cat(1,MeanVectorAP_12,MeanVectorAP_13,MeanVectorAP_14);
    SDVectorAP = cat(1,SDVectorAP_12,SDVectorAP_13,SDVectorAP_14);
    NParticlesAP = cat(1,NParticlesAP_12,NParticlesAP_13,NParticlesAP_14);
    FractionON = cat(1,FractionON_Instant_12,FractionON_Instant_13,FractionON_Instant_14);
    %MeanVectorAP_ForSum = cat(1,MeanVectorAP_12*0.25,MeanVectorAP_13*0.5,MeanVectorAP_14);
elseif NC==13
    MeanVectorAP = cat(1,MeanVectorAP_13,MeanVectorAP_14);
    SDVectorAP = cat(1,SDVectorAP_13,SDVectorAP_14);
    NParticlesAP = cat(1,NParticlesAP_13,NParticlesAP_14);
    FractionON = cat(1,FractionON_Instant_13,FractionON_Instant_14);
    %MeanVectorAP_ForSum = cat(1,MeanVectorAP_13*0.5,MeanVectorAP_14);
elseif NC==14
    MeanVectorAP = MeanVectorAP_14;
    SDVectorAP = SDVectorAP_14;
    NParticlesAP = NParticlesAP_14;
    FractionON = FractionON_Instant_14;
    %MeanVectorAP_ForSum = cat(1,MeanVectorAP_13*0.5,MeanVectorAP_14);
else 
    warning('This part is left as an option. You can designate earlier cycles by editing this code.')

end

%% Plot to check before the averaging (Save the MeanVectorAP, etc. from individual embryos for future plots)
% AP = 13; % You can change this
% hold on
% for i=1:numEmbryos
%     errorbar(ElapsedTime,MeanVectorAP(:,AP,i),SDVectorAP(:,AP,i))
% end
% 
% errorbar(ElapsedTime, MeanVectorAP_combined(:,AP), SDVectorAP(:,AP))
%% Define synchronized fields from individual embryos
MeanVectorAP_individual = MeanVectorAP;
SDVectorAP_individual = SDVectorAP;
NParticlesAP_individual = NParticlesAP;
FractionON_individual = FractionON;
%MeanVectorAP_ForSum_individual = MeanVectorAP_ForSum;

%% Convert NaNs to zeros
% This is for a convenient averaging
MeanVectorAP(isnan(MeanVectorAP)) = 0;
SDVectorAP(isnan(SDVectorAP)) = 0;
NParticlesAP(isnan(NParticlesAP)) = 0;
FractionON(isnan(FractionON)) = 0;

%% Option
% For r3, I have negative MeanVectorAP for some APbins...
% I'll assign zeros for negative MeanVectorAPs
MeanVectorAP(MeanVectorAP<0) = 0;
NParticlesAP(MeanVectorAP<0) = 0;

%% Averaging - weighted sum
sumMean = zeros(NewFrameLength,numAPBins);
sumSD = zeros(NewFrameLength,numAPBins);
sumNParticles = zeros(NewFrameLength,numAPBins);
%sumMean_ForSum = zeros(NewFrameLength,numAPBins);

for i=1:numEmbryos
    sumMean = sumMean + squeeze(MeanVectorAP(:,:,i).*NParticlesAP(:,:,i));
    sumSD = sumSD + squeeze(SDVectorAP(:,:,i).^2.*NParticlesAP(:,:,i));
    sumNParticles = sumNParticles + squeeze(NParticlesAP(:,:,i));
    % MeanVectorAP considering the doubling of nuclei in each cycle
    %sumMean_ForSum = sumMean_ForSum + squeeze(MeanVectorAP_ForSum(:,:,i).*NParticlesAP(:,:,i));
end
    
MeanVectorAP_combined = sumMean./sumNParticles;
SDVectorAP_combined = sqrt(sumSD./sumNParticles);
NParticlesAP_combined = sumNParticles;
%MeanVectorAP_ForSum_combined = sumMean_ForSum./sumNParticles;

% For FractionON, we will use nanmean, to ignore the APbin(of embryos) that
% doesn't have any values.
FractionON_temp = FractionON;
FractionON_temp(FractionON_temp==0) = nan;
FractionON_combined = nanmean(FractionON_temp,3);
FractionON_combined(isnan(FractionON_combined)) = 0;

%% Plot averaged MS2 traces along with individual embryos
% EmbryosLegend = {'Embryo1','Embryo2','Embryo3','Embryo4','Embryo5'};
% EmbryosLegend = EmbryosLegend(1:numEmbryos)
% for AP = 1:41
%     APbin = (AP-1)*2.5;
%     clf
%     hold on
%     for i=1:numEmbryos
%         errorbar(ElapsedTime,MeanVectorAP_forPlot(:,AP,i),SDVectorAP_forPlot(:,AP,i))
%     end
% 
%     errorbar(ElapsedTime,MeanVectorAP(:,AP),SEVectorAP(:,AP))
% 
%     title({'Averaged MS2 spot fluorescence';['@ APbin = ',num2str(APbin),'%']})
%     xlabel('Time (min)')
%     ylabel('MS2 Spot Fluorescence (AU)')
%     legend(EmbryosLegend,'Mean')
%     yMax = max(max(MeanVectorAP_forPlot(:,AP,:))) + max(max(SDVectorAP_forPlot(:,AP,:)));
%     if yMax==0
%         yMax=1;
%     end
%     ylim([0 yMax])
%     
%     standardizeFigure_YJK(gca,legend,'yellow','cyan','magenta','lightblue','red')
%     saveas(gcf,['E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\Figures-OpposingGradients\hbP2-r0123-AveragedMS2Traces\',...
%             DataType,'AveragedTraces_IndividualEmbryos-AP=',num2str(APbin),'%'],'tif');
% end

%% Define the fields that needs to be saved
%% Nuclear cycle (all AP bins are synchronized)
if NC==12
    nc12 = 1;
    nc13 = nc12 + L12;
    nc14 = nc13 + L13;
elseif NC==13
    nc12 = 0;
    nc13 = 1;
    nc14 = nc13 + L13;
elseif NC==14
    nc12 = 0;
    nc13 = 0;
    nc14 = 1;
end
%% Define the fields to be saved
MeanVectorAP = MeanVectorAP_combined;
SDVectorAP = SDVectorAP_combined;
NParticlesAP = NParticlesAP_combined;
SEVectorAP = SDVectorAP./sqrt(NParticlesAP); % Standard error of mean (SD / sqrt(number of observation, here, particles)
ElapsedTime = ElapsedTime;
FractionON = FractionON_combined;

%% Convert zeros to Nans
% Convert the zeros in APbins to Nans, since those were not measured.
MeanVectorAP(MeanVectorAP ==0) = nan;
SDVectorAP(SDVectorAP ==0) = nan;
SEVectorAP(SEVectorAP ==0) = nan; % Standard error of mean (SD / sqrt(number of observation)
NParticlesAP(NParticlesAP==0) = nan;

% fields from individual embryos
MeanVectorAP_individual(MeanVectorAP_individual ==0) = nan;
SDVectorAP_individual(SDVectorAP_individual ==0) = nan;
NParticlesAP_individual(NParticlesAP_individual ==0) = nan;


%% Save the fields in .mat file
    if ~isempty(savePath)
        save([savePath,filesep,DataType,'.mat'],...
            'MeanVectorAP','SDVectorAP','SEVectorAP','NParticlesAP','FractionON',...
            'ElapsedTime','nc12', 'nc13', 'nc14',...
            'MeanVectorAP_individual','SDVectorAP_individual',...
            'NParticlesAP_individual','FractionON_individual')

    else
        warning('Define the File Path in the argument above')
    end
end