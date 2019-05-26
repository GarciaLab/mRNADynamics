%% Averaging multiple datasets
function AverageDatasets(DataType,varargin)
% Author : Yang Joon Kim (yjkim90@berkeley.edu)
% This is edited from Meghan's CombineMultipleEmbryos.m script

% DESCRIPTION
% This function has input of datatype in DataStatus.xls, grabs all datasets
% in that tab,and calculates 
% 1) Averaged MS2 spot fluorescence (weighted sum) ,Standard
% Deviation, and the total number of MS2 spots from multiple embryos in
% nc12, nc13 and nc14.
% 2) Fraction ON : This can be calculated for each embryo, then averaged in
% multiple ways.
% 3) Accumulated mRNA (Accumulated fluorescence over nc13 and nc14) per
% nucleus, since it's using MeanVectorAP
% 4) Accumulated mRNA (considering Fraction ON) :
% This will be more appropriate for deciding the boundary profiles, such as
% boundary position, sharpness, and width.

% Note. This is assuming that you're interested in nc12, nc13 and nc14.
% If you don't have the whole nc12, you might need to comment that dataset
% out in the DataStatus.xlsx
%
% OPTIONS
% 'NC', N : designate the nuclear cycle to begin avearging, this can be 12,13 or
% 14.

% PARAMETERS
% DataType: This is a string that is identical to the name of the tab in
% dataStatus.xlsx that you wish to analyze.
%
% OUTPUT
% Variables for plotting, or more detailed analysis with the Averaged spot
% fluorescence over time. Save as 'Name of the DataType'.mat file
% (nc12, nc13, nc14, NParticlesAP,MeanVectorsAP, SDVectorAP, ElapsedTime,
%  AccumulatedmRNA, FractionON, AccumulatedmRNA_FractionON  ) 
% corresponding to the embryos combined

[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;

if ischar(DataType)
    Data = LoadMS2Sets(DataType);
else
    Data = DataType;
end

% Save path option
savePath = '';

for embryo=1:length(varargin)
    if strcmpi(varargin{embryo}, 'savePath')
        savePath = varargin{embryo+1};
    end
end

NC = 13; % Default
for embryo = 1:length(varargin)
   if strcmpi(varargin{embryo}, 'NC')
        NC = varargin{embryo+1};  
   end
end

numEmbryos=length(Data);

%Find the total number of frames for each embryo
numFrames = zeros(1, numEmbryos);
maxAPIndex = ones(1, numEmbryos)*41; %initialize with the value 41 just in case this isn't set by compileparticles
maxTime = zeros(1, numEmbryos);
for embryo = 1:numEmbryos
    if NC~=12
        numFrames(embryo) = size(Data(embryo).ElapsedTime, 2);
        nc13(embryo) = Data(embryo).nc13;
        length_total_NC(embryo) = numFrames(embryo) - nc13(embryo)+1; % length of frames from nc13 to the last frame
        if ~isempty(Data(embryo).MaxAPIndex)
            maxAPIndex(embryo) = Data(embryo).MaxAPIndex; % Usually 41, in 2.5% binning
        end
        maxTime(embryo) = Data(embryo).ElapsedTime(numFrames(embryo));
    elseif NC==12
        numFrames(embryo) = size(Data(embryo).ElapsedTime, 2);
        nc12(embryo) = Data(embryo).nc12;
        length_total_NC(embryo) = numFrames(embryo) - nc12(embryo)+1; % length of frames from nc13 to the last frame
        if ~isempty(Data(embryo).MaxAPIndex)
            maxAPIndex(embryo) = Data(embryo).MaxAPIndex; % Usually 41, in 2.5% binning
        end
        maxTime(embryo) = Data(embryo).ElapsedTime(numFrames(embryo));
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
% empty values can be plugged with Nans.

% This ElapsedTime variable has evenly spaced time points estimated from diff(ElapsedTime)
% (This is assumption that we took the data with negligible time between serieses, which is pretty fair)

%% For all AP bins, define the new time vector
if NC==12
    Length12 = zeros(numEmbryos,numAPBins);
end
Length13 = zeros(numEmbryos,numAPBins);
Length14 = zeros(numEmbryos,numAPBins);

for APBin=1:numAPBins
    % For all embryos, go through to find the longest nuclear cycle (number
    % of frames)
    for embryo=1:numEmbryos
        % Define the nuclear cycle (In case we start from nc 12)
        if NC==12
            nc12(embryo,APBin) = Data(embryo).APDivision(12,APBin);
        end
            nc13(embryo,APBin) = Data(embryo).APDivision(13,APBin);
            nc14(embryo,APBin) = Data(embryo).APDivision(14,APBin);
            
        % Calculate the number of frames (during nc)
        if NC==12 && nc12(embryo,APBin)~=0
            Length12 (embryo,APBin) = nc13(embryo,APBin) - nc12(embryo,APBin);
            Length13 (embryo,APBin) = nc14(embryo,APBin) - nc13(embryo,APBin);
            Length14 (embryo,APBin) = numFrames(embryo) - nc14(embryo,APBin);
        elseif NC==13 && nc13(embryo,APBin)~=0
            Length13 (embryo,APBin) = nc14(embryo,APBin) - nc13(embryo,APBin);
            Length14 (embryo,APBin) = numFrames(embryo) - nc14(embryo,APBin);
        end
    end
    % Find the maximum length for each cycle
        numFrames13(APBin) = max(Length13(:,APBin));
        numFrames14(APBin) = max(Length14(:,APBin));
        TotalFrames(APBin) = numFrames13(APBin) + numFrames14(APBin);
    if NC==12
        numFrames12(APBin) = max(Length12(:,APBin));
        TotalFrames(APBin) = numFrames12(APBin) + numFrames13(APBin) + numFrames14(APBin);
    end
   
end
 
% NewFrameLength = max(TotalFrames);

%% Define empty matrices (filled with Nans)

if NC==12
    L12 = max(max(Length12));
else
    L12 = 0;
end
L13 = max(max(Length13));
L14 = max(max(Length14)); % Get the minimum for now, we can fix this better later.

if NC==12
    MeanVectorAP_12 = NaN(L12,numAPBins,numEmbryos);
    SDVectorAP_12 = NaN(L12,numAPBins,numEmbryos);
    NParticlesAP_12 = NaN(L12,numAPBins,numEmbryos);
    FractionON_Instant_12 = NaN(L12,numAPBins,numEmbryos);
end
MeanVectorAP_13 = NaN(L13,numAPBins,numEmbryos);
SDVectorAP_13 = NaN(L13,numAPBins,numEmbryos);
NParticlesAP_13 = NaN(L13,numAPBins,numEmbryos);
FractionON_Instant_13 = NaN(L13,numAPBins,numEmbryos);

MeanVectorAP_14 = NaN(L14,numAPBins,numEmbryos);
SDVectorAP_14 = NaN(L14,numAPBins,numEmbryos);
NParticlesAP_14 = NaN(L14,numAPBins,numEmbryos);
FractionON_Instant_14 = NaN(L14,numAPBins,numEmbryos);

% Total matrices
NewFrameLength = L12+L13+L14;

MeanVectorAP = NaN(NewFrameLength,numAPBins,numEmbryos);
SDVectorAP = NaN(NewFrameLength,numAPBins,numEmbryos);
NParticlesAP = NaN(NewFrameLength,numAPBins,numEmbryos);
FractionON_Instant = NaN(NewFrameLength,numAPBins,numEmbryos);

% Synchornize all fields as all of them starts from nc 13 (or 12)
% This synchronization should be done for each AP bin, since they might
% have different anaphase time point.

for APBin=1:numAPBins       
    for embryo=1:numEmbryos
        % First, 
        Data(embryo).NEllipsesAP(Data(embryo).NEllipsesAP==0)=inf;
        %assuming just one spot channel
        if iscell(Data(embryo).MeanVectorAP)
            Data(embryo).MeanVectorAP = Data(embryo).MeanVectorAP{1};
            Data(embryo).SDVectorAP = Data(embryo).SDVectorAP{1};
            Data(embryo).NParticlesAP = Data(embryo).NParticlesAP{1};
        end
        if NC==12 && nc12(embryo,APBin)~=0
            % Sync the fields for each nc
            % NC12
            MeanVectorAP_12(1:L12,APBin,embryo) = Data(embryo).MeanVectorAP(nc12(embryo,APBin):nc12(embryo,APBin)+L12-1,APBin);
            SDVectorAP_12(1:L12,APBin,embryo) = Data(embryo).SDVectorAP(nc12(embryo,APBin):nc12(embryo,APBin)+L12-1,APBin);
            NParticlesAP_12(1:L12,APBin,embryo) = Data(embryo).NParticlesAP(nc12(embryo,APBin):nc12(embryo,APBin)+L12-1,APBin);
            FractionON_Instant_12(1:L12,APBin,embryo) = Data(embryo).NParticlesAP(nc12(embryo,APBin):nc12(embryo,APBin)+L12-1,APBin)./...
                                                Data(embryo).NEllipsesAP(nc12(embryo,APBin):nc12(embryo,APBin)+L12-1,APBin);
            % NC13
            MeanVectorAP_13(1:L13,APBin,embryo) = Data(embryo).MeanVectorAP(nc13(embryo,APBin):nc13(embryo,APBin)+L13-1,APBin);
            SDVectorAP_13(1:L13,APBin,embryo) = Data(embryo).SDVectorAP(nc13(embryo,APBin):nc13(embryo,APBin)+L13-1,APBin);
            NParticlesAP_13(1:L13,APBin,embryo) = Data(embryo).NParticlesAP(nc13(embryo,APBin):nc13(embryo,APBin)+L13-1,APBin);
            FractionON_Instant_13(1:L13,APBin,embryo) = Data(embryo).NParticlesAP(nc13(embryo,APBin):nc13(embryo,APBin)+L13-1,APBin)./...
                                                Data(embryo).NEllipsesAP(nc13(embryo,APBin):nc13(embryo,APBin)+L13-1,APBin);
            % NC14                          
            MeanVectorAP_14(1:numFrames(embryo)-nc14(embryo,APBin),APBin,embryo) = Data(embryo).MeanVectorAP(nc14(embryo,APBin):numFrames(embryo)-1,APBin);
            SDVectorAP_14(1:numFrames(embryo)-nc14(embryo,APBin),APBin,embryo) = Data(embryo).SDVectorAP(nc14(embryo,APBin):numFrames(embryo)-1,APBin);
            NParticlesAP_14(1:numFrames(embryo)-nc14(embryo,APBin),APBin,embryo) = Data(embryo).NParticlesAP(nc14(embryo,APBin):numFrames(embryo)-1,APBin);
            FractionON_Instant_14(1:numFrames(embryo)-nc14(embryo,APBin),APBin,embryo) = Data(embryo).NParticlesAP(nc14(embryo,APBin):numFrames(embryo)-1,APBin)./...
                                                Data(embryo).NEllipsesAP(nc14(embryo,APBin):numFrames(embryo)-1,APBin);
%             SDVectorAP(1:numFrames(i)-nc12(i)+1,:,i) = Data(i).SDVectorAP(nc12(i):numFrames(i),:);
%             NParticlesAP(1:numFrames(i)-nc12(i)+1,:,i) = Data(i).NParticlesAP(nc12(i):numFrames(i),:);
%             FractionOn_Instant(1:numFrames(i)-nc13(i)+1,:,i) = Data(i).NParticlesAP(nc13(i):numFrames(i),:)./Data(i).NEllipsesAP(nc13(i):numFrames(i),:)

        %elseif nc13(i)==0
        %    error('Check the Movie if it really does not start from nc13, then you should edit this code or make that dataset as an exception')
        elseif nc13(embryo,APBin)~=0 && nc14(embryo,APBin)~=0
%             MeanVectorAP(1:numFrames(i)-nc13(i)+1,:,i) = Data(i).MeanVectorAP(nc13(i):numFrames(i),:);
%             SDVectorAP(1:numFrames(i)-nc13(i)+1,:,i) = Data(i).SDVectorAP(nc13(i):numFrames(i),:);
%             NParticlesAP(1:numFrames(i)-nc13(i)+1,:,i) = Data(i).NParticlesAP(nc13(i):numFrames(i),:);
%             FractionOn_Instant(1:numFrames(i)-nc13(i)+1,:,i) = Data(i).NParticlesAP(nc13(i):numFrames(i),:)./Data(i).NEllipsesAP(nc13(i):numFrames(i),:);

            % NC13
            MeanVectorAP_13(1:L13,APBin,embryo) = Data(embryo).MeanVectorAP(nc13(embryo,APBin):nc13(embryo,APBin)+L13-1,APBin);
            SDVectorAP_13(1:L13,APBin,embryo) = Data(embryo).SDVectorAP(nc13(embryo,APBin):nc13(embryo,APBin)+L13-1,APBin);
            NParticlesAP_13(1:L13,APBin,embryo) = Data(embryo).NParticlesAP(nc13(embryo,APBin):nc13(embryo,APBin)+L13-1,APBin);
            FractionON_Instant_13(1:L13,APBin,embryo) = Data(embryo).NParticlesAP(nc13(embryo,APBin):nc13(embryo,APBin)+L13-1,APBin)./...
                                                Data(embryo).NEllipsesAP(nc13(embryo,APBin):nc13(embryo,APBin)+L13-1,APBin);
            % NC14                          
            MeanVectorAP_14(1:numFrames(embryo)-nc14(embryo,APBin),APBin,embryo) = Data(embryo).MeanVectorAP(nc14(embryo,APBin):numFrames(embryo)-1,APBin);
            SDVectorAP_14(1:numFrames(embryo)-nc14(embryo,APBin),APBin,embryo) = Data(embryo).SDVectorAP(nc14(embryo,APBin):numFrames(embryo)-1,APBin);
            NParticlesAP_14(1:numFrames(embryo)-nc14(embryo,APBin),APBin,embryo) = Data(embryo).NParticlesAP(nc14(embryo,APBin):numFrames(embryo)-1,APBin);
            FractionON_Instant_14(1:numFrames(embryo)-nc14(embryo,APBin),APBin,embryo) = Data(embryo).NParticlesAP(nc14(embryo,APBin):numFrames(embryo)-1,APBin)./...
                                                Data(embryo).NEllipsesAP(nc14(embryo,APBin):numFrames(embryo)-1,APBin);
        end
    end
end

% % Take the most frequent value of dT from the ElapsedTime. It's because the
% % dT can be different in case we stopped and restarted the movie.
deltaT = mode(diff(Data(1).ElapsedTime)); 
ElapsedTime = deltaT*(0:NewFrameLength-1);

%% Average all fields at each time point
%% Make Nans as zeros
% MeanVectorAP(isnan(MeanVectorAP)) = 0;
% SDVectorAP(isnan(SDVectorAP)) = 0;
% NParticlesAP(isnan(NParticlesAP)) = 0;

if NC==12
    MeanVectorAP_12(isnan(MeanVectorAP_12)) = 0;
    SDVectorAP_12(isnan(SDVectorAP_12)) = 0;
    NParticlesAP_12(isnan(NParticlesAP_12)) = 0;
    %FractionON_Instant_12(isnan(FractionON_Instant_12)) = 0;

end

MeanVectorAP_13(isnan(MeanVectorAP_13)) = 0;
SDVectorAP_13(isnan(SDVectorAP_13)) = 0;
NParticlesAP_13(isnan(NParticlesAP_13)) = 0;
%FractionON_Instant_13(isnan(FractionON_Instant_13)) = 0;

MeanVectorAP_14(isnan(MeanVectorAP_14)) = 0;
SDVectorAP_14(isnan(SDVectorAP_14)) = 0;
NParticlesAP_14(isnan(NParticlesAP_14)) = 0;
%FractionON_Instant_14(isnan(FractionON_Instant_14)) = 0;

%% Concatenate the MeanVectorAP_12,13,and 14 into MeanVectorAP 
% (same for SD, NParticles, and FractionON_Instant
if NC==12
    MeanVectorAP = cat(1,MeanVectorAP_12,MeanVectorAP_13,MeanVectorAP_14);
    SDVectorAP = cat(1,SDVectorAP_12,SDVectorAP_13,SDVectorAP_14);
    NParticlesAP = cat(1,NParticlesAP_12,NParticlesAP_13,NParticlesAP_14);
    FractionON = cat(1,FractionON_Instant_12,FractionON_Instant_13,FractionON_Instant_14);
    MeanVectorAP_ForSum = cat(1,MeanVectorAP_12*0.25,MeanVectorAP_13*0.5,MeanVectorAP_14);
elseif NC==13
    MeanVectorAP = cat(1,MeanVectorAP_13,MeanVectorAP_14);
    SDVectorAP = cat(1,SDVectorAP_13,SDVectorAP_14);
    NParticlesAP = cat(1,NParticlesAP_13,NParticlesAP_14);
    FractionON = cat(1,FractionON_Instant_13,FractionON_Instant_14);
    MeanVectorAP_ForSum = cat(1,MeanVectorAP_13*0.5,MeanVectorAP_14);
else 
    warning('This part is left as an option. You can designate earlier cycles by editing this code.')

end

%% Plot to check before the averaging (Save the MeanVectorAP, etc. from individual embryos for future plots)
%AP = 10; % You can change this
%hold on
% for i=1:numEmbryos
%     errorbar(ElapsedTime,MeanVectorAP(:,AP,i),SDVectorAP(:,AP,i))
% end
MeanVectorAP_forPlot = MeanVectorAP;
SDVectorAP_forPlot = SDVectorAP;
NParticlesAP_forPlot = NParticlesAP;

%% Averaging
sumMean = zeros(NewFrameLength,numAPBins);
sumSD = zeros(NewFrameLength,numAPBins);
sumNParticles = zeros(NewFrameLength,numAPBins);
sumMean_ForSum = zeros(NewFrameLength,numAPBins);

for embryo=1:numEmbryos
    sumMean = sumMean + squeeze(MeanVectorAP(:,:,embryo).*NParticlesAP(:,:,embryo));
    sumSD = sumSD + squeeze(SDVectorAP(:,:,embryo).^2.*NParticlesAP(:,:,embryo));
    sumNParticles = sumNParticles + squeeze(NParticlesAP(:,:,embryo));
    % MeanVectorAP considering the doubling of nuclei in each cycle
    sumMean_ForSum = sumMean_ForSum + squeeze(MeanVectorAP_ForSum(:,:,embryo).*NParticlesAP(:,:,embryo));
end
    
MeanVectorAPTemp = sumMean./sumNParticles;
SDVectorAPTemp = sqrt(sumSD./sumNParticles);
NParticlesAPTemp = sumNParticles;
MeanVectorAP_ForSumTemp = sumMean_ForSum./sumNParticles;
% For FractionON, we will use nanmean, to ignore the APbin(of embryos) that
% doesn't have any values.
FractionONTemp = nanmean(FractionON,3);

MeanVectorAP = MeanVectorAPTemp;
SDVectorAP = SDVectorAPTemp;
SEVectorAP = SDVectorAP/sqrt(numEmbryos); % Standard error of mean (SD / sqrt(number of observation)
NParticlesAP = NParticlesAPTemp;
ElapsedTime = ElapsedTime;
FractionON = FractionONTemp;
MeanVectorAP_ForSum = MeanVectorAP_ForSumTemp;

%% Accumulate mRNA over time (This can be made as an optional)
% I will calculate the Integrated mRNA from the MeanVectorAP
NFrames = length(ElapsedTime);
nAPbins = max(maxAPIndex);

AccumulatedmRNA = zeros(NFrames,nAPbins);
AccumulatedmRNA_SD =  zeros(NFrames,nAPbins);
% Consider the FractionON
% For MeanVectorAP, we can just multiply the FractonOn (instantaneous)
% But, for error estimation, we need to think hard. For example, we can use
% bootstrapping, or 
AccumulatedmRNA_FractionON = zeros(NFrames,nAPbins);
AccumulatedmRNA_FractionON_SD =  zeros(NFrames,nAPbins);

MeanVectorAP(isnan(MeanVectorAP))=0;
SDVectorAP(isnan(SDVectorAP))=0;
FractionON(isnan(FractionON))=0;
NParticlesAP(isnan(NParticlesAP))=0;
MeanVectorAP_ForSum(isnan(MeanVectorAP_ForSum))=0;

for embryo=1:maxAPIndex
    for APBin=2:length(ElapsedTime)
        AccumulatedmRNA(APBin,embryo) = trapz(ElapsedTime(1:APBin),MeanVectorAP(1:APBin,embryo));
        AccumulatedmRNA_SD(APBin,embryo) = sqrt(trapz(ElapsedTime(1:APBin),SDVectorAP(1:APBin,embryo).^2));
        
        AccumulatedmRNA_FractionON(APBin,embryo) = trapz(ElapsedTime(1:APBin),MeanVectorAP_ForSum(1:APBin,embryo).*FractionON(1:APBin,embryo));
        %AccumulatedmRNA_SD(j,i) = sqrt(trapz(ElapsedTime(1:j),SDVectorAP(1:j,i).^2));
    end
end

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
elseif NC==13
    nc12 = 0;
    nc13 = 1;
end
    %nc13 = nc12 + L12;
    nc14 = nc13 + L13;
    
%% Fraction ON 
%(as defined in Garcia, 2013, couting the fraction of nuclei 
% if they are ever turned on in one nuclear cycle)

%(1) Averaging the Fraction ON from each embryo.
FractionON_individual = zeros(41,3,numEmbryos);
for embryo=1:numEmbryos
    FractionON_individual(:,:,embryo) = Data(embryo).EllipsesOnAP./Data(embryo).TotalEllipsesAP;
end
FractionON_Average = nanmean(FractionON_individual,3);
FractionON_Average_Error = nanstd(FractionON_individual,[],3) / sqrt(numEmbryos);
%  %% Check by plotting (FractionON)
% NC = 12;
% hold on
% plot(0:0.025:1,FractionON_individual(:,NC-11,1))
% plot(0:0.025:1,FractionON_individual(:,NC-11,2))
% 
% errorbar(0:0.025:1,FractionON_Average(:,NC-11),FractionON_Average_Error(:,NC-11))
% 
% legend('embryo1','embryo2','Average')

%% (2) Summing up the number of ON nuclei over all embryos, then divide
%by the total number of nuclei
TotalONNuclei = zeros(41,3);
TotalNuclei = zeros(41,3);
    
for embryo=1:numEmbryos
    TotalONNuclei = TotalONNuclei + Data(embryo).EllipsesOnAP;
    TotalNuclei = TotalNuclei + Data(embryo).TotalEllipsesAP;
end

FractionONTemp = TotalONNuclei ./ TotalNuclei;
FractionON_Global = FractionONTemp;
%% Save the fields in .mat file
    if ~isempty(savePath)
        save([savePath,filesep,DataType,'_FromNC12.mat'],...
            'MeanVectorAP','SDVectorAP','SEVectorAP','NParticlesAP','ElapsedTime',...
            'FractionON', 'FractionON_Global',...
            'FractionON_Average','FractionON_Average_Error',...
            'nc12', 'nc13', 'nc14',...
            'AccumulatedmRNA','AccumulatedmRNA_SD', 'AccumulatedmRNA_FractionON')
    else
        warning('Define the File Path in the argument above')
    end
end