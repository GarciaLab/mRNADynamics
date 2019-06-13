%% Averaging multiple datasets
function AverageDatasets(DataType,varargin)
% Author : Yang Joon Kim (yjkim90@berkeley.edu)
% This is edited from Meghan's CombineMultipleEmbryos.m script
% Last Updated : 6/2/2019

% DESCRIPTION
% This function has input of datatype in DataStatus.xls, grabs all datasets
% in that tab,and calculates 
% 1) Averaged MS2 spot fluorescence (weighted sum) ,Standard
% Deviation, and the total number of MS2 spots from multiple embryos in
% nc12, nc13 and nc14.
% 2) Fraction ON : This can be calculated for each embryo, then averaged in
% multiple ways.
% 3) Accumulated mRNA (Accumulated fluorescence over nc13 and nc14) among ON nuclei,
% since it's using MeanVectorAP
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

Data = LoadMS2Sets(DataType);
%Prefix = DataType;
%Data = load(['E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\',Prefix,'\CompiledParticles.mat']);

% Save path option
savePath = '';

for i=1:length(varargin)
    if strcmpi(varargin{i}, 'savePath')
        savePath = varargin{i+1};
    end
end

NC = 13; % Default
for i = 1:length(varargin)
   if strcmpi(varargin{i}, 'NC')
        NC = varargin{i+1};  
   end
end

numEmbryos=length(Data);

%Find the total number of frames for each embryo
numFrames = NaN(1, numEmbryos);
maxAPIndex = NaN(1, numEmbryos);
maxTime = NaN(1, numEmbryos);

for i = 1:numEmbryos
    if NC~=12
        numFrames(i) = size(Data(i).ElapsedTime, 2);
        nc13(i) = Data(i).nc13;
        length_total_NC(i) = numFrames(i) - nc13(i)+1; % length of frames from nc13 to the last frame
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
end
Length13 = NaN(numEmbryos,numAPBins);
Length14 = NaN(numEmbryos,numAPBins);

for j=1:numAPBins
    % For all embryos, go through to find the longest nuclear cycle (number
    % of frames)
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
        end
    end
    % Find the maximum length for each cycle
        numFrames13(j) = max(Length13(:,j));
        numFrames14(j) = max(Length14(:,j));
        TotalFrames(j) = numFrames13(j) + numFrames14(j);
    if NC==12
        numFrames12(j) = max(Length12(:,j));
        TotalFrames(j) = numFrames12(j) + numFrames13(j) + numFrames14(j);
    end
   
end
 
% NewFrameLength = max(TotalFrames);

%% Define empty matrices (filled with zeros)
% These zeros that are left all the way to the end should be converted to
% zeross.

if NC==12
    L12 = max(max(Length12));
else
    L12 = 0;
end
L13 = max(max(Length13));
L14 = max(max(Length14)); % Get the minimum for now, we can fix this better later.

if NC==12
    MeanVectorAP_12 = zeros(L12,numAPBins,numEmbryos);
    SDVectorAP_12 = zeros(L12,numAPBins,numEmbryos);
    NParticlesAP_12 = zeros(L12,numAPBins,numEmbryos);
    FractionON_Instant_12 = zeros(L12,numAPBins,numEmbryos);
end
MeanVectorAP_13 = zeros(L13,numAPBins,numEmbryos);
SDVectorAP_13 = zeros(L13,numAPBins,numEmbryos);
NParticlesAP_13 = zeros(L13,numAPBins,numEmbryos);
FractionON_Instant_13 = zeros(L13,numAPBins,numEmbryos);

MeanVectorAP_14 = zeros(L14,numAPBins,numEmbryos);
SDVectorAP_14 = zeros(L14,numAPBins,numEmbryos);
NParticlesAP_14 = zeros(L14,numAPBins,numEmbryos);
FractionON_Instant_14 = zeros(L14,numAPBins,numEmbryos);

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
        NEllipsesAPTemp =Data(i).NEllipsesAP;  %Data(i).NEllipsesAP;
        
        if NC==12 && nc12(i,j)~=0
            % Sync the fields for each nc
            % NC12
            MeanVectorAP_12(1:L12,j,i) = MeanVectorAPTemp(nc12(i,j):nc12(i,j)+L12-1,j);
            SDVectorAP_12(1:L12,j,i) = SDVectorAPTemp(nc12(i,j):nc12(i,j)+L12-1,j);
            NParticlesAP_12(1:L12,j,i) = NParticlesAPTemp(nc12(i,j):nc12(i,j)+L12-1,j);
            FractionON_Instant_12(1:L12,j,i) = NParticlesAPTemp(nc12(i,j):nc12(i,j)+L12-1,j)./...
                                                NEllipsesAPT
emp(nc12(i,j):nc12(i,j)+L12-1,j);
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
%             SDVectorAP(1:numFrames(i)-nc12(i)+1,:,i) = Data(i).SDVectorAP(nc12(i):numFrames(i),:);
%             NParticlesAP(1:numFrames(i)-nc12(i)+1,:,i) = Data(i).NParticlesAP(nc12(i):numFrames(i),:);
%             FractionOn_Instant(1:numFrames(i)-nc13(i)+1,:,i) = Data(i).NParticlesAP(nc13(i):numFrames(i),:)./Data(i).NEllipsesAP(nc13(i):numFrames(i),:)

        %elseif nc13(i)==0  
        %    error('Check the Movie if it really does not start from nc13, then you should edit this code or make that dataset as an exception')
        elseif nc13(i,j)~=0 && nc14(i,j)~=0
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
        end
    end
end

% % Take the most frequent value of dT from the ElapsedTime. It's because the
% % dT can be different in case we stopped and restarted the movie.
deltaT = mode(diff(Data(1).ElapsedTime)); 
ElapsedTime = deltaT*(0:NewFrameLength-1);

%% Average all fields at each time point
% %% Convert NaNs to zeros
% 
% % % Note that this should be done only for the APbins that has values.
% % % If we convert all NaNs to zeros for APbins that we didn't measure, it's
% % % misleading.
% 
% % MeanVectorAP(iszeros(MeanVectorAP)) = 0;
% % SDVectorAP(iszeros(SDVectorAP)) = 0;
% % NParticlesAP(iszeros(NParticlesAP)) = 0;
% 
%     
% if NC==12
%     MeanVectorAP_12(iszeros(MeanVectorAP_12)) = 0;
%     SDVectorAP_12(iszeros(SDVectorAP_12)) = 0;
%     NParticlesAP_12(iszeros(NParticlesAP_12)) = 0;
%     %FractionON_Instant_12(iszeros(FractionON_Instant_12)) = 0;
% 
% end
% 
% MeanVectorAP_13(iszeros(MeanVectorAP_13)) = 0;
% SDVectorAP_13(iszeros(SDVectorAP_13)) = 0;
% NParticlesAP_13(iszeros(NParticlesAP_13)) = 0;
% %FractionON_Instant_13(iszeros(FractionON_Instant_13)) = 0;
% 
% MeanVectorAP_14(iszeros(MeanVectorAP_14)) = 0;
% SDVectorAP_14(iszeros(SDVectorAP_14)) = 0;
% NParticlesAP_14(iszeros(NParticlesAP_14)) = 0;
% %FractionON_Instant_14(iszeros(FractionON_Instant_14)) = 0;

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
FractionON__individual = FractionON;
MeanVectorAP_ForSum_individual = MeanVectorAP_ForSum;

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

%% Averaging - weigthed sum
sumMean = zeros(NewFrameLength,numAPBins);
sumSD = zeros(NewFrameLength,numAPBins);
sumNParticles = zeros(NewFrameLength,numAPBins);
sumMean_ForSum = zeros(NewFrameLength,numAPBins);

for i=1:numEmbryos
    sumMean = sumMean + squeeze(MeanVectorAP(:,:,i).*NParticlesAP(:,:,i));
    sumSD = sumSD + squeeze(SDVectorAP(:,:,i).^2.*NParticlesAP(:,:,i));
    sumNParticles = sumNParticles + squeeze(NParticlesAP(:,:,i));
    % MeanVectorAP considering the doubling of nuclei in each cycle
    sumMean_ForSum = sumMean_ForSum + squeeze(MeanVectorAP_ForSum(:,:,i).*NParticlesAP(:,:,i));
end
    
MeanVectorAP_combined = sumMean./sumNParticles;
SDVectorAP_combined = sqrt(sumSD./sumNParticles);
NParticlesAP_combined = sumNParticles;
MeanVectorAP_ForSum_combined = sumMean_ForSum./sumNParticles;

% For FractionON, we will use nanmean, to ignore the APbin(of embryos) that
% doesn't have any values.
FractionON_temp = FractionON;
FractionON_temp(FractionON_temp==0) = nan;
FractionON_combined = nanmean(FractionON_temp,3);
FractionON_combined(isnan(FractionON_combined)) = 0;

% Convert Nans to Zeros for (Mean)VectorAP_combined
MeanVectorAP_combined(isnan(MeanVectorAP_combined)) = 0;  
SDVectorAP_combined(isnan(SDVectorAP_combined)) = 0;
NParticlesAP_combined(isnan(NParticlesAP_combined)) = 0;
MeanVectorAP_ForSum_combined(isnan(MeanVectorAP_ForSum_combined)) = 0;

%% Accumulate mRNA over time (This can be made as an option, or a separate function)
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
AccumulatedmRNA_FractionON_SE =  zeros(NFrames,nAPbins);

for i=1:maxAPIndex
    for j=2:length(ElapsedTime)
        % This is assuming that all nuclei are on, then thinking about how
        % many mRNA molecules are made in 
        AccumulatedmRNA(j,i) = trapz(ElapsedTime(1:j),MeanVectorAP_combined(1:j,i)); 
        AccumulatedmRNA_SD(j,i) = sqrt(trapz(ElapsedTime(1:j),SDVectorAP_combined(1:j,i).^2));
        
        AccumulatedmRNA_FractionON(j,i) = trapz(ElapsedTime(1:j),MeanVectorAP_combined(1:j,i).*FractionON_combined(1:j,i));
        % Define the error for the Accumulated mRNA, this should be
        % revisited later for with more justification.
        AccumulatedmRNA_FractionON_SD(j,i) = sqrt(trapz(ElapsedTime(1:j),SDVectorAP(1:j,i).^2.*FractionON_combined(1:j,i)));
        AccumulatedmRNA_FractionON_SE(j,i) = sqrt(trapz(ElapsedTime(1:j),SDVectorAP(1:j,i).^2.*FractionON_combined(1:j,i)./numEmbryos));
    end
end

 %% Quick plot to check - Accumulated mRNA
 
% APaxis = 0:0.025:1;
% 
% % hold on
% % for i=1:length(ElapsedTime)
% %     errorbar(APaxis, AccumulatedmRNA(i,:), AccumulatedmRNA_SD(i,:))
% %     pause
% % end
% 
% hold on
% for i=1:length(ElapsedTime)
%     errorbar(APaxis, AccumulatedmRNA_FractionON(i,:), AccumulatedmRNA_FractionON_SD(i,:))
%     pause
% end
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
% FractionON_individual = zeros(41,3,numEmbryos);
% for i=1:numEmbryos
%     FractionON_individual(:,:,i) = Data(i).EllipsesOnAP./Data(i).TotalEllipsesAP;
% end
% FractionON_Average = zerosmean(FractionON_individual,3);
% FractionON_Average_Error = zerosstd(FractionON_individual,[],3) / sqrt(numEmbryos);
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
% TotalONNuclei = zeros(41,3);
% TotalNuclei = zeros(41,3);
%     
% for i=1:numEmbryos
%     TotalONNuclei = TotalONNuclei + Data(i).EllipsesOnAP;
%     TotalNuclei = TotalNuclei + Data(i).TotalEllipsesAP;
% end
% 
% FractionONTemp = TotalONNuclei ./ TotalNuclei;
% FractionON_Global = FractionONTemp;

%% Define the fields to be saved
MeanVectorAP = MeanVectorAP_combined;
SDVectorAP = SDVectorAP_combined;
SEVectorAP = SDVectorAP_combined/sqrt(numEmbryos); % Standard error of mean (SD / sqrt(number of observation)
NParticlesAP = NParticlesAP_combined;
ElapsedTime = ElapsedTime;
FractionON = FractionON_combined;
MeanVectorAP_ForSum = MeanVectorAP_ForSum_combined;

% AccumulatedmRNA 
% AccumulatedmRNA_SD 
% 
% AccumulatedmRNA_FractionON 
% AccumulatedmRNA_FractionON_SD 
% AccumulatedmRNA_FractionON_SE 
%% Convert zeros to Nans
% Convert the zeros in APbins to Nans, since those were not measured.
MeanVectorAP(MeanVectorAP ==0) = nan;
SDVectorAP(SDVectorAP ==0) = nan;
SEVectorAP(SEVectorAP ==0) = nan; % Standard error of mean (SD / sqrt(number of observation)
NParticlesAP(NParticlesAP==0) = nan;
MeanVectorAP_ForSum(MeanVectorAP_ForSum==0) = nan;

AccumulatedmRNA(AccumulatedmRNA ==0) = nan; 
AccumulatedmRNA_SD(AccumulatedmRNA_SD ==0) = nan;

AccumulatedmRNA_FractionON(AccumulatedmRNA_FractionON==0) = nan;
AccumulatedmRNA_FractionON_SD(AccumulatedmRNA_FractionON_SD==0) = nan; 
AccumulatedmRNA_FractionON_SE(AccumulatedmRNA_FractionON_SE ==0) = nan; 
%% Save the fields in .mat file
    if ~isempty(savePath)
        save([savePath,filesep,DataType,'.mat'],...
            'MeanVectorAP','SDVectorAP','SEVectorAP','NParticlesAP','ElapsedTime',...
                        'nc12', 'nc13', 'nc14',...
            'AccumulatedmRNA','AccumulatedmRNA_SD', 'AccumulatedmRNA_FractionON',...
            'AccumulatedmRNA_FractionON_SD',...
            'MeanVectorAP_individual','SDVectorAP_individual','NParticlesAP_individual')
            %'FractionON', 'FractionON_Global',...
            %'FractionON_Average','FractionON_Average_Error',...

    else
        warning('Define the File Path in the argument above')
    end
end