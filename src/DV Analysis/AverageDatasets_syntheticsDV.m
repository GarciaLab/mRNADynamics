%% Averaging multiple datasets
function AverageDatasets_syntheticsDV(DataType,varargin)
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

if ischar(DataType)
    Data= LoadMS2Sets(DataType);
else   
    Data = DataType;
end

% Save path option
savePath = '';

for i=1:length(varargin)
    if strcmpi(varargin{i}, 'savePath')
        savePath = varargin{i+1};
    end
end

NC = 12; % Default
for i = 1:length(varargin)
   if strcmpi(varargin{i}, 'NC')
        NC = varargin{i+1};  
   end
end

numEmbryos=length(Data);

%Find the total number of frames for each embryo
numFrames = NaN(1, numEmbryos);
maxDVIndex = NaN(1, numEmbryos);
maxTime = NaN(1, numEmbryos);

for i = 1:numEmbryos
    if NC~=12
        numFrames(i) = size(Data(i).Particles.ElapsedTime, 2);
        nc13(i) = Data(i).Particles.nc13;
        length_total_NC(i) = numFrames(i) - nc13(i)+1; % length of frames from nc13 to the last frame
        maxDVIndex(i) = size(Data(i).Particles.DVbinID,2);
        maxTime(i) = Data(i).Particles.ElapsedTime(numFrames(i));
    elseif NC==12
        numFrames(i) = size(Data(i).Particles.ElapsedTime, 2);
        nc12(i) = Data(i).Particles.nc12;
        length_total_NC(i) = numFrames(i) - nc12(i)+1; % length of frames from nc13 to the last frame
        %maxAPIndex(i) = Data(i).Particles.MaxAPIndex; % Usually 41, in 2.5% binning
        maxDVIndex(i) = size(Data(i).Particles.DVbinID,2);
        maxTime(i) = Data(i).Particles.ElapsedTime(numFrames(i));
    end
end

%Store the number of AP bins (this should always be 41).
%numAPBins = maxAPIndex(1);
numAPBins = 41;
numDVBins = maxDVIndex(1);

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
    Length12 = zeros(numEmbryos,numDVBins);
end
Length13 = zeros(numEmbryos,numDVBins);
Length14 = zeros(numEmbryos,numDVBins);

for j=1:numDVBins
    % For all embryos, go through to find the longest nuclear cycle (number
    % of frames)
    for i=1:numEmbryos
        % Define the nuclear cycle (In case we start from nc 12)
        if NC==12
            %nc12(i,j) = Data(i).Particles.APDivision(12,j);
            nc12(i,j) = max(Data(i).Particles.APDivision(12,:));
        end
            %nc13(i,j) = Data(i).Particles.APDivision(13,j);
            nc13(i,j) = max(Data(i).Particles.APDivision(13,:));
            %nc14(i,j) = Data(i).Particles.APDivision(14,j);
            nc14(i,j) = max(Data(i).Particles.APDivision(14,:));
            
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
 
%NewFrameLength = max(TotalFrames);

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
    MeanVectorDV_12 = zeros(L12,numDVBins,numEmbryos);
    SDVectorDV_12 = zeros(L12,numDVBins,numEmbryos);
    NParticlesDV_12 = zeros(L12,numDVBins,numEmbryos);
    FractionON_Instant_12 = zeros(L12,numDVBins,numEmbryos);
    NEllipsesDV_12 = zeros(L12,numDVBins,numEmbryos);
    
    MeanVectorAll_12 = zeros(L12,numEmbryos);
    SDVectorAll_12 = zeros(L12,numEmbryos);
    NParticlesAll_12 = zeros(L12,numEmbryos);
    FractionONAll_Instant_12 = zeros(L12,numEmbryos);
    NEllipsesAll_12 = zeros(L12,numEmbryos);
end
MeanVectorDV_13 = zeros(L13,numDVBins,numEmbryos);
SDVectorDV_13 = zeros(L13,numDVBins,numEmbryos);
NParticlesDV_13 = zeros(L13,numDVBins,numEmbryos);
FractionON_Instant_13 = zeros(L13,numDVBins,numEmbryos);
NEllipsesDV_13 = zeros(L13,numDVBins,numEmbryos);

MeanVectorAll_13 = zeros(L13,numEmbryos);
SDVectorAll_13 = zeros(L13,numEmbryos);
NParticlesAll_13 = zeros(L13,numEmbryos);
FractionONAll_Instant_13 = zeros(L13,numEmbryos);
NEllipsesAll_13 = zeros(L13,numEmbryos);

MeanVectorDV_14 = zeros(L14,numDVBins,numEmbryos);
SDVectorDV_14 = zeros(L14,numDVBins,numEmbryos);
NParticlesDV_14 = zeros(L14,numDVBins,numEmbryos);
FractionON_Instant_14 = zeros(L14,numDVBins,numEmbryos);
NEllipsesDV_14 = zeros(L14,numDVBins,numEmbryos);

MeanVectorAll_14 = zeros(L14,numEmbryos);
SDVectorAll_14 = zeros(L14,numEmbryos);
NParticlesAll_14 = zeros(L14,numEmbryos);
FractionONAll_Instant_14 = zeros(L14,numEmbryos);
NEllipsesAll_14 = zeros(L14,numEmbryos);

% Total matrices
NewFrameLength = L12+L13+L14;

MeanVectorAll = zeros(NewFrameLength,numEmbryos);
MeanVectorDV = zeros(NewFrameLength,numDVBins,numEmbryos);
SDVectorAll = zeros(NewFrameLength,numEmbryos);
SDVectorDV = zeros(NewFrameLength,numDVBins,numEmbryos);
NParticlesDV = zeros(NewFrameLength,numDVBins,numEmbryos);
NParticlesAll = zeros(NewFrameLength,numEmbryos);
FractionON_Instant = zeros(NewFrameLength,numDVBins,numEmbryos);
FractionON_InstantAll = zeros(NewFrameLength,numEmbryos);

% Synchornize all fields as all of them starts from nc 13 (or 12)
% This synchronization should be done for each AP bin, since they might
% have different anaphase time point.

% First, calculate for all particles

% Then calculate for each bins
for i=1:numEmbryos
    % First, 
    MeanVectorAllTemp = cell2mat(Data(i).Particles.MeanVectorAll); %Data(i).MeanVectorAll;
    SDVectorAllTemp = cell2mat(Data(i).Particles.SDVectorAll); %Data(i).SDVectorAll;
    NParticlesAllTemp = cell2mat(Data(i).Particles.NParticlesAll); %Data(i).NParticlesAll;
    NEllipsesDVTemp =Data(i).Particles.NEllipsesDV;  %Data(i).NEllipsesAll;
    for k = 1:size(NEllipsesDVTemp,1)
        NEllipsesAllTemp(k) = sum(NEllipsesDVTemp(k,~isinf(NEllipsesDVTemp(k,:))));
    end
    nc12all = Data(i).Particles.nc12;
    nc13all = Data(i).Particles.nc13;
    nc14all = Data(i).Particles.nc14;
    
    
    if NC==12 && nc12all~=0
        % Sync the fields for each nc
        % NC12
        MeanVectorAll_12(1:L12,i) = MeanVectorAllTemp(nc12all:nc12all+L12-1);
        SDVectorAll_12(1:L12,i) = SDVectorAllTemp(nc12all:nc12all+L12-1);
        NParticlesAll_12(1:L12,i) = NParticlesAllTemp(nc12all:nc12all+L12-1);
        FractionONAll_Instant_12(1:L12,i) = NParticlesAllTemp(nc12all:nc12all+L12-1)./...
                                            NEllipsesAllTemp(nc12all:nc12all+L12-1);
        NEllipsesAll_12(1:L12,i) = NEllipsesAllTemp(nc12all:nc12all+L12-1);
        % NC13
        MeanVectorAll_13(1:min(L13,nc14all-nc13all),i) = MeanVectorAllTemp(nc13all:nc13all+min(L13,nc14all-nc13all)-1);
        SDVectorAll_13(1:min(L13,nc14all-nc13all),i) = SDVectorAllTemp(nc13all:nc13all+min(L13,nc14all-nc13all)-1);
        NParticlesAll_13(1:min(L13,nc14all-nc13all),i) = NParticlesAllTemp(nc13all:nc13all+min(L13,nc14all-nc13all)-1);
        FractionONAll_Instant_13(1:min(L13,nc14all-nc13all),i) = NParticlesAllTemp(nc13all:nc13all+min(L13,nc14all-nc13all)-1)./...
                                            NEllipsesAllTemp(nc13all:nc13all+min(L13,nc14all-nc13all)-1);
        NEllipsesAll_13(1:min(L13,nc14all-nc13all),i) = NEllipsesAllTemp(nc13all:nc13all+min(L13,nc14all-nc13all)-1);
        % NC14                          
        MeanVectorAll_14(1:min(L14,numFrames(i)-nc14all),i) = MeanVectorAllTemp(nc14all:nc14all+min(L14,numFrames(i)-nc14all)-1);
        SDVectorAll_14(1:min(L14,numFrames(i)-nc14all),i) = SDVectorAllTemp(nc14all:nc14all+min(L14,numFrames(i)-nc14all)-1);
        NParticlesAll_14(1:min(L14,numFrames(i)-nc14all),i) = NParticlesAllTemp(nc14all:nc14all+min(L14,numFrames(i)-nc14all)-1);
        FractionONAll_Instant_14(1:min(L14,numFrames(i)-nc14all),i) = NParticlesAll(nc14all:nc14all+min(L14,numFrames(i)-nc14all)-1)./...
                                            NEllipsesAllTemp(nc14all:nc14all+min(L14,numFrames(i)-nc14all)-1);
        NEllipsesAll_14(1:min(L14,numFrames(i)-nc14all),i) = NEllipsesAllTemp(nc14all:nc14all+min(L14,numFrames(i)-nc14all)-1);
%             SDVectorAP(1:numFrames(i)-nc12(i)+1,:,i) = Data(i).SDVectorAP(nc12(i):numFrames(i),:);
%             NParticlesAP(1:numFrames(i)-nc12(i)+1,:,i) = Data(i).NParticlesAP(nc12(i):numFrames(i),:);
%             FractionOn_Instant(1:numFrames(i)-nc13(i)+1,:,i) = Data(i).NParticlesAP(nc13(i):numFrames(i),:)./Data(i).NEllipsesAP(nc13(i):numFrames(i),:)
    %elseif nc13(i)==0  
    %    error('Check the Movie if it really does not start from nc13, then you should edit this code or make that dataset as an exception')
    %{
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
%}    
    end
end

for j=1:numDVBins       
    for i=1:numEmbryos
        % First, 
        Data(i).Particles.NEllipsesDV(Data(i).Particles.NEllipsesDV==0)=inf;
        MeanVectorDVTemp = cell2mat(Data(i).Particles.MeanVectorDV); %Data(i).MeanVectorAP;
        SDVectorDVTemp = cell2mat(Data(i).Particles.SDVectorDV); %Data(i).SDVectorAP;
        NParticlesDVTemp = cell2mat(Data(i).Particles.NParticlesDV); %Data(i).NParticlesAP;
        NEllipsesDVTemp =Data(i).Particles.NEllipsesDV;  %Data(i).NEllipsesAP;
        
        if NC==12 && nc12(i,j)~=0
            % Sync the fields for each nc
            % NC12
            MeanVectorDV_12(1:L12,j,i) = MeanVectorDVTemp(nc12(i,j):nc12(i,j)+L12-1,j);
            SDVectorDV_12(1:L12,j,i) = SDVectorDVTemp(nc12(i,j):nc12(i,j)+L12-1,j);
            NParticlesDV_12(1:L12,j,i) = NParticlesDVTemp(nc12(i,j):nc12(i,j)+L12-1,j);
            FractionON_Instant_12(1:L12,j,i) = NParticlesDVTemp(nc12(i,j):nc12(i,j)+L12-1,j)./...
                                                NEllipsesDVTemp(nc12(i,j):nc12(i,j)+L12-1,j);
            NEllipsesDV_12(1:L12,j,i) = NEllipsesDVTemp(nc12(i,j):nc12(i,j)+L12-1,j);
            % NC13
            MeanVectorDV_13(1:min(L13,nc14(i,j)-nc13(i,j)),j,i) = MeanVectorDVTemp(nc13(i,j):nc13(i,j)+min(L13,nc14(i,j)-nc13(i,j))-1,j);
            SDVectorDV_13(1:min(L13,nc14(i,j)-nc13(i,j)),j,i) = SDVectorDVTemp(nc13(i,j):nc13(i,j)+min(L13,nc14(i,j)-nc13(i,j))-1,j);
            NParticlesDV_13(1:min(L13,nc14(i,j)-nc13(i,j)),j,i) = NParticlesDVTemp(nc13(i,j):nc13(i,j)+min(L13,nc14(i,j)-nc13(i,j))-1,j);
            FractionON_Instant_13(1:min(L13,nc14(i,j)-nc13(i,j)),j,i) = NParticlesDVTemp(nc13(i,j):nc13(i,j)+min(L13,nc14(i,j)-nc13(i,j))-1,j)./...
                                                NEllipsesDVTemp(nc13(i,j):nc13(i,j)+min(L13,nc14(i,j)-nc13(i,j))-1,j);
            NEllipsesDV_13(1:min(L13,nc14(i,j)-nc13(i,j)),j,i) = NEllipsesDVTemp(nc13(i,j):nc13(i,j)+min(L13,nc14(i,j)-nc13(i,j))-1,j);
            % NC14                          
            MeanVectorDV_14(1:numFrames(i)-nc14(i,j),j,i) = MeanVectorDVTemp(nc14(i,j):numFrames(i)-1,j);
            SDVectorDV_14(1:numFrames(i)-nc14(i,j),j,i) = SDVectorDVTemp(nc14(i,j):numFrames(i)-1,j);
            NParticlesDV_14(1:numFrames(i)-nc14(i,j),j,i) = NParticlesDVTemp(nc14(i,j):numFrames(i)-1,j);
            FractionON_Instant_14(1:numFrames(i)-nc14(i,j),j,i) = NParticlesDVTemp(nc14(i,j):numFrames(i)-1,j)./...
                                                NEllipsesDVTemp(nc14(i,j):numFrames(i)-1,j);
            NEllipsesDV_14(1:numFrames(i)-nc14(i,j),j,i) = NEllipsesDVTemp(nc14(i,j):numFrames(i)-1,j);
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
            MeanVectorDV_13(1:L13,j,i) = MeanVectorDVTemp(nc13(i,j):nc13(i,j)+L13-1,j);
            SDVectorDV_13(1:L13,j,i) = SDVectorDVTemp(nc13(i,j):nc13(i,j)+L13-1,j);
            NParticlesDV_13(1:L13,j,i) = NParticlesDVTemp(nc13(i,j):nc13(i,j)+L13-1,j);
            FractionON_Instant_13(1:L13,j,i) = NParticlesDVTemp(nc13(i,j):nc13(i,j)+L13-1,j)./...
                                               NEllipsesDVTemp(nc13(i,j):nc13(i,j)+L13-1,j);
            NEllipsesDV_13(1:min(L13,nc14(i,j)-nc13(i,j)),j,i) = NEllipsesDVTemp(nc13(i,j):nc13(i,j)+min(L13,nc14(i,j)-nc13(i,j))-1,j);
            % NC14                          
            MeanVectorDV_14(1:numFrames(i)-nc14(i,j),j,i) = MeanVectorDVTemp(nc14(i,j):numFrames(i)-1,j);
            SDVectorDV_14(1:numFrames(i)-nc14(i,j),j,i) = SDVectorDVTemp(nc14(i,j):numFrames(i)-1,j);
            NParticlesDV_14(1:numFrames(i)-nc14(i,j),j,i) = NParticlesDVTemp(nc14(i,j):numFrames(i)-1,j);
            FractionON_Instant_14(1:numFrames(i)-nc14(i,j),j,i) = NParticlesDVTemp(nc14(i,j):numFrames(i)-1,j)./...
                                               NEllipsesDVTemp(nc14(i,j):numFrames(i)-1,j);
            NEllipsesDV_14(1:numFrames(i)-nc14(i,j),j,i) = NEllipsesDVTemp(nc14(i,j):numFrames(i)-1,j);
        end
    end
end

% % Take the most frequent value of dT from the ElapsedTime. It's because the
% % dT can be different in case we stopped and restarted the movie.
deltaT = mode(diff(Data(1).Particles.ElapsedTime)); 
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
    MeanVectorDV = cat(1,MeanVectorDV_12,MeanVectorDV_13,MeanVectorDV_14);
    SDVectorDV = cat(1,SDVectorDV_12,SDVectorDV_13,SDVectorDV_14);
    NParticlesDV = cat(1,NParticlesDV_12,NParticlesDV_13,NParticlesDV_14);
    FractionON = cat(1,FractionON_Instant_12,FractionON_Instant_13,FractionON_Instant_14);
    MeanVectorDV_ForSum = cat(1,MeanVectorDV_12*0.25,MeanVectorDV_13*0.5,MeanVectorDV_14);
    NEllipsesDV = cat(1,NEllipsesDV_12,NEllipsesDV_13,NEllipsesDV_14);
    % For all
    MeanVectorAll = cat(1,MeanVectorAll_12,MeanVectorAll_13,MeanVectorAll_14);
    SDVectorAll = cat(1,SDVectorAll_12,SDVectorAll_13,SDVectorAll_14);
    NParticlesAll = cat(1,NParticlesAll_12,NParticlesAll_13,NParticlesAll_14);
    FractionONAll = cat(1,FractionONAll_Instant_12,FractionONAll_Instant_13,FractionONAll_Instant_14);
    MeanVectorAll_ForSum = cat(1,MeanVectorAll_12*0.25,MeanVectorAll_13*0.5,MeanVectorAll_14);
    NEllipsesAll = cat(1,NEllipsesAll_12,NEllipsesAll_13,NEllipsesAll_14);
elseif NC==13
    MeanVectorDV = cat(1,MeanVectorDV_13,MeanVectorDV_14);
    SDVectorDV = cat(1,SDVectorDV_13,SDVectorDV_14);
    NParticlesDV = cat(1,NParticlesDV_13,NParticlesDV_14);
    FractionON = cat(1,FractionON_Instant_13,FractionON_Instant_14);
    MeanVectorDV_ForSum = cat(1,MeanVectorDV_13*0.5,MeanVectorDV_14);
    NEllipsesDV = cat(1,NEllipsesDV_13,NEllipsesDV_14);
    % For all
    MeanVectorAll = cat(1,MeanVectorAll_13,MeanVectorAll_14);
    SDVectorAll = cat(1,SDVectorAll_13,SDVectorAll_14);
    NParticlesAll = cat(1,NParticlesAll_13,NParticlesAll_14);
    FractionONAll = cat(1,FractionONAll_Instant_13,FractionONAll_Instant_14);
    MeanVectorAll_ForSum = cat(1,MeanVectorAll_13*0.5,MeanVectorAll_14);
    NEllipsesAll = cat(1,NEllipsesAll_13,NEllipsesAll_14);
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
MeanVectorDV_individual = MeanVectorDV;
SDVectorDV_individual = SDVectorDV;
NParticlesDV_individual = NParticlesDV;
FractionON__individual = FractionON;
MeanVectorDV_ForSum_individual = MeanVectorDV_ForSum;
NEllipsesDV_individual = NEllipsesDV;


MeanVectorAll_individual = MeanVectorAll;
SDVectorAll_individual = SDVectorAll;
NParticlesAll_individual = NParticlesAll;
FractionONAll__individual = FractionONAll;
MeanVectorAll_ForSum_individual = MeanVectorAll_ForSum;
NEllipsesAll_individual = NEllipsesAll;

%% Convert NaNs to zeros
% This is for a convenient averaging
MeanVectorDV(isnan(MeanVectorDV)) = 0;
SDVectorDV(isnan(SDVectorDV)) = 0;
NParticlesDV(isnan(NParticlesDV)) = 0;
NEllipsesDV(isinf(NEllipsesDV)) = 0;
FractionON(isnan(FractionON)) = 0;

MeanVectorAll(isnan(MeanVectorAll)) = 0;
SDVectorAll(isnan(SDVectorAll)) = 0;
NParticlesAll(isnan(NParticlesAll)) = 0;
NEllipsesAll(isnan(NEllipsesAll)) = 0;
FractionONAll(isnan(FractionONAll)) = 0;

%% Option
% For r3, I have negative MeanVectorAP for some APbins...
% I'll assign zeros for negative MeanVectorAPs
MeanVectorDV(MeanVectorDV<0) = 0;
NParticlesDV(MeanVectorDV<0) = 0;

%% Averaging - weigthed sum
sumMean = zeros(NewFrameLength,numDVBins);
sumSD = zeros(NewFrameLength,numDVBins);
sumNParticles = zeros(NewFrameLength,numDVBins);
sumMean_ForSum = zeros(NewFrameLength,numDVBins);
sumNEllipses = zeros(NewFrameLength,numDVBins);

sumMeanAll = zeros(NewFrameLength,1);
sumSDAll = zeros(NewFrameLength,1);
sumNParticlesAll = zeros(NewFrameLength,1);
sumMeanAll_ForSum = zeros(NewFrameLength,1);
sumNEllipsesAll = zeros(NewFrameLength,1);

for i=1:numEmbryos
    sumMean = sumMean + squeeze(MeanVectorDV(:,:,i).*NParticlesDV(:,:,i));
    sumSD = sumSD + squeeze(SDVectorDV(:,:,i).^2.*NParticlesDV(:,:,i));
    sumNParticles = sumNParticles + squeeze(NParticlesDV(:,:,i));
    % MeanVectorAP considering the doubling of nuclei in each cycle
    sumMean_ForSum = sumMean_ForSum + squeeze(MeanVectorDV_ForSum(:,:,i).*NParticlesDV(:,:,i));
    sumNEllipses = sumNEllipses + squeeze(NEllipsesDV(:,:,i));
    
    sumMeanAll = sumMeanAll + squeeze(MeanVectorAll(:,i).*NParticlesAll(:,i));
    sumSDAll = sumSDAll + squeeze(SDVectorAll(:,i).^2.*NParticlesAll(:,i));
    sumNParticlesAll = sumNParticlesAll + squeeze(NParticlesAll(:,i));
    % MeanVectorAP considering the doubling of nuclei in each cycle
    sumMeanAll_ForSum = sumMeanAll_ForSum + squeeze(MeanVectorAll_ForSum(:,i).*NParticlesAll(:,i));
    sumNEllipsesAll = sumNEllipsesAll + squeeze(NEllipsesAll(:,i));
end
    
MeanVectorDV_combined = sumMean./sumNParticles;
SDVectorDV_combined = sqrt(sumSD./sumNParticles);
NParticlesDV_combined = sumNParticles;
NEllipsesDV_combined = sumNEllipses;
MeanVectorDV_ForSum_combined = sumMean_ForSum./sumNParticles;
FractionON_combined = sumNParticles./sumNEllipses;

MeanVectorAll_combined = sumMeanAll./sumNParticlesAll;
SDVectorAll_combined = sqrt(sumSDAll./sumNParticlesAll);
NParticlesAll_combined = sumNParticles;
NEllipsesAll_combined = sumNEllipsesAll;
MeanVectorAll_ForSum_combined = sumMeanAll_ForSum./sumNParticlesAll;
FractionONAll_combined = sumNParticlesAll./sumNEllipsesAll;

%{
% For FractionON, we will use nanmean, to ignore the APbin(of embryos) that
% doesn't have any values.
FractionON_temp = FractionON;
FractionON_temp(FractionON_temp==0) = nan;
FractionON_combined = nanmean(FractionON_temp,3);
FractionON_combined(isnan(FractionON_combined)) = 0;

FractionONAll_temp = FractionONAll;
FractionONAll_temp(FractionONAll_temp==0) = nan;
FractionONAll_combined = mean(FractionONAll_temp,2);
FractionONAll_combined(isnan(FractionONAll_combined)) = 0;
%}

% Convert Nans to Zeros for (Mean)VectorAP_combined
MeanVectorDV_combined(isnan(MeanVectorDV_combined)) = 0;  
SDVectorDV_combined(isnan(SDVectorDV_combined)) = 0;
NParticlesDV_combined(isnan(NParticlesDV_combined)) = 0;
MeanVectorDV_ForSum_combined(isnan(MeanVectorDV_ForSum_combined)) = 0;
NEllipsesDV_combined(isnan(NEllipsesDV_combined)) = 0;

MeanVectorAll_combined(isnan(MeanVectorAll_combined)) = 0;  
SDVectorAll_combined(isnan(SDVectorAll_combined)) = 0;
NParticlesAll_combined(isnan(NParticlesAll_combined)) = 0;
MeanVectorAll_ForSum_combined(isnan(MeanVectorAll_ForSum_combined)) = 0;
NEllipsesAll_combined(isnan(NEllipsesAll_combined)) = 0;

%% Accumulate mRNA over time (This can be made as an option, or a separate function)
% I will calculate the Integrated mRNA from the MeanVectorAP
NFrames = length(ElapsedTime);
nDVbins = max(maxDVIndex);

AccumulatedmRNA = zeros(NFrames,nDVbins);
AccumulatedmRNA_SD =  zeros(NFrames,nDVbins);

AccumulatedmRNAAll = zeros(NFrames,1);
AccumulatedmRNAAll_SD = zeros(NFrames,1);

% Consider the FractionON
% For MeanVectorAP, we can just multiply the FractonOn (instantaneous)
% But, for error estimation, we need to think hard. For example, we can use
% bootstrapping, or 
AccumulatedmRNA_FractionON = zeros(NFrames,nDVbins);
AccumulatedmRNA_FractionON_SD =  zeros(NFrames,nDVbins);
AccumulatedmRNA_FractionON_SE =  zeros(NFrames,nDVbins);

AccumulatedmRNAAll_FractionON = zeros(NFrames,1);
AccumulatedmRNAAll_FractionON_SD =  zeros(NFrames,1);
AccumulatedmRNAAll_FractionON_SE =  zeros(NFrames,1);

for i=1:maxDVIndex
    for j=2:length(ElapsedTime)
        % This is assuming that all nuclei are on, then thinking about how
        % many mRNA molecules are made in 
        AccumulatedmRNA(j,i) = trapz(ElapsedTime(1:j),MeanVectorDV_combined(1:j,i)); 
        AccumulatedmRNA_SD(j,i) = sqrt(trapz(ElapsedTime(1:j),SDVectorDV_combined(1:j,i).^2));
        
        AccumulatedmRNA_FractionON(j,i) = trapz(ElapsedTime(1:j),MeanVectorDV_combined(1:j,i).*FractionON_combined(1:j,i));
        % Define the error for the Accumulated mRNA, this should be
        % revisited later for with more justification.
        AccumulatedmRNA_FractionON_SD(j,i) = sqrt(trapz(ElapsedTime(1:j),SDVectorDV(1:j,i).^2.*FractionON_combined(1:j,i)));
        AccumulatedmRNA_FractionON_SE(j,i) = sqrt(trapz(ElapsedTime(1:j),SDVectorDV(1:j,i).^2.*FractionON_combined(1:j,i)./numEmbryos));
    end
end

for i=2:length(ElapsedTime)
    % This is assuming that all nuclei are on, then thinking about how
    % many mRNA molecules are made in 
    AccumulatedmRNAAll(i) = trapz(ElapsedTime(1:i),MeanVectorAll_combined(1:i)); 
    AccumulatedmRNAAll_SD(i) = sqrt(trapz(ElapsedTime(1:i),SDVectorAll_combined(1:i).^2));

    AccumulatedmRNAAll_FractionON(i) = trapz(ElapsedTime(1:i),MeanVectorAll_combined(1:i).*FractionONAll_combined(1:i));
    % Define the error for the Accumulated mRNA, this should be
    % revisited later for with more justification.
    AccumulatedmRNAAll_FractionON_SD(i) = sqrt(trapz(ElapsedTime(1:i),SDVectorAll_combined(1:i).^2.*FractionONAll_combined(1:i)));
    AccumulatedmRNAAll_FractionON_SE(i) = sqrt(trapz(ElapsedTime(1:i),SDVectorAll_combined(1:i).^2.*FractionONAll_combined(1:i)./numEmbryos));
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
MeanVectorDV = MeanVectorDV_combined;
SDVectorDV = SDVectorDV_combined;
SEVectorDV = SDVectorDV_combined/sqrt(numEmbryos); % Standard error of mean (SD / sqrt(number of observation)
NParticlesDV = NParticlesDV_combined;
NEllipsesDV = NEllipsesDV_combined;
ElapsedTime = ElapsedTime;
FractionON = FractionON_combined;
MeanVectorDV_ForSum = MeanVectorDV_ForSum_combined;

MeanVectorAll = MeanVectorAll_combined;
SDVectorAll = SDVectorAll_combined;
SEVectorAll = SDVectorAll_combined/sqrt(numEmbryos); % Standard error of mean (SD / sqrt(number of observation)
NParticlesAll = NParticlesAll_combined;
NEllipsesAll = NEllipsesAll_combined;
FractionONAll = FractionONAll_combined;
MeanVectorAll_ForSum = MeanVectorAll_ForSum_combined;

% AccumulatedmRNA 
% AccumulatedmRNA_SD 
% 
% AccumulatedmRNA_FractionON 
% AccumulatedmRNA_FractionON_SD 
% AccumulatedmRNA_FractionON_SE 
%% Convert zeros to Nans
% Convert the zeros in APbins to Nans, since those were not measured.
MeanVectorDV(MeanVectorDV ==0) = nan;
SDVectorDV(SDVectorDV ==0) = nan;
SEVectorDV(SEVectorDV ==0) = nan; % Standard error of mean (SD / sqrt(number of observation)
NParticlesDV(NParticlesDV==0) = nan;
NEllipsesDV(NEllipsesDV==0) = nan;
MeanVectorDV_ForSum(MeanVectorDV_ForSum==0) = nan;

AccumulatedmRNA(AccumulatedmRNA ==0) = nan; 
AccumulatedmRNA_SD(AccumulatedmRNA_SD ==0) = nan;

AccumulatedmRNA_FractionON(AccumulatedmRNA_FractionON==0) = nan;
AccumulatedmRNA_FractionON_SD(AccumulatedmRNA_FractionON_SD==0) = nan; 
AccumulatedmRNA_FractionON_SE(AccumulatedmRNA_FractionON_SE ==0) = nan; 

MeanVectorAll(MeanVectorAll ==0) = nan;
SDVectorAll(SDVectorAll ==0) = nan;
SEVectorAll(SEVectorAll ==0) = nan; % Standard error of mean (SD / sqrt(number of observation)
NParticlesAll(NParticlesAll==0) = nan;
NEllipsesAll(NEllipsesAll==0) = nan;
MeanVectorAll_ForSum(MeanVectorAll_ForSum==0) = nan;

AccumulatedmRNAAll(AccumulatedmRNAAll ==0) = nan; 
AccumulatedmRNAAll_SD(AccumulatedmRNAAll_SD ==0) = nan;

AccumulatedmRNAAll_FractionON(AccumulatedmRNAAll_FractionON==0) = nan;
AccumulatedmRNAAll_FractionON_SD(AccumulatedmRNAAll_FractionON_SD==0) = nan; 
AccumulatedmRNAAll_FractionON_SE(AccumulatedmRNAAll_FractionON_SE ==0) = nan; 
%% Save the fields in .mat file
    if ~isempty(savePath)
        save([savePath,filesep,DataType,'.mat'],...
            'MeanVectorDV','SDVectorDV','SEVectorDV','NParticlesDV','ElapsedTime',...
                        'nc12', 'nc13', 'nc14',...
            'AccumulatedmRNA','AccumulatedmRNA_SD', 'AccumulatedmRNA_FractionON',...
            'AccumulatedmRNA_FractionON_SD',...
            'MeanVectorDV_individual','SDVectorDV_individual','NParticlesDV_individual',...
            'MeanVectorAll','SDVectorAll','SEVectorAll','NParticlesAll','ElapsedTime',...
            'AccumulatedmRNAAll','AccumulatedmRNAAll_SD', 'AccumulatedmRNAAll_FractionON',...
            'AccumulatedmRNAAll_FractionON_SD',...
            'MeanVectorAll_individual','SDVectorAll_individual','NParticlesAll_individual',...
            'FractionON','FractionONAll','NEllipsesDV','NEllipsesAll');
            
            %'FractionON', 'FractionON_Global',...
            %'FractionON_Average','FractionON_Average_Error',...

    else
        warning('Define the File Path in the argument above')
    end
end