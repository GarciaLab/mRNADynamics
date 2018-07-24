% Function to plot the movement of the gene expression line as a function of
% time
function [LineOfRepression,ExpressionLineWithoutNoise,YPositionOfAllParticlesInFrame,FluorescenceOfAllParticlesInFrame, er, NumberOfRuns]=...
    FindEctogenicBoundary(CompiledParticles,NumberOfFrames,Direction,StartFrame,LinesPerFrame, PixelsPerLine,cf)
%% Initialisation
close all;
% We first get all the required variables

% In order to process bidirectional repression, we split the image into
% two. This functionality could be extended, I guess.
% The number of runs is two if bidirectional, one otherwise
NumberOfRuns=1;
if Direction=='b'
    NumberOfRuns=2;
end
%ErrorOfAllParticlesInFrame
% Note that the image is inverted in matlab, so 0 is the top in the y-axis
LineOfRepression=zeros(NumberOfRuns,1);
for RunIndex=1:NumberOfRuns
    if NumberOfRuns==1
        UpperBoundary=LinesPerFrame;
        LowerBoundary=0;
    elseif NumberOfRuns==2
        if RunIndex==1
            Direction='u';
            UpperBoundary=LinesPerFrame/2;
            LowerBoundary=0;
        elseif RunIndex==2
            Direction='d';
            UpperBoundary=LinesPerFrame;
            LowerBoundary=LinesPerFrame/2;
        end
    end
    
    %% Clean up the negative values
    % This shouldn't be much of an issue. Negative values come up due to
    % background correction, and are almost always very small values
    for i=1:length(CompiledParticles)
        for j=1:length(CompiledParticles(i).Fluo)
            a=CompiledParticles(i).Fluo;
            if sum(a<0)>0
                CompiledParticles(i).Fluo(a<0)=0;
            end
        end
    end
    
    %% Find the maximum particle y posn
    if Direction == 'u'
    LimitOfExpressionDetected=LowerBoundary*ones(1,NumberOfFrames); %This is the array of minimum y posn
    end
    % We stipulate that this line represents the line of repression
    if (Direction=='d')
        LimitOfExpressionDetected=ones(1,NumberOfFrames)*UpperBoundary;
    end
    
    NumberOfExpressingSitesInFrame=zeros(1,NumberOfFrames);
    
    h=waitbar(0,'Finding Boundary and Characterising Expression');
    for i=1:length(CompiledParticles)
        waitbar(i/length(CompiledParticles));
        FramesOfCurrentParticle=CompiledParticles(i).Frame;
        % The frames in which this particle appears
        for j=StartFrame:NumberOfFrames
            % Check if each frame is in the array of frames for the particle
            IndexOfFrame=find(FramesOfCurrentParticle==j);
            if (length(IndexOfFrame)==1)
                NumberOfExpressingSitesInFrame(j)=NumberOfExpressingSitesInFrame(j)+1;
                if (Direction=='u')
                    if (CompiledParticles(i).yPos(IndexOfFrame)>LimitOfExpressionDetected(j)) && (CompiledParticles(i).yPos(IndexOfFrame)>LowerBoundary) ...
                            && (CompiledParticles(i).yPos(IndexOfFrame)<UpperBoundary)
                        LimitOfExpressionDetected(j)=CompiledParticles(i).yPos(IndexOfFrame);
                    end
                elseif (Direction=='d')
                    if (CompiledParticles(i).yPos(IndexOfFrame)<LimitOfExpressionDetected(j)) && (CompiledParticles(i).yPos(IndexOfFrame)<UpperBoundary) ...
                            && (CompiledParticles(i).yPos(IndexOfFrame)>LowerBoundary)
                        LimitOfExpressionDetected(j)=CompiledParticles(i).yPos(IndexOfFrame);
                    end
                end
            end
        end
    end
    close(h);
    
    %% Minimum Expression Noise Plots
    % Plot with noise
    figure('Name','Noisy')
    if Direction == 'u'
        LimitOfExpressionDetected(LimitOfExpressionDetected==UpperBoundary)=NaN;
    elseif Direction== 'd'
        LimitOfExpressionDetected(LimitOfExpressionDetected==LowerBoundary)=NaN;
    end
    %if(inverted=='u')
    plot(StartFrame:NumberOfFrames, (LimitOfExpressionDetected(StartFrame:NumberOfFrames)));
    %elseif(inverted=='d')
    %    plot(nc13:NumberOfFrames, ones(-nc13+NumberOfFrames+1,1)'-expline(nc13:NumberOfFrames)/ht);
    %end
    title('Line of minimum expression (with noise)')
    xlabel('Frame (Only nc13 and nc14 are included)')
    ylabel('Minimum y-coordinate of gene expression')
    saveas(gcf,['Line of minimum expression (with noise)_',num2str(RunIndex),'.png'])
    
    % Plot without noise
    figure('Name','Noiseless')
    ExpressionLineWithoutNoise=LimitOfExpressionDetected;
    ExpressionLineWithoutNoise(NumberOfExpressingSitesInFrame<3)=NaN;
    %if(inverted=='u')
        plot(StartFrame:NumberOfFrames, (ExpressionLineWithoutNoise(StartFrame:NumberOfFrames)));
    %elseif(inverted=='d')
     %   plot(nc13:NumberOfFrames, ones(-nc13+NumberOfFrames+1,1)'-allpartsy2(nc13:NumberOfFrames)/ht);
    %end
    title('Line of minimum expression (frames with at least three particles)')
    xlabel('Frame (Only nc13 and nc14 are included)')
    ylabel('Minimum y-coordinate of gene expression')
    saveas(gcf,['Line of minimum expression (without noise=3)_',num2str(RunIndex),'.png'])
    
    %% Compile Traces
    FluorescenceOfAllParticlesInFrame=cell(NumberOfFrames,1);
    YPositionOfAllParticlesInFrame=cell(NumberOfFrames,1);
    XPositionOfAllParticlesInFrame=cell(NumberOfFrames,1);
    er=cell(NumberOfFrames,1);
    
    for i=1:NumberOfFrames
        FluorescenceOfAllParticlesInFrame{i,1}=[];
        YPositionOfAllParticlesInFrame{i,1}=[];
        XPositionOfAllParticlesInFrame{i,1}=[];
        er{i,1}=[];
    end
    
    for i=1:length(CompiledParticles)
        FramesOfCurrentParticle=CompiledParticles(i).Frame;
        % The frames in which this particle appears
        for j=StartFrame:NumberOfFrames
            % Check if each frame is in the array of frames for the particle
            IndexOfFrame=find(FramesOfCurrentParticle==j);
            if length(IndexOfFrame)>1
                 fprintf('Something is wrong, there are %i frames with the same number %i, check particle %i \n',length(IndexOfFrame),j,i)
            end
            if (length(IndexOfFrame)==1)
                YPositionOfAllParticlesInFrame{j,1}=[YPositionOfAllParticlesInFrame{j,1} CompiledParticles(i).yPos(IndexOfFrame)];
                XPositionOfAllParticlesInFrame{j,1}=[XPositionOfAllParticlesInFrame{j,1} CompiledParticles(i).xPos(IndexOfFrame)];
                FluorescenceOfAllParticlesInFrame{j,1}=[FluorescenceOfAllParticlesInFrame{j,1} CompiledParticles(i).Fluo(IndexOfFrame)];
                er{j,1}=[er{j,1}, ones(size(CompiledParticles(1).Fluo))*...
                                CompiledParticles(1).FluoError];
            end
        end
    end
    
    
    %% Compute Statistics
    % Define constant ectogenic boundary based on when the expression line
    % plateaus (doesn't move significantly). The longest such plateau is
    % defined as the boundary
    FinalEctogenicBoundary=0;
    Plateau=0;
    LengthOfCurrentPlateau=0;
    for i=2:cf
        if ~isnan(LimitOfExpressionDetected(i)) && ~isnan(LimitOfExpressionDetected(i-1))
            if LimitOfExpressionDetected(i)-LimitOfExpressionDetected(i-1)<10
                LengthOfCurrentPlateau=LengthOfCurrentPlateau+1;
                if LengthOfCurrentPlateau>Plateau
                    Plateau=LengthOfCurrentPlateau;
                    FinalEctogenicBoundary=mean(LimitOfExpressionDetected(i-LengthOfCurrentPlateau+1:i));
                end
            else
                LengthOfCurrentPlateau=0;
            end
        end
    end
    LineOfRepression(RunIndex)=FinalEctogenicBoundary;
end   
%     means1=zeros(NumberOfFrames,1);
%     means2=zeros(NumberOfFrames,1);
%     
%     min1=zeros(NumberOfFrames,1);
%     min2=zeros(NumberOfFrames,1);
%     
%     max1=zeros(NumberOfFrames,1);
%     max2=zeros(NumberOfFrames,1);
%     
%     std1=zeros(NumberOfFrames,1);
%     std2=zeros(NumberOfFrames,1);
%     
%     count1=zeros(NumberOfFrames,1);
%     count2=zeros(NumberOfFrames,1);
%     
%     sum1=zeros(NumberOfFrames,1);
%     sum2=zeros(NumberOfFrames,1);
%     
%     for i=nc13:NumberOfFrames
%         if (inverted=='u')
%             traces1=find(alltracesy{i,1}(:)>lrep & alltracesy{i,1}(:)<lrep+100);
%             traces2=find(alltracesy{i,1}(:)<lrep & alltracesy{i,1}(:)>lrep-100);
%         elseif (inverted=='d')
%             traces1=find(alltracesy{i,1}(:)<lrep & alltracesy{i,1}(:)>lrep-100);
%             traces2=find(alltracesy{i,1}(:)>lrep & alltracesy{i,1}(:)<lrep+100);
%         end
%         if ~isempty(traces1)
%             min1(i)=min(alltraces{i,1}(traces1));
%             if min1(i)<0
%                 fprintf('Check frame %i with intensity %d\n',i,min1(i));
%             end
%             max1(i)=max(alltraces{i,1}(traces1));
%             std1(i)=std(alltraces{i,1}(traces1));
%             count1(i)=length(traces1);
%             sum1(i)=sum(alltraces{i,1}(traces1));
%         end
%         
%         if ~isempty(traces2)
%             min2(i)=min(alltraces{i,1}(traces2));
%             if min1(i)<0
%                 fprintf('Check frame %i with intensity %d\n',i,min1(i));
%             end
%             max2(i)=max(alltraces{i,1}(traces2));
%             std2(i)=std(alltraces{i,1}(traces2));
%             count2(i)=length(alltraces{i,1}(traces2));
%             sum2(i)=sum(alltraces{i,1}(traces2));
%         end
%         
%         means1(i)=mean(alltraces{i,1}(traces1));
%         means2(i)=mean(alltraces{i,1}(traces2));
%     end
%     
%     %% Clean up the statistics
%     for i=nc13:NumberOfFrames
%         bd=isnan(max1); % Finds the nans (bad data)
%         gd=find(~bd); % Indices of actual points (good data)
%         bd([1:(min(gd)-1) (max(gd)+1):end])=0; % Fills extremes with 0
%         max1(bd)=interp1(gd,max1(gd),find(bd)); % Interpolates data
%         min1(bd)=interp1(gd,min1(gd),find(bd)); % Interpolates data
%         std1(bd)=interp1(gd,std1(gd),find(bd)); % Interpolates data
%         means1(bd)=0; % Interpolates data
%         
%         bd=isnan(max2); % Finds the nans (bad data)
%         gd=find(~bd); % Indices of actual points (good data)
%         bd([1:(min(gd)-1) (max(gd)+1):end])=0; % Fills extremes with 0
%         max2(bd)=interp1(gd,max2(gd),find(bd)); % Interpolates data
%         min2(bd)=interp1(gd,min2(gd),find(bd)); % Interpolates data
%         std2(bd)=interp1(gd,std2(gd),find(bd)); % Interpolates data
%         means2(bd)=0; % Interpolates data
%     end
%     
%     %% Plot the statistics
%     figure('Name','Mean Intensity')
%     plot(nc13:NumberOfFrames, means1(nc13:NumberOfFrames),'g')
%     hold on
%     plot(nc13:NumberOfFrames, means2(nc13:NumberOfFrames),'b')
%     xlabel('Frame (Only nc13 and nc14 are included)')
%     ylabel('Mean Intensity (arb. units)')
%     legend('Mesoderm','Neurogenic Ectoderm');
%     saveas(gcf,'Mean Intensity Plots.png');
%     
%     figure('Name','Minimum Intensity')
%     plot(nc13:NumberOfFrames, min1(nc13:NumberOfFrames),'g')
%     hold on
%     plot(nc13:NumberOfFrames, min2(nc13:NumberOfFrames),'b')
%     xlabel('Frame (Only nc13 and nc14 are included)');
%     ylabel('Minimum Intensity (arb. units)');
%     legend('Mesoderm','Neurogenic Ectoderm');
%     saveas(gcf, 'Minimum Intensity Plots.png');
%     
%     figure('Name','Maximum Intensity Plots')
%     plot(nc13:NumberOfFrames, max1(nc13:NumberOfFrames),'g')
%     hold on
%     plot(nc13:NumberOfFrames, max2(nc13:NumberOfFrames),'b')
%     xlabel('Frame (Only nc13 and nc14 are included)')
%     ylabel('Maximum Intensity (arb. units)');
%     legend('Mesoderm','Neurogenic Ectoderm');
%     saveas(gcf, 'Maximum Intensity Plots.png');
%     
%     figure('Name','Standard Deviation of Intensity')
%     plot(nc13:NumberOfFrames, std1(nc13:NumberOfFrames),'g')
%     hold on
%     plot(nc13:NumberOfFrames, std2(nc13:NumberOfFrames),'b')
%     xlabel('Frame (Only nc13 and nc14 are included)')
%     ylabel('Standard Deviation of Intensity (arb. units)');
%     legend('Mesoderm','Neurogenic Ectoderm');
%     saveas(gcf, 'Standard Deviation Plots.png');
%     
%     figure('Name','Density of Particles')
%     plot(nc13:NumberOfFrames, count1(nc13:NumberOfFrames),'g')
%     hold on
%     plot(nc13:NumberOfFrames, count2(nc13:NumberOfFrames),'b')
%     xlabel('Frame (Only nc13 and nc14 are included)')
%     ylabel('Spatial Density of Expressing Particles (arb. units)');
%     legend('Mesoderm','Neurogenic Ectoderm');
%     saveas(gcf, 'Spatial Density Plots.png');
%     
%     figure('Name','Sum intensity of Particles')
%     plot(nc13:NumberOfFrames, sum1(nc13:NumberOfFrames),'g')
%     hold on
%     plot(nc13:NumberOfFrames, sum2(nc13:NumberOfFrames),'b')
%     xlabel('Frame (Only nc13 and nc14 are included)')
%     ylabel('Total Intensity (arb. units)');
%     legend('Mesoderm','Neurogenic Ectoderm');
%     saveas(gcf, 'Sum Intensity Plots.png');