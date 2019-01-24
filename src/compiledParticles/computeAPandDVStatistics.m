function [AllTracesVector, AllTracesAP, AllTracesDV, MeanVectorAP_ROI, ...
    SDVectorAP_ROI, NParticlesAP_ROI, MeanVectorAP_nonROI, SDVectorAP_nonROI, ...
    NParticlesAP_nonROI, MeanVectorAP, SDVectorAP, NParticlesAP, MeanVectorDV_ROI, ...
    SDVectorDV_ROI, NParticlesDV_ROI, MeanVectorDV_nonROI, SDVectorDV_nonROI, ...
    NParticlesDV_nonROI, MeanVectorDV, SDVectorDV, NParticlesDV, ...
    MeanVectorAnterior, MeanVectorAll, SDVectorAll, NParticlesAll, MaxFrame] =...
    computeAPandDVStatistics(NChannels, CompiledParticles, FrameInfo, ExperimentAxis, ...
    APFilter, ROI, CompiledParticles_ROI, CompiledParticles_nonROI, ...
    APFilter_ROI, APFilter_nonROI, NewCyclePos, MaxFrame, ...
    MeanVectorAP_ROI, SDVectorAP_ROI, NParticlesAP_ROI, MeanVectorAP_nonROI, ...
    SDVectorAP_nonROI, NParticlesAP_nonROI, MeanVectorAP, SDVectorAP, ...
    NParticlesAP, MeanVectorDV_ROI, SDVectorDV_ROI, NParticlesDV_ROI, ...
    MeanVectorDV_nonROI, SDVectorDV_nonROI, NParticlesDV_nonROI, ...
    MeanVectorDV, SDVectorDV, NParticlesDV, MeanVectorAnterior, DVFilter_ROI, ...
    DVFilter_nonROI, DVFilter)
%computeAPandDVStatistics Summary of this function goes here
%   Detailed explanation goes here

    for ChN=1:NChannels
        if ~isempty(CompiledParticles{ChN})

            %Get the data for the individual particles in a matrix that has the frame
            %number and the particle number as dimensions. Also, get a vector that
            %reports the mean AP position.

            [AllTracesVector{ChN},AllTracesAP{ChN},AllTracesDV{ChN}]=AllTraces(FrameInfo,CompiledParticles{ChN},'NoAP');

            if strcmpi(ExperimentAxis,'AP') || strcmpi(ExperimentAxis,'DV')
                %Mean plot for different AP positions

                %Figure out the AP range to use
                MinAPIndex=1;%min(find(sum(APFilter)));
                MaxAPIndex=size(APFilter{ChN},2);%max(find(sum(APFilter)));

                if ROI

                    %Get the corresponding mean information (ROI, CompiledParticles_ROI)
                    k=1;
                    for i=MinAPIndex:MaxAPIndex
                        [MeanVectorAPTemp_ROI,SDVectorAPTemp_ROI,NParticlesAPTemp_ROI]=AverageTraces(FrameInfo,...
                            CompiledParticles_ROI{ChN}(APFilter_ROI{ChN}(:,i)));
                        MeanVectorAPCell_ROI{k}=MeanVectorAPTemp_ROI';
                        SDVectorAPCell_ROI{k}=SDVectorAPTemp_ROI';
                        NParticlesAPCell_ROI{k}=NParticlesAPTemp_ROI';
                        k=k+1;
                    end
                    MeanVectorAP_ROI{ChN}=cell2mat(MeanVectorAPCell_ROI);
                    SDVectorAP_ROI{ChN}=cell2mat(SDVectorAPCell_ROI);
                    NParticlesAP_ROI{ChN}=cell2mat(NParticlesAPCell_ROI);

                    % Get the corresponding mean information
                    %(nonROI, CompiledParticles_nonROI -> save all in MeanVectorAP_nonROI
                    k=1;
                    for i=MinAPIndex:MaxAPIndex
                        [MeanVectorAPTemp_nonROI,SDVectorAPTemp_nonROI,NParticlesAPTemp_nonROI]=AverageTraces(FrameInfo,...
                            CompiledParticles_nonROI{ChN}(APFilter_nonROI{ChN}(:,i)));
                        MeanVectorAPCell_nonROI{k}=MeanVectorAPTemp_nonROI';
                        SDVectorAPCell_nonROI{k}=SDVectorAPTemp_nonROI';
                        NParticlesAPCell_nonROI{k}=NParticlesAPTemp_nonROI';
                        k=k+1;
                    end
                    MeanVectorAP_nonROI{ChN}=cell2mat(MeanVectorAPCell_nonROI);
                    SDVectorAP_nonROI{ChN}=cell2mat(SDVectorAPCell_nonROI);
                    NParticlesAP_nonROI{ChN}=cell2mat(NParticlesAPCell_nonROI);
    %%                
                    %AP Means

                    % Get the mean information for all of the CompiledParticles
                    % (Save this in MeanVectorAP)
                    k=1;
                    for i=MinAPIndex:MaxAPIndex
                        [MeanVectorAPTemp,SDVectorAPTemp,NParticlesAPTemp]=AverageTraces(FrameInfo,...
                            CompiledParticles{ChN}(APFilter{ChN}(:,i)));
                        MeanVectorAPCell{k}=MeanVectorAPTemp';
                        SDVectorAPCell{k}=SDVectorAPTemp';
                        NParticlesAPCell{k}=NParticlesAPTemp';
                        k=k+1;
                    end
                    MeanVectorAP{ChN}=cell2mat(MeanVectorAPCell);
                    SDVectorAP{ChN}=cell2mat(SDVectorAPCell);
                    NParticlesAP{ChN}=cell2mat(NParticlesAPCell);

                    %Calculate the mean for only anterior particles
                    try
                        MeanVectorAPAnterior{ChN} = MeanVectorAP{ChN}(:,5:15); %Only average particles within window of 10% to 35% w/ 2.5% AP resolution. P2P expression is relatively flat here.
                        MeanVectorAnterior{ChN} = nanmean(MeanVectorAPAnterior{ChN},2);
                    catch
                        %That didn't work
                    end



                else % This is the case which we don't use ROI option
                    %Get the corresponding mean information
                    k=1;
                    for i=MinAPIndex:MaxAPIndex
                        [MeanVectorAPTemp,SDVectorAPTemp,NParticlesAPTemp]=AverageTraces(FrameInfo,...
                            CompiledParticles{ChN}(APFilter{ChN}(:,i)));
                        MeanVectorAPCell{k}=MeanVectorAPTemp';
                        SDVectorAPCell{k}=SDVectorAPTemp';
                        NParticlesAPCell{k}=NParticlesAPTemp';
                        k=k+1;
                    end
                    MeanVectorAP{ChN}=cell2mat(MeanVectorAPCell);
                    SDVectorAP{ChN}=cell2mat(SDVectorAPCell);
                    NParticlesAP{ChN}=cell2mat(NParticlesAPCell);

                    %Calculate the mean for only anterior particles
                    try
                        MeanVectorAPAnterior{ChN} = MeanVectorAP{ChN}(:,5:15); %Only average particles within window of 10% to 35% w/ 2.5% AP resolution. P2P expression is relatively flat here.
                        MeanVectorAnterior{ChN} = nanmean(MeanVectorAPAnterior{ChN},2);
                    catch
                        %That didn't work
                    end

                end
            end
    %%        
            %DV ROI and non-ROI means

            if strcmpi(ExperimentAxis,'DV') 
                %Mean plot for different  positions

                %Figure out the DV range to use
                MinDVIndex=1;
                MaxDVIndex=size(DVFilter{ChN},2);

                if ROI

                    %Get the corresponding mean information (ROI, CompiledParticles_ROI)
                    k=1;
                    for i=MinDVIndex:MaxDVIndex
                        [MeanVectorDVTemp_ROI,SDVectorDVTemp_ROI,NParticlesDVTemp_ROI]=AverageTraces(FrameInfo,...
                            CompiledParticles_ROI{ChN}(DVFilter_ROI{ChN}(:,i)));
                        MeanVectorDVCell_ROI{k}=MeanVectorDVTemp_ROI';
                        SDVectorDVCell_ROI{k}=SDVectorDVTemp_ROI';
                        NParticlesDVCell_ROI{k}=NParticlesDVTemp_ROI';
                        k=k+1;
                    end
                    MeanVectorDV_ROI{ChN}=cell2mat(MeanVectorDVCell_ROI);
                    SDVectorDV_ROI{ChN}=cell2mat(SDVectorDVCell_ROI);
                    NParticlesDV_ROI{ChN}=cell2mat(NParticlesDVCell_ROI);

                    % Get the corresponding mean information 
                    %(nonROI, CompiledParticles_nonROI -> save all in MeanVectorDV_nonROI
                    k=1;
                    for i=MinDVIndex:MaxDVIndex
                        [MeanVectorDVTemp_nonROI,SDVectorDVTemp_nonROI,NParticlesDVTemp_nonROI]=AverageTraces(FrameInfo,...
                            CompiledParticles_nonROI{ChN}(DVFilter_nonROI{ChN}(:,i)));
                        MeanVectorDVCell_nonROI{k}=MeanVectorDVTemp_nonROI';
                        SDVectorDVCell_nonROI{k}=SDVectorDVTemp_nonROI';
                        NParticlesDVCell_nonROI{k}=NParticlesDVTemp_nonROI';
                        k=k+1;
                    end
                    MeanVectorDV_nonROI=cell2mat(MeanVectorDVCell_nonROI);
                    SDVectorDV_nonROI=cell2mat(SDVectorDVCell_nonROI);
                    NParticlesDV_nonROI=cell2mat(NParticlesDVCell_nonROI);

                    % Get the mean information for all of the CompiledParticles
                    % (Save this in MeanVectorDV)
                    k=1;
                    for i=MinDVIndex:MaxDVIndex
                        [MeanVectorDVTemp,SDVectorDVTemp,NParticlesDVTemp]=AverageTraces(FrameInfo,...
                            CompiledParticles{ChN}(DVFilter{ChN}(:,i)));
                        MeanVectorDVCell{k}=MeanVectorDVTemp';
                        SDVectorDVCell{k}=SDVectorDVTemp';
                        NParticlesDVCell{k}=NParticlesDVTemp';
                        k=k+1;
                    end
                    MeanVectorDV{ChN}=cell2mat(MeanVectorDVCell);
                    SDVectorDV{ChN}=cell2mat(SDVectorDVCell);
                    NParticlesDV{ChN}=cell2mat(NParticlesDVCell);


                else % This is the case which we don't use ROI option
                %Get the corresponding mean information
                k=1;
                for i=MinDVIndex:MaxDVIndex
                    [MeanVectorDVTemp,SDVectorDVTemp,NParticlesDVTemp]=AverageTraces(FrameInfo,...
                        CompiledParticles{ChN}(DVFilter{ChN}(:,i)));
                    MeanVectorDVCell{k}=MeanVectorDVTemp';
                    SDVectorDVCell{k}=SDVectorDVTemp';
                    NParticlesDVCell{k}=NParticlesDVTemp';
                    k=k+1;
                end
                MeanVectorDV{ChN}=cell2mat(MeanVectorDVCell);
                SDVectorDV{ChN}=cell2mat(SDVectorDVCell);
                NParticlesDV{ChN}=cell2mat(NParticlesDVCell);

                end
            end

    %%       
            %Calculate the mean for all AP bins
            [MeanVectorAll{ChN},SDVectorAll{ChN},NParticlesAll{ChN}]=AverageTraces(FrameInfo,CompiledParticles{ChN});

            %Now find the different maxima in each nc

            MaxFrame{ChN}=[];
            for i=1:length(NewCyclePos)
                if i==1
                    [~,MaxIndex]=max(MeanVectorAll{ChN}(1:NewCyclePos(1)));
                    MaxFrame{ChN}=[MaxFrame{ChN},MaxIndex];
                elseif i<=length(NewCyclePos)
                    [~,MaxIndex]=max(MeanVectorAll{ChN}(NewCyclePos(i-1):NewCyclePos(i)));
                    MaxFrame{ChN}=[MaxFrame{ChN},NewCyclePos(i-1)+MaxIndex-1];
                end
            end
            [~,MaxIndex]=max(MeanVectorAll{ChN}(NewCyclePos(i):end));
            if ~isempty(NewCyclePos)        %Why is this empty sometimes?
                %I think this only occurs with suboptimal
                %data
                MaxFrame{ChN}=[MaxFrame{ChN},NewCyclePos(i)+MaxIndex-1];
            end
        end
    end
end

