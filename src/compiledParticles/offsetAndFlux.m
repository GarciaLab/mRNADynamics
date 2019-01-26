function [MeanOffsetVector, SDOffsetVector, NOffsetParticles] = offsetAndFlux(NChannels, ...
    SkipFluctuations, ncFilter, ElapsedTime, CompiledParticles, DropboxFolder, ...
    Prefix, ExperimentAxis, intArea, MeanVectorAll, SDVectorAll, MaxFrame, numFrames)
%OFFSETANDFLUX Summary of this function goes here
%   Detailed explanation goes here
if NChannels==1
    
    %Is there a correlation between the fluctuations coming from the offset and
    %those observed in the traces? In order to figure this out I'll fit the
    %nc13 and nc14 intensity data with splines and compute the deviations with
    %respect to them. I'll look into different versions of the data such as
    %with and without the offset subtracted
    
    
    if ~SkipFluctuations & ~isempty(ncFilter)
        
        FilteredParticles=find(ncFilter(:,end)|ncFilter(:,end-1));
        
        OffsetFluct=[];
        DataRawFluct=[];
        %DataOldFluct=[];
        DataSplineFluct=[];
        
        for j=1:length(FilteredParticles)
            
            try
                %Deviation from offset with respect to spline
                optFit = adaptiveSplineFit(double([ElapsedTime(CompiledParticles{1}(FilteredParticles(j)).Frame)]),...
                    double([CompiledParticles{1}(FilteredParticles(j)).Off*intArea]),5);
                SplineValues=ppval(optFit,double([ElapsedTime(CompiledParticles{1}(FilteredParticles(j)).Frame)]));
                
                
                
                %Deviation of the raw data, without background subtraction, with
                %respect to a spline.
                DataFitRaw = adaptiveSplineFit(double([ElapsedTime(CompiledParticles{1}(FilteredParticles(j)).Frame)]),...
                    double(CompiledParticles{1}(FilteredParticles(j)).FluoRaw),10);
                DataFitRawValues=ppval(DataFitRaw,double([ElapsedTime(CompiledParticles{1}(FilteredParticles(j)).Frame)]));
                
                
                %                 %Deviation of the raw data minus the actual offset with respect to a
                %                 %spline
                %                 DataFitOld = adaptiveSplineFit(double([ElapsedTime(CompiledParticles{1}(FilteredParticles(j)).Frame)]),...
                %                     double(CompiledParticles{1}(FilteredParticles(j)).FluoOld),10);
                %                 DataSplineValuesOld=ppval(DataFitOld,double([ElapsedTime(CompiledParticles{1}(FilteredParticles(j)).Frame)]));
                
                
                
                %Deviation of the raw data minues the spline offset
                DataFit = adaptiveSplineFit(double([ElapsedTime(CompiledParticles{1}(FilteredParticles(j)).Frame)]),...
                    double(CompiledParticles{1}(FilteredParticles(j)).Fluo),10);
                DataSplineValues=ppval(DataFit,double([ElapsedTime(CompiledParticles{1}(FilteredParticles(j)).Frame)]));
                
                
                
                %Put all the data together for the plot
                OffsetFluct=[OffsetFluct,CompiledParticles{1}(FilteredParticles(j)).Off*intArea-SplineValues];
                DataRawFluct=[DataRawFluct,double(CompiledParticles{1}(FilteredParticles(j)).FluoRaw)-DataFitRawValues];
                %DataOldFluct=[DataOldFluct,double(CompiledParticles{1}(FilteredParticles(j)).FluoOld)-DataSplineValuesOld];
                DataSplineFluct=[DataSplineFluct,double(CompiledParticles{1}(FilteredParticles(j)).Fluo)-DataSplineValues];
            end
        end
        
        lim = 4500; %AR 1/12/18 why 4500?
        xRange=linspace(-lim,lim);
        
        figure(4)
        plot(OffsetFluct,DataRawFluct,'.k')
        xlabel('Offset fluctuation')
        ylabel('Fluctuations in raw data')
        axis square
        xlim([-lim,lim])
        ylim([-lim,lim])
        R = corrcoef(OffsetFluct,DataRawFluct);
        title(['Correlation: ',num2str(R(2,1))])
        saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'TracesFluctuations',filesep,'Fluct-OffsetVsRawData.tif'])
        
        %         figure(5)
        %         plot(OffsetFluct,DataOldFluct,'.k')
        %         xlabel('Offset fluctuation')
        %         ylabel('Fluctuations with instantaneous offset subtraction')
        %         axis square
        %         xlim([-lim,lim])
        %         ylim([-lim,lim])
        %         R = corrcoef(OffsetFluct,DataOldFluct)
        %         title(['Correlation: ',num2str(R(2,1))])
        %         saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'TracesFluctuations',filesep,'Fluct-OffsetVsInstData.tif'])
        
        
        
        figure(6)
        plot(OffsetFluct,DataSplineFluct,'.k')
        xlabel('Offset fluctuation')
        ylabel('Fluctuations with spline offset subtraction')
        axis square
        xlim([-lim,lim])
        ylim([-lim,lim])
        R = corrcoef(OffsetFluct,DataSplineFluct);
        title(['Correlation: ',num2str(R(2,1))])
        saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'TracesFluctuations',filesep,'Fluct-OffsetVsSplineData.tif'])
    end
    
    
    
    %Look at the offset of each particle. Do they all look the same? This is
    %only for AP for now
    if ~isempty(CompiledParticles{1})
        if strcmpi(ExperimentAxis,'AP')
            figure(7)
            subplot(1,2,1)
            MaxAP=0;
            MinAP=inf;
            hold all
            if ~isempty(MaxFrame{1})
                for i=1:length(CompiledParticles{1})
                    if sum(CompiledParticles{1}(i).Frame==MaxFrame{1}(end-1))
                        MaxAP=max([CompiledParticles{1}(i).MeanAP,MaxAP]);
                        MinAP=min([CompiledParticles{1}(i).MeanAP,MinAP]);
                        FramePos=find(CompiledParticles{1}(i).Frame==MaxFrame{1}(end-1));
                        plot(CompiledParticles{1}(i).MeanAP,CompiledParticles{1}(i).Off(FramePos),'.k')
                    end   
                end
            else
                warning('MaxFrame is empty. Unable to check offset of particles')
            end 
            hold off
            if MinAP<MaxAP
                xlim([MinAP*0.8,MaxAP*1.2])
            end
            xlabel('AP position')
            ylabel('Offset fluorescence')
            title('Offset at maximum in nc13')
            axis square
            
            subplot(1,2,2)
            MaxAP=0;
            MinAP=inf;
            hold all
            if ~isempty(MaxFrame{1})
                for i=1:length(CompiledParticles)
                    if sum(CompiledParticles{1}(i).Frame==MaxFrame{1}(end))
                        MaxAP=max([CompiledParticles{1}(i).MeanAP,MaxAP]);
                        MinAP=min([CompiledParticles{1}(i).MeanAP,MinAP]);
                        FramePos=find(CompiledParticles{1}(i).Frame==MaxFrame{1}(end));
                        plot(CompiledParticles{1}(i).MeanAP,CompiledParticles{1}(i).Off(FramePos),'.k')
                    end
                end
            else
                warning('MaxFrame is empty. Unable to check offset of particles')
            end
            hold off
            if MinAP < Inf && MaxAP > 0
                % ES 2014-09-12: This change is for cases in which no spots were
                % detected during the final time point.
                xlim([MinAP*0.8,MaxAP*1.2]);
            end
            xlabel('AP position')
            ylabel('Offset fluorescence')
            title('Offset at maximum in nc14')
            axis square
            saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'Offset',filesep,'OffsetVsAP.tif'])
        end
        
        %Average over all time points
        
        %Average and SD over each time point. In order to do this we'll generate a
        %cell array with all the values for a given time point
        
        OffsetCell=cell(numFrames,1);
        
        
        for i=1:length(CompiledParticles{1})
            for j=1:length(CompiledParticles{1}(i).Frame)
                OffsetCell{CompiledParticles{1}(i).Frame(j)}=[OffsetCell{CompiledParticles{1}(i).Frame(j)},...
                    CompiledParticles{1}(i).Off(j)];
            end
        end
        
        
        MeanOffsetTrace=cellfun(@nanmean,OffsetCell,'UniformOutput',false);
        SDOffsetTrace=cellfun(@nanstd,OffsetCell,'UniformOutput',false);
        NParticlesOffsetTrace=cellfun(@length,OffsetCell,'UniformOutput',false);
        
        
        MeanOffsetVector=[MeanOffsetTrace{:}];
        SDOffsetVector=[SDOffsetTrace{:}];
        NOffsetParticles=[NParticlesOffsetTrace{:}];
        
        
        if strcmpi(ExperimentAxis,'AP')
            figure(8)
            intArea=109;  %109 pixels is the default area when the pixels are assumed to be 212nm x 212 nm AR 9/3/18
            errorbar(1:length(MeanOffsetVector),MeanOffsetVector*intArea,...
                SDOffsetVector*intArea,'.-r')
            hold on
            errorbar(1:length(MeanVectorAll{1}),MeanVectorAll{1},...
                SDVectorAll{1},'.-k')
            hold off
            xlabel('Frame')
            ylabel('Fluorescence (au)')
            xlim([0,length(MeanOffsetVector)*1.1])
            legend('Offset','Spots','Location','Best')
            saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'Offset',filesep,'OffsetAndFluoTime.tif'])
            
            
            figure(9)
            errorbar(1:length(MeanOffsetVector),MeanOffsetVector*intArea-min(MeanOffsetVector*intArea)+...
                min(MeanVectorAll{1}),...
                SDOffsetVector*intArea,'.-r')
            hold on
            errorbar(1:length(MeanVectorAll{1}),MeanVectorAll{1},...
                SDVectorAll{1},'.-k')
            hold off
            xlabel('Frame')
            ylabel('Fluorescence (au)')
            xlim([0,length(MeanOffsetVector)*1.1])
            legend('Offset (displaced)','Spots','Location','Best')
            axis square
            saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'Offset',filesep,'OffsetAndFluoTime-Displaced.tif'])
        end
    else
        MeanOffsetVector=[];
        SDOffsetVector=[];
        NOffsetParticles=[];
    end
else
    %This is just to make 2spot2color work by skipping the calculations for
    %this function (JL 1/23/19).
    MeanOffsetVector=[];
    SDOffsetVector=[];
    NOffsetParticles=[];
end
end

