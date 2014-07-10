function [TimeNC14,MeanNC14,SENC14,MeanOnRatioNC14,SEOnRatioNC14,...
    IntegralNC14,SEIntegralNC14,IntegralONNC14,SEIntegralONNC14,...
    TimeNC13,MeanNC13,SENC13,MeanOnRatioNC13,SEOnRatioNC13,...
    IntegralNC13,SEIntegralNC13,...
    IntegralDegNC14,SEIntegralDegNC14,IntegralDegONNC14,SEIntegralDegONNC14,...
    IntegralDegNC13,SEIntegralDegNC13,IntegralDegONNC13,SEIntegralDegONNC13]=...
    PlotMeanAPMovies(Data,Label)

%This function plots the mean AP profile of a data set for the maximum in
%nc13 and nc14 and generates a movie as well.

%IntegralNC14 is the total amount produced per ALL nuclei
%IntegralONNC14 is the total amount produced per ON nuclei 

%% Parameters and data sets

%Parameters:
MinParticles=3;     %Minimum number of particles necessary in a bin for
                    %it to be considered
MinTimePoints=5;    %Minimum number of time points for the interpolation
MinEmbryos=2;       %Minimum number of embryos
           
ElongationRate=1.54;   %In kb/min


if findstr(lower(Data(1).SetName),'eve')
    HalfLifeForPlot=7;  %Half-life to extract from the calculated mRNA accumulation
                        %by Jacques
    %What experiment are we dealing with? This is useful to figure out the
    %elongation time, for example.
    ElongationTime=6.443/ElongationRate;
elseif findstr(lower(Data(1).SetName),'hbbac')
    HalfLifeForPlot=1000;    %Infinite half-life for now
    ElongationTime=6.443/ElongationRate;
else
    error('Add information for this construct')
end



                    
%Some labels to use for plots
Labels='k.r.g.b.y.c.m.ksrsgsbsyscsmskorogoboyocomok^r^g^b^y^c^m^';
                    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,SchnitzcellsFolder]=...
    DetermineLocalFolders(Data(1).SetName(11:end-1));






%Create the folders we need and empty them if necessary

%Folder for report figures
mkdir([DropboxFolder,'\Reports'])
mkdir([DropboxFolder,'\Reports\APPlots\APMovies'])
mkdir([DropboxFolder,'\Reports\APPlots\APMovies\APMovie',Label])
mkdir([DropboxFolder,'\Reports\APPlots\APMovies\APMovie',Label,'SE'])
mkdir([DropboxFolder,'\Reports\APPlots\APMovies\APMovieMean',Label])


%Make a movie of all the different embryos as a function of time and AP.
%We'll set a new zero of time at the beginning of each nc.

%Determine the maximum fluorescence value for the whole movie
MaxValue=max(cellfun(@(x) max(max(x)),{Data.MeanVectorAP}));

%Different movies will have different time resolutions. We'll pick the
%smallest time intervals and interpolate the values of all the data sets
%accordingly. We are only doing this for nc13 and nc14.
for i=1:length(Data)
    nc13Resolution(i)=mean(diff(Data(i).ElapsedTime(Data(i).nc13:Data(i).nc14)));
    nc14Resolution(i)=mean(diff(Data(i).ElapsedTime(Data(i).nc14:end)));
end
nc13Resolution=min(nc13Resolution);
nc14Resolution=min(nc14Resolution);

%Figure out the time window nc14
for i=1:length(Data)
    nc14Length(i)=Data(i).ElapsedTime(end)-Data(i).ElapsedTime(Data(i).nc14);
end
nc14Length=min(nc14Length);
%All time points for nc14
TimeNC14=0:nc14Resolution:nc14Length;


%Figure out the time window nc13
for i=1:length(Data)
    nc13Length(i)=Data(i).ElapsedTime(Data(i).nc14)-Data(i).ElapsedTime(Data(i).nc13);
end
nc13Length=min(nc13Length);
%All time points for nc13
TimeNC13=0:nc13Resolution:nc13Length;



for j=1:length(Data)

    %Get the AP vector related to each nc
    
    %Filter the bins that have at least the minimum number the particles
    NParticlesAPFilter=Data(j).NParticlesAP>=MinParticles;
    NParticlesAPFilterNC13=NParticlesAPFilter(Data(j).nc13:Data(j).nc14,:);
    NParticlesAPFilterNC14=NParticlesAPFilter(Data(j).nc14:end,:);
    
    
    %Mean levels
    MeanVectorAPNC13{j}=Data(j).MeanVectorAP(Data(j).nc13:Data(j).nc14,:);
    MeanVectorAPNC14{j}=Data(j).MeanVectorAP(Data(j).nc14:end,:);
    
    SDVectorAPNC13{j}=Data(j).SDVectorAP(Data(j).nc13:Data(j).nc14,:);
    SDVectorAPNC14{j}=Data(j).SDVectorAP(Data(j).nc14:end,:);
    
    SEVectorAPNC13{j}=Data(j).SDVectorAP(Data(j).nc13:Data(j).nc14,:)./...
        sqrt(Data(j).NParticlesAP(Data(j).nc13:Data(j).nc14,:));
    SEVectorAPNC14{j}=Data(j).SDVectorAP(Data(j).nc14:end,:)./...
        sqrt(Data(j).NParticlesAP(Data(j).nc14:end,:));
       
    MeanVectorAPNC13{j}(~NParticlesAPFilterNC13)=nan;
    MeanVectorAPNC14{j}(~NParticlesAPFilterNC14)=nan;
    
    SDVectorAPNC13{j}(~NParticlesAPFilterNC13)=nan;
    SDVectorAPNC14{j}(~NParticlesAPFilterNC14)=nan;
    
    SEVectorAPNC13{j}(~NParticlesAPFilterNC13)=nan;
    SEVectorAPNC14{j}(~NParticlesAPFilterNC14)=nan;
    
    
    %New version of integration. This is using the trapezoidal integration
    %now
    %nc14
    m=1;
    for i=Data(j).nc14:length(Data(j).ElapsedTime)
        [TotalVectorAPNC14{j}(m,:),Dummy,TotalVectorONAPNC14{j}(m,:),Dummy]=...
            IntegrateIndivSetTimeWindow(Data(j),MinParticles,...
            14,[0,Data(j).ElapsedTime(i)-Data(j).ElapsedTime(Data(j).nc14)]);
        m=m+1;
    end
    
    
    %nc13
    m=1;
    for i=Data(j).nc13:Data(j).nc14
        [TotalVectorAPNC13{j}(m,:),Dummy,TotalVectorONAPNC13{j}(m,:),Dummy]=...
            IntegrateIndivSetTimeWindow(Data(j),MinParticles,...
            13,[0,Data(j).ElapsedTime(i)-Data(j).ElapsedTime(Data(j).nc13)]);
        m=m+1;
    end
    
    
    %If the mRNA accumulation accounting for life times exists then load it
	%and process it
    if isfield(Data,'AccumulationData')
        
        %Get the grames
        FramesNC13=Data(j).nc13:Data(j).nc14;
        FramesNC14=Data(j).nc14:length(Data(j).ElapsedTime);
        
        %Initialize the vectors
        AccumNC13{j}=nan(length(FramesNC13),length(Data(j).APbinID));
        AccumNC13ON{j}=nan(length(FramesNC13),length(Data(j).APbinID));
        AccumNC14{j}=nan(length(FramesNC13),length(Data(j).APbinID));
        AccumNC14ON{j}=nan(length(FramesNC13),length(Data(j).APbinID));
        
            
        %Find the right accumulation given the life time
        IndexHalfLife=find([Data(j).AccumulationData.halflife]==HalfLifeForPlot);
        
        if isempty(IndexHalfLife)
            error(['Half life not found for set ',Data(j).SetName])
        end
        
        FilterNC13=Data(j).ncFilter(:,end-1);       
        FilterNC14=Data(j).ncFilter(:,end);  
        
        %Extract the accumulated amount per ALL nuclei
        %nc13
        for k=1:length(FramesNC13)
            for n=1:length(Data(j).APbinID)
                
                %Which particles are in this NC and AP bin?
                ParticlesToCheck=find(FilterNC13&(Data(j).APFilter(:,n)));
                
                %Get the particles in this frame and AP bin
                IntegralTemp=[];
                
                for m=1:length(ParticlesToCheck)
                    if sum(Data(j).AccumulationData(IndexHalfLife).CompiledParticles(ParticlesToCheck(m)).Frames==FramesNC13(k))
                        FrameIndex=...
                            find(Data(j).AccumulationData(IndexHalfLife).CompiledParticles(ParticlesToCheck(m)).Frames==...
                            FramesNC13(k));
                        IntegralTemp=...
                            [IntegralTemp,Data(j).AccumulationData(IndexHalfLife).CompiledParticles(m).mRNAacum(FrameIndex)];
                    end
                end
                
                
                
                if length(IntegralTemp)>=MinParticles
                    
                    AccumNC13{j}(k,n)=mean(IntegralTemp)*Data(j).EllipsesOnAP(n,end-1)/Data(j).TotalEllipsesAP(n,end-1);
                    AccumNC13ON{j}(k,n)=mean(IntegralTemp);
                end
                
            end
            
            
            
        end
        
        
        %nc14
        for k=1:length(FramesNC14)
            for n=1:length(Data(j).APbinID)
                
                %Which particles are in this NC and AP bin?
                ParticlesToCheck=find(FilterNC14&(Data(j).APFilter(:,n)));
                
                
                %Get the particles in this frame and AP bin
                IntegralTemp=[];
                for m=1:length(ParticlesToCheck)
                    if sum(Data(j).AccumulationData(IndexHalfLife).CompiledParticles(ParticlesToCheck(m)).Frames==FramesNC14(k))
                        FrameIndex=...
                            find(Data(j).AccumulationData(IndexHalfLife).CompiledParticles(ParticlesToCheck(m)).Frames==...
                            FramesNC14(k));
                        IntegralTemp=...
                            [IntegralTemp,Data(j).AccumulationData(IndexHalfLife).CompiledParticles(ParticlesToCheck(m)).mRNAacum(FrameIndex)];
                    end
                end
                
                if length(IntegralTemp)>=MinParticles
                    AccumNC14{j}(k,n)=mean(IntegralTemp)*Data(j).EllipsesOnAP(n,end)/Data(j).TotalEllipsesAP(n,end);
                    AccumNC14ON{j}(k,n)=mean(IntegralTemp);
                end
                
            end
            
            
            
        end
    end   
          

    %Fraction of active nuclei vs. time
    OnRatioAPNC14{j}=Data(j).OnRatioAP(Data(j).nc14:end,:);
    OnRatioAPNC13{j}=Data(j).OnRatioAP(Data(j).nc13:Data(j).nc14,:);
    
    
    %Interpolate each AP bin as a function of time
    MeanVectorAPNC14Interp{j}=nan(length(TimeNC14),length(Data(j).APbinID));
    MeanOnRatioAPNC14Interp{j}=nan(length(TimeNC14),length(Data(j).APbinID));
    TotalVectorAPNC14Interp{j}=nan(length(TimeNC14),length(Data(j).APbinID));
    TotalVectorONAPNC14Interp{j}=nan(length(TimeNC14),length(Data(j).APbinID));
    AccumVectorAPNC14Interp{j}=nan(length(TimeNC14),length(Data(j).APbinID));
    AccumVectorONAPNC14Interp{j}=nan(length(TimeNC14),length(Data(j).APbinID));
    
    MeanVectorAPNC13Interp{j}=nan(length(TimeNC13),length(Data(j).APbinID));
    MeanOnRatioAPNC13Interp{j}=nan(length(TimeNC13),length(Data(j).APbinID));
    TotalVectorAPNC13Interp{j}=nan(length(TimeNC13),length(Data(j).APbinID));
    TotalVectorONAPNC13Interp{j}=nan(length(TimeNC13),length(Data(j).APbinID));
    AccumVectorAPNC13Interp{j}=nan(length(TimeNC13),length(Data(j).APbinID));
    AccumVectorONAPNC13Interp{j}=nan(length(TimeNC13),length(Data(j).APbinID));
    

    for i=1:length(Data(j).APbinID)
        
        if sum(~isnan(MeanVectorAPNC13{j}(:,i)))>1

            TimeWindow=Data(j).ElapsedTime(Data(j).nc13:Data(j).nc14)-...
                Data(j).ElapsedTime(Data(j).nc13);
            
            %Use only data points where we had MinParticles to calculate
            %the mean
            FilterTemp=NParticlesAPFilter(Data(j).nc13:Data(j).nc14,i);
            
            if sum(FilterTemp)>=MinTimePoints
            
                MeanVectorAPNC13Interp{j}(:,i)=pchip(TimeWindow(FilterTemp),...
                    MeanVectorAPNC13{j}(FilterTemp,i)',...
                    TimeNC13);

                
                %Check if there are enough statistics
                if sum(~isnan(OnRatioAPNC13{j}(FilterTemp,i)))>=MinParticles
                    MeanOnRatioAPNC13Interp{j}(:,i)=pchip(TimeWindow(FilterTemp),...
                        OnRatioAPNC13{j}(FilterTemp,i)',...
                        TimeNC13);
                else
                    MeanOnRatioAPNC13Interp{j}(:,i)=nan;
                end

                if sum(~isnan(TotalVectorAPNC13{j}(FilterTemp,i)))>=MinParticles
                    TotalVectorAPNC13Interp{j}(:,i)=pchip(TimeWindow(FilterTemp),...
                        TotalVectorAPNC13{j}(FilterTemp,i)',...
                        TimeNC13);
                else
                    TotalVectorAPNC13Interp{j}(:,i)=nan;
                end
                
                if sum(~isnan(TotalVectorONAPNC13{j}(FilterTemp,i)))>=MinParticles
                    TotalVectorONAPNC13Interp{j}(:,i)=pchip(TimeWindow(FilterTemp),...
                        TotalVectorONAPNC13{j}(FilterTemp,i)',...
                        TimeNC13);
                else
                    TotalVectorONAPNC13Interp{j}(:,i)=nan;
                end
                
                if sum(~isnan(AccumNC13{j}(FilterTemp,i)))>=MinParticles
                    AccumVectorAPNC13Interp{j}(:,i)=pchip(TimeWindow(FilterTemp),...
                        AccumNC13{j}(FilterTemp,i)',...
                        TimeNC13);
                else
                    AccumVectorAPNC13Interp{j}(:,i)=nan;
                end
            
                
                if sum(~isnan(AccumNC13ON{j}(FilterTemp,i)))>=MinParticles
                    AccumVectorONAPNC13Interp{j}(:,i)=pchip(TimeWindow(FilterTemp),...
                        AccumNC13ON{j}(FilterTemp,i)',...
                        TimeNC13);
                else
                    AccumVectorONAPNC13Interp{j}(:,i)=nan;
                end
                
                
                MeanVectorAPNC13Interp{j}(TimeNC13<min(TimeWindow(FilterTemp)),i)=nan;
                MeanVectorAPNC13Interp{j}(TimeNC13>max(TimeWindow(FilterTemp)),i)=nan;
                
                MeanOnRatioAPNC13Interp{j}(TimeNC13<min(TimeWindow(FilterTemp)),i)=nan;
                MeanOnRatioAPNC13Interp{j}(TimeNC13>max(TimeWindow(FilterTemp)),i)=nan;
                
                %For the total amount produced I need 0 before the time
                %window and the maximum after the time window
                TotalVectorAPNC13Interp{j}(TimeNC13<min(TimeWindow(FilterTemp)),i)=0;
                TotalVectorAPNC13Interp{j}(TimeNC13>max(TimeWindow(FilterTemp)),i)=...
                    max(TotalVectorAPNC13Interp{j}(TimeNC13<max(TimeWindow(FilterTemp)),i));
                
                TotalVectorONAPNC13Interp{j}(TimeNC13<min(TimeWindow(FilterTemp)),i)=0;
                TotalVectorONAPNC13Interp{j}(TimeNC13>max(TimeWindow(FilterTemp)),i)=...
                    max(TotalVectorONAPNC13Interp{j}(TimeNC13<max(TimeWindow(FilterTemp)),i));
                
                AccumVectorAPNC13Interp{j}(TimeNC13<min(TimeWindow(FilterTemp)-ElongationTime),i)=0;
                AccumVectorAPNC13Interp{j}(TimeNC13>max(TimeWindow(FilterTemp)),i)=...
                    max(AccumVectorAPNC13Interp{j}(TimeNC13<max(TimeWindow(FilterTemp)),i));
                
                AccumVectorONAPNC13Interp{j}(TimeNC13<min(TimeWindow(FilterTemp)-ElongationTime),i)=0;
                AccumVectorONAPNC13Interp{j}(TimeNC13>max(TimeWindow(FilterTemp)),i)=...
                    max(AccumVectorONAPNC13Interp{j}(TimeNC13<max(TimeWindow(FilterTemp)),i));
                
            else
                MeanVectorAPNC13Interp{j}(:,i)=nan;
                MeanOnRatioAPNC13Interp{j}(:,i)=nan;
                TotalVectorAPNC13Interp{j}(:,i)=nan;
                TotalVectorONAPNC13Interp{j}(:,i)=nan;
                AccumVectorAPNC13Interp{j}(:,i)=nan;
                AccumVectorONAPNC13Interp{j}(:,i)=nan;
            end
        end
            
       
        
        if sum(~isnan(MeanVectorAPNC14{j}(:,i)))>1
            
            TimeWindow=Data(j).ElapsedTime(Data(j).nc14:end)-...
                Data(j).ElapsedTime(Data(j).nc14);
            
            %Use only data points where we had MinParticles to calculate
            %the mean
            FilterTemp=NParticlesAPFilter(Data(j).nc14:end,i);
            
            if sum(FilterTemp)>=MinTimePoints
            
                MeanVectorAPNC14Interp{j}(:,i)=pchip(TimeWindow(FilterTemp),...
                    MeanVectorAPNC14{j}(FilterTemp,i)',...
                    TimeNC14);
                
                MeanOnRatioAPNC14Interp{j}(:,i)=pchip(TimeWindow(FilterTemp),...
                    OnRatioAPNC14{j}(FilterTemp,i)',...
                    TimeNC14);

                TotalVectorAPNC14Interp{j}(:,i)=pchip(TimeWindow(FilterTemp),...
                    TotalVectorAPNC14{j}(FilterTemp,i)',...
                    TimeNC14);
                
                
                TotalVectorONAPNC14Interp{j}(:,i)=pchip(TimeWindow(FilterTemp),...
                    TotalVectorONAPNC14{j}(FilterTemp,i)',...
                    TimeNC14);
                
                AccumVectorAPNC14Interp{j}(:,i)=pchip(TimeWindow(FilterTemp),...
                    AccumNC14{j}(FilterTemp,i)',TimeNC14);
                AccumVectorONAPNC14Interp{j}(:,i)=pchip(TimeWindow(FilterTemp),...
                    AccumNC14ON{j}(FilterTemp,i)',TimeNC14);
                
                
                
            
                MeanVectorAPNC14Interp{j}(TimeNC14<min(TimeWindow(FilterTemp)),i)=nan;
                MeanVectorAPNC14Interp{j}(TimeNC14>max(TimeWindow(FilterTemp)),i)=nan;
                
                MeanOnRatioAPNC14Interp{j}(TimeNC14<min(TimeWindow(FilterTemp)),i)=nan;
                MeanOnRatioAPNC14Interp{j}(TimeNC14>max(TimeWindow(FilterTemp)),i)=nan;
                
                %For the total amount produced I need 0 before the time
                %window and the maximum after the time window
                TotalVectorAPNC14Interp{j}(TimeNC14<min(TimeWindow(FilterTemp)),i)=0;
                TotalVectorAPNC14Interp{j}(TimeNC14>max(TimeWindow(FilterTemp)),i)=...
                    max(TotalVectorAPNC14Interp{j}(TimeNC14<max(TimeWindow(FilterTemp)),i));
                
                TotalVectorONAPNC14Interp{j}(TimeNC14<min(TimeWindow(FilterTemp)),i)=0;
                TotalVectorONAPNC14Interp{j}(TimeNC14>max(TimeWindow(FilterTemp)),i)=...
                    max(TotalVectorONAPNC14Interp{j}(TimeNC14<max(TimeWindow(FilterTemp)),i));
                
                AccumVectorAPNC14Interp{j}(TimeNC14<min(TimeWindow(FilterTemp)+ElongationTime),i)=0;
                AccumVectorAPNC14Interp{j}(TimeNC14>max(TimeWindow(FilterTemp)),i)=...
                    max(AccumVectorAPNC14Interp{j}(TimeNC14<max(TimeWindow(FilterTemp)),i));
                
                AccumVectorONAPNC14Interp{j}(TimeNC14<min(TimeWindow(FilterTemp)+ElongationTime),i)=0;
                AccumVectorONAPNC14Interp{j}(TimeNC14>max(TimeWindow(FilterTemp)),i)=...
                    max(AccumVectorONAPNC14Interp{j}(TimeNC14<max(TimeWindow(FilterTemp)),i));

            else
                MeanVectorAPNC14Interp{j}(:,i)=nan;
                MeanOnRatioAPNC14Interp{j}(:,i)=nan;
                TotalVectorAPNC14Interp{j}(:,i)=nan;
                TotalVectorONAPNC14Interp{j}(:,i)=nan;
                AccumVectorAPNC14Interp{j}(:,i)=nan;
                AccumVectorONAPNC14Interp{j}(:,i)=nan;
            end
        end
    end
end



MovieFigure=figure;

%% Movie with SD
% 
% %Plot the resulting traces for nc13 - Using Standard Deviation
% [MaxTimeIndex,MaxTimeSet]=max(cellfun(@length,TimeNC13));
% 
% for i=1:MaxTimeIndex
%     LegendLabel={};
%     k=1;
%     clf
%     hold on
%     PlotHandle=[];
%     figure(MovieFigure)
%     for j=1:length(Data)
%         if i<=size(MeanVectorAPNC13{j},1)
%             LegendLabel{k}=num2str(j);
%             k=k+1;
%             PlotHandle=[PlotHandle,errorbar(Data(j).APbinID,MeanVectorAPNC13{j}(i,:),...
%                 SDVectorAPNC13{j}(i,:),...
%                 [Labels((j-1)*2+1:(j-1)*2+2),'-'],'MarkerFaceColor',Labels((j-1)*2+1))];
%         end
%     end
%     hold off
%     title(['nc13, ',num2str(TimeNC13{MaxTimeSet}(i)-TimeNC13{MaxTimeSet}(1)),' min'])
%     ylim([0,MaxValue])
%     xlim([0.1,0.8])
%     set(gca,'XTick',[0.1:0.1:0.8])
%     legend(PlotHandle,LegendLabel,'Location','NorthEast')
%     box on
%     %drawnow
%     StandardFigure(PlotHandle,gca)
%     saveas(gcf,[DropboxFolder,'\Reports\APPlots\APMovies\APMovie',Label,'\nc13-',iIndex(i,3),'.tif'])
% end
%     
% 
% %Plot the resulting traces for nc14  - Using Standard Deviation
% [MaxTimeIndex,MaxTimeSet]=max(cellfun(@length,TimeNC14));
% 
% for i=1:MaxTimeIndex
%     LegendLabel={};
%     k=1;
%     figure(MovieFigure)
%     clf
%     hold on
%     PlotHandle=[];
%     for j=1:length(Data)
%         if i<=size(MeanVectorAPNC14{j},1)
%             LegendLabel{k}=num2str(j);
%             k=k+1;
%             PlotHandle=[PlotHandle,errorbar(Data(j).APbinID,MeanVectorAPNC14{j}(i,:),...
%                 SDVectorAPNC14{j}(i,:),...
%                 [Labels((j-1)*2+1:(j-1)*2+2),'-'],'MarkerFaceColor',Labels((j-1)*2+1))];
%         end
%     end
%     hold off
%     title(['nc14, ',num2str(TimeNC14{MaxTimeSet}(i)-TimeNC14{MaxTimeSet}(1)),' min'])
%     ylim([0,MaxValue])
%     xlim([0.1,0.8])
%     set(gca,'XTick',[0.1:0.1:0.8])
%     legend(PlotHandle,LegendLabel,'Location','NorthEast')
%     box on
%     %drawnow
%     StandardFigure(PlotHandle,gca)
%     saveas(gcf,[DropboxFolder,'\Reports\APPlots\APMovies\APMovie',Label,'\nc14-',iIndex(i,3),'.tif'])
% end
%     


%% Movie with SE

% 
% %Plot the resulting traces for nc13 - Using Standard Error
% [MaxTimeIndex,MaxTimeSet]=max(cellfun(@length,TimeNC13));
% 
% for i=1:MaxTimeIndex
%     LegendLabel={};
%     k=1;
%     figure(MovieFigure)
%     clf
%     hold on
%     PlotHandle=[];
%     for j=1:length(Data)
%         if i<=size(MeanVectorAPNC13{j},1)
%             LegendLabel{k}=num2str(j);
%             k=k+1;
%             PlotHandle=[PlotHandle,errorbar(Data(j).APbinID,MeanVectorAPNC13{j}(i,:),...
%                 SEVectorAPNC13{j}(i,:),...
%                 [Labels((j-1)*2+1:(j-1)*2+2),'-'],'MarkerFaceColor',Labels((j-1)*2+1))];
%         end
%     end
%     hold off
%     title(['nc13, ',num2str(TimeNC13{MaxTimeSet}(i)-TimeNC13{MaxTimeSet}(1)),' min'])
%     ylim([0,MaxValue])
%     xlim([0.1,0.8])
%     set(gca,'XTick',[0.1:0.1:0.8])
%     legend(PlotHandle,LegendLabel,'Location','NorthEast')
%     box on
%     %drawnow
%     StandardFigure(PlotHandle,gca)
%     saveas(gcf,[DropboxFolder,'\Reports\APPlots\APMovies\APMovie',Label,'SE\nc13-',iIndex(i,3),'.tif'])
% end
% 
% 
% %Plot the resulting traces for nc14  - Using Standard Error
% [MaxTimeIndex,MaxTimeSet]=max(cellfun(@length,TimeNC14));
% 
% for i=1:MaxTimeIndex
%     LegendLabel={};
%     k=1;
%     figure(MovieFigure)
%     clf
%     hold on
%     PlotHandle=[];
%     for j=1:length(Data)
%         if i<=size(MeanVectorAPNC14{j},1)
%             LegendLabel{k}=num2str(j);
%             k=k+1;
%             PlotHandle=[PlotHandle,errorbar(Data(j).APbinID,MeanVectorAPNC14{j}(i,:),...
%                 SEVectorAPNC14{j}(i,:),...
%                 [Labels((j-1)*2+1:(j-1)*2+2),'-'],'MarkerFaceColor',Labels((j-1)*2+1))];
%         end
%     end
%     hold off
%     title(['nc14, ',num2str(TimeNC14{MaxTimeSet}(i)-TimeNC14{MaxTimeSet}(1)),' min'])
%     ylim([0,MaxValue])
%     xlim([0.1,0.8])
%     set(gca,'XTick',[0.1:0.1:0.8])
%     legend(PlotHandle,LegendLabel,'Location','NorthEast')
%     box on
%     %drawnow
%     StandardFigure(PlotHandle,gca)
%     saveas(gcf,[DropboxFolder,'\Reports\APPlots\APMovies\APMovie',Label,'SE\nc14-',iIndex(i,3),'.tif'])
% end
%     

%% Movie of averages

%nc13:

MaxTimeIndex=length(TimeNC13);

%Average fluorescence
for i=1:MaxTimeIndex

    figure(MovieFigure)
    clf
    
    
    %Populate a matrix with the values as a function of AP for this time
    %point.
    MeanTemp=[];
    for j=1:length(Data)
        if i<=size(MeanVectorAPNC13Interp{j},1)
            MeanTemp=[MeanTemp;MeanVectorAPNC13Interp{j}(i,:)];
        end
    end
    
    if ~isempty(MeanTemp)
        %Now calculale the mean and SD. We need to be careful with the Nans!
        for j=1:size(MeanTemp,2)
            if sum(~isnan(MeanTemp(:,j)))>=MinEmbryos
                MeanNC13(i,j)=mean(MeanTemp(~isnan(MeanTemp(:,j)),j));
                SDNC13(i,j)=std(MeanTemp(~isnan(MeanTemp(:,j)),j));
                SENC13(i,j)=std(MeanTemp(~isnan(MeanTemp(:,j)),j))/sqrt(sum(~isnan(MeanTemp(:,j))));
            else
                MeanNC13(i,j)=nan;
                SDNC13(i,j)=nan;
                SENC13(i,j)=nan;
            end
        end


       PlotHandle=errorbar(Data(1).APbinID,MeanNC13(i,:),...
            SENC13(i,:),'.-k');
        errorbar_tick(PlotHandle,0);
        hold off
        %title(['nc14, ',num2str(TimeNC14{MaxTimeSet}(i)-TimeNC14{MaxTimeSet}(1)),' min'])
        ylim([0,1.5E4])
        xlim([0.3,0.5])
        set(gca,'XTick',[0.1:0.1:0.8])
        box on
        pause(0.1)
        drawnow
        %StandardFigure(PlotHandle,gca)
        %saveas(gcf,[DropboxFolder,'\Reports\APPlots\APMovies\APMovieMean',Label,'\nc14-',iIndex(i,3),'.tif'])
    end
end


%Ratio of on nuclei
for i=1:MaxTimeIndex

    figure(MovieFigure)
    clf
    
    
    %Populate a matrix with the values as a function of AP for this time
    %point.
    MeanTemp=[];
    for j=1:length(Data)
        if i<=size(OnRatioAPNC13{j},1)
            MeanTemp=[MeanTemp;OnRatioAPNC13{j}(i,:)];
        end
    end
    
    if ~isempty(MeanTemp)
        %Now calculale the mean and SD. We need to be careful with the Nans!
        for j=1:size(MeanTemp,2)
            if sum(~isnan(MeanTemp(:,j)))>=MinEmbryos
                MeanOnRatioNC13(i,j)=mean(MeanTemp(~isnan(MeanTemp(:,j)),j));
                SEOnRatioNC13(i,j)=std(MeanTemp(~isnan(MeanTemp(:,j)),j))/sqrt(sum(~isnan(MeanTemp(:,j))));
            else
                MeanOnRatioNC13(i,j)=nan;
                SEOnRatioNC13(i,j)=nan;
            end
        end


       PlotHandle=errorbar(Data(1).APbinID,MeanOnRatioNC13(i,:),...
            SEOnRatioNC13(i,:),'.-k');
        errorbar_tick(PlotHandle,0);
        hold off
        %title(['nc14, ',num2str(TimeNC14{MaxTimeSet}(i)-TimeNC14{MaxTimeSet}(1)),' min'])
        ylim([0,1])
        xlim([0.3,0.5])
        set(gca,'XTick',[0.1:0.1:0.8])
        box on
        pause(0.1)
        drawnow
        %StandardFigure(PlotHandle,gca)
        %saveas(gcf,[DropboxFolder,'\Reports\APPlots\APMovies\APMovieMean',Label,'\nc14-',iIndex(i,3),'.tif'])
    end
end


%Plot the accumulation vs. time. This is the accumulation per ALL nuclei
MaxValue=max(cellfun(@(x) max(max(x)),TotalVectorAPNC13Interp));
for i=1:MaxTimeIndex
    figure(MovieFigure)
    clf
    
    
    MeanTemp=[];
    for j=1:length(Data)
        if i<=size(TotalVectorAPNC13Interp{j},1)
            MeanTemp=[MeanTemp;TotalVectorAPNC13Interp{j}(i,:)];
        end
    end
    
    if ~isempty(MeanTemp)
        %Now calculale the mean and SD. We need to be careful with the Nans!
        for j=1:size(MeanTemp,2)
            if sum(~isnan(MeanTemp(:,j)))>=MinEmbryos
                IntegralNC13(i,j)=mean(MeanTemp(~isnan(MeanTemp(:,j)),j));
                SDIntegralNC13(i,j)=std(MeanTemp(~isnan(MeanTemp(:,j)),j));
                SEIntegralNC13(i,j)=std(MeanTemp(~isnan(MeanTemp(:,j)),j))/sqrt(sum(~isnan(MeanTemp(:,j))));
            else
                IntegralNC13(i,j)=nan;
                SDIntegralNC13(i,j)=nan;
                SEIntegralNC13(i,j)=nan;
            end
        end


       PlotHandle=errorbar(Data(1).APbinID,IntegralNC13(i,:),...
            SEIntegralNC13(i,:),'.-k');
        errorbar_tick(PlotHandle,0);
        hold off
        %title(['nc14, ',num2str(TimeNC14{MaxTimeSet}(i)-TimeNC14{MaxTimeSet}(1)),' min'])
        ylim([0,MaxValue])
        xlim([0.3,0.5])
        
        set(gca,'XTick',[0.1:0.1:0.8])
        box on
        pause(0.1)
        drawnow
        %StandardFigure(PlotHandle,gca)
        %saveas(gcf,[DropboxFolder,'\Reports\APPlots\APMovies\APMovieMean',Label,'\nc14-',iIndex(i,3),'.tif'])
    end
end





%nc14:
MaxTimeIndex=length(TimeNC14);

%Average fluorescence
for i=1:MaxTimeIndex

    figure(MovieFigure)
    clf
    
    
    %Populate a matrix with the values as a function of AP for this time
    %point.
    MeanTemp=[];
    for j=1:length(Data)
        if i<=size(MeanVectorAPNC14Interp{j},1)
            MeanTemp=[MeanTemp;MeanVectorAPNC14Interp{j}(i,:)];
        end
    end
    
    if ~isempty(MeanTemp)
        %Now calculale the mean and SD. We need to be careful with the Nans!
        for j=1:size(MeanTemp,2)
            if sum(~isnan(MeanTemp(:,j)))>=MinEmbryos
                MeanNC14(i,j)=mean(MeanTemp(~isnan(MeanTemp(:,j)),j));
                SDNC14(i,j)=std(MeanTemp(~isnan(MeanTemp(:,j)),j));
                SENC14(i,j)=std(MeanTemp(~isnan(MeanTemp(:,j)),j))/sqrt(sum(~isnan(MeanTemp(:,j))));
            else
                MeanNC14(i,j)=nan;
                SDNC14(i,j)=nan;
                SENC14(i,j)=nan;
            end
        end


       PlotHandle=errorbar(Data(1).APbinID,MeanNC14(i,:),...
            SENC14(i,:),'.-k');
        errorbar_tick(PlotHandle,0);
        hold off
        %title(['nc14, ',num2str(TimeNC14{MaxTimeSet}(i)-TimeNC14{MaxTimeSet}(1)),' min'])
        ylim([0,1.5E4])
        xlim([0.3,0.5])
        set(gca,'XTick',[0.1:0.1:0.8])
        box on
        pause(0.1)
        drawnow
        %StandardFigure(PlotHandle,gca)
        %saveas(gcf,[DropboxFolder,'\Reports\APPlots\APMovies\APMovieMean',Label,'\nc14-',iIndex(i,3),'.tif'])
    end
end
    

%Ratio of on nuclei
for i=1:MaxTimeIndex

    figure(MovieFigure)
    clf
    
    
    %Populate a matrix with the values as a function of AP for this time
    %point.
    MeanTemp=[];
    for j=1:length(Data)
        if i<=size(OnRatioAPNC14{j},1)
            MeanTemp=[MeanTemp;OnRatioAPNC14{j}(i,:)];
        end
    end
    
    if ~isempty(MeanTemp)
        %Now calculale the mean and SD. We need to be careful with the Nans!
        for j=1:size(MeanTemp,2)
            if sum(~isnan(MeanTemp(:,j)))>=MinEmbryos
                MeanOnRatioNC14(i,j)=mean(MeanTemp(~isnan(MeanTemp(:,j)),j));
                SEOnRatioNC14(i,j)=std(MeanTemp(~isnan(MeanTemp(:,j)),j))/sqrt(sum(~isnan(MeanTemp(:,j))));
            else
                MeanOnRatioNC14(i,j)=nan;
                SEOnRatioNC14(i,j)=nan;
            end
        end


       PlotHandle=errorbar(Data(1).APbinID,MeanOnRatioNC14(i,:),...
            SEOnRatioNC14(i,:),'.-k');
        errorbar_tick(PlotHandle,0);
        hold off
        %title(['nc14, ',num2str(TimeNC14{MaxTimeSet}(i)-TimeNC14{MaxTimeSet}(1)),' min'])
        ylim([0,1])
        xlim([0.3,0.5])
        set(gca,'XTick',[0.1:0.1:0.8])
        box on
        pause(0.1)
        drawnow
        %StandardFigure(PlotHandle,gca)
        %saveas(gcf,[DropboxFolder,'\Reports\APPlots\APMovies\APMovieMean',Label,'\nc14-',iIndex(i,3),'.tif'])
    end
end
    

%Plot the accumulation vs. time. This is the accumulation per ALL nuclei
MaxValue=max(cellfun(@(x) max(max(x)),TotalVectorAPNC14Interp));
for i=1:MaxTimeIndex
    figure(MovieFigure)
    clf
    
    
    MeanTemp=[];
    for j=1:length(Data)
        if i<=size(TotalVectorAPNC14Interp{j},1)
            MeanTemp=[MeanTemp;TotalVectorAPNC14Interp{j}(i,:)];
        end
    end
    
    if ~isempty(MeanTemp)
        %Now calculale the mean and SD. We need to be careful with the Nans!
        for j=1:size(MeanTemp,2)
            if sum(~isnan(MeanTemp(:,j)))>=MinEmbryos
                IntegralNC14(i,j)=mean(MeanTemp(~isnan(MeanTemp(:,j)),j));
                SDIntegralNC14(i,j)=std(MeanTemp(~isnan(MeanTemp(:,j)),j));
                SEIntegralNC14(i,j)=std(MeanTemp(~isnan(MeanTemp(:,j)),j))/sqrt(sum(~isnan(MeanTemp(:,j))));
            else
                IntegralNC14(i,j)=nan;
                SDIntegralNC14(i,j)=nan;
                SEIntegralNC14(i,j)=nan;
            end
        end


       PlotHandle=errorbar(Data(1).APbinID,IntegralNC14(i,:),...
            SEIntegralNC14(i,:),'.-k');
        errorbar_tick(PlotHandle,0);
        hold off
        %title(['nc14, ',num2str(TimeNC14{MaxTimeSet}(i)-TimeNC14{MaxTimeSet}(1)),' min'])
        ylim([0,MaxValue])
        xlim([0.3,0.5])
        
        set(gca,'XTick',[0.1:0.1:0.8])
        box on
        pause(0.1)
        drawnow
        %StandardFigure(PlotHandle,gca)
        %saveas(gcf,[DropboxFolder,'\Reports\APPlots\APMovies\APMovieMean',Label,'\nc14-',iIndex(i,3),'.tif'])
    end
end


MovieFigure2=figure;

%Plot the accumulation vs. time. This is the accumulation per ON nuclei
MaxValue=max(cellfun(@(x) max(max(x)),TotalVectorONAPNC14Interp));
for i=1:MaxTimeIndex
    figure(MovieFigure2)
    clf
    
    
    MeanTemp=[];
    for j=1:length(Data)
        if i<=size(TotalVectorONAPNC14Interp{j},1)
            MeanTemp=[MeanTemp;TotalVectorONAPNC14Interp{j}(i,:)];
        end
    end
    
    if ~isempty(MeanTemp)
        %Now calculale the mean and SD. We need to be careful with the Nans!
        for j=1:size(MeanTemp,2)
            if sum(~isnan(MeanTemp(:,j)))>=MinEmbryos
                IntegralONNC14(i,j)=mean(MeanTemp(~isnan(MeanTemp(:,j)),j));
                SDIntegralONNC14(i,j)=std(MeanTemp(~isnan(MeanTemp(:,j)),j));
                SEIntegralONNC14(i,j)=std(MeanTemp(~isnan(MeanTemp(:,j)),j))/sqrt(sum(~isnan(MeanTemp(:,j))));
            else
                IntegralONNC14(i,j)=nan;
                SDIntegralONNC14(i,j)=nan;
                SEIntegralONNC14(i,j)=nan;
            end
        end


       PlotHandle=errorbar(Data(1).APbinID,IntegralONNC14(i,:),...
            SEIntegralONNC14(i,:),'.-k');
        errorbar_tick(PlotHandle,0);
        hold off
        %title(['nc14, ',num2str(TimeNC14{MaxTimeSet}(i)-TimeNC14{MaxTimeSet}(1)),' min'])
        ylim([0,MaxValue])
        xlim([0.3,0.5])
        
        set(gca,'XTick',[0.1:0.1:0.8])
        box on
        pause(0.1)
        drawnow
        %StandardFigure(PlotHandle,gca)
        %saveas(gcf,[DropboxFolder,'\Reports\APPlots\APMovies\APMovieMean',Label,'\nc14-',iIndex(i,3),'.tif'])
    end
end


%Accumulation accounting for degradation

%Plot the accumulation+degradation vs. time. This is the accumulation per ALL nuclei
MaxValue=max(cellfun(@(x) max(max(x)),AccumVectorAPNC14Interp));
for i=1:MaxTimeIndex
    figure(MovieFigure)
    clf
    
    
    MeanTemp=[];
    for j=1:length(Data)
        if i<=size(AccumVectorAPNC14Interp{j},1)
            MeanTemp=[MeanTemp;AccumVectorAPNC14Interp{j}(i,:)];
        end
    end
    
    if ~isempty(MeanTemp)
        %Now calculale the mean and SD. We need to be careful with the Nans!
        for j=1:size(MeanTemp,2)
            if sum(~isnan(MeanTemp(:,j)))>=MinEmbryos
                IntegralDegNC14(i,j)=mean(MeanTemp(~isnan(MeanTemp(:,j)),j));
                SDIntegralDegNC14(i,j)=std(MeanTemp(~isnan(MeanTemp(:,j)),j));
                SEIntegralDegNC14(i,j)=std(MeanTemp(~isnan(MeanTemp(:,j)),j))/sqrt(sum(~isnan(MeanTemp(:,j))));
            else
                IntegralDegNC14(i,j)=nan;
                SDIntegralDegNC14(i,j)=nan;
                SEIntegralDegNC14(i,j)=nan;
            end
        end


       PlotHandle=errorbar(Data(1).APbinID,IntegralDegNC14(i,:),...
            SEIntegralDegNC14(i,:),'.-k');
        errorbar_tick(PlotHandle,0);
        hold off
        %title(['nc14, ',num2str(TimeNC14{MaxTimeSet}(i)-TimeNC14{MaxTimeSet}(1)),' min'])
        ylim([0,MaxValue])
        xlim([0.3,0.5])
        
        title('Integral+Degradation ALL nuclei, nc14')
        
        set(gca,'XTick',[0.1:0.1:0.8])
        box on
        pause(0.1)
        drawnow
        %StandardFigure(PlotHandle,gca)
        %saveas(gcf,[DropboxFolder,'\Reports\APPlots\APMovies\APMovieMean',Label,'\nc14-',iIndex(i,3),'.tif'])
    end
end


MovieFigure2=figure;

%Plot the accumulation vs. time. This is the accumulation per ON nuclei
MaxValue=max(cellfun(@(x) max(max(x)),AccumVectorONAPNC14Interp));
for i=1:MaxTimeIndex
    figure(MovieFigure2)
    clf
    
    
    MeanTemp=[];
    for j=1:length(Data)
        if i<=size(AccumVectorONAPNC14Interp{j},1)
            MeanTemp=[MeanTemp;AccumVectorONAPNC14Interp{j}(i,:)];
        end
    end
    
    if ~isempty(MeanTemp)
        %Now calculale the mean and SD. We need to be careful with the Nans!
        for j=1:size(MeanTemp,2)
            if sum(~isnan(MeanTemp(:,j)))>=MinEmbryos
                IntegralDegONNC14(i,j)=mean(MeanTemp(~isnan(MeanTemp(:,j)),j));
                SDIntegralDegONNC14(i,j)=std(MeanTemp(~isnan(MeanTemp(:,j)),j));
                SEIntegralDegONNC14(i,j)=std(MeanTemp(~isnan(MeanTemp(:,j)),j))/sqrt(sum(~isnan(MeanTemp(:,j))));
            else
                IntegralDegONNC14(i,j)=nan;
                SDIntegralDegONNC14(i,j)=nan;
                SEIntegralDegONNC14(i,j)=nan;
            end
        end


       PlotHandle=errorbar(Data(1).APbinID,IntegralDegONNC14(i,:),...
            SEIntegralDegONNC14(i,:),'.-k');
        errorbar_tick(PlotHandle,0);
        hold off
        %title(['nc14, ',num2str(TimeNC14{MaxTimeSet}(i)-TimeNC14{MaxTimeSet}(1)),' min'])
        ylim([0,MaxValue])
        xlim([0.3,0.5])
        
        title('Integral+Degradation ON nuclei, nc14')
        
        set(gca,'XTick',[0.1:0.1:0.8])
        box on
        pause(0.1)
        drawnow
        %StandardFigure(PlotHandle,gca)
        %saveas(gcf,[DropboxFolder,'\Reports\APPlots\APMovies\APMovieMean',Label,'\nc14-',iIndex(i,3),'.tif'])
    end
end



%Plot the accumulation+degradation vs. time. This is the accumulation per ALL nuclei
MaxValue=max(cellfun(@(x) max(max(x)),AccumVectorAPNC13Interp));
for i=1:MaxTimeIndex
    figure(MovieFigure)
    clf
    
    
    MeanTemp=[];
    for j=1:length(Data)
        if i<=size(AccumVectorAPNC13Interp{j},1)
            MeanTemp=[MeanTemp;AccumVectorAPNC13Interp{j}(i,:)];
        end
    end
    
    if ~isempty(MeanTemp)
        %Now calculale the mean and SD. We need to be careful with the Nans!
        for j=1:size(MeanTemp,2)
            if sum(~isnan(MeanTemp(:,j)))>=MinEmbryos
                IntegralDegNC13(i,j)=mean(MeanTemp(~isnan(MeanTemp(:,j)),j));
                SDIntegralDegNC13(i,j)=std(MeanTemp(~isnan(MeanTemp(:,j)),j));
                SEIntegralDegNC13(i,j)=std(MeanTemp(~isnan(MeanTemp(:,j)),j))/sqrt(sum(~isnan(MeanTemp(:,j))));
            else
                IntegralDegNC13(i,j)=nan;
                SDIntegralDegNC13(i,j)=nan;
                SEIntegralDegNC13(i,j)=nan;
            end
        end


       PlotHandle=errorbar(Data(1).APbinID,IntegralDegNC13(i,:),...
            SEIntegralDegNC13(i,:),'.-k');
        errorbar_tick(PlotHandle,0);
        hold off
        %title(['NC13, ',num2str(TimeNC13{MaxTimeSet}(i)-TimeNC13{MaxTimeSet}(1)),' min'])
        ylim([0,MaxValue])
        xlim([0.3,0.5])
        
        title('Integral+Degradation ALL nuclei, NC13')
        
        set(gca,'XTick',[0.1:0.1:0.8])
        box on
        pause(0.1)
        drawnow
        %StandardFigure(PlotHandle,gca)
        %saveas(gcf,[DropboxFolder,'\Reports\APPlots\APMovies\APMovieMean',Label,'\NC13-',iIndex(i,3),'.tif'])
    end
end


MovieFigure2=figure;

%Plot the accumulation vs. time. This is the accumulation per ON nuclei
MaxValue=max(cellfun(@(x) max(max(x)),AccumVectorONAPNC13Interp));
for i=1:MaxTimeIndex
    figure(MovieFigure2)
    clf
    
    
    MeanTemp=[];
    for j=1:length(Data)
        if i<=size(AccumVectorONAPNC13Interp{j},1)
            MeanTemp=[MeanTemp;AccumVectorONAPNC13Interp{j}(i,:)];
        end
    end
    
    if ~isempty(MeanTemp)
        %Now calculale the mean and SD. We need to be careful with the Nans!
        for j=1:size(MeanTemp,2)
            if sum(~isnan(MeanTemp(:,j)))>=MinEmbryos
                IntegralDegONNC13(i,j)=mean(MeanTemp(~isnan(MeanTemp(:,j)),j));
                SDIntegralDegONNC13(i,j)=std(MeanTemp(~isnan(MeanTemp(:,j)),j));
                SEIntegralDegONNC13(i,j)=std(MeanTemp(~isnan(MeanTemp(:,j)),j))/sqrt(sum(~isnan(MeanTemp(:,j))));
            else
                IntegralDegONNC13(i,j)=nan;
                SDIntegralDegONNC13(i,j)=nan;
                SEIntegralDegONNC13(i,j)=nan;
            end
        end


       PlotHandle=errorbar(Data(1).APbinID,IntegralDegONNC13(i,:),...
            SEIntegralDegONNC13(i,:),'.-k');
        errorbar_tick(PlotHandle,0);
        hold off
        %title(['NC13, ',num2str(TimeNC13{MaxTimeSet}(i)-TimeNC13{MaxTimeSet}(1)),' min'])
        ylim([0,MaxValue])
        xlim([0.3,0.5])
        
        title('Integral+Degradation ON nuclei, NC13')
        
        set(gca,'XTick',[0.1:0.1:0.8])
        box on
        pause(0.1)
        drawnow
        %StandardFigure(PlotHandle,gca)
        %saveas(gcf,[DropboxFolder,'\Reports\APPlots\APMovies\APMovieMean',Label,'\NC13-',iIndex(i,3),'.tif'])
    end
end