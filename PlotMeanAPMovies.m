function [TimeNC14,MeanNC14,SENC14,MeanOnRatioNC14,SEOnRatioNC14]=...
    PlotMeanAPMovies(Data,Label)

%This function plots the mean AP profile of a data set for the maximum in
%nc13 and nc14 and generates a movie as well.


%% Parameters and data sets

%Parameters:
MinParticles=3;     %Minimum number of particles necessary in a bin for
                    %it to be considered
MinTimePoints=5;    %Minimum number of time points for the interpolation
MinEmbryos=2;       %Minimum number of embryos
                    
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

%Figure out the time window
for i=1:length(Data)
    nc14Length(i)=Data(i).ElapsedTime(end)-Data(i).ElapsedTime(Data(i).nc14);
end
nc14Length=min(nc14Length);

%All time points for nc14
TimeNC14=0:nc14Resolution:nc14Length;


for j=1:length(Data)

    %Create the corresponding time vector for this set
    TimeNC13{j}=Data(j).ElapsedTime(Data(j).nc13):...
        nc13Resolution:Data(j).ElapsedTime(Data(j).nc14);
%     TimeNC14{j}=Data(j).ElapsedTime(Data(j).nc14):...
%         nc14Resolution:Data(j).ElapsedTime(end);
    
    
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
    
    
    
    %Fraction of active nuclei vs. time
    OnRatioAPNC14{j}=Data(j).OnRatioAP(Data(j).nc14:end,:);
    
    
    %Interpolate each AP bin as a function of time
    MeanVectorAPNC14Interp{j}=nan(length(TimeNC14),length(Data(j).APbinID));
    MeanOnRatioAPNC14Interp{j}=nan(length(TimeNC14),length(Data(j).APbinID));
    for i=1:length(Data(j).APbinID)
        
        if sum(~isnan(MeanVectorAPNC13{j}(:,i)))>1
            MeanVectorAPNC13Interp{j}(:,i)=pchip(Data(j).ElapsedTime(Data(j).nc13:Data(j).nc14),...
                MeanVectorAPNC13{j}(:,i)',...
                TimeNC13{j});
            SDVectorAPNC13Interp{j}(:,i)=pchip(Data(j).ElapsedTime(Data(j).nc13:Data(j).nc14),...
                SDVectorAPNC13{j}(:,i)',...
                TimeNC13{j});
            SEVectorAPNC13Interp{j}(:,i)=pchip(Data(j).ElapsedTime(Data(j).nc13:Data(j).nc14),...
                SEVectorAPNC13{j}(:,i)',...
                TimeNC13{j});
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
            
                %Make Nan everything outside the range of the original data
                MeanVectorAPNC14Interp{j}(1:(min(find(~isnan(MeanVectorAPNC14{j}(:,i))))-1),i)=nan;
                MeanVectorAPNC14Interp{j}((max(find(~isnan(MeanVectorAPNC14{j}(:,i))))+1):end,i)=nan;
                
                MeanOnRatioAPNC14Interp{j}(1:(min(find(~isnan(OnRatioAPNC14{j}(:,i))))-1),i)=nan;
                MeanOnRatioAPNC14Interp{j}((max(find(~isnan(OnRatioAPNC14{j}(:,i))))+1):end,i)=nan;
            else
                MeanVectorAPNC14Interp{j}(:,i)=nan;
                MeanOnRatioAPNC14Interp{j}(:,i)=nan;
            end
                
                
            %Compare the interpolation to the actual data
%             figure(10)
%             plot(Data(j).ElapsedTime(Data(j).nc14:end)-...
%                 Data(j).ElapsedTime(Data(j).nc14),...
%                 MeanVectorAPNC14{j}(:,i)','.-k')
%             hold on
%             plot(TimeNC14,MeanVectorAPNC14Interp{j}(:,i)','-r')
%             hold off
%             pause
            
            
%             SDVectorAPNC14Interp{j}(:,i)=pchip(Data(j).ElapsedTime(Data(j).nc14:end),...
%                 SDVectorAPNC14{j}(:,i)',...
%                 TimeNC14{j});
%             SEVectorAPNC14Interp{j}(:,i)=pchip(Data(j).ElapsedTime(Data(j).nc14:end),...
%                 SEVectorAPNC14{j}(:,i)',...
%                 TimeNC14{j});
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

%% Movie of average

%I removed this because the interpolation wasn't working well


%Plot the resulting traces for nc13 - Using Standard Deviation
% [MaxTimeIndex,MaxTimeSet]=max(cellfun(@length,TimeNC13));
% 
% MovieFigure=figure;
% for i=1:MaxTimeIndex
% 
%     clf
%     PlotHandle=[];
%     figure(MovieFigure)
%     
%     %Populate a matrix with the values as a function of AP for this time
%     %point.
%     MeanTemp=[];
%     for j=1:length(Data)
%         if i<=size(MeanVectorAPNC13Interp{j},1)
%             MeanTemp=[MeanTemp;MeanVectorAPNC13Interp{j}(i,:)];
%         end
%     end
%     
%     %Now calculate the mean and SD. We need to be careful with the Nans!
%     for j=1:size(MeanTemp,2)
%         MeanNC13(i,j)=mean(MeanTemp(~isnan(MeanTemp(:,j)),j));
%         SDNC13(i,j)=std(MeanTemp(~isnan(MeanTemp(:,j)),j));
%     end
%     
%     PlotHandle=errorbar(Data(1).APbinID,MeanNC13(i,:),...
%         SDNC13(i,:),'.-k');
% 
%     title(['nc13, ',num2str(TimeNC13{MaxTimeSet}(i)-TimeNC13{MaxTimeSet}(1)),' min'])
%     ylim([0,MaxValue])
%     xlim([0.1,0.8])
%     set(gca,'XTick',[0.1:0.1:0.8])
%     box on
%     %drawnow
%     StandardFigure(PlotHandle,gca)
%     saveas(gcf,[DropboxFolder,'\Reports\APPlots\APMovies\APMovieMean',Label,'\nc13-',iIndex(i,3),'.tif'])
% end
%     

%Plot the resulting traces for nc14  - Using Standard Error
%[MaxTimeIndex,MaxTimeSet]=max(cellfun(@length,TimeNC14));

MaxTimeIndex=length(TimeNC14);

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
    
