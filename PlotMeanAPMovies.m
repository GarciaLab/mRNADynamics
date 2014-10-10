function [TimeNC14,MeanNC14,SENC14,MeanAllNC14,SEAllNC14...
    MeanOnRatioNC14,SEOnRatioNC14,...
    IntegralNC14,SEIntegralNC14,IntegralONNC14,SEIntegralONNC14,...
    TimeNC13,MeanNC13,SENC13,MeanAllNC13,SEAllNC13,MeanOnRatioNC13,SEOnRatioNC13,...
    IntegralNC13,SEIntegralNC13,...
    IntegralDegNC14,SEIntegralDegNC14,IntegralDegONNC14,SEIntegralDegONNC14,...
    IntegralDegNC13,SEIntegralDegNC13,IntegralDegONNC13,SEIntegralDegONNC13]=...
    PlotMeanAPMovies(Data,Label,MinEmbryos)

%This function plots the mean AP profile of a data set for the maximum in
%nc13 and nc14 and generates a movie as well.

%IntegralNC14 is the total amount produced per ALL nuclei
%IntegralONNC14 is the total amount produced per ON nuclei 

%% Parameters and data sets

%Parameters:
MinParticles=3;     %Minimum number of particles necessary in a bin for
                    %it to be considered
MinTimePoints=5;    %Minimum number of time points for the interpolation
if ~exist('MinEmbryos')
    MinEmbryos=2;       %Minimum number of embryos
end
           
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
    warning('Confirm size of construct!')
elseif findstr(lower(Data(1).SetName),'knibac')
    HalfLifeForPlot=1000;    %Infinite half-life for now
    ElongationTime=6.443/ElongationRate;
    warning('Confirm size of construct!')
else
    error('Add information for this construct')
end



                    
%Some labels to use for plots
Labels='k.r.g.b.y.c.m.ksrsgsbsyscsmskorogoboyocomok^r^g^b^y^c^m^';
                    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,SchnitzcellsFolder]=...
    DetermineLocalFolders(Data(1).SetName(11:end-1));






%Create the folders we need and empty them if necessary

%Folder for report figures
mkdir([DropboxFolder,filesep,'Reports'])
mkdir([DropboxFolder,filesep,'Reports',filesep,'APPlots',filesep,'APMovies'])
mkdir([DropboxFolder,filesep,'Reports',filesep,'APPlots',filesep,'APMovies',filesep,'APMovie',Label])
mkdir([DropboxFolder,filesep,'Reports',filesep,'APPlots',filesep,'APMovies',filesep,'APMovie',Label,'SE'])
mkdir([DropboxFolder,filesep,'Reports',filesep,'APPlots',filesep,'APMovies',filesep,'APMovieMean',Label])

%Make a movie of all the different embryos as a function of time and AP.
%We'll set a new zero of time at the beginning of each nc.

%Determine the maximum fluorescence value for the whole movie
MaxValue=max(cellfun(@(x) max(max(x)),{Data.MeanVectorAP}));

%Different movies will have different time resolutions. We'll pick the
%smallest time intervals and interpolate the values of all the data sets
%accordingly. We are only doing this for nc13 and nc14.
nc13Resolution=[];
nc14Resolution=[];
for i=1:length(Data)
    if Data(i).nc13>0
        nc13Resolution=[nc13Resolution,mean(diff(Data(i).ElapsedTime(Data(i).nc13:Data(i).nc14)))];
    end
    nc14Resolution=[nc14Resolution,mean(diff(Data(i).ElapsedTime(Data(i).nc14:end)))];
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
nc13Length=[];
for i=1:length(Data)
    if Data(i).nc13>0
        nc13Length=[nc13Length,...
            Data(i).ElapsedTime(Data(i).nc14)-Data(i).ElapsedTime(Data(i).nc13)];
    end
end
nc13Length=min(nc13Length);
%All time points for nc13
TimeNC13=0:nc13Resolution:nc13Length;

%I'm building up all the variables for nc13 in case there are data sets
%that did not catch the mitosis entering into that cycle.
MeanVectorAPNC13=[];
SDVectorAPNC13=[];
SEVectorAPNC13=[];

MeanVectorAllAPNC13=[];
SDVectorAllAPNC13=[];
SEVectorAllAPNC13=[];

TotalVectorAPNC13=[];
TotalVectorONAPNC13=[];
OnRatioAPNC13=[];

MeanVectorAPNC13Interp=[];
MeanVectorAllAPNC13Interp=[];
MeanOnRatioAPNC13Interp=[];
TotalVectorAPNC13Interp=[];
TotalVectorONAPNC13Interp=[];

for j=1:length(Data)
    %Get the AP vector related to each nc
    
    %Filter the bins that have at least the minimum number the particles
    NParticlesAPFilter=Data(j).NParticlesAP>=MinParticles;
    
    %nc13:
    if Data(j).nc13>0
        NParticlesAPFilterNC13=NParticlesAPFilter(Data(j).nc13:Data(j).nc14,:);
        
        %Mean levels of ON nuclei
        MeanVectorAPNC13{end+1}=Data(j).MeanVectorAP(Data(j).nc13:Data(j).nc14,:);
        SDVectorAPNC13{end+1}=Data(j).SDVectorAP(Data(j).nc13:Data(j).nc14,:);
        SEVectorAPNC13{end+1}=Data(j).SDVectorAP(Data(j).nc13:Data(j).nc14,:)./...
            sqrt(Data(j).NParticlesAP(Data(j).nc13:Data(j).nc14,:));
    
        %Mean levels of ALL nuclei
        MeanVectorAllAPNC13{end+1}=Data(j).MeanVectorAP(Data(j).nc13:Data(j).nc14,:).*...
            Data(j).OnRatioAP(Data(j).nc13:Data(j).nc14,:);
        SDVectorAllAPNC13{end+1}=Data(j).SDVectorAP(Data(j).nc13:Data(j).nc14,:).*...
            Data(j).OnRatioAP(Data(j).nc13:Data(j).nc14,:);
        SEVectorAllAPNC13{end+1}=Data(j).SDVectorAP(Data(j).nc13:Data(j).nc14,:).*...
            Data(j).OnRatioAP(Data(j).nc13:Data(j).nc14,:)./...
            sqrt(Data(j).NParticlesAP(Data(j).nc13:Data(j).nc14,:)./...
            Data(j).OnRatioAP(Data(j).nc13:Data(j).nc14,:));
        
        MeanVectorAPNC13{end}(~NParticlesAPFilterNC13)=nan;
        SDVectorAPNC13{end}(~NParticlesAPFilterNC13)=nan;
        SEVectorAPNC13{end}(~NParticlesAPFilterNC13)=nan;
  
        MeanVectorAllAPNC13{end}(~NParticlesAPFilterNC13)=nan;
        SDVectorAllAPNC13{end}(~NParticlesAPFilterNC13)=nan;
        SEVectorAllAPNC13{end}(~NParticlesAPFilterNC13)=nan;
        
        
        %New version of integration. This is using the trapezoidal integration
        %now
        m=1;
        for i=Data(j).nc13:Data(j).nc14
            [TotalVectorAPNC13{end+1}(m,:),Dummy,TotalVectorONAPNC13{end+1}(m,:),Dummy]=...
                IntegrateIndivSetTimeWindow(Data(j),MinParticles,...
                13,[0,Data(j).ElapsedTime(i)-Data(j).ElapsedTime(Data(j).nc13)]);
            m=m+1;
        end
        
        %Fraction of active nuclei vs. time
        OnRatioAPNC13{end+1}=Data(j).OnRatioAP(Data(j).nc13:Data(j).nc14,:);
        
        %Interpolate each AP bin as a function of time
        MeanVectorAPNC13Interp{end+1}=nan(length(TimeNC13),length(Data(j).APbinID));
        MeanVectorAllAPNC13Interp{end+1}=nan(length(TimeNC13),length(Data(j).APbinID));
        MeanOnRatioAPNC13Interp{end+1}=nan(length(TimeNC13),length(Data(j).APbinID));
        TotalVectorAPNC13Interp{end+1}=nan(length(TimeNC13),length(Data(j).APbinID));
        TotalVectorONAPNC13Interp{end+1}=nan(length(TimeNC13),length(Data(j).APbinID));
        
        for i=1:length(Data(j).APbinID)
        
            if sum(~isnan(MeanVectorAPNC13{end}(:,i)))>1

                TimeWindow=Data(j).ElapsedTime(Data(j).nc13:Data(j).nc14)-...
                    Data(j).ElapsedTime(Data(j).nc13);

                %Use only data points where we had MinParticles to calculate
                %the mean
                FilterTemp=NParticlesAPFilter(Data(j).nc13:Data(j).nc14,i);

                if sum(FilterTemp)>=MinTimePoints

                    MeanVectorAPNC13Interp{end}(:,i)=pchip(TimeWindow(FilterTemp),...
                        MeanVectorAPNC13{end}(FilterTemp,i)',...
                        TimeNC13);
                    

                    %Check if there are enough statistics
                    if sum(~isnan(MeanVectorAllAPNC13{end}(FilterTemp,i)))>=MinParticles
                        MeanVectorAllAPNC13Interp{end}(:,i)=pchip(TimeWindow(FilterTemp),...
                            MeanVectorAllAPNC13{end}(FilterTemp,i)',...
                            TimeNC13);
                    else
                        MeanVectorAllAPNC13Interp{end}(:,i)=nan;
                    end
                    
                    if sum(~isnan(OnRatioAPNC13{end}(FilterTemp,i)))>=MinParticles
                        MeanOnRatioAPNC13Interp{end}(:,i)=pchip(TimeWindow(FilterTemp),...
                            OnRatioAPNC13{end}(FilterTemp,i)',...
                            TimeNC13);
                    else
                        MeanOnRatioAPNC13Interp{end}(:,i)=nan;
                    end

                    if sum(~isnan(TotalVectorAPNC13{end}(FilterTemp,i)))>=MinParticles
                        TotalVectorAPNC13Interp{end}(:,i)=pchip(TimeWindow(FilterTemp),...
                            TotalVectorAPNC13{end}(FilterTemp,i)',...
                            TimeNC13);
                    else
                        TotalVectorAPNC13Interp{end}(:,i)=nan;
                    end

                    if sum(~isnan(TotalVectorONAPNC13{end}(FilterTemp,i)))>=MinParticles
                        TotalVectorONAPNC13Interp{end}(:,i)=pchip(TimeWindow(FilterTemp),...
                            TotalVectorONAPNC13{end}(FilterTemp,i)',...
                            TimeNC13);
                    else
                        TotalVectorONAPNC13Interp{end}(:,i)=nan;
                    end



                    MeanVectorAPNC13Interp{end}(TimeNC13<min(TimeWindow(FilterTemp)),i)=nan;
                    MeanVectorAPNC13Interp{end}(TimeNC13>max(TimeWindow(FilterTemp)),i)=nan;
                    
                    MeanVectorAllAPNC13Interp{end}(TimeNC13<min(TimeWindow(FilterTemp)),i)=nan;
                    MeanVectorAllAPNC13Interp{end}(TimeNC13>max(TimeWindow(FilterTemp)),i)=nan;

                    MeanOnRatioAPNC13Interp{end}(TimeNC13<min(TimeWindow(FilterTemp)),i)=nan;
                    MeanOnRatioAPNC13Interp{end}(TimeNC13>max(TimeWindow(FilterTemp)),i)=nan;

                    %For the total amount produced I need 0 before the time
                    %window and the maximum after the time window
                    TotalVectorAPNC13Interp{end}(TimeNC13<min(TimeWindow(FilterTemp)),i)=0;
                    TotalVectorAPNC13Interp{end}(TimeNC13>max(TimeWindow(FilterTemp)),i)=...
                        max(TotalVectorAPNC13Interp{end}(TimeNC13<max(TimeWindow(FilterTemp)),i));

                    TotalVectorONAPNC13Interp{end}(TimeNC13<min(TimeWindow(FilterTemp)),i)=0;
                    TotalVectorONAPNC13Interp{end}(TimeNC13>max(TimeWindow(FilterTemp)),i)=...
                        max(TotalVectorONAPNC13Interp{end}(TimeNC13<max(TimeWindow(FilterTemp)),i));

                else
                    MeanVectorAPNC13Interp{end}(:,i)=nan;
                    MeanVectorAllAPNC13Interp{end}(:,i)=nan;
                    MeanOnRatioAPNC13Interp{end}(:,i)=nan;
                    TotalVectorAPNC13Interp{end}(:,i)=nan;
                    TotalVectorONAPNC13Interp{end}(:,i)=nan;
                end
            end
        end
        
        
    end
    
    
    
    
    
    
    %nc14:
    NParticlesAPFilterNC14=NParticlesAPFilter(Data(j).nc14:end,:);
    
   
    %Mean levels
    MeanVectorAPNC14{j}=Data(j).MeanVectorAP(Data(j).nc14:end,:);
    SDVectorAPNC14{j}=Data(j).SDVectorAP(Data(j).nc14:end,:);
    SEVectorAPNC14{j}=Data(j).SDVectorAP(Data(j).nc14:end,:)./...
        sqrt(Data(j).NParticlesAP(Data(j).nc14:end,:));
    
    %Mean levels of ALL nuclei
    MeanVectorAllAPNC14{j}=Data(j).MeanVectorAP(Data(j).nc14:end,:,:).*...
        Data(j).OnRatioAP(Data(j).nc14:end,:,:);
    SDVectorAllAPNC14{j}=Data(j).SDVectorAP(Data(j).nc14:end,:,:).*...
        Data(j).OnRatioAP(Data(j).nc14:end,:,:);
    SEVectorAllAPNC14{j}=Data(j).SDVectorAP(Data(j).nc14:end,:,:).*...
        Data(j).OnRatioAP(Data(j).nc14:end,:,:)./...
        sqrt(Data(j).NParticlesAP(Data(j).nc14:end,:,:)./...
        Data(j).OnRatioAP(Data(j).nc14:end,:,:));

    
    MeanVectorAPNC14{j}(~NParticlesAPFilterNC14)=nan;
    SDVectorAPNC14{j}(~NParticlesAPFilterNC14)=nan;
    SEVectorAPNC14{j}(~NParticlesAPFilterNC14)=nan;
    
    MeanVectorAllAPNC14{j}(~NParticlesAPFilterNC14)=nan;
    SDVectorAllAPNC14{j}(~NParticlesAPFilterNC14)=nan;
    SEVectorAllAPNC14{j}(~NParticlesAPFilterNC14)=nan;
    
    
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
        
        
        %Interpolate each AP bin as a function of time
        AccumVectorAPNC14Interp{j}=nan(length(TimeNC14),length(Data(j).APbinID));
        AccumVectorONAPNC14Interp{j}=nan(length(TimeNC14),length(Data(j).APbinID));

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

                    AccumVectorAPNC13Interp{j}(TimeNC13<min(TimeWindow(FilterTemp)-ElongationTime),i)=0;
                    AccumVectorAPNC13Interp{j}(TimeNC13>max(TimeWindow(FilterTemp)),i)=...
                        max(AccumVectorAPNC13Interp{j}(TimeNC13<max(TimeWindow(FilterTemp)),i));

                    AccumVectorONAPNC13Interp{j}(TimeNC13<min(TimeWindow(FilterTemp)-ElongationTime),i)=0;
                    AccumVectorONAPNC13Interp{j}(TimeNC13>max(TimeWindow(FilterTemp)),i)=...
                        max(AccumVectorONAPNC13Interp{j}(TimeNC13<max(TimeWindow(FilterTemp)),i));

                else
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

                    AccumVectorAPNC14Interp{j}(:,i)=pchip(TimeWindow(FilterTemp),...
                        AccumNC14{j}(FilterTemp,i)',TimeNC14);
                    AccumVectorONAPNC14Interp{j}(:,i)=pchip(TimeWindow(FilterTemp),...
                        AccumNC14ON{j}(FilterTemp,i)',TimeNC14);

                    %For the total amount produced I need 0 before the time
                    %window and the maximum after the time window
                    AccumVectorAPNC14Interp{j}(TimeNC14<min(TimeWindow(FilterTemp)+ElongationTime),i)=0;
                    AccumVectorAPNC14Interp{j}(TimeNC14>max(TimeWindow(FilterTemp)),i)=...
                        max(AccumVectorAPNC14Interp{j}(TimeNC14<max(TimeWindow(FilterTemp)),i));

                    AccumVectorONAPNC14Interp{j}(TimeNC14<min(TimeWindow(FilterTemp)+ElongationTime),i)=0;
                    AccumVectorONAPNC14Interp{j}(TimeNC14>max(TimeWindow(FilterTemp)),i)=...
                        max(AccumVectorONAPNC14Interp{j}(TimeNC14<max(TimeWindow(FilterTemp)),i));

                else
                    AccumVectorAPNC14Interp{j}(:,i)=nan;
                    AccumVectorONAPNC14Interp{j}(:,i)=nan;
                end
            end
        end
        
    end   
          

    %Fraction of active nuclei vs. time
    OnRatioAPNC14{j}=Data(j).OnRatioAP(Data(j).nc14:end,:);
   
    
    %Interpolate each AP bin as a function of time
    MeanVectorAPNC14Interp{j}=nan(length(TimeNC14),length(Data(j).APbinID));
    MeanVectorAllAPNC14Interp{j}=nan(length(TimeNC14),length(Data(j).APbinID));
    MeanOnRatioAPNC14Interp{j}=nan(length(TimeNC14),length(Data(j).APbinID));
    TotalVectorAPNC14Interp{j}=nan(length(TimeNC14),length(Data(j).APbinID));
    TotalVectorONAPNC14Interp{j}=nan(length(TimeNC14),length(Data(j).APbinID));
    

    for i=1:length(Data(j).APbinID)
        
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
                
                if sum(~isnan(MeanVectorAllAPNC14{j}(FilterTemp,i)))>=MinParticles
                    MeanVectorAllAPNC14Interp{j}(:,i)=pchip(TimeWindow(FilterTemp),...
                        MeanVectorAllAPNC14{j}(FilterTemp,i)',...
                        TimeNC14);
                else
                    MeanVectorAllAPNC14Interp{j}(:,i)=nan;
                end

                if sum(~isnan(OnRatioAPNC14{j}(FilterTemp,i)))>=MinParticles
                    MeanOnRatioAPNC14Interp{j}(:,i)=pchip(TimeWindow(FilterTemp),...
                        OnRatioAPNC14{j}(FilterTemp,i)',...
                        TimeNC14);
                else
                    MeanOnRatioAPNC14Interp{j}(:,i)=nan;
                end

                if sum(~isnan(TotalVectorAPNC14{j}(FilterTemp,i)))>=MinParticles
                    TotalVectorAPNC14Interp{j}(:,i)=pchip(TimeWindow(FilterTemp),...
                        TotalVectorAPNC14{j}(FilterTemp,i)',...
                        TimeNC14);
                else
                    TotalVectorAPNC14Interp{j}(:,i)=nan;
                end

                if sum(~isnan(TotalVectorONAPNC14{j}(FilterTemp,i)))>=MinParticles
                    TotalVectorONAPNC14Interp{j}(:,i)=pchip(TimeWindow(FilterTemp),...
                        TotalVectorONAPNC14{j}(FilterTemp,i)',...
                        TimeNC14);
                else
                    TotalVectorONAPNC14Interp{j}(:,i)=nan;
                end
 
            
                MeanVectorAPNC14Interp{j}(TimeNC14<min(TimeWindow(FilterTemp)),i)=nan;
                MeanVectorAPNC14Interp{j}(TimeNC14>max(TimeWindow(FilterTemp)),i)=nan;
                
                MeanVectorAllAPNC14Interp{j}(TimeNC14<min(TimeWindow(FilterTemp)),i)=nan;
                MeanVectorAllAPNC14Interp{j}(TimeNC14>max(TimeWindow(FilterTemp)),i)=nan;
                
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
                

            else
                MeanVectorAPNC14Interp{j}(:,i)=nan;
                MeanVectorAllAPNC14Interp{j}(:,i)=nan;
                MeanOnRatioAPNC14Interp{j}(:,i)=nan;
                TotalVectorAPNC14Interp{j}(:,i)=nan;
                TotalVectorONAPNC14Interp{j}(:,i)=nan;
            end
        end
    end
end





%% Movie of averages

MovieFigure=figure;

%nc13:

MaxTimeIndex=length(TimeNC13);

%Average fluorescence per ON nuclei
for i=1:MaxTimeIndex

    figure(MovieFigure)
    clf
    
    %Populate a matrix with the values as a function of AP for this time
    %point.
    MeanTemp=[];
    for j=1:length(MeanVectorAPNC13Interp)
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
        ylim([0,1.5E4])
        xlim([0,1])
        set(gca,'XTick',[0.1:0.1:0.8])
        box on
        pause(0.1)
        xlabel('AP position (x/L)')
        ylabel('Mean fluorescence per ON nuclei (au)')
        title('nc13')
        drawnow
        %StandardFigure(PlotHandle,gca)
        %saveas(gcf,[DropboxFolder,'\Reports\APPlots\APMovies\APMovieMean',Label,'\nc14-',iIndex(i,3),'.tif'])
    end
end


%Average fluorescence per ALL nuclei
for i=1:MaxTimeIndex

    figure(MovieFigure)
    clf
    
    %Populate a matrix with the values as a function of AP for this time
    %point.
    MeanTemp=[];
    for j=1:length(MeanVectorAllAPNC13Interp)
        if i<=size(MeanVectorAllAPNC13Interp{j},1)
            MeanTemp=[MeanTemp;MeanVectorAllAPNC13Interp{j}(i,:)];
        end
    end
    
    if ~isempty(MeanTemp)
        %Now calculale the mean and SD. We need to be careful with the Nans!
        for j=1:size(MeanTemp,2)
            if sum(~isnan(MeanTemp(:,j)))>=MinEmbryos
                MeanAllNC13(i,j)=mean(MeanTemp(~isnan(MeanTemp(:,j)),j));
                SDAllNC13(i,j)=std(MeanTemp(~isnan(MeanTemp(:,j)),j));
                SEAllNC13(i,j)=std(MeanTemp(~isnan(MeanTemp(:,j)),j))/sqrt(sum(~isnan(MeanTemp(:,j))));
            else
                MeanAllNC13(i,j)=nan;
                SDAllNC13(i,j)=nan;
                SEAllNC13(i,j)=nan;
            end
        end

       PlotHandle=errorbar(Data(1).APbinID,MeanAllNC13(i,:),...
            SEAllNC13(i,:),'.-k');
        errorbar_tick(PlotHandle,0);
        hold off
        ylim([0,1.5E4])
        xlim([0,1])
        set(gca,'XTick',[0.1:0.1:0.8])
        box on
        pause(0.1)
        xlabel('AP position (x/L)')
        ylabel('Mean fluorescence per ALL nuclei (au)')
        title('nc13')
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
    for j=1:length(OnRatioAPNC13)
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
MaxValues=[];
for j=1:length(TotalVectorAPNC13Interp)
    if ~isnan(max(max(TotalVectorAPNC13Interp{j})))
        MaxValues=[MaxValues,max(max(TotalVectorAPNC13Interp{j}))];
    end
end
MaxValue=max(MaxValues);

if ~isempty(MaxValue)
    for i=1:MaxTimeIndex
        figure(MovieFigure)
        clf


        MeanTemp=[];
        for j=1:length(TotalVectorAPNC13Interp)
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
            %ylim([0,MaxValue])
            xlim([0.3,0.5])

            set(gca,'XTick',[0.1:0.1:0.8])
            box on
            pause(0.1)
            drawnow
            %StandardFigure(PlotHandle,gca)
            %saveas(gcf,[DropboxFolder,'\Reports\APPlots\APMovies\APMovieMean',Label,'\nc14-',iIndex(i,3),'.tif'])
        end
    end
else
    IntegralNC13=nan;
    SDIntegralNC13=nan;
    SEIntegralNC13=nan;
end





%nc14:
MaxTimeIndex=length(TimeNC14);

%Average fluorescence of ON nuclei
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
    

%nc14:


%Average fluorescence of ALL nuclei
for i=1:MaxTimeIndex

    figure(MovieFigure)
    clf
    
    
    %Populate a matrix with the values as a function of AP for this time
    %point.
    MeanTemp=[];
    for j=1:length(Data)
        if i<=size(MeanVectorAllAPNC14Interp{j},1)
            MeanTemp=[MeanTemp;MeanVectorAllAPNC14Interp{j}(i,:)];
        end
    end
    
    if ~isempty(MeanTemp)
        %Now calculale the mean and SD. We need to be careful with the Nans!
        for j=1:size(MeanTemp,2)
            if sum(~isnan(MeanTemp(:,j)))>=MinEmbryos
                MeanAllNC14(i,j)=mean(MeanTemp(~isnan(MeanTemp(:,j)),j));
                SDAllNC14(i,j)=std(MeanTemp(~isnan(MeanTemp(:,j)),j));
                SEAllNC14(i,j)=std(MeanTemp(~isnan(MeanTemp(:,j)),j))/sqrt(sum(~isnan(MeanTemp(:,j))));
            else
                MeanAllNC14(i,j)=nan;
                SDAllNC14(i,j)=nan;
                SEAllNC14(i,j)=nan;
            end
        end


       PlotHandle=errorbar(Data(1).APbinID,MeanAllNC14(i,:),...
            SEAllNC14(i,:),'.-k');
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

if exist('AccumVectorAPNC14Interp')
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
else
    IntegralDegNC14=nan;
    SEIntegralDegNC14=nan;
    IntegralDegONNC14=nan;
    SEIntegralDegONNC14=nan;
    IntegralDegNC13=nan;
    SEIntegralDegNC13=nan;
    IntegralDegONNC13=nan;
    SEIntegralDegONNC13=nan;
end