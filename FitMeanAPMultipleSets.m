function FitMeanAPMultipleSets(DataType, minParticles, minTimePoints)

%This function performs fits to the mean fluorescence as a function of time
%of a particular dataset.

%2013/08/18: Changed this so it can automatically detect whether we are
%dealing with a 5' or 3' data set

%Fitting:
%a,z: On time
%s,x: Off time
%d,c: Rate, fine
%D,C: Rate, coarse

%Moving around:
%, .: Move in AP
%n,m: Move in nc
%k,l: Change fit range from the right
%h,j: Change fit range from the left


%Parameters:
MinParticles=minParticles;     %Minimum number of particles in an AP bin
MinTimePoints=minTimePoints;    %Minimum number of time points where we'll have at least
                    %the minimum number of particles.
ElongationRate=1.54;    %In kb/minutes.
GeneLength5=5.296;      %Distance from the first MS2 site to the end of the
                        %TUB3'UTR in kb for the 5' constrcut.
GeneLength3=1.941;      %Distance from the first MS2 site to the end of the
                        %TUB3'UTR in kb for the 3' constrcut.                        
                                   
close all

%Find out which computer this is. That will determine the folder structure.
%Information about about folders

% ES 2013-10-29: Required for multiple users to be able to analyze data on
% one computer
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;
% 
% 
% if ~isempty(varargin)
%     Prefix=varargin;
               

FolderTemp=uigetdir(DropboxFolder,'Choose folder with files to analyze. If multiple sets, pick one.');
Dashes=strfind(FolderTemp,'\');
Prefix=FolderTemp((Dashes(end)+1):end);

all_data = LoadMS2Sets(DataType);
MeanVectorAPTot = [];
ElapsedTimeTot = [];
NParticlesAPTot = [];
SDVectorAPTot = [];
OnRatioAPTot = [];
FramesTot = [];
APDivisionTot = 0*all_data(1).APDivision;
t12 = [];
t13 = [];
t14 = [];
for k=1:length(all_data)
    APDivisionTot = APDivisionTot + all_data(k).APDivision+all_data(k).APDivision;
    MeanVectorAPTot = vertcat(MeanVectorAPTot, all_data(k).MeanVectorAP);
    ElapsedTimeTot = horzcat(ElapsedTimeTot, all_data(k).ElapsedTime);
    NParticlesAPTot = vertcat(NParticlesAPTot,all_data(k).NParticlesAP);
    SDVectorAPTot = vertcat(SDVectorAPTot,all_data(k).SDVectorAP);
    OnRatioAPTot = vertcat(OnRatioAPTot,all_data(k).OnRatioAP);
    t12 = [t12, all_data(k).ElapsedTime(all_data(k).nc12)];
    t13 = [t13, all_data(k).ElapsedTime(all_data(k).nc13)];
    t14 = [t14, all_data(k).ElapsedTime(all_data(k).nc14)];
end
tmins = [min(t12), min(t13), min(t14)];
tminsindices = [find(ElapsedTimeTot == tmins(1)), find(tmins(2) == ElapsedTimeTot), find(tmins(3) == ElapsedTimeTot)];

%Rough frame window to consider in the fits

for k=1:length(all_data)
    
    %Some data sets won't ha ve nc12
    %% 'hu
    if all_data(k).nc12>0
        FrameWindow12{k}=[all_data(k).nc12:all_data(k).nc13];
    else
        FrameWindow12{k}=[];
    end
    %Some data sets won't have nc13 either
    if all_data(k).nc13>0
        FrameWindow13{k}=[all_data(k).nc13:all_data(k).nc14];
    else
        FrameWindow13{k}=[];
    end



    FrameWindow14{k}=[all_data(k).nc14:length(all_data(k).ElapsedTime)];      


    %Detect what type of data set we're dealing with so we can set the delay
    if (~isempty(findstr(Prefix,'X1')))|(~isempty(findstr(Prefix,'P2P')))|...
            (~isempty(findstr(Prefix,'evePr')))
        Delay=GeneLength5/ElongationRate;    %Minutes for PolII to fall off after reaching
                                            %the first MS2 site.
        display('Treating data set as 5''')

        %Initial parameters for fits
        Rate012=4E3;     %Rate per minute
        TimeStart012=3;
        TimeEnd012=7;

        Rate013=4E3;     %Rate per minute
        TimeStart013=5;
        TimeEnd013=10;

        Rate014=4E3;     %Rate per minute
        TimeStart014=5;
        TimeEnd014=1000;  
    elseif ~isempty(findstr(Prefix,'X2'))
        Delay=GeneLength3/ElongationRate;
        display('Treating data set as 3''')

        Rate012=4E3;     %Rate per minute
        TimeStart012=3;
        TimeEnd012=7;

        Rate013=4E3;     %Rate per minute
        TimeStart013=7.5;
        TimeEnd013=10;

        Rate014=4E3;     %Rate per minute
        TimeStart014=7.5;
        TimeEnd014=1000; 
    else
        % ES 2013-10-14: I don't use X1 or X2 prefixes
        if exist('StemLoopEnd', 'var') && strcmp(StemLoopEnd, '5''')
            Delay=GeneLength5/ElongationRate;    %Minutes for PolII to fall off after reaching
            %the first MS2 site.
            display('Treating data set as 5''')

            Rate012=4E3;     %Rate per minute
            TimeStart012=3;
            TimeEnd012=7;

            Rate013=4E3;     %Rate per minute
            TimeStart013=5;
            TimeEnd013=10;

            Rate014=4E3;     %Rate per minute
            TimeStart014=5;
            TimeEnd014=1000;
        elseif exist('StemLoopEnd', 'var') && strcmp(StemLoopEnd, '3''')
            Delay=GeneLength3/ElongationRate;
            display('Treating data set as 3''')

            Rate012=4E3;     %Rate per minute
            TimeStart012=3;
            TimeEnd012=7;

            Rate013=4E3;     %Rate per minute
            TimeStart013=7.5;
            TimeEnd013=10;

            Rate014=4E3;     %Rate per minute
            TimeStart014=7.5;
            TimeEnd014=1000;
        else
            error('Could not recognize data type from the Prefix or from the value of StemLoopEnd in MovieDatabase.')
        end
    end
end


%Set the first guess for the parameters for each AP bin and also
%disapproved the ones that did not have enough data points. FitResults
%is a structure with the fits corresponding to each AP position at nc12, nc13
%or nc14

%Going to assume all data sets have the same AP binning

FitResults(length(all_data(1).APbinID),3).Rate0=[];

%Set default starting values for nc 13 and nc14
%nc12
for i=1:length(all_data(1).APbinID)
    if isempty(FitResults(i,1).Rate0)
        FitResults(i,1).Rate0=Rate012;    
        FitResults(i,1).TimeStart0=TimeStart012;
        FitResults(i,1).TimeEnd0=TimeEnd012;
        FitResults(i,1).FrameFilter=[];
        FitResults(i,1).FitFrameRange=[];
        p = 0*all_data(1).NParticlesAP(FrameWindow12{1}, i);
        for k=1:length(all_data)
            p = SumArrays( p ,all_data(k).NParticlesAP(FrameWindow12{k}, i));
        end
        if sum(p>=MinParticles)>=MinTimePoints
            FitResults(i,1).Approved=0;
        else
            FitResults(i,1).Approved=-1;
        end
    end
end

%nc13
for i=1:length(all_data(1).APbinID)
    if isempty(FitResults(i,2).Rate0)
        FitResults(i,2).Rate0=Rate013;    
        FitResults(i,2).TimeStart0=TimeStart013;
        FitResults(i,2).TimeEnd0=TimeEnd013;
        FitResults(i,2).FrameFilter=[];
        FitResults(i,2).FitFrameRange=[];
        p = 0*all_data(1).NParticlesAP(FrameWindow13{1}, i);
        for k=1:length(all_data)
            p = SumArrays( p ,all_data(k).NParticlesAP(FrameWindow13{k}, i));
        end
        if sum(p>=MinParticles)>=MinTimePoints
            FitResults(i,2).Approved=0;
        else
            FitResults(i,2).Approved=-1;
        end

    end
end
%nc14
for i=1:length(all_data(k).APbinID)
    if isempty(FitResults(i,3).Rate0)
        FitResults(i,3).Rate0=Rate014;    
        FitResults(i,3).TimeStart0=TimeStart014;
        FitResults(i,3).TimeEnd0=[];
        FitResults(i,3).FrameFilter=[];
        FitResults(i,3).FitFrameRange=[];        
        p = 0*all_data(1).NParticlesAP(FrameWindow14{1}, i);
        for k=1:length(all_data)
            p = SumArrays( p ,all_data(k).NParticlesAP(FrameWindow14{k}, i));
        end
        if sum(p>=MinParticles)>=MinTimePoints
            FitResults(i,3).Approved=0;
        else
            FitResults(i,3).Approved=-1;
        end
    end
end


%Figure out which nc we can use
if any(all_data.nc12)
    CurrentNC=12;
    MinNC=12;       %We'll use this to keep track of the minimum nc
elseif any(all_data.nc13)
    all_data(k).CurrentNC=13;
    MinNC=13;
elseif any(all_data.nc14)
    all_data(k).CurrentNC=14;
    MinNC=14;
else
    error('There is a problem. Are the ncs defined?')
end



%Go through each AP bin
FitFigure=figure;
i=min(find(sum(NParticlesAPTot)));
cc=1;

while (cc~=13)
    
    figure(FitFigure)
    clf
    
    if FitResults(i,CurrentNC-11).Approved==-1
        set(gcf,'Color','r')
    elseif FitResults(i,CurrentNC-11).Approved==1
        set(gcf,'Color','g')
    else
        set(gcf,'Color','default')
    end
    
    
%WILL NEED TO UPDATE THIS TO GET FRAME INFO FROM APDIVISIONS
    if APDivisionTot(CurrentNC,i)
        if CurrentNC~=14
            FrameWindow = tminsindices(CurrentNC-11):tminsindices(CurrentNC-11+1);
        else
            FrameWindow=tminsindices(CurrentNC-11):length(ElapsedTimeTot);
        end

        %Check that we have the minimum number of particles for a minimum
        %amount of time
        if sum(NParticlesAPTot(FrameWindow, i)>=MinParticles)>=MinTimePoints

            %Extract the data for this range of frames
            FluoData=MeanVectorAPTot(FrameWindow,i);
            SDFluoData=SDVectorAPTot(FrameWindow,i);
            NData=NParticlesAPTot(FrameWindow,i);
            TimeData=ElapsedTimeTot(FrameWindow);
            OnRatioData=OnRatioAPTot(FrameWindow,i);

            %Now filter them according the number of particles
            FrameFilter=NData>=MinParticles;


            %As an initial guess, use FrameFilter to determine the range of the
            %fit
            if isempty(FitResults(i,CurrentNC-11).FitFrameRange)
                FitFrameRange=FrameWindow(FrameFilter);
                if CurrentNC==14
                    %NEED TO CHANGE tminindices to something that depends
                    %on apdivisions
                    FitFrameRange=FitFrameRange((ElapsedTimeTot(FitFrameRange)-ElapsedTimeTot(tminsindices(CurrentNC-11)))<12);
                end
                FitResults(i,CurrentNC-11).FitFrameRange=FitFrameRange;
            else
                FitFrameRange=FitResults(i,CurrentNC-11).FitFrameRange;
            end

            %Filter the frames according to FitFrameRange
            FitFrameFilter=ismember(FrameWindow,FitFrameRange);
            
            OnRatioDataForFit=OnRatioData(FitFrameFilter);
            MaxOnRatioForFit=max(OnRatioData);
            OnRatioDataForFit=OnRatioDataForFit/MaxOnRatioForFit;

            FluoDataForFit=FluoData(FitFrameFilter).*OnRatioDataForFit;
            SDFluoDataForFit=SDFluoData(FitFrameFilter).*OnRatioDataForFit;
            NDataForFit=NData(FitFrameFilter);
            TimeDataForFit=TimeData(FitFrameFilter);
            

            %These is the maximum range of data for the fit
            OnRatioData=OnRatioData(FrameFilter);
            MaxOnRatio=max(OnRatioData);
            OnRatioData=OnRatioData/MaxOnRatio;

            FluoData=FluoData(FrameFilter).*OnRatioData;
            SDFluoData=SDFluoData(FrameFilter).*OnRatioData;
            NData=NData(FrameFilter);
            TimeData=TimeData(FrameFilter);

            if CurrentNC~=14
                %Do the fit
                x0=[FitResults(i,CurrentNC-11).TimeStart0,FitResults(i,CurrentNC-11).TimeEnd0,FitResults(i,CurrentNC-11).Rate0];


                
                %Get rid of any NaN in the data
                NanFilter=~isnan(FluoDataForFit);

                if ~isempty(TimeData(NanFilter))

                    [xFit,resnorm,residual,exitflag,output,lambda,jacobian]=...
                        lsqnonlin(@(x) lsqnonlinFitFluorescenceCurve(TimeDataForFit(NanFilter)-ElapsedTimeTot(FrameWindow(1)),...
                        FluoDataForFit(NanFilter),Delay,...
                        ElapsedTimeTot(FrameWindow(end))-ElapsedTimeTot(FrameWindow(1)),x),x0);

                    FitResults(i,CurrentNC-11).TimeStart=xFit(1);
                    FitResults(i,CurrentNC-11).TimeEnd=xFit(2);
                    FitResults(i,CurrentNC-11).RateFit=xFit(3);

                    %Estimate an error bar out of the confidence intervals
                    FitResults(i,CurrentNC-11).CI=nlparci(xFit,residual,'jacobian',jacobian);

                    FitResults(i,CurrentNC-11).SDTimeStart=(FitResults(i,CurrentNC-11).CI(1,2)-FitResults(i,CurrentNC-11).CI(1,1))/2;
                    FitResults(i,CurrentNC-11).SDTimeEnd=(FitResults(i,CurrentNC-11).CI(2,2)-FitResults(i,CurrentNC-11).CI(2,1))/2;
                    FitResults(i,CurrentNC-11).SDRateFit=(FitResults(i,CurrentNC-11).CI(3,2)-FitResults(i,CurrentNC-11).CI(3,1))/2;

                    %Plot the results
                    %Get the corresponding fitted curve
                    [TimeFit,FluoFit]=FluorescenceCurve(ElapsedTimeTot(FrameWindow(end))-...
                        ElapsedTimeTot(FrameWindow(1)),...
                        xFit(1),xFit(2),xFit(3),Delay);
                    %Plot all the data
                    PlotHandle=errorbar(ElapsedTimeTot(FrameWindow)-ElapsedTimeTot(FrameWindow(1)),...
                        MeanVectorAP(FrameWindow,i).*OnRatioAP(FrameWindow,i)/MaxOnRatio,...
                        SDVectorAP(FrameWindow,i)./...
                        sqrt(NParticlesAP(FrameWindow,i)).*OnRatioAP(FrameWindow,i)/MaxOnRatio,'.-k');
                    hold on
                    %Plot the data that could be used for the fit
                    PlotHandle(end+1)=plot(ElapsedTimeTot(FrameWindow(FrameFilter))-ElapsedTimeTot(FrameWindow(1)),...
                        FluoData,'or');
                    %Plot the data that was actually used for the fit
                    PlotHandle(end+1)=plot(ElapsedTimeTot(FitFrameRange)-ElapsedTimeTot(FrameWindow(1)),...
                        FluoData(ismember(FrameWindow(FrameFilter),FitFrameRange)),'or','MarkerFaceColor','r');
                    
                    %Plot the fit
                    PlotHandle(end+1)=plot(TimeFit,FluoFit,'-r');
                    hold off
                    %StandardFigure(PlotHandle,gca)
                    ylabel('Mean fluorescence nucleus')
                    xlabel('Time into nc (min)')
                    
                    try
                        ylim([0,max(MeanVectorAPTot(FrameWindow,i).*OnRatioAPTot(FrameWindow,i)/MaxOnRatio+...
                            SDVectorAP(FrameWindow,i)./...
                            sqrt(NParticlesAPTot(FrameWindow,i)).*OnRatioAPTot(FrameWindow,i)/MaxOnRatio)])
                    catch
                        display('Error in displaying the plot')
                    end

                    legend(['tON=',num2str(FitResults(i,CurrentNC-11).TimeStart),' \pm ',num2str(FitResults(i,CurrentNC-11).SDTimeStart)],...
                        ['tOFF=',num2str(FitResults(i,CurrentNC-11).TimeEnd),' \pm ',num2str(FitResults(i,CurrentNC-11).SDTimeEnd)],...
                        ['Rate=',num2str(FitResults(i,CurrentNC-11).RateFit),' \pm ',num2str(FitResults(i,CurrentNC-11).SDRateFit)],...
                        'Location','SouthOutside')
                end
            elseif CurrentNC==14
                
                %Do the fit
                x0=[FitResults(i,CurrentNC-11).TimeStart0,FitResults(i,CurrentNC-11).Rate0];


                
                %Get rid of any NaN in the data
                NanFilter=~isnan(FluoDataForFit);

                if ~isempty(TimeData(NanFilter))

                    [xFit,resnorm,residual,exitflag,output,lambda,jacobian]=...
                        lsqnonlin(@(x) lsqnonlinFitFluorescenceCurveNC14(TimeDataForFit(NanFilter)-...
                        ElapsedTimeTot(FrameWindow(1)),...
                        FluoDataForFit(NanFilter),Delay,...
                        ElapsedTimeTot(FrameWindow(end))-ElapsedTimeTot(FrameWindow(1)),x),x0);

                    FitResults(i,CurrentNC-11).TimeStart=xFit(1);
                    FitResults(i,CurrentNC-11).RateFit=xFit(2);

                    %Estimate an error bar out of the confidence intervals
                    FitResults(i,CurrentNC-11).CI=nlparci(xFit,residual,'jacobian',jacobian);

                    FitResults(i,CurrentNC-11).SDTimeStart=(FitResults(i,CurrentNC-11).CI(1,2)-FitResults(i,CurrentNC-11).CI(1,1))/2;
                    FitResults(i,CurrentNC-11).SDRateFit=(FitResults(i,CurrentNC-11).CI(2,2)-FitResults(i,CurrentNC-11).CI(2,1))/2;




                    %Plot the results
                    %Get the corresponding fitted curve
                    [TimeFit,FluoFit]=FluorescenceCurve(ElapsedTimeTot(FrameWindow(end))-...
                        ElapsedTimeTot(FrameWindow(1)),...
                        xFit(1),1000,xFit(2),Delay);
                    %Plot all the data
                    PlotHandle=errorbar(ElapsedTimeTot(FrameWindow)-ElapsedTimeTot(FrameWindow(1)),...
                        MeanVectorAPTot(FrameWindow,i).*OnRatioAPTot(FrameWindow,i)/MaxOnRatio,...
                        SDVectorAPTot(FrameWindow,i)./...
                        sqrt(NParticlesAPTot(FrameWindow,i)).*OnRatioAPTot(FrameWindow,i)/MaxOnRatio,'.-k');
                    hold on
                    %Plot the data that could be used for the fit
                    PlotHandle(end+1)=plot(ElapsedTimeTot(FrameWindow(FrameFilter))-ElapsedTimeTot(FrameWindow(1)),...
                        FluoData,'or');
                    %Plot the data that was actually used for the fit
                    PlotHandle(end+1)=plot(ElapsedTimeTot(FitFrameRange)-ElapsedTimeTot(FrameWindow(1)),...
                        FluoData(ismember(FrameWindow(FrameFilter),FitFrameRange)),'or','MarkerFaceColor','r');
                    
                    %Plot the fit
                    PlotHandle(end+1)=plot(TimeFit,FluoFit,'-r');
                    hold off
                    ylabel('Mean fluorescence nucleus')
                    xlabel('Time into nc (min)')
                    
                    
                    try
                        ylim([0,max(MeanVectorAPTot(FrameWindow,i).*OnRatioAPTot(FrameWindow,i)/MaxOnRatio+...
                            SDVectorAPTot(FrameWindow,i)./...
                            sqrt(NParticlesAPTot(FrameWindow,i)).*OnRatioAPTot(FrameWindow,i)/MaxOnRatio)])
                    catch
                        display('Error in displaying the plot')
                    end

                    legend(['tON=',num2str(FitResults(i,CurrentNC-11).TimeStart),' \pm ',num2str(FitResults(i,CurrentNC-11).SDTimeStart)],...
                        ['Rate=',num2str(FitResults(i,CurrentNC-11).RateFit),' \pm ',num2str(FitResults(i,CurrentNC-11).SDRateFit)],...
                        'Location','SouthOutside')
                end
            end

        end
    end
    
    title([num2str(all_data(1).APbinID(i)),' AP, TimeStart0=',num2str(FitResults(i,CurrentNC-11).TimeStart0),...
        ', TimeEnd0=',num2str(FitResults(i,CurrentNC-11).TimeEnd0),', Rate=',num2str(FitResults(i,CurrentNC-11).Rate0),...
        ', nc',num2str(CurrentNC)])
    
    
    %Set the limits on the x-axis
    if CurrentNC==14
        %xlim([0,ElapsedTime(end)])
        xlim([0,60])
%     else
%         xlim([0,ElapsedTime(eval(['nc',num2str(nc+1)]))-...
%             ElapsedTime(eval(['nc',num2str(nc)]))])
    end
    
    
    
    figure(FitFigure)
    ct=waitforbuttonpress;
    cc=get(FitFigure,'currentcharacter');
    cm=get(gca,'CurrentPoint');
    
    %Move between AP positions
    if (ct~=0)&(cc=='.')&(i<length(all_data(1).APbinID))
        i=i+1;
    elseif (ct~=0)&(cc==',')&(i>1)
        i=i-1;
    
    %Approve, disapprove fit
    elseif (ct~=0)&(cc=='q')
        if FitResults(i,CurrentNC-11).Approved==0
            FitResults(i,CurrentNC-11).Approved=1;
        elseif FitResults(i,CurrentNC-11).Approved==1
            FitResults(i,CurrentNC-11).Approved=0;
        end

    
    %Disapprove, disapprove fit
    elseif (ct~=0)&(cc=='w')
        if FitResults(i,CurrentNC-11).Approved==0
            FitResults(i,CurrentNC-11).Approved=-1;
        elseif FitResults(i,CurrentNC-11).Approved==-1
            FitResults(i,CurrentNC-11).Approved=0;
        end
  
    
        
    %Move right range of fit
    elseif (ct~=0)&(cc=='k')&(length(FitResults(i,CurrentNC-11).FitFrameRange)>2)
        FitResults(i,CurrentNC-11).FitFrameRange=FitResults(i,CurrentNC-11).FitFrameRange(1:end-1);
    elseif (ct~=0)&(cc=='l')
        if ~isempty(find(~ismember(FrameWindow(FrameFilter),FitResults(i,CurrentNC-11).FitFrameRange)))
            FilteredFramesTemp=FrameWindow(FrameFilter);
            FitResults(i,CurrentNC-11).FitFrameRange(end+1)=...
                FilteredFramesTemp(min(find(~ismember(FilteredFramesTemp,FitResults(i,CurrentNC-11).FitFrameRange))));
        end
    %Move left range of fit
    elseif (ct~=0)&(cc=='j')&(length(FitResults(i,CurrentNC-11).FitFrameRange)>2)
        FitResults(i,CurrentNC-11).FitFrameRange=FitResults(i,CurrentNC-11).FitFrameRange(2:end);
    elseif (ct~=0)&(cc=='h')
        if ~isempty(find(~ismember(FrameWindow(FrameFilter),FitResults(i,CurrentNC-11).FitFrameRange)))
            FilteredFramesTemp=FrameWindow(FrameFilter);
            FitResults(i,CurrentNC-11).FitFrameRange=...
                [FilteredFramesTemp(max(find(~ismember(FilteredFramesTemp,FitResults(i,CurrentNC-11).FitFrameRange)))),...
                FitResults(i,CurrentNC-11).FitFrameRange];
        end
    %Reset frame fit range
     elseif (ct~=0)&(cc=='r')   
        FitResults(i,CurrentNC-11).FitFrameRange=FrameWindow(FrameFilter);

        
        
    %Change the initial parameters
    %TimeStart
    elseif (ct~=0)&(cc=='a')&((CurrentNC==14)|(CurrentNC~=14&FitResults(i,CurrentNC-11).TimeStart0<FitResults(i,CurrentNC-11).TimeEnd0))
        FitResults(i,CurrentNC-11).TimeStart0=FitResults(i,CurrentNC-11).TimeStart0+1;
    elseif (ct~=0)&(cc=='z')&(FitResults(i,CurrentNC-11).TimeStart0>1)
        FitResults(i,CurrentNC-11).TimeStart0=FitResults(i,CurrentNC-11).TimeStart0-1;
    %TimeEnd
    elseif (ct~=0)&(cc=='s')&(FitResults(i,CurrentNC-11).TimeEnd0<ElapsedTimeTot(FrameWindow(end))-ElapsedTimeTot(FrameWindow(1)))
        FitResults(i,CurrentNC-11).TimeEnd0=FitResults(i,CurrentNC-11).TimeEnd0+1;
    elseif (ct~=0)&(cc=='x')&(FitResults(i,CurrentNC-11).TimeEnd0>FitResults(i,CurrentNC-11).TimeStart0)
        FitResults(i,CurrentNC-11).TimeEnd0=FitResults(i,CurrentNC-11).TimeEnd0-1;
    %Rate, fine
    elseif (ct~=0)&(cc=='c')&(FitResults(i,CurrentNC-11).Rate0>100)
        FitResults(i,CurrentNC-11).Rate0=FitResults(i,CurrentNC-11).Rate0-100;
    elseif (ct~=0)&(cc=='d')
        FitResults(i,CurrentNC-11).Rate0=FitResults(i,CurrentNC-11).Rate0+100;    
    %Rate, coarse
    elseif (ct~=0)&(cc=='C')&(FitResults(i,CurrentNC-11).Rate0>100)
        FitResults(i,CurrentNC-11).Rate0=FitResults(i,CurrentNC-11).Rate0-500;
    elseif (ct~=0)&(cc=='D')
        FitResults(i,CurrentNC-11).Rate0=FitResults(i,CurrentNC-11).Rate0+500; 
    
    %Switch NCs
    elseif (ct~=0)&(cc=='m')&CurrentNC<14
        CurrentNC=CurrentNC+1;
    elseif (ct~=0)&(cc=='n')&CurrentNC>MinNC
        CurrentNC=CurrentNC-1;
        
        
    %Save
    elseif (ct~=0)&(cc=='v')
        save([DropboxFolder,filesep,Prefix,filesep,'MeanFits.mat'],...
        'FitResults')
    display('MeanFits.mat saved')
        
    %Debug mode
    elseif (ct~=0)&(cc=='9')
        keyboard
    
    end
    
end


%Save the information
save([DropboxFolder,filesep,Prefix,filesep,'MeanFits.mat'],...
    'FitResults')
display('MeanFits.mat saved')        
        
close(FitFigure)