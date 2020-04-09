function FitMeanAPTiltedTrapezoid(varargin)

% Edited from FitMeanAPAsymmetric by Yang Joon Kim (yjkim90@berkeley.edu)
%V3, separate the rate of turning on, in the middle, and rate of turning off. 
%I found that the MS2 fluorescence profile does not have a plateau, but
%tilted in the interphase, meaning that there are three phases, initiation
%with rate r1, then another rate of RNAP loading, r2, then stop of RNAP
%loading.
%This code will grab results from V1 and convert them to V2 in the first
%round. The mat files will not be backwards compatible.
%I also modified the errors to report the 68% CI.

%This function performs fits to the mean fluorescence as a function of time
%of a particular dataset.
%OUTPUT: MeanFitsV3.mat
%It gives you n columns each representing a nuclear cycle and m rows each
%representing a bin number

%Fitting:
%a,z: On time
%s,x: Off time
%d,c: Rate, fine
%D,C: Rate, coarse
%q,w: approve or disapprove a fit
%R : manual fitting, clicking 4 points, as described below in detail.
%e: save
%Moving around:
%, .: Move in AP
%n,m: Move in nc
%k,l: Change fit range from the right
%h,j: Change fit range from the left


%Parameters:
MinParticles=2;     %Minimum number of particles in an AP bin
MinTimePoints=5;    %Minimum number of time points where we'll have at least
                    %the minimum number of particles.
ElongationRate=1.54;    %In kb/minutes.
GeneLength=5.296;       %Distance from the first MS2 site to the end of the
                        %TUB3'UTR in kb.
Delay=GeneLength/ElongationRate;    %Minutes for PolII to fall off after reaching
                                    %the first MS2 site.

                                    
close all


%Get the default folders
[SourcePath,FISHPath,DropboxFolder,MS2CodePath]=...
    DetermineLocalFolders;

if ~isempty(varargin)
    Prefix=varargin{1};
               
else
    FolderTemp=uigetdir(DropboxFolder,'Choose folder with files to analyze');
    Dashes=strfind(FolderTemp,'\');
    Prefix=FolderTemp((Dashes(end)+1):end);
end

%Get the relevant folders now:
[SourcePath,FISHPath,DropboxFolder,MS2CodePath]=...
    DetermineLocalFolders(Prefix);
        


%Load the complied particles and the division information                                    
load([DropboxFolder,filesep,Prefix,'\CompiledParticles.mat'])

if exist([DropboxFolder,filesep,Prefix,'\APDivision.mat'])
    load([DropboxFolder,filesep,Prefix,'\APDivision.mat'])
else
    error('Could not load APDivision.mat. Make sure to have done the manual check of division.')
end

% YJK, 2018-12-00. I noticed that some datasets have MeanVectorAP,
% SDVectorAP, and NParticlesAP in a cell (1x1). For those, we should grab
% them out.

if iscell(MeanVectorAP) || iscell(SDVectorAP) || iscell(NParticlesAP)
    MeanVectorAP = cell2mat(MeanVectorAP);
    SDVectorAP = cell2mat(SDVectorAP);
    NParticlesAP = cell2mat(NParticlesAP);
end

%Initial parameters for fits. We will estimate the maximum rate based on
%the elongation time and the maximum average fluorescence of the data set.
MaxRate=max(max(MeanVectorAP))/Delay;

Rate012=MaxRate;     %Rate per minute
Rate012_1 = MaxRate/2; %Rate per minute
TimeStart012=3;
TimeEnd012=7;

Rate013=MaxRate;     %Rate per minute
Rate013_1 = MaxRate/2; %Rate per minute
TimeStart013=5;
TimeEnd013=10;

Rate014=MaxRate;     %Rate per minute
Rate014_1 = MaxRate/2; %Rate per minute
TimeStart014=5;
TimeEnd014=1000;                                    
                           


 
%Rough frame window to consider in the fits

%Some data sets won't have nc12
if nc12>0
    FrameWindow12=[nc12:nc13];
else
    FrameWindow12=[];
end
FrameWindow13=[nc13:nc14];
FrameWindow14=[nc14:length(ElapsedTime)];      


         

%Set the first guess for the parameters for each AP bin and also
%dissaproved the ones that did not have enough data points. Fit results has
%is a structure with the fits corresponding to each AP position and nc13
%or nc14
if exist([DropboxFolder,filesep,Prefix,filesep,'MeanFitsV3.mat'])
    load([DropboxFolder,filesep,Prefix,filesep,'MeanFitsV3.mat']);
    if isempty(FitResults)
        FitResults(length(APbinID),3).Rate0=[];
    end
else
    FitResults(length(APbinID),3).Rate0=[];
end





%Set default starting values for nc 13 and nc14
%nc12
for i=1:length(APbinID)
    if isempty(FitResults(i,1).Rate0)
        FitResults(i,1).Rate0=Rate012;
        FitResults(i,1).Rate1=Rate012_1;
        FitResults(i,1).TimeStart0=TimeStart012;
        FitResults(i,1).TimeEnd0=TimeEnd012;
        FitResults(i,1).RateOff0=-Rate012;
        FitResults(i,1).FrameFilter=[];
        FitResults(i,1).FitFrameRange=[];
        if sum(NParticlesAP(FrameWindow12,i)>=MinParticles)>=MinTimePoints
            FitResults(i,1).Approved=0;
        else
            FitResults(i,1).Approved=-1;
        end
    end
end
%nc13
for i=1:length(APbinID)
    if isempty(FitResults(i,2).Rate0)
        FitResults(i,2).Rate0=Rate013;    
        FitResults(i,2).Rate1=Rate013_1;
        FitResults(i,2).TimeStart0=TimeStart013;
        FitResults(i,2).TimeEnd0=TimeEnd013;
        FitResults(i,2).RateOff0=-Rate013;
        FitResults(i,2).FrameFilter=[];
        FitResults(i,2).FitFrameRange=[];
        if sum(NParticlesAP(FrameWindow13,i)>=MinParticles)>=MinTimePoints
            FitResults(i,2).Approved=0;
        else
            FitResults(i,2).Approved=-1;
        end
    end
end
%nc14
for i=1:length(APbinID)
    if isempty(FitResults(i,3).Rate0)
        FitResults(i,3).Rate0=Rate014;  
        FitResults(i,3).Rate1=Rate014_1;
        FitResults(i,3).TimeStart0=TimeStart014;
        FitResults(i,3).TimeEnd0=[];
        FitResults(i,3).RateOff0=[];
        FitResults(i,3).FrameFilter=[];
        FitResults(i,3).FitFrameRange=[];        
        if sum(NParticlesAP(FrameWindow14,i)>=MinParticles)>=MinTimePoints
            FitResults(i,3).Approved=0;
        else
            FitResults(i,3).Approved=-1;
        end
    end
end





         



%Go through each AP bin

FitFigure=figure;
CurrentNC=12;
i=min(find(sum(NParticlesAP)));
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
    
    
    if APDivision(CurrentNC,i)
        if CurrentNC~=14
            FrameWindow=APDivision(CurrentNC,i):APDivision(CurrentNC+1,i);
        else
            FrameWindow=APDivision(CurrentNC,i):length(ElapsedTime);
        end

        %Check that we have the minimum number of particles for a minimum
        %amount of time
        if (sum(NParticlesAP(FrameWindow,i)>=MinParticles)>=MinTimePoints)

            %Extract the data for this range of frames
            FluoData=MeanVectorAP(FrameWindow,i);
            SDFluoData=SDVectorAP(FrameWindow,i);
            NData=NParticlesAP(FrameWindow,i);
            TimeData=ElapsedTime(FrameWindow);
            %OnRatioData=OnRatioAP(FrameWindow,i);

            %Now filter them according the number of particles
            FrameFilter=NData>=MinParticles;


            %As an initial guess, use FrameFilter to determine the range of the
            %fit
            if isempty(FitResults(i,CurrentNC-11).FitFrameRange)
                FitFrameRange=FrameWindow(FrameFilter);
                if CurrentNC==14
                    FitFrameRange=FitFrameRange((ElapsedTime(FitFrameRange)-ElapsedTime(APDivision(CurrentNC,i)))<12);
                end
                FitResults(i,CurrentNC-11).FitFrameRange=FitFrameRange;
            else
                FitFrameRange=FitResults(i,CurrentNC-11).FitFrameRange;
            end

            %Filter the frames according to FitFrameRange
            FitFrameFilter=ismember(FrameWindow,FitFrameRange);
            
            
%             %HG: Removing the OnRatioFit stuff here. I'm not too happy
%             %about it.
%             OnRatioDataForFit=OnRatioData(FitFrameFilter);
%             MaxOnRatioForFit=max(OnRatioData);
%             OnRatioDataForFit=OnRatioDataForFit/MaxOnRatioForFit;

            FluoDataForFit=FluoData(FitFrameFilter);%.*OnRatioDataForFit;
            SDFluoDataForFit=SDFluoData(FitFrameFilter);%.*OnRatioDataForFit;
            NDataForFit=NData(FitFrameFilter);
            TimeDataForFit=TimeData(FitFrameFilter);
            

%             %These is the maximum range of data for the fit
%             OnRatioData=OnRatioData(FrameFilter);
%             MaxOnRatio=max(OnRatioData);
%             OnRatioData=OnRatioData/MaxOnRatio;

            FluoData=FluoData(FrameFilter)%.*OnRatioData;
            SDFluoData=SDFluoData(FrameFilter)%.*OnRatioData;
            NData=NData(FrameFilter);
            TimeData=TimeData(FrameFilter);

            if CurrentNC~=14
                %Do the fit
                x0=[FitResults(i,CurrentNC-11).TimeStart0,...
                    FitResults(i,CurrentNC-11).TimeEnd0,...
                    FitResults(i,CurrentNC-11).Rate0,...
                    FitResults(i,CurrentNC-11).RateOff0,...
                    FitResults(i,CurrentNC-11).Rate1];


                
                %Get rid of any NaN in the data
                NanFilter=~isnan(FluoDataForFit);

                if ~isempty(TimeData(NanFilter))

                    [xFit,resnorm,residual,exitflag,output,lambda,jacobian]=...
                        lsqnonlin(@(x) lsqnonlinFitFluorescenceCurveV3(TimeDataForFit(NanFilter)-...
                        ElapsedTime(FrameWindow(1)),...
                        FluoDataForFit(NanFilter),Delay,...
                        ElapsedTime(FrameWindow(end))-ElapsedTime(FrameWindow(1)),x),x0);

                    FitResults(i,CurrentNC-11).TimeStart=xFit(1);
                    FitResults(i,CurrentNC-11).TimeEnd=xFit(2);
                    FitResults(i,CurrentNC-11).RateFit=xFit(3);
                    FitResults(i,CurrentNC-11).RateOffFit=xFit(4);

                    %Estimate an error bar out of the confidence intervals
                    FitResults(i,CurrentNC-11).CI=nlparci(xFit,residual,'jacobian',jacobian,'alpha',0.68);

                    FitResults(i,CurrentNC-11).SDTimeStart=(FitResults(i,CurrentNC-11).CI(1,2)-FitResults(i,CurrentNC-11).CI(1,1))/2;
                    FitResults(i,CurrentNC-11).SDTimeEnd=(FitResults(i,CurrentNC-11).CI(2,2)-FitResults(i,CurrentNC-11).CI(2,1))/2;
                    FitResults(i,CurrentNC-11).SDRateFit=(FitResults(i,CurrentNC-11).CI(3,2)-FitResults(i,CurrentNC-11).CI(3,1))/2;
                    FitResults(i,CurrentNC-11).SDRateOffFit=(FitResults(i,CurrentNC-11).CI(4,2)-FitResults(i,CurrentNC-11).CI(4,1))/2;


                    %Plot the results
                    %Get the corresponding fitted curve
                    [TimeFit,FluoFit]=FluorescenceCurveV2(ElapsedTime(FrameWindow(end))-...
                        ElapsedTime(FrameWindow(1)),...
                        xFit(1),xFit(2),xFit(3),xFit(4),Delay);
                    %Plot all the data
                    PlotHandle=errorbar(ElapsedTime(FrameWindow)-ElapsedTime(FrameWindow(1)),...
                        MeanVectorAP(FrameWindow,i),...%.*OnRatioAP(FrameWindow,i)/MaxOnRatio,...
                        SDVectorAP(FrameWindow,i)./...
                        sqrt(NParticlesAP(FrameWindow,i)),'.-k');%.*OnRatioAP(FrameWindow,i)/MaxOnRatio,'.-k');
                    hold on
                    %Plot the data that could be used for the fit
                    PlotHandle(end+1)=plot(ElapsedTime(FrameWindow(FrameFilter))-ElapsedTime(FrameWindow(1)),...
                        FluoData,'or');
                    %Plot the data that was actually used for the fit
                    PlotHandle(end+1)=plot(ElapsedTime(FitFrameRange)-ElapsedTime(FrameWindow(1)),...
                        FluoData(ismember(FrameWindow(FrameFilter),FitFrameRange)),'or','MarkerFaceColor','r');
                    
                    %Plot the fit
                    PlotHandle(end+1)=plot(TimeFit,FluoFit,'-r');
                    hold off
                    %StandardFigure(PlotHandle,gca)
                    ylabel('Mean fluorescence nucleus')
                    xlabel('Time into nc (min)')
                    
                    try
                        ylim([0,max(MeanVectorAP(FrameWindow,i)+...%.*OnRatioAP(FrameWindow,i)/MaxOnRatio+...
                            SDVectorAP(FrameWindow,i)./...
                            sqrt(NParticlesAP(FrameWindow,i)))])%.*OnRatioAP(FrameWindow,i)/MaxOnRatio)])
                    catch
                        display('Error in displaying the plot')
                    end

                    legend(PlotHandle,['tON=',num2str(FitResults(i,CurrentNC-11).TimeStart),' \pm ',num2str(FitResults(i,CurrentNC-11).SDTimeStart)],...
                        ['tOFF=',num2str(FitResults(i,CurrentNC-11).TimeEnd),' \pm ',num2str(FitResults(i,CurrentNC-11).SDTimeEnd)],...
                        ['Rate=',num2str(FitResults(i,CurrentNC-11).RateFit),' \pm ',num2str(FitResults(i,CurrentNC-11).SDRateFit)],...
                        ['RateOff=',num2str(FitResults(i,CurrentNC-11).RateOffFit),' \pm ',num2str(FitResults(i,CurrentNC-11).SDRateOffFit)],...
                        'Location','Best')
                    
  
                end
            elseif CurrentNC==14
                
                %Do the fit
                x0=[FitResults(i,CurrentNC-11).TimeStart0,FitResults(i,CurrentNC-11).Rate0];


                
                %Get rid of any NaN in the data
                NanFilter=~isnan(FluoDataForFit);

                if ~isempty(TimeData(NanFilter))

                    [xFit,resnorm,residual,exitflag,output,lambda,jacobian]=...
                        lsqnonlin(@(x) lsqnonlinFitFluorescenceCurveNC14(TimeDataForFit(NanFilter)-...
                        ElapsedTime(FrameWindow(1)),...
                        FluoDataForFit(NanFilter),Delay,...
                        ElapsedTime(FrameWindow(end))-ElapsedTime(FrameWindow(1)),x),x0);

                    FitResults(i,CurrentNC-11).TimeStart=xFit(1);
                    FitResults(i,CurrentNC-11).RateFit=xFit(2);

                    %Estimate an error bar out of the confidence intervals
                    FitResults(i,CurrentNC-11).CI=nlparci(xFit,residual,'jacobian',jacobian);

                    FitResults(i,CurrentNC-11).SDTimeStart=(FitResults(i,CurrentNC-11).CI(1,2)-FitResults(i,CurrentNC-11).CI(1,1))/2;
                    FitResults(i,CurrentNC-11).SDRateFit=(FitResults(i,CurrentNC-11).CI(2,2)-FitResults(i,CurrentNC-11).CI(2,1))/2;




                    %Plot the results
                    %Get the corresponding fitted curve
                    [TimeFit,FluoFit]=FluorescenceCurve(ElapsedTime(FrameWindow(end))-...
                        ElapsedTime(FrameWindow(1)),...
                        xFit(1),1000,xFit(2),Delay);
                    %Plot all the data
                    PlotHandle=errorbar(ElapsedTime(FrameWindow)-ElapsedTime(FrameWindow(1)),...
                        MeanVectorAP(FrameWindow,i),...%.*OnRatioAP(FrameWindow,i)/MaxOnRatio,...
                        SDVectorAP(FrameWindow,i)./...
                        sqrt(NParticlesAP(FrameWindow,i)),'.-k');%.*OnRatioAP(FrameWindow,i)/MaxOnRatio,'.-k');
                    hold on
                    %Plot the data that could be used for the fit
                    PlotHandle(end+1)=plot(ElapsedTime(FrameWindow(FrameFilter))-ElapsedTime(FrameWindow(1)),...
                        FluoData,'or');
                    %Plot the data that was actually used for the fit
                    PlotHandle(end+1)=plot(ElapsedTime(FitFrameRange)-ElapsedTime(FrameWindow(1)),...
                        FluoData(ismember(FrameWindow(FrameFilter),FitFrameRange)),'or','MarkerFaceColor','r');
                    
                    %Plot the fit
                    PlotHandle(end+1)=plot(TimeFit,FluoFit,'-r');
                    hold off
                    ylabel('Mean fluorescence nucleus')
                    xlabel('Time into nc (min)')
                    
                    
                    try
                        ylim([0,max(MeanVectorAP(FrameWindow,i)+...%.*OnRatioAP(FrameWindow,i)/MaxOnRatio+...
                            SDVectorAP(FrameWindow,i)./...
                            sqrt(NParticlesAP(FrameWindow,i)))])%.*OnRatioAP(FrameWindow,i)/MaxOnRatio)])
                    catch
                        display('Error in displaying the plot')
                    end

                    legend(PlotHandle,['tON=',num2str(FitResults(i,CurrentNC-11).TimeStart),' \pm ',num2str(FitResults(i,CurrentNC-11).SDTimeStart)],...
                        ['Rate=',num2str(FitResults(i,CurrentNC-11).RateFit),' \pm ',num2str(FitResults(i,CurrentNC-11).SDRateFit)],...
                        ['RateOff=',num2str(FitResults(i,CurrentNC-11).RateOffFit),' \pm ',num2str(FitResults(i,CurrentNC-11).SDRateOffFit)],...
                        'Location','Best')
                end
            end

        end
    end
    
    title([num2str(APbinID(i)),' AP, tStart0=',num2str(FitResults(i,CurrentNC-11).TimeStart0),...
        ', tEnd0=',num2str(FitResults(i,CurrentNC-11).TimeEnd0),', Rate=',num2str(FitResults(i,CurrentNC-11).Rate0),...
        ', RateOff=',num2str(FitResults(i,CurrentNC-11).RateOff0),', nc',num2str(CurrentNC)])
    
    
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
    if (ct~=0)&(cc=='.')&(i<length(APbinID))
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
    elseif (ct~=0)&(cc=='s')&(FitResults(i,CurrentNC-11).TimeEnd0<ElapsedTime(FrameWindow(end))-ElapsedTime(FrameWindow(1)))
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
    %RateOff, fine
    elseif (ct~=0)&(cc=='v')&(FitResults(i,CurrentNC-11).RateOff0<-100)
        FitResults(i,CurrentNC-11).RateOff0=FitResults(i,CurrentNC-11).RateOff0-100;
    elseif (ct~=0)&(cc=='f')
        FitResults(i,CurrentNC-11).RateOff0=FitResults(i,CurrentNC-11).RateOff0+100;    
    %RateOff, coarse
    elseif (ct~=0)&(cc=='V')&(FitResults(i,CurrentNC-11).RateOff0<-100)
        FitResults(i,CurrentNC-11).RateOff0=FitResults(i,CurrentNC-11).RateOff0-500;
    elseif (ct~=0)&(cc=='F')
        FitResults(i,CurrentNC-11).RateOff0=FitResults(i,CurrentNC-11).RateOff0+500;  
        
    % Options that YJK added for manual fitting
    elseif (ct~=0)&(cc=='R')
        % click 4 points, inital rise, peak, 2nd pick, and down slope.
        [X,Y] = ginput(4)
        FitResults(i,CurrentNC-11).Rate0 = (Y(2)-Y(1))/(X(2)-X(1));
        FitResults(i,CurrentNC-11).TimeStart0 = (Y(1)*(X(1)-X(2)) + X(1)*(Y(2)-Y(1))) / (Y(2) - Y(1));
        FitResults(i,CurrentNC-11).TimeEnd0 = X(3);
        FitResults(i,CurrentNC-11).RateOff0 = (Y(4)-Y(3))/(X(4)-X(3));
        FitResults(i,CurrentNC-11).Rate1 = (Y(3)-Y(2))/(X(3)-X(2)); %- (Y(4)-Y(3))/(X(4)-X(3));
        FitResults(i,CurrentNC-11).TimeSwitch1 = X(2); % This is the time point when the fluo starts to drop.
    
    %Switch NCs
    elseif (ct~=0)&(cc=='m')&CurrentNC<14
        CurrentNC=CurrentNC+1;
    elseif (ct~=0)&(cc=='n')&CurrentNC>12
        CurrentNC=CurrentNC-1;
        
        
    %Save
    elseif (ct~=0)&(cc=='e')
        save([DropboxFolder,filesep,Prefix,filesep,'MeanFitsV3.mat'],...
        'FitResults')
    display('MeanFitsV3.mat saved')
        
    %Debug mode
    elseif (ct~=0)&(cc=='9')
        keyboard
    
    end
    
end


%Save the information
save([DropboxFolder,filesep,Prefix,filesep,'MeanFitsV3.mat'],...
    'FitResults')
display('MeanFitsV3.mat saved')        
        
close(FitFigure)