function FitResults = FitMeanAll(varargin)
%function FitResults = FitMeanAll(varargin)
%
%DESCRIPTION
%This function performs fits to the mean fluorescence as a function of time
%of a particular dataset.
%It performs the fit by averaging over all traces. This is useful for the
%NoBcd controls
%
%ARGUMENTS
%varargin: This should be a string containing the prefix for this dataset.
%
%OPTIONS
%No options
%
%CONTROLS
%%Fitting:
%a,z: On time
%s,x: Off time
%d,c: Rate, fine
%D,C: Rate, coarse
%
%Moving around:
%n,m: Move in nc
%k,l: Change fit range from the right
%h,j: Change fit range from the left
%
%OUTPUT
%FitResults: A structure containing fitted rates, on times, and off times for
%average traces. This is saved in MeanFits.mat
%
%Author (contact): Hernan Garcia (hgarcia@berkeley.edu)
%Created: Unknown
%Last Updated: Unknown
%
%
%Some things I (HG) need to do:
%-Align the different traces according to APdivisions and not to the regular
% nc times specified in the XLS spreadsheet
%-Do something to account for the on ratio as I do in the case of the
% single AP bins


%Parameters:
MinParticles=2;     %Minimum number of particles in an AP bin
MinTimePoints=5;    %Minimum number of time points where we'll have at least
                    %the minimum number of particles.
ElongationRate=1.54;    %In kb/minutes.
GeneLength=5.296;       %Distance from the first MS2 site to the end of the
                        %TUB3'UTR in kb.
Delay=GeneLength/ElongationRate;    %Minutes for PolII to fall off after reaching
                                    %the first MS2 site.

                                    

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
                                    
                                    
                                    
                                    
                                    
close all

%Find out which computer this is. That will determine the folder structure.
%Information about about folders
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,SchnitzcellsFolder]=...
    DetermineLocalFolders;    

if ~isempty(varargin)
    Prefix=varargin{1};
               
else
    FolderTemp=uigetdir(DropboxFolder,'Choose folder with files to analyze');
    Dashes=strfind(FolderTemp,'\');
    Prefix=FolderTemp((Dashes(end)+1):end);
end


[SourcePath,FISHPath,DropboxFolder,MS2CodePath,SchnitzcellsFolder]=...
    DetermineLocalFolders(Prefix);    
                                    
%Load the complied particles and the division information                                    
load([DropboxFolder,filesep,Prefix,'\CompiledParticles.mat'])


%Rough frame window to consider in the fits

%Some data sets won't have nc12
if nc12>0
    FrameWindow12=[nc12:nc13];
else
    FrameWindow12=[];
end
if nc13>0
    FrameWindow13=[nc13:nc14];
else
    FrameWindow13=[];
end
FrameWindow14=[nc14:length(ElapsedTime)];      


         

%Set the first guess for the parameters for each AP bin and also
%dissaproved the ones that did not have enough data points. Fit results has
%is a structure with the fits corresponding to each AP position and nc13
%or nc14
if exist([DropboxFolder,filesep,Prefix,filesep,'MeanFits.mat'])
    load([DropboxFolder,filesep,Prefix,filesep,'MeanFits.mat']);
    if isempty(FitResults)
        FitResults(3).Rate0=[];
    end
else
    FitResults(3).Rate0=[];
end





%Set default starting values for nc 13 and nc14
%nc12
if isempty(FitResults(1).Rate0)
    FitResults(1).Rate0=Rate012;    
    FitResults(1).TimeStart0=TimeStart012;
    FitResults(1).TimeEnd0=TimeEnd012;
    FitResults(1).FrameFilter=[];
    FitResults(1).FitFrameRange=[];
    if sum(NParticlesAll(FrameWindow12)>=MinParticles)>=MinTimePoints
        FitResults(1).Approved=0;
    else
        FitResults(1).Approved=-1;
    end
end

%nc13
if isempty(FitResults(2).Rate0)
    FitResults(2).Rate0=Rate013;    
    FitResults(2).TimeStart0=TimeStart013;
    FitResults(2).TimeEnd0=TimeEnd013;
    FitResults(2).FrameFilter=[];
    FitResults(2).FitFrameRange=[];
    if sum(NParticlesAll(FrameWindow13)>=MinParticles)>=MinTimePoints
        FitResults(2).Approved=0;
    else
        FitResults(2).Approved=-1;
    end
end

%nc14
if isempty(FitResults(3).Rate0)
    FitResults(3).Rate0=Rate014;    
    FitResults(3).TimeStart0=TimeStart014;
    FitResults(3).TimeEnd0=[];;
    FitResults(3).FrameFilter=[];
    FitResults(3).FitFrameRange=[];        
    if sum(NParticlesAll(FrameWindow14)>=MinParticles)>=MinTimePoints
        FitResults(3).Approved=0;
    else
        FitResults(3).Approved=-1;
    end
end





%Go through each AP bin

FitFigure=figure;
CurrentNC=12;
cc=1;

while (cc~=13)
    
    figure(FitFigure)
    clf
    
    if FitResults(CurrentNC-11).Approved==-1
        set(gcf,'Color','r')
    elseif FitResults(CurrentNC-11).Approved==1
        set(gcf,'Color','g')
    else
        set(gcf,'Color','default')
    end
    
    
    if eval(['nc',num2str(CurrentNC)])
        if CurrentNC~=14
            FrameWindow=eval(['nc',num2str(CurrentNC)]):...
                eval(['nc',num2str(CurrentNC+1)]);
        else
            FrameWindow=eval(['nc',num2str(CurrentNC)]):length(ElapsedTime);
        end

        %Check that we have the minimum number of particles for a minimum
        %amount of time
        if (sum(NParticlesAll(FrameWindow)>=MinParticles)>=MinTimePoints)

            %Extract the data for this range of frames
            FluoData=MeanVectorAll(FrameWindow);
            SDFluoData=SDVectorAll(FrameWindow);
            NData=NParticlesAll(FrameWindow);
            TimeData=ElapsedTime(FrameWindow);
            OnRatioData=ones(size(FrameWindow));%OnRatioAP(FrameWindow,:);

            %Now filter them according the number of particles
            FrameFilter=NData>=MinParticles;


            %As an initial guess, use FrameFilter to determine the range of the
            %fit
            if isempty(FitResults(CurrentNC-11).FitFrameRange)
                FitFrameRange=FrameWindow(FrameFilter);
                if CurrentNC==14
                    FitFrameRange=FitFrameRange((ElapsedTime(FitFrameRange)-ElapsedTime(eval(['nc',num2str(CurrentNC)]))<12));
                end
                FitResults(CurrentNC-11).FitFrameRange=FitFrameRange;
            else
                FitFrameRange=FitResults(CurrentNC-11).FitFrameRange;
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
                x0=[FitResults(CurrentNC-11).TimeStart0,FitResults(CurrentNC-11).TimeEnd0,FitResults(CurrentNC-11).Rate0];


                
                %Get rid of any NaN in the data
                NanFilter=~isnan(FluoDataForFit);

                if ~isempty(TimeData(NanFilter))

                    [xFit,resnorm,residual,exitflag,output,lambda,jacobian]=...
                        lsqnonlin(@(x) lsqnonlinFitFluorescenceCurve(TimeDataForFit(NanFilter)-...
                        ElapsedTime(FrameWindow(1)),...
                        FluoDataForFit(NanFilter)',Delay,...
                        ElapsedTime(FrameWindow(end))-ElapsedTime(FrameWindow(1)),x),x0);

                    FitResults(CurrentNC-11).TimeStart=xFit(1);
                    FitResults(CurrentNC-11).TimeEnd=xFit(2);
                    FitResults(CurrentNC-11).RateFit=xFit(3);

                    %Estimate an error bar out of the confidence intervals
                    FitResults(CurrentNC-11).CI=nlparci(xFit,residual,'jacobian',jacobian);

                    FitResults(CurrentNC-11).SDTimeStart=(FitResults(CurrentNC-11).CI(1,2)-FitResults(CurrentNC-11).CI(1,1))/2;
                    FitResults(CurrentNC-11).SDTimeEnd=(FitResults(CurrentNC-11).CI(2,2)-FitResults(CurrentNC-11).CI(2,1))/2;
                    FitResults(CurrentNC-11).SDRateFit=(FitResults(CurrentNC-11).CI(3,2)-FitResults(CurrentNC-11).CI(3,1))/2;




                    %Plot the results
                    %Get the corresponding fitted curve
                    [TimeFit,FluoFit]=FluorescenceCurve(ElapsedTime(FrameWindow(end))-...
                        ElapsedTime(FrameWindow(1)),...
                        xFit(1),xFit(2),xFit(3),Delay);
                    %Plot all the data
%                     PlotHandle=errorbar(ElapsedTime(FrameWindow)-ElapsedTime(FrameWindow(1)),...
%                         MeanVectorAP(FrameWindow,i).*OnRatioAP(FrameWindow,i)/MaxOnRatio,...
%                         SDVectorAP(FrameWindow,i)./...
%                         sqrt(NParticlesAP(FrameWindow,i)).*OnRatioAP(FrameWindow,i)/MaxOnRatio,'.-k');
                    PlotHandle=errorbar(ElapsedTime(FrameWindow)-ElapsedTime(FrameWindow(1)),...
                        MeanVectorAll(FrameWindow),...
                        SDVectorAll(FrameWindow)./...
                        sqrt(NParticlesAll(FrameWindow)),'.-k');
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
%                         ylim([0,max(MeanVectorAP(FrameWindow,i).*OnRatioAP(FrameWindow,i)/MaxOnRatio+...
%                             SDVectorAP(FrameWindow,i)./...
%                             sqrt(NParticlesAP(FrameWindow,i)).*OnRatioAP(FrameWindow,i)/MaxOnRatio)])
                        ylim([0,max(MeanVectorAll(FrameWindow)+...
                            SDVectorAll(FrameWindow)./...
                            sqrt(NParticlesAll(FrameWindow)))])
                    catch
                        display('Error in displaying the plot')
                    end

                    legend(['tON=',num2str(FitResults(CurrentNC-11).TimeStart),' \pm ',num2str(FitResults(CurrentNC-11).SDTimeStart)],...
                        ['tOFF=',num2str(FitResults(CurrentNC-11).TimeEnd),' \pm ',num2str(FitResults(CurrentNC-11).SDTimeEnd)],...
                        ['Rate=',num2str(FitResults(CurrentNC-11).RateFit),' \pm ',num2str(FitResults(CurrentNC-11).SDRateFit)],...
                        'Location','Best')
                end
            elseif CurrentNC==14
                
                %Do the fit
                x0=[FitResults(CurrentNC-11).TimeStart0,FitResults(CurrentNC-11).Rate0];


                
                %Get rid of any NaN in the data
                NanFilter=~isnan(FluoDataForFit);

                if ~isempty(TimeData(NanFilter))

                    [xFit,resnorm,residual,exitflag,output,lambda,jacobian]=...
                        lsqnonlin(@(x) lsqnonlinFitFluorescenceCurveNC14(TimeDataForFit(NanFilter)-...
                        ElapsedTime(FrameWindow(1)),...
                        FluoDataForFit(NanFilter)',Delay,...
                        ElapsedTime(FrameWindow(end))-ElapsedTime(FrameWindow(1)),x),x0);

                    FitResults(CurrentNC-11).TimeStart=xFit(1);
                    FitResults(CurrentNC-11).RateFit=xFit(2);

                    %Estimate an error bar out of the confidence intervals
                    FitResults(CurrentNC-11).CI=nlparci(xFit,residual,'jacobian',jacobian);

                    FitResults(CurrentNC-11).SDTimeStart=(FitResults(CurrentNC-11).CI(1,2)-FitResults(CurrentNC-11).CI(1,1))/2;
                    FitResults(CurrentNC-11).SDRateFit=(FitResults(CurrentNC-11).CI(2,2)-FitResults(CurrentNC-11).CI(2,1))/2;




                    %Plot the results
                    %Get the corresponding fitted curve
                    [TimeFit,FluoFit]=FluorescenceCurve(ElapsedTime(FrameWindow(end))-...
                        ElapsedTime(FrameWindow(1)),...
                        xFit(1),1000,xFit(2),Delay);
                    %Plot all the data
%                     PlotHandle=errorbar(ElapsedTime(FrameWindow)-ElapsedTime(FrameWindow(1)),...
%                         MeanVectorAP(FrameWindow,i).*OnRatioAP(FrameWindow,i)/MaxOnRatio,...
%                         SDVectorAP(FrameWindow,i)./...
%                         sqrt(NParticlesAP(FrameWindow,i)).*OnRatioAP(FrameWindow,i)/MaxOnRatio,'.-k');
                    PlotHandle=errorbar(ElapsedTime(FrameWindow)-ElapsedTime(FrameWindow(1)),...
                        MeanVectorAll(FrameWindow),...
                        SDVectorAll(FrameWindow)./...
                        sqrt(NParticlesAll(FrameWindow)),'.-k');
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
%                         ylim([0,max(MeanVectorAP(FrameWindow,i).*OnRatioAP(FrameWindow,i)/MaxOnRatio+...
%                             SDVectorAP(FrameWindow,i)./...
%                             sqrt(NParticlesAP(FrameWindow,i)).*OnRatioAP(FrameWindow,i)/MaxOnRatio)])
                        ylim([0,max(MeanVectorAll(FrameWindow)+...
                            SDVectorAll(FrameWindow)./...
                            sqrt(NParticlesAll(FrameWindow)))])
                    catch
                        display('Error in displaying the plot')
                    end

                    legend(['tON=',num2str(FitResults(CurrentNC-11).TimeStart),' \pm ',num2str(FitResults(CurrentNC-11).SDTimeStart)],...
                        ['Rate=',num2str(FitResults(CurrentNC-11).RateFit),' \pm ',num2str(FitResults(CurrentNC-11).SDRateFit)],...
                        'Location','Best')
                end
            end

        end
    end
    
    title(['TimeStart0=',num2str(FitResults(CurrentNC-11).TimeStart0),...
        ', TimeEnd0=',num2str(FitResults(CurrentNC-11).TimeEnd0),', Rate=',num2str(FitResults(CurrentNC-11).Rate0),...
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
    
    
    %Approve, disapprove fit
    if (ct~=0)&(cc=='q')
        if FitResults(CurrentNC-11).Approved==0
            FitResults(CurrentNC-11).Approved=1;
        elseif FitResults(CurrentNC-11).Approved==1
            FitResults(CurrentNC-11).Approved=0;
        end

    
    %Disapprove, disapprove fit
    elseif (ct~=0)&(cc=='w')
        if FitResults(CurrentNC-11).Approved==0
            FitResults(CurrentNC-11).Approved=-1;
        elseif FitResults(CurrentNC-11).Approved==-1
            FitResults(CurrentNC-11).Approved=0;
        end
  
    
        
    %Move right range of fit
    elseif (ct~=0)&(cc=='k')&(length(FitResults(CurrentNC-11).FitFrameRange)>2)
        FitResults(CurrentNC-11).FitFrameRange=FitResults(CurrentNC-11).FitFrameRange(1:end-1);
    elseif (ct~=0)&(cc=='l')
        if ~isempty(find(~ismember(FrameWindow(FrameFilter),FitResults(CurrentNC-11).FitFrameRange)))
            FilteredFramesTemp=FrameWindow(FrameFilter);
            FitResults(CurrentNC-11).FitFrameRange(end+1)=...
                FilteredFramesTemp(min(find(~ismember(FilteredFramesTemp,FitResults(CurrentNC-11).FitFrameRange))));
        end
    %Move left range of fit
    elseif (ct~=0)&(cc=='j')&(length(FitResults(CurrentNC-11).FitFrameRange)>2)
        FitResults(CurrentNC-11).FitFrameRange=FitResults(CurrentNC-11).FitFrameRange(2:end);
    elseif (ct~=0)&(cc=='h')
        if ~isempty(find(~ismember(FrameWindow(FrameFilter),FitResults(CurrentNC-11).FitFrameRange)))
            FilteredFramesTemp=FrameWindow(FrameFilter);
            FitResults(CurrentNC-11).FitFrameRange=...
                [FilteredFramesTemp(max(find(~ismember(FilteredFramesTemp,FitResults(CurrentNC-11).FitFrameRange)))),...
                FitResults(CurrentNC-11).FitFrameRange];
        end
    %Reset frame fit range
     elseif (ct~=0)&(cc=='r')   
        FitResults(CurrentNC-11).FitFrameRange=FrameWindow(FrameFilter);

        
        
    %Change the initial parameters
    %TimeStart
    elseif (ct~=0)&(cc=='a')&((CurrentNC==14)|(CurrentNC~=14&FitResults(CurrentNC-11).TimeStart0<FitResults(CurrentNC-11).TimeEnd0))
        FitResults(CurrentNC-11).TimeStart0=FitResults(CurrentNC-11).TimeStart0+1;
    elseif (ct~=0)&(cc=='z')&(FitResults(CurrentNC-11).TimeStart0>1)
        FitResults(CurrentNC-11).TimeStart0=FitResults(CurrentNC-11).TimeStart0-1;
    %TimeEnd
    elseif (ct~=0)&(cc=='s')&(FitResults(CurrentNC-11).TimeEnd0<ElapsedTime(FrameWindow(end))-ElapsedTime(FrameWindow(1)))
        FitResults(CurrentNC-11).TimeEnd0=FitResults(CurrentNC-11).TimeEnd0+1;
    elseif (ct~=0)&(cc=='x')&(FitResults(CurrentNC-11).TimeEnd0>FitResults(CurrentNC-11).TimeStart0)
        FitResults(CurrentNC-11).TimeEnd0=FitResults(CurrentNC-11).TimeEnd0-1;
    %Rate, fine
    elseif (ct~=0)&(cc=='c')&(FitResults(CurrentNC-11).Rate0>100)
        FitResults(CurrentNC-11).Rate0=FitResults(CurrentNC-11).Rate0-100;
    elseif (ct~=0)&(cc=='d')
        FitResults(CurrentNC-11).Rate0=FitResults(CurrentNC-11).Rate0+100;    
    %Rate, coarse
    elseif (ct~=0)&(cc=='C')&(FitResults(CurrentNC-11).Rate0>100)
        FitResults(CurrentNC-11).Rate0=FitResults(CurrentNC-11).Rate0-500;
    elseif (ct~=0)&(cc=='D')
        FitResults(CurrentNC-11).Rate0=FitResults(CurrentNC-11).Rate0+500; 
    
    %Switch NCs
    elseif (ct~=0)&(cc=='m')&CurrentNC<14
        CurrentNC=CurrentNC+1;
    elseif (ct~=0)&(cc=='n')&CurrentNC>12
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