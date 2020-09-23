function FitMeanAP_LinearSlope_ExponentialDecay_NC14(varargin)

%DESCRIPTION
% This code is for generating FigS in Garcia 2013 paper.
% Basically, we want to fit a linear slope in the initial phase, then fit
% an exponential decay curve on the "decay regime" after the peak (or
% plauteu if there's any?), assuming that the decay is happening in an
% exponential scale. Then, the decay half-life (tau) can be used to
% calculate the "turn-off time".

% probably from HG : I also modified the errors to report the 68% CI.

% This function performs fits to the mean fluorescence as a function of time
% of a particular dataset.

% INPUT : 
%  Prefix
% 'NCWindow', [nc14start, nc14end]: Start and end
%                           times of fitting for cycle 14.
%                           Pass this in as a 1x2 array of doubles.
%                           Read as [time after nc14 start, time after nc14 start].
%                           Default values are [1, 25].

% OUTPUT: MeanFits_LinearRise_ExponentDecay.mat

% It gives you 1 column representing a nuclear cycle, nc14 and m rows each
% representing a bin number

% Fitting: (This needs editing later)
%a,z: On time
%s,x: Off time
%d,c: Rate, fine
%D,C: Rate, coarse
%q,w: approve or disapprove a fit
%S : save the fitted plot
%e: save
%Moving around:
%, .: Move in AP
%n,m: Move in nc
%k,l: Change fit range from the right
%h,j: Change fit range from the left


%% Variable inputs

%Default settings
LoadPrefix = true; %User selection of Prefix by default.
ncwindow = [5, 30]; %Time of fitting for nc14, [min]
KeepPool = false; %By default, shut down parallel pool after running code

%User specified settings
UseLocalMovieDatabase = false; %By default, don't use a local MovieDatabase
for i=1:length(varargin)
    if strcmpi(varargin{i},'Prefix')
        Prefix = varargin{i+1};
        LoadPrefix = false;
    end
    if strcmpi(varargin{i},'NCWindow')
        ncwindow = varargin{i+1};
    end
    if strcmpi(varargin{i},'LocalMovieDatabase')
        UseLocalMovieDatabase = true;
        LocalDropboxFolderString = varargin{i+1};
    end
end

%Extract nuclear cycle fit windows, [minutes]
firstnc14time = ncwindow(1,1);
lastnc14time = ncwindow(1,2);


%% Parameters:
MinParticles=1;     %Minimum number of particles in an AP bin % default is 2
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
    if length(varargin)>1
        DropboxFolder = varargin{2};
        existDropboxFolder = 1;
    else
        existDropboxFolder = 0;
    end
else
    FolderTemp=uigetdir(DropboxFolder,'Choose folder with files to analyze');
    Dashes=strfind(FolderTemp,'\');
    Prefix=FolderTemp((Dashes(end)+1):end);
end

%Get the relevant folders now:
if ~existDropboxFolder
    [SourcePath,FISHPath,DropboxFolder,MS2CodePath]=...
        DetermineLocalFolders(Prefix);
end
        


%Load the complied particles and the division information                                    
load([DropboxFolder,filesep,Prefix,'\CompiledParticles.mat'])

if exist([DropboxFolder,filesep,Prefix,'\APDivision.mat'], 'file')
    load([DropboxFolder,filesep,Prefix,'\APDivision.mat'], 'APDivision')
else
    error('Could not load APDivision.mat. Make sure to have done the manual check of division.')
end


% Extract the fields from the cell structure (This is for fields like MeanVectorAP
% that are saved inside {}.
channel = 1;

if iscell(MeanVectorAP)
    MeanVectorAP = MeanVectorAP{channel};
    SDVectorAP = SDVectorAP{channel};
% For instantaneous Fraction ON considered 
% (basically averaged spot fluo over ALL nuclei)
%     try
%         OnRatioAP = OnRatioAP{channel};
%     catch
%         error('No instantaneous fraction on. Check if it has Ellipses.mat')
%     end
    OnRatioAP = ones(size(MeanVectorAP));
%     try 
%         OnRatioAP = OnRatioAP{channel};
%     catch
%         OnRatioAP = ones(size(MeanVectorAP));
%     end
    NParticlesAP = NParticlesAP{channel};
end

%Initial parameters for fits. We will estimate the maximum rate based on
%the elongation time and the maximum average fluorescence of the data set.
MaxRate=max(max(MeanVectorAP))/Delay;

Rate012=MaxRate;     %Rate per minute
TimeStart012=3;
TimeEnd012=7;

Rate013=MaxRate;     %Rate per minute
TimeStart013=5;
TimeEnd013=10;

Rate014=MaxRate;     %Rate per minute
TimeStart014=5;
TimeEnd014=1000;                                    
                           


 
%Rough frame window to consider in the fits

%Some data sets won't have nc12 or even nc13
% if  nc13>0 && nc12>0 % both nc12 and nc13 are captured
%     FrameWindow13=[nc13:nc14];
%     FrameWindow12=[nc12:nc13];
% elseif nc13>0 && nc12 ==0 % only nc13 is captured
%     FrameWindow13=[nc13:nc14];
%     FrameWindow12=[];
% elseif nc13==0 % neither nc13 nor nc12 are captured
%     FrameWindow13=[];
%     FrameWindow12=[];
% end

FrameWindow14=[nc14:length(ElapsedTime)];      


         

%Set the first guess for the parameters for each AP bin and also
%dissaproved the ones that did not have enough data points. Fit results has
% a structure with the fits corresponding to each AP position and nc13
%or nc14
if exist([DropboxFolder,filesep,Prefix,filesep,'MeanFits_LinearRise_ExponentDecay.mat'])
    load([DropboxFolder,filesep,Prefix,filesep,'MeanFits_LinearRise_ExponentDecay.mat']);
    if isempty(FitResults)
        FitResults(length(APbinID),1).Rate0=[];
    end
% elseif exist([DropboxFolder,filesep,Prefix,filesep,'MeanFitsAsymmetric.mat'])
%     load([DropboxFolder,filesep,Prefix,filesep,'MeanFitsAsymmetric.mat']);
%     if isempty(FitResults)
%         FitResults(length(APbinID),3).Rate0=[];
%     end
else
    FitResults(length(APbinID),1).Rate0=[];
end




% Set default starting values for nc14
% nc14
for i=1:length(APbinID)
    if isempty(FitResults(i,1).Rate0)
        FitResults(i,1).Rate0 = Rate014; % 5min
        FitResults(i,1).TimeStart0 = TimeStart014;
        FitResults(i,1).TimePeak0= 10; %[min], 10min
        FitResults(i,1).TimeDown0= 12; %[min], 12min
        FitResults(i,1).Tau0= 30; % [min], half life of the decay curve, 20min
        FitResults(i,1).Fluo_basal = 50; % [AU], basal fluo
        FitResults(i,1).FrameFilter = [];
        FitResults(i,1).FitFrameRange = [];        
        if sum(NParticlesAP(FrameWindow14,i)>=MinParticles)>=MinTimePoints
            FitResults(i,1).Approved = 0;
        else
            FitResults(i,1).Approved = -1;
        end
    end
end







% Go through each AP bin
% Note that the FitResult is a structure with dimension of (APbins x 1)

FitFigure=figure;
CurrentNC=14;
% find the most anterior APbin which has a non-zero particle number
i=min(find(sum(NParticlesAP))); % i : APbin index
cc=1;

while (cc~=13)
    
    figure(FitFigure)
    clf
    
    % Change the background color depending on approval/disapproval
    % of the fit
    if FitResults(i,1).Approved==-1
        set(gcf,'Color','r')
    elseif FitResults(i,1).Approved==1
        set(gcf,'Color','g')
    else
        set(gcf,'Color','default')
    end
    
    
    if APDivision(CurrentNC,i)
        
        FrameWindow=APDivision(CurrentNC,i):length(ElapsedTime);


        %Check that we have the minimum number of particles for a minimum
        %amount of time
        if (sum(NParticlesAP(FrameWindow,i)>=MinParticles)>=MinTimePoints)

            %Extract the data for this range of frames
            FluoData=MeanVectorAP(FrameWindow,i);
            SDFluoData=SDVectorAP(FrameWindow,i);
            NData=NParticlesAP(FrameWindow,i);
            TimeData=ElapsedTime(FrameWindow);
            OnRatioData=OnRatioAP(FrameWindow,i);

            %Now filter them according the number of particles
            FrameFilter=NData>=MinParticles;


            %As an initial guess, use FrameFilter to determine the range of the
            %fit
            if isempty(FitResults(i,1).FitFrameRange)
                FitFrameRange=FrameWindow(FrameFilter);
                % Limit the frame range for fitting 
                % Default is [3,25] min in nc14, but this should be modular
                % using the optional input keys.
                tFitStart = 3; %[min]
                tFitEnd = 25; % [min]
                frameRange1 = (ElapsedTime(FitFrameRange)-ElapsedTime(APDivision(CurrentNC,i)))>tFitStart;
                frameRange2 = (ElapsedTime(FitFrameRange)-ElapsedTime(APDivision(CurrentNC,i)))<tFitEnd;
                FitFrameRange=FitFrameRange(frameRange1 & frameRange2);
                FitResults(i,1).FitFrameRange=FitFrameRange;
            else
                FitFrameRange=FitResults(i,1).FitFrameRange;
            end

            %Filter the frames according to FitFrameRange
            FitFrameFilter=ismember(FrameWindow,FitFrameRange);
            
            
            %HG: Removing the OnRatioFit stuff here. I'm not too happy
            %about it.
%             OnRatioDataForFit=OnRatioData(FitFrameFilter);
%             MaxOnRatioForFit=max(OnRatioData);
%             OnRatioDataForFit=OnRatioDataForFit/MaxOnRatioForFit;
% 
%             FluoDataForFit=FluoData(FitFrameFilter).*OnRatioDataForFit;
%             SDFluoDataForFit=SDFluoData(FitFrameFilter).*OnRatioDataForFit;
%             NDataForFit=NData(FitFrameFilter);
%             TimeDataForFit=TimeData(FitFrameFilter);

            % There's no OnRaioAP data anymore used for calculating the
            % MeanVectorAP. So, we're only thinking about the ON nuclei
            % behavior.
%             OnRatioDataForFit=OnRatioData(FitFrameFilter);
%             MaxOnRatioForFit=max(OnRatioData);
%             OnRatioDataForFit=OnRatioDataForFit/MaxOnRatioForFit;
            % filter out the Fluo, error, number of particles, Time with
            % the frame filter.
            FluoDataForFit=FluoData(FitFrameFilter);
            SDFluoDataForFit=SDFluoData(FitFrameFilter);
            NDataForFit=NData(FitFrameFilter);
            TimeDataForFit=TimeData(FitFrameFilter);
            

            %These is the maximum range of data for the fit
%             OnRatioData=OnRatioData(FrameFilter);
%             MaxOnRatio=max(OnRatioData);
%             OnRatioData=OnRatioData/MaxOnRatio;

%             FluoData=FluoData(FrameFilter).*OnRatioData;
%             SDFluoData=SDFluoData(FrameFilter).*OnRatioData;
%             NData=NData(FrameFilter);
%             TimeData=TimeData(FrameFilter);

            % Do the same filtering for the raw data (for plotting the raw
            % fluo traces, etc.)
            FluoData=FluoData(FrameFilter);
            SDFluoData=SDFluoData(FrameFilter);
            NData=NData(FrameFilter);
            TimeData=TimeData(FrameFilter);
                
            % Do the fit
            % initial condition : 
%             TimeStart0=x0(1); % turn-on time
%             Rate0=x0(2); % initial rate
%             TimePeak0 = x0(3) % peak time (time at which the fluorescence peaks)
%             TimeDown0 = x0(4); % time point when the fluo starts to go down
%             Tau0 =x0(5); % half-life
%             Fluo_basal = x0(6); % basal fluorescence level (basal expression level given that the fluorescence doesn't go to zero even in late NC14)
            
            x0=[FitResults(i,1).TimeStart0, FitResults(i,1).Rate0,...
                FitResults(i,1).TimePeak0,  FitResults(i,1).TimeDown0,...
                 FitResults(i,1).Tau0,  FitResults(i,1).Fluo_basal];


                
                %Get rid of any NaN in the data
                NanFilter=~isnan(FluoDataForFit);

                if ~isempty(TimeData(NanFilter))

                    [xFit,resnorm,residual,exitflag,output,lambda,jacobian]=...
                        lsqnonlin(@(x) lsqnonlinFitFluorescenceCurveNC14LinearSlopeExpDecay(TimeDataForFit(NanFilter)-...
                        ElapsedTime(FrameWindow(1)),...
                        FluoDataForFit(NanFilter),x),x0);

                    FitResults(i,1).TimeStart=xFit(1);
                    FitResults(i,1).RateFit=xFit(2);
                    FitResults(i,1).TimePeak=xFit(3);
                    FitResults(i,1).TimeDown=xFit(4);
                    FitResults(i,1).Tau=xFit(5);
                    FitResults(i,1).FluoBasal=xFit(6);
%                     FitResults(i,1).RateOffFit=nan; % just for syntax: YJK
%                     FitResults(i,1).SDRateOffFit=nan; % just for syntax: YJK

                    %Estimate an error bar out of the confidence intervals
                    FitResults(i,1).CI=nlparci(xFit,residual,'jacobian',jacobian);

                    FitResults(i,1).SDTimeStart=(FitResults(i,1).CI(1,2)-FitResults(i,1).CI(1,1))/2;
                    FitResults(i,1).SDRateFit=(FitResults(i,1).CI(2,2)-FitResults(i,1).CI(2,1))/2;
                    
                    FitResults(i,1).SDTau=(FitResults(i,1).CI(5,2)-FitResults(i,1).CI(5,1))/2;




                    %Plot the results
                    %Get the corresponding fitted curve
                    [TimeFit,FluoFit]=FluorescenceCurveLinearSlopeExpDecay(ElapsedTime(FrameWindow(end))-...
                        ElapsedTime(FrameWindow(1)), xFit);
                    %Plot all the data
                    PlotHandle=errorbar(ElapsedTime(FrameWindow)-ElapsedTime(FrameWindow(1)),...
                        MeanVectorAP(FrameWindow,i),...
                        SDVectorAP(FrameWindow,i)./...
                        sqrt(NParticlesAP(FrameWindow,i)),'.-k');
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
                        ylim([0,max(MeanVectorAP(FrameWindow,i)+...
                            SDVectorAP(FrameWindow,i)./...
                            sqrt(NParticlesAP(FrameWindow,i)))])
                    catch
                        display('Error in displaying the plot')
                    end

                    legend(PlotHandle,['tON=',num2str(FitResults(i,1).TimeStart),' \pm ',num2str(FitResults(i,1).SDTimeStart)],...
                        ['Rate=',num2str(FitResults(i,1).RateFit),' \pm ',num2str(FitResults(i,1).SDRateFit)],...
                        ['Tau=',num2str(FitResults(i,1).Tau),' \pm ',num2str(FitResults(i,1).SDTau)],...
                        'Location','Best')
                end
            end

        end
    
    title([num2str(APbinID(i)),' AP, tStart0=',num2str(FitResults(i,1).TimeStart0),...
        ', tPeak=',num2str(FitResults(i,1).TimePeak),', Rate=',num2str(FitResults(i,1).Rate0),...
        ', Tau=',num2str(FitResults(i,1).Tau),', nc',num2str(CurrentNC)])
    
    
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
        if FitResults(i,1).Approved==0
            FitResults(i,1).Approved=1;
        elseif FitResults(i,1).Approved==1
            FitResults(i,1).Approved=0;
        end

    
    %Disapprove, disapprove fit
    elseif (ct~=0)&(cc=='w')
        if FitResults(i,1).Approved==0
            FitResults(i,1).Approved=-1;
        elseif FitResults(i,1).Approved==-1
            FitResults(i,1).Approved=0;
        end
  
    
        
    %Move right range of fit
    elseif (ct~=0)&(cc=='k')&(length(FitResults(i,1).FitFrameRange)>2)
        FitResults(i,1).FitFrameRange=FitResults(i,1).FitFrameRange(1:end-1);
    elseif (ct~=0)&(cc=='l')
        if ~isempty(find(~ismember(FrameWindow(FrameFilter),FitResults(i,1).FitFrameRange)))
            FilteredFramesTemp=FrameWindow(FrameFilter);
            FitResults(i,1).FitFrameRange(end+1)=...
                FilteredFramesTemp(min(find(~ismember(FilteredFramesTemp,FitResults(i,1).FitFrameRange))));
        end
    %Move left range of fit
    elseif (ct~=0)&(cc=='j')&(length(FitResults(i,1).FitFrameRange)>2)
        FitResults(i,1).FitFrameRange=FitResults(i,1).FitFrameRange(2:end);
    elseif (ct~=0)&(cc=='h')
        if ~isempty(find(~ismember(FrameWindow(FrameFilter),FitResults(i,1).FitFrameRange)))
            FilteredFramesTemp=FrameWindow(FrameFilter);
            FitResults(i,1).FitFrameRange=...
                [FilteredFramesTemp(max(find(~ismember(FilteredFramesTemp,FitResults(i,1).FitFrameRange)))),...
                FitResults(i,1).FitFrameRange];
        end
    %Reset frame fit range
     elseif (ct~=0)&(cc=='r')   
        FitResults(i,1).FitFrameRange=FrameWindow(FrameFilter);

        
        
    %Change the initial parameters
    %TimeStart
    elseif (ct~=0)&(cc=='a')&((CurrentNC==14)|(CurrentNC~=14&FitResults(i,1).TimeStart0<FitResults(i,1).TimeEnd0))
        FitResults(i,1).TimeStart0=FitResults(i,1).TimeStart0+1;
    elseif (ct~=0)&(cc=='z')&(FitResults(i,1).TimeStart0>1)
        FitResults(i,1).TimeStart0=FitResults(i,1).TimeStart0-1;
    %TimeEnd
%     elseif (ct~=0)&(cc=='s')&(FitResults(i,1).TimeEnd0<ElapsedTime(FrameWindow(end))-ElapsedTime(FrameWindow(1)))
%         FitResults(i,1).TimeEnd0=FitResults(i,1).TimeEnd0+1;
%     elseif (ct~=0)&(cc=='x')&(FitResults(i,1).TimeEnd0>FitResults(i,1).TimeStart0)
%         FitResults(i,1).TimeEnd0=FitResults(i,1).TimeEnd0-1;
    %Rate, fine
    elseif (ct~=0)&(cc=='c')&(FitResults(i,1).Rate0>100)
        FitResults(i,1).Rate0=FitResults(i,1).Rate0-100;
    elseif (ct~=0)&(cc=='d')
        FitResults(i,1).Rate0=FitResults(i,1).Rate0+100;    
    %Rate, coarse
    elseif (ct~=0)&(cc=='C')&(FitResults(i,1).Rate0>100)
        FitResults(i,1).Rate0=FitResults(i,1).Rate0-500;
    elseif (ct~=0)&(cc=='D')
        FitResults(i,1).Rate0=FitResults(i,1).Rate0+500; 
    %RateOff, fine
    elseif (ct~=0)&(cc=='v')&(FitResults(i,1).RateOff0<-100)
        FitResults(i,1).RateOff0=FitResults(i,1).RateOff0-100;
    elseif (ct~=0)&(cc=='f')
        FitResults(i,1).RateOff0=FitResults(i,1).RateOff0+100;    
    %RateOff, coarse
    elseif (ct~=0)&(cc=='V')&(FitResults(i,1).RateOff0<-100)
        FitResults(i,1).RateOff0=FitResults(i,1).RateOff0-500;
    elseif (ct~=0)&(cc=='F')
        FitResults(i,1).RateOff0=FitResults(i,1).RateOff0+500;     
        
    
    %Switch NCs
%     elseif (ct~=0)&(cc=='m')&CurrentNC<14
%         CurrentNC=CurrentNC+1;
%     elseif (ct~=0)&(cc=='n')&CurrentNC>12
%         CurrentNC=CurrentNC-1;
    
        % Save the fitted slope & raw data plot
    elseif (ct~=0)&(cc=='S')
        if ~exist([DropboxFolder,filesep,Prefix,filesep,'AsymmetricFit_snapshots'])
            mkdir([DropboxFolder,filesep,Prefix,filesep,'AsymmetricFit_snapshots'])
        end
        FigPath = [DropboxFolder,filesep,Prefix,filesep,'AsymmetricFit_snapshots'];
        % Save the figures as .tif and .pdf
        StandardFigure(FitFigure,FitFigure.CurrentAxes)
        saveas(FitFigure,[FigPath,filesep, 'LinearSlope_ExpDecay_AP=', num2str(APbinID(i)*100),'%' , '_NC',num2str(CurrentNC) , '.tif']); 
        saveas(FitFigure,[FigPath,filesep, 'LinearSlope_ExpDecay_AP=', num2str(APbinID(i)*100),'%' , '_NC',num2str(CurrentNC) , '.pdf']); 
        display('Plot for the fitted slope and raw data is saved')
        
    %Save
    elseif (ct~=0)&(cc=='e')
        save([DropboxFolder,filesep,Prefix,filesep,'MeanFits_LinearRise_ExponentDecay.mat'],...
        'FitResults')
    display('MeanFits_LinearRise_ExponentDecay.mat saved')
        
    %Debug mode
    elseif (ct~=0)&(cc=='9')
        keyboard
    
    end
    
end


%Save the information
save([DropboxFolder,filesep,Prefix,filesep,'MeanFits_LinearRise_ExponentDecay.mat'],...
    'FitResults')
display('MeanFits_LinearRise_ExponentDecay.mat saved')        
        
close(FitFigure)