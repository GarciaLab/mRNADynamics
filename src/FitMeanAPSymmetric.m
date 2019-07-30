function FitMeanAPSymmetric(varargin)

%This function performs fits to the mean fluorescence as a function of time
%of a particular dataset.

%MT, 2016-11-0: This function now provides the option to fit a single 
%dataset or to combine multiple datasets and fit the resulting dataset.
%Multiple datasets are combined by calling CombineMultipleEmbryos and the
%resulting dataset is fit by this script just like a single dataset.
%OPTIONAL INPUT: varargin{1} = Prefix
%                varargin{2} = 'multiple'
%                               Include string, 'multiple', if you need to 
%                               combine multiple embryos before fitting the
%                               data with mean. If you are fitting a single
%                               embryo dataset, do not include.
%                varargin{3} =  DataType
%                               If you include the 'multiple' parameter,
%                               you need to include the DataType
%                               of the embryo datasets you wish to analyze.
%                               DataType is the name of the tab in the
%                               DataStatus Excel file, which should include
%                               all the information for the datasets you 
%                               want to combine.
%NB: If you do not enter a Prefix, but enter 'multiple' and a DataType, the
%code is still able to handle this type of input.

%2013/08/18: Changed this so it can automatically detect whether we are
%dealing with a 5' or 3' data set
%OUTPUT:MeanFits.mat
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

%Approve/Reject:
% You need to approve or reject fits.
% w: reject; q: approve. Enter saves MeanFits.mat.

%SAVE: v

%Parameters:
MinParticles=1;     %Minimum number of particles in an AP bin
MinTimePoints=5;    %Minimum number of time points where we'll have at least
                    %the minimum number of particles.
ElongationRate=1.54;    %In kb/minutes.
GeneLength5=5.296;      %Distance from the first MS2 site to the end of the
                        %TUB3'UTR in kb for the 5' constrcut.
GeneLength3=1.941;      %Distance from the first MS2 site to the end of the
                        %TUB3'UTR in kb for the 3' constrcut.                        

MultipleEmbryos = 0;    %Keeps track of whether or not multiple embryos 
                        %need to be combined before performing the fits
                                                                   
                                    
close all

%Find out which computer this is. That will determine the folder structure.
%Information about about folders
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;


if ~isempty(varargin)
    %MT, 2016-11-4
    %Determine if multiple embryos need to be combined
    
    %Case where 'multiple' is indicated, but neither Prefix nor DataType is
    %provided
    if length(varargin) == 1
        if strcmp(varargin{1},'multiple')
            FolderTemp=uigetdir(DropboxFolder,'Choose folder with files to analyze');
            Dashes=strfind(FolderTemp,filesep);
            Prefix=FolderTemp((Dashes(end)+1):end);
            
            MultipleEmbryos = 1;
            
            prompt = 'You''ve indicated that you have multiple embryos you need to combine,\nbut have not included a DataType.\nPlease review documentation and enter the appropriate string:';
            userInput = input(prompt, 's');
            if isempty(userInput)
                error('You did not provide a DataType, which is required to combine multiple embryos. Please review documentation and re-run the function with the appropriate input parameters.');
            else
                DataType = strtrim(userInput);
            end
        %Case where Prefix is provided for a single embryo
        else
            Prefix = varargin{1};
        end
    elseif length(varargin) == 2
        %Case where Prefix is not provided, but 'multiple' and DataType are
        if strcmp(varargin{1},'multiple')
            FolderTemp=uigetdir(DropboxFolder,'Choose folder with files to analyze');
            Dashes=strfind(FolderTemp,filesep);
            Prefix=FolderTemp((Dashes(end)+1):end);
            MultipleEmbryos = 1;
        %Case where Prefix is provided and 'multiple' is indicated, but no
        %DataType is provided
        elseif strcmp(varargin{2},'multiple')
            Prefix = varargin{1};
            MultipleEmbryos = 1;
            
            prompt = 'You''ve indicated that you have multiple embryos you need to combine, \nbut have not included a DataType. \nPlease review documentation and enter the appropriate string:';
            userInput = input(prompt, 's');
            if isempty(userInput)
                error('You did not provide a DataType, which is required to combine multiple embryos. Please review documentation and re-run the function with the appropriate input parameters.');
            else
                DataType = strtrim(userInput);
            end
        end
    elseif length(varargin) == 3
        %Case where Prefix is provided, 'multiple' is indicated, and 
        %DataType is provided
        if strcmp(varargin{2},'multiple')
            Prefix = varargin{1};
            MultipleEmbryos = 1;
            DataType = varargin{3};
        end
    end          
else
    FolderTemp=uigetdir(DropboxFolder,'Choose folder with files to analyze');
    Dashes=strfind(FolderTemp,filesep);
    Prefix=FolderTemp((Dashes(end)+1):end);
end


[~,~,DropboxFolder,~,~]=...
    DetermineLocalFolders(Prefix);
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

%Approve/Reject:
% You need to approve or reject fits.
% w: reject; q: approve. Enter saves MeanFits.mat.

%SAVE: v

buttonFig = figure('units', 'normalized', 'Position',[.3 .3 .3 .3]);
pb = uicontrol(buttonFig,'Style','text','String',['Controls-',...
'a,z: On time',...
's,x: Off time',...
'd,c: Rate, fine',...
'D,C: Rate, coarse',...
...
'Moving around:',...
', .: Move in AP',...
'n,m: Move in nc',...
'k,l: Change fit range from the right',...
'h,j: Change fit range from the left',...
...
'Approve/Reject:',...
'You need to approve or reject fits',...
'w: reject; q: approve. Enter saves MeanFits.mat.'],...
'units', 'normalized', 'Position',[0 0 1 1]);

%Run CombineMultipleEmbryos if required
if MultipleEmbryos
    saveFolder = CombineMultipleEmbryos(DataType);
    disp('Embryos successfully combined')
end

                                    
%Load the complied particles and the division information
if MultipleEmbryos
    load([saveFolder,filesep,DataType,'_Combined_CompiledParticles.mat'])

    if exist([saveFolder,filesep,DataType,'_Combined_APDivision.mat'], 'file')
        load([saveFolder,filesep,DataType,'_Combined_APDivision.mat'])
    else
        error('Could not load Combined_APDivision.mat. Make sure to have done the manual check of division.')
    end
else
    load([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat'])
    
    if exist([DropboxFolder,filesep,Prefix,filesep,'APDivision.mat'], 'file')
        load([DropboxFolder,filesep,Prefix,filesep,'APDivision.mat'], 'APDivision')
    else
        error('Could not load APDivision.mat. Make sure to have done the manual check of division.')
    end
end
                                  
 
%Rough frame window to consider in the fits

%Some data sets won't have nc12
if nc12>0
    FrameWindow12=nc12:nc13;
else
    FrameWindow12=[];
end

%Some data sets won't have nc13 either
if nc13>0
    FrameWindow13=nc13:nc14;
else
    FrameWindow13=[];
end



FrameWindow14=nc14:length(ElapsedTime);      

prime = 5;

if (~isempty(findstr(Prefix,'X1')))|(~isempty(findstr(Prefix,'P2P')))|...
        (~isempty(findstr(Prefix,'evePr')))...
        | (exist('StemLoopEnd', 'var') && strcmp(StemLoopEnd, '5'''))
    prime = 5;
elseif ~isempty(findstr(Prefix,'X2')) | (exist('StemLoopEnd', 'var') && strcmp(StemLoopEnd, '3'''))
    prime = 3;
end

if prime == 5
    disp('Treating data set as 5''')
    Delay=GeneLength5/ElongationRate;
elseif prime == 3
    Delay=GeneLength3/ElongationRate;
end

MaxRate=max(max(MeanVectorAP))/Delay;

%Initial parameters for fits

Rate012=MaxRate;     %Rate per minute
TimeStart012=3;
TimeEnd012=7;

Rate013=MaxRate;     %Rate per minute
TimeStart013=5; 
TimeEnd013=10;

Rate014=MaxRate;     %Rate per minute
TimeStart014=5;
TimeEnd014=1000;  

if prime == 3
    TimeStart013 = TimeStart013 + 2.5;
    TimeStart014 = TimeStart013 + 2.5;
end

%Set the first guess for the parameters for each AP bin and also
%disapproved the ones that did not have enough data points. Fit results has
%is a structure with the fits corresponding to each AP position and nc13
%or nc14
if MultipleEmbryos
    if exist([DropboxFolder,filesep,DataType,'_Combined_MeanFits.mat'], 'file')
        load([DropboxFolder,filesep,DataType,'_Combined_MeanFits.mat']);
        if isempty(FitResults)
            FitResults(length(APbinID),3).Rate0=[];
        end
    else
        FitResults(length(APbinID),3).Rate0=[];
    end
else
    if exist([DropboxFolder,filesep,Prefix,filesep,'MeanFits.mat'], 'file')
        load([DropboxFolder,filesep,Prefix,filesep,'MeanFits.mat']);
        if isempty(FitResults)
            FitResults(length(APbinID),3).Rate0=[];
        end
    else
        FitResults(length(APbinID),3).Rate0=[];
    end
end






%Set default starting values for nc 13 and nc14
%nc12
for APBin=1:length(APbinID)
    if isempty(FitResults(APBin,1).Rate0)
        FitResults(APBin,1).Rate0=Rate012;    
        FitResults(APBin,1).TimeStart0=TimeStart012;
        FitResults(APBin,1).TimeEnd0=TimeEnd012;
        FitResults(APBin,1).FrameFilter=[];
        FitResults(APBin,1).FitFrameRange=[];
        if sum(NParticlesAP(FrameWindow12,APBin)>=MinParticles)>=MinTimePoints
            FitResults(APBin,1).Approved=0;
        else
            FitResults(APBin,1).Approved=-1;
        end
    end
end

%nc13


for APBin=1:length(APbinID)
    if isempty(FitResults(APBin,2).Rate0)
        FitResults(APBin,2).Rate0=Rate013;    
        FitResults(APBin,2).TimeStart0=TimeStart013;
        FitResults(APBin,2).TimeEnd0=TimeEnd013;
        FitResults(APBin,2).FrameFilter=[];
        FitResults(APBin,2).FitFrameRange=[];
        
        if sum(NParticlesAP(FrameWindow13,APBin)>=MinParticles)>=MinTimePoints
            FitResults(APBin,2).Approved=0;
        else
            FitResults(APBin,2).Approved=-1;
        end

    end
end

%nc14

for APBin=1:length(APbinID)
    if isempty(FitResults(APBin,3).Rate0)
        FitResults(APBin,3).Rate0=Rate014;    
        FitResults(APBin,3).TimeStart0=TimeStart014;
        FitResults(APBin,3).TimeEnd0=[];
        FitResults(APBin,3).FrameFilter=[];
        FitResults(APBin,3).FitFrameRange=[];        
        if sum(NParticlesAP(FrameWindow14,APBin)>=MinParticles)>=MinTimePoints
            FitResults(APBin,3).Approved=0;
        else
            FitResults(APBin,3).Approved=-1;
        end
    end
end





%Figure out which nc we can use
if nc12
    CurrentNC=12;
    MinNC=12;       %We'll use this to keep track of the minimum nc
elseif nc13
    CurrentNC=13;
    MinNC=13;
elseif nc14
    CurrentNC=14;
    MinNC=14;
else
    error('There is a problem. Are the ncs defined?')
end


%Go through each AP bin
FitFigure=figure;
APBin=find(sum(NParticlesAP), 1 ); %index of the first AP bin that has a non-zero number of particles
cc=1;

 lsqOptions=optimset('Display','none');

minRate = 0;
maxRate = Inf;

while (cc~=13)
    
    figure(FitFigure)
    clf
    
    if FitResults(APBin,CurrentNC-11).Approved==-1
        set(gcf,'Color','r')
    elseif FitResults(APBin,CurrentNC-11).Approved==1
        set(gcf,'Color','g')
    else
        set(gcf,'Color','default')
    end
    
    
    if APDivision(CurrentNC,APBin)
        if CurrentNC~=14
            FrameWindow=APDivision(CurrentNC,APBin):APDivision(CurrentNC+1,APBin);
        else
            FrameWindow=APDivision(CurrentNC,APBin):length(ElapsedTime);
        end
        
        if isempty(FrameWindow)
            error('There might be a problem with the division times. Check in CheckDivisionTimes.m')
        end

        %Check that we have the minimum number of particles for a minimum
        %amount of time
        if (sum(NParticlesAP(FrameWindow,APBin)>=MinParticles)>=MinTimePoints)

            %Extract the data for this range of frames
            FluoData=MeanVectorAP(FrameWindow,APBin);
            SDFluoData=SDVectorAP(FrameWindow,APBin);
            NData=NParticlesAP(FrameWindow,APBin);
            TimeData=ElapsedTime(FrameWindow);
          
            %Now filter them according the number of particles
            FrameFilter=NData>=MinParticles;


            %As an initial guess, use FrameFilter to determine the range of the
            %fit
            if isempty(FitResults(APBin,CurrentNC-11).FitFrameRange)
                FitFrameRange=FrameWindow(FrameFilter);
                if CurrentNC==14
                    FitFrameRange=FitFrameRange((ElapsedTime(FitFrameRange)-ElapsedTime(APDivision(CurrentNC,APBin)))<12);
                end
                FitResults(APBin,CurrentNC-11).FitFrameRange=FitFrameRange;
            else
                FitFrameRange=FitResults(APBin,CurrentNC-11).FitFrameRange;
            end

            %Filter the frames according to FitFrameRange
            FitFrameFilter=ismember(FrameWindow,FitFrameRange);
            
           

            FluoDataForFit=FluoData(FitFrameFilter);
            SDFluoDataForFit=SDFluoData(FitFrameFilter);
            NDataForFit=NData(FitFrameFilter);
            TimeDataForFit=TimeData(FitFrameFilter);
            


            FluoData=FluoData(FrameFilter);
            SDFluoData=SDFluoData(FrameFilter);
            NData=NData(FrameFilter);
            TimeData=TimeData(FrameFilter);

            if CurrentNC~=14
                %Do the fit
                x0=[FitResults(APBin,CurrentNC-11).TimeStart0,...
                    FitResults(APBin,CurrentNC-11).TimeEnd0,...
                    FitResults(APBin,CurrentNC-11).Rate0];

                
                %Get rid of any NaN in the data
                NanFilter=~isnan(FluoDataForFit);

                if ~isempty(TimeDataForFit(NanFilter))
                    
                 
                   if CurrentNC == 12
                        lb = [3, 4, minRate];
                        ub = [12 12 maxRate];
                   elseif CurrentNC == 13
                        lb = [3, 4, maxRate];
                        ub = [20 20 minRate];
                   end

                    [xFit,resnorm,residual,exitflag,output,lambda,jacobian]=...
                        lsqnonlin(@(x) lsqnonlinFitFluorescenceCurve(TimeDataForFit(NanFilter)-...
                        ElapsedTime(FrameWindow(1)),...
                        FluoDataForFit(NanFilter),Delay,...
                        ElapsedTime(FrameWindow(end))-ElapsedTime(FrameWindow(1)),x),x0, lb,ub, lsqOptions);

                    FitResults(APBin,CurrentNC-11).TimeStart=xFit(1);
                    FitResults(APBin,CurrentNC-11).TimeEnd=xFit(2);
                    FitResults(APBin,CurrentNC-11).RateFit=xFit(3);

                    %Estimate an error bar out of the confidence intervals
                    FitResults(APBin,CurrentNC-11).CI=nlparci(xFit,residual,'jacobian',jacobian);

                    FitResults(APBin,CurrentNC-11).SDTimeStart=(FitResults(APBin,CurrentNC-11).CI(1,2)-FitResults(APBin,CurrentNC-11).CI(1,1))/2;
                    FitResults(APBin,CurrentNC-11).SDTimeEnd=(FitResults(APBin,CurrentNC-11).CI(2,2)-FitResults(APBin,CurrentNC-11).CI(2,1))/2;
                    FitResults(APBin,CurrentNC-11).SDRateFit=(FitResults(APBin,CurrentNC-11).CI(3,2)-FitResults(APBin,CurrentNC-11).CI(3,1))/2;




                    %Plot the results
                    %Get the corresponding fitted curve
                    [TimeFit,FluoFit]=FluorescenceCurve(ElapsedTime(FrameWindow(end))-...
                        ElapsedTime(FrameWindow(1)),...
                        xFit(1),xFit(2),xFit(3),Delay);
                    %Plot all the data
                    PlotHandle=errorbar(ElapsedTime(FrameWindow)-ElapsedTime(FrameWindow(1)),...
                        MeanVectorAP(FrameWindow,APBin),...
                        SDVectorAP(FrameWindow,APBin)./...
                        sqrt(NParticlesAP(FrameWindow,APBin)),'.-k');
                    hold on
                    %Plot the data that could be used for the fit
                    PlotHandle(end+1)=plot(ElapsedTime(FrameWindow(FrameFilter))-ElapsedTime(FrameWindow(1)),...
                        FluoData,'or');
                    %Plot the data that was actually used for the fit
                    PlotHandle(end+1)=plot(ElapsedTime(FitFrameRange)-ElapsedTime(FrameWindow(1)),...
                        FluoDataForFit(ismember(FrameWindow(FrameFilter),FitFrameRange)),'or','MarkerFaceColor','r');
                    
                    %Plot the fit
                    PlotHandle(end+1)=plot(TimeFit,FluoFit,'-r');
                    hold off
                    %StandardFigure(PlotHandle,gca)
                    ylabel('Mean fluorescence nucleus')
                    xlabel('Time into nc (min)')
                    
                    try
                        ylim([0,max(MeanVectorAP(FrameWindow,APBin)+...
                            SDVectorAP(FrameWindow,APBin)./...
                            sqrt(NParticlesAP(FrameWindow,APBin)))]);
                    catch
                        disp('Error in displaying the plot')
                    end

                    legend(['tON=',num2str(FitResults(APBin,CurrentNC-11).TimeStart),' \pm ',num2str(FitResults(APBin,CurrentNC-11).SDTimeStart)],...
                        ['tOFF=',num2str(FitResults(APBin,CurrentNC-11).TimeEnd),' \pm ',num2str(FitResults(APBin,CurrentNC-11).SDTimeEnd)],...
                        ['Rate=',num2str(FitResults(APBin,CurrentNC-11).RateFit),' \pm ',num2str(FitResults(APBin,CurrentNC-11).SDRateFit)],...
                        'Location','SouthOutside')
                end
            elseif CurrentNC==14
                
                %Do the fit
                x0=[FitResults(APBin,CurrentNC-11).TimeStart0,FitResults(APBin,CurrentNC-11).Rate0];


                
                %Get rid of any NaN in the data
                NanFilter=~isnan(FluoDataForFit);

                if ~isempty(TimeData(NanFilter))

                    [xFit,resnorm,residual,exitflag,output,lambda,jacobian]=...
                        lsqnonlin(@(x) lsqnonlinFitFluorescenceCurveNC14(TimeDataForFit(NanFilter)-...
                        ElapsedTime(FrameWindow(1)),...
                        FluoDataForFit(NanFilter),Delay,...
                        ElapsedTime(FrameWindow(end))-ElapsedTime(FrameWindow(1)),x),x0);

                    FitResults(APBin,CurrentNC-11).TimeStart=xFit(1);
                    FitResults(APBin,CurrentNC-11).RateFit=xFit(2);

                    %Estimate an error bar out of the confidence intervals
                    FitResults(APBin,CurrentNC-11).CI=nlparci(xFit,residual,'jacobian',jacobian);

                    FitResults(APBin,CurrentNC-11).SDTimeStart=(FitResults(APBin,CurrentNC-11).CI(1,2)-FitResults(APBin,CurrentNC-11).CI(1,1))/2;
                    FitResults(APBin,CurrentNC-11).SDRateFit=(FitResults(APBin,CurrentNC-11).CI(2,2)-FitResults(APBin,CurrentNC-11).CI(2,1))/2;




                    %Plot the results
                    %Get the corresponding fitted curve
                    [TimeFit,FluoFit]=FluorescenceCurve(ElapsedTime(FrameWindow(end))-...
                        ElapsedTime(FrameWindow(1)),...
                        xFit(1),1000,xFit(2),Delay);
                    %Plot all the data
                    PlotHandle=errorbar(ElapsedTime(FrameWindow)-ElapsedTime(FrameWindow(1)),...
                        MeanVectorAP(FrameWindow,APBin),...
                        SDVectorAP(FrameWindow,APBin)./...
                        sqrt(NParticlesAP(FrameWindow,APBin)),'.-k');
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
                        ylim([0,max(MeanVectorAP(FrameWindow,APBin)+...
                            SDVectorAP(FrameWindow,APBin)./...
                            sqrt(NParticlesAP(FrameWindow,APBin)))])
                    catch
                        disp('Error in displaying the plot')
                    end

                    legend(['tON=',num2str(FitResults(APBin,CurrentNC-11).TimeStart),' \pm ',num2str(FitResults(APBin,CurrentNC-11).SDTimeStart)],...
                        ['Rate=',num2str(FitResults(APBin,CurrentNC-11).RateFit),' \pm ',num2str(FitResults(APBin,CurrentNC-11).SDRateFit)],...
                        'Location','SouthOutside')
                end
            end

        end
    end
    
    title([num2str(APbinID(APBin)),' AP, TimeStart0=',num2str(FitResults(APBin,CurrentNC-11).TimeStart0),...
        ', TimeEnd0=',num2str(FitResults(APBin,CurrentNC-11).TimeEnd0),', Rate=',num2str(FitResults(APBin,CurrentNC-11).Rate0),...
        ', nc',num2str(CurrentNC)])
    
    
    %Set the limits on the x-axis
    if CurrentNC==14
        xlim([0,60]) %min
    end
    
    
    
    figure(FitFigure)
    try
        ct=waitforbuttonpress;
    catch
        error('Fits not saved.');
    end
    cc=get(FitFigure,'currentcharacter');
    cm=get(gca,'CurrentPoint');
    
    %Move between AP positions
    if (ct~=0)&(cc=='.')&(APBin<length(APbinID))
        APBin=APBin+1;
    elseif (ct~=0)&(cc==',')&(APBin>1)
        APBin=APBin-1;
    
    %Approve, disapprove fit
    elseif (ct~=0)&(cc=='q')
        if FitResults(APBin,CurrentNC-11).Approved~=1
            FitResults(APBin,CurrentNC-11).Approved=1;
        else
            FitResults(APBin,CurrentNC-11).Approved=0;
        end

    
    %Disapprove, disapprove fit
    elseif (ct~=0)&(cc=='w')
        if FitResults(APBin,CurrentNC-11).Approved~=-1
            FitResults(APBin,CurrentNC-11).Approved=-1;
        else
            FitResults(APBin,CurrentNC-11).Approved=0;
        end
  
    
        
    %Move right range of fit
    elseif (ct~=0)&(cc=='k')&(length(FitResults(APBin,CurrentNC-11).FitFrameRange)>2)
        FitResults(APBin,CurrentNC-11).FitFrameRange=FitResults(APBin,CurrentNC-11).FitFrameRange(1:end-1);
    elseif (ct~=0)&(cc=='K')&(length(FitResults(APBin,CurrentNC-11).FitFrameRange)>6)
        FitResults(APBin,CurrentNC-11).FitFrameRange=FitResults(APBin,CurrentNC-11).FitFrameRange(1:end-5);
    elseif (ct~=0)&(cc=='l')
        if ~isempty(find(~ismember(FrameWindow(FrameFilter),FitResults(APBin,CurrentNC-11).FitFrameRange), 1))
            FilteredFramesTemp=FrameWindow(FrameFilter);
            %HG added
            FilteredFramesTemp=FilteredFramesTemp(~ismember(FilteredFramesTemp,FitResults(APBin,CurrentNC-11).FitFrameRange));
            FitResults(APBin,CurrentNC-11).FitFrameRange(end+1)=...
                min(FilteredFramesTemp(FilteredFramesTemp>max(FitResults(APBin,CurrentNC-11).FitFrameRange)));
        end
    %Move left range of fit
    elseif (ct~=0)&(cc=='j')&(length(FitResults(APBin,CurrentNC-11).FitFrameRange)>2)
        FitResults(APBin,CurrentNC-11).FitFrameRange=FitResults(APBin,CurrentNC-11).FitFrameRange(2:end);
    elseif (ct~=0)&(cc=='J')&(length(FitResults(APBin,CurrentNC-11).FitFrameRange)>6)
        FitResults(APBin,CurrentNC-11).FitFrameRange=FitResults(APBin,CurrentNC-11).FitFrameRange(6:end);
    elseif (ct~=0)&(cc=='h')
        if ~isempty(find(~ismember(FrameWindow(FrameFilter),FitResults(APBin,CurrentNC-11).FitFrameRange)))
            FilteredFramesTemp=FrameWindow(FrameFilter);
            %Modified HG
            FitResults(APBin,CurrentNC-11).FitFrameRange=...
                [max(FilteredFramesTemp(FilteredFramesTemp<min(FitResults(APBin,CurrentNC-11).FitFrameRange))),...
                FitResults(APBin,CurrentNC-11).FitFrameRange];
        end
    %Reset frame fit range
     elseif (ct~=0)&(cc=='r')   
        FitResults(APBin,CurrentNC-11).FitFrameRange=FrameWindow(FrameFilter);

        
        
    %Change the initial parameters
    %TimeStart
    elseif (ct~=0)&(cc=='a')&((CurrentNC==14)|(CurrentNC~=14&FitResults(APBin,CurrentNC-11).TimeStart0<FitResults(APBin,CurrentNC-11).TimeEnd0))
        FitResults(APBin,CurrentNC-11).TimeStart0=FitResults(APBin,CurrentNC-11).TimeStart0+1;
    elseif (ct~=0)&(cc=='z')&(FitResults(APBin,CurrentNC-11).TimeStart0>1)
        FitResults(APBin,CurrentNC-11).TimeStart0=FitResults(APBin,CurrentNC-11).TimeStart0-1;
    %TimeEnd
    elseif (ct~=0)&(cc=='s')&(FitResults(APBin,CurrentNC-11).TimeEnd0<ElapsedTime(FrameWindow(end))-ElapsedTime(FrameWindow(1)))
        FitResults(APBin,CurrentNC-11).TimeEnd0=FitResults(APBin,CurrentNC-11).TimeEnd0+1;
    elseif (ct~=0)&(cc=='x')&(FitResults(APBin,CurrentNC-11).TimeEnd0>FitResults(APBin,CurrentNC-11).TimeStart0)
        FitResults(APBin,CurrentNC-11).TimeEnd0=FitResults(APBin,CurrentNC-11).TimeEnd0-1;
    %Rate, fine
    elseif (ct~=0)&(cc=='c')&(FitResults(APBin,CurrentNC-11).Rate0>100)
        FitResults(APBin,CurrentNC-11).Rate0=FitResults(APBin,CurrentNC-11).Rate0-10;
    elseif (ct~=0)&(cc=='d')
        FitResults(APBin,CurrentNC-11).Rate0=FitResults(APBin,CurrentNC-11).Rate0+10;    
    %Rate, coarse
    elseif (ct~=0)&(cc=='C')&(FitResults(APBin,CurrentNC-11).Rate0>100)
        FitResults(APBin,CurrentNC-11).Rate0=FitResults(APBin,CurrentNC-11).Rate0-100;
    elseif (ct~=0)&(cc=='D')
        FitResults(APBin,CurrentNC-11).Rate0=FitResults(APBin,CurrentNC-11).Rate0+100; 
    
    %Switch NCs
    elseif (ct~=0)&(cc=='m')&CurrentNC<14
        CurrentNC=CurrentNC+1;
    elseif (ct~=0)&(cc=='n')&CurrentNC>MinNC
        CurrentNC=CurrentNC-1;
        
    %Set lower and upper bounds for fitting the initiation rate
    elseif (ct~=0)&(cc=='g')
        
        inRates = inputdlg({'Select min rate:', ...
                'Select max rate'});
        minRate = str2double(inRates{1});
        maxRate = str2double(inRates{2}); 
        
    %Save
    elseif ct~=0 & cc=='v'
        if MultipleEmbryos
            save([saveFolder,filesep,DataType,'_Combined_MeanFits.mat'],...
                'FitResults')
            disp('DataType_Combined_MeanFits.mat saved')
        else
            save([DropboxFolder,filesep,Prefix,filesep,'MeanFits.mat'],...
                'FitResults')
            disp('MeanFits.mat saved')
        end
        
    %Debug mode
    elseif (ct~=0)&(cc=='9')
        keyboard
    
    end
    
end


%Save the information
if MultipleEmbryos
    save([saveFolder,filesep,DataType,'_Combined_MeanFits.mat'],...
        'FitResults')
    disp('DataType_Combined_MeanFits.mat saved')
else
    save([DropboxFolder,filesep,Prefix,filesep,'MeanFits.mat'],...
        'FitResults')
    disp('MeanFits.mat saved')        
end

close(pb);
close(FitFigure);