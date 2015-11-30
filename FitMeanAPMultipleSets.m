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
%<,>: Move between data sets
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

%AR 5/31/15: THIS WILL BE FOR EXTRACTING PREFIX FROM DATASTATUS.XLSX. (TO
%DO)
% PausingXLSName='DataStatus.xlsx';
% [StatusNum,StatusTxt]=xlsread([DropboxFolder,filesep,PausingXLSName],DataType);
% CompileRow=find(strcmp(StatusTxt(:,1),'AnalyzeLiveData Compile Particles'));
% CompiledSets=find(strcmp(StatusTxt(CompileRow,:),'READY')|strcmp(StatusTxt(CompileRow,:),'ApproveAll'));
% for i=1:length(CompiledSets)
%     SetName=StatusTxt{6,CompiledSets(i)};
%     Quotes=strfind(SetName,'''');
%     Prefix=SetName((Quotes(1)+1):(Quotes(end)-1));

FolderTemp=uigetdir(DropboxFolder,'Choose folder with files to analyze. If multiple sets, pick one.');
Dashes=strfind(FolderTemp,'\');
Prefix=FolderTemp((Dashes(end)+1):end);
Dashes2 = strfind(Prefix, '-');
[XLSNum,XLSTxt,XLSRaw]=xlsread([DropboxFolder,filesep,'MovieDatabase.xlsx']);
DataFolderColumn=find(strcmp(XLSTxt(1,:),'DataFolder'));
PrefixRow = find(strcmp(XLSTxt(:, DataFolderColumn),...
    [Prefix(1:Dashes2(3)-1), '\', Prefix(Dashes2(3)+1:end)]));
APResolutionColumn = find(strcmp(XLSRaw(1,:),'APResolution'));
APResolution = XLSRaw{PrefixRow,APResolutionColumn};



                                    
%Load the complied particles and the division information                                    
load([DropboxFolder,filesep,Prefix,'\CompiledParticles.mat'])

data = LoadMS2Sets(DataType);
CP = {};
for i=1:length(data)
    CP{i} = data(i).CompiledParticles;
end
NParticles = 0;
for i = 1:length(CP)
    NParticles = NParticles + length(CP{i});
end
Fluos = {};
for j = 1:length(CP)
    for i = 1:length(CP{j})
        Fluos{j,i} = CP{j}(i).Fluo;
    end
end
FluoErrors = {};
for j = 1:length(CP)
    for i = 1:length(CP{j})
        FluoErrors{j,i} = CP{j}(i).FluoError;
    end
end
Frames = {};
for j = 1:length(CP)
    for i = 1:length(CP{j})
        Frames{j,i} = CP{j}(i).Frame;
    end
end
MeanAPs = {};
for j = 1:length(CP)
    for i = 1:length(CP{j})
        MeanAPs{j,i} = CP{j}(i).MeanAP;
    end
end

%Rough frame window to consider in the fits

%Some data sets won't have nc12
for i = 1:length(data)
    
    if data(i).nc12>0
        FrameWindow12{i}=[data(i).nc12:data(i).nc13];
    else
        FrameWindow12{i}=[];
    end
    %Some data sets won't have nc13 either
    if data(i).nc13>0
        FrameWindow13{i}=[data(i).nc13:data(i).nc14];
    else
        FrameWindow13{i}=[];
    end
    
    FrameWindow14{i}=[data(i).nc14:length(data(i).ElapsedTime)]; 

end


%Detect what type of data set we're dealing with so we can set the delay
if (~isempty(findstr(Prefix,'X1')))|(~isempty(findstr(Prefix,'P2P')))|...
        (~isempty(findstr(Prefix,'evePr')))
    Delay=GeneLength5/ElongationRate;    %Minutes for PolII to fall off after reaching
                                        %the first MS2 site.
    display('Treating data set as 5''')
    
    %Initial parameters for fits
    Rate012=400;     %Rate per minute
    TimeStart012=3;
    TimeEnd012=7;

    Rate013=400;     %Rate per minute
    TimeStart013=5;
    TimeEnd013=10;

    Rate014=400;     %Rate per minute
    TimeStart014=5;
    TimeEnd014=1000;  
elseif ~isempty(findstr(Prefix,'X2'))
    Delay=GeneLength3/ElongationRate;
    display('Treating data set as 3''')
    
    Rate012=400;     %Rate per minute
    TimeStart012=3;
    TimeEnd012=7;

    Rate013=400;     %Rate per minute
    TimeStart013=7.5;
    TimeEnd013=10;

    Rate014=400;     %Rate per minute
    TimeStart014=7.5;
    TimeEnd014=1000; 
else
    % ES 2013-10-14: I don't use X1 or X2 prefixes
    if exist('StemLoopEnd', 'var') && strcmp(StemLoopEnd, '5''')
        Delay=GeneLength5/ElongationRate;    %Minutes for PolII to fall off after reaching
        %the first MS2 site.
        display('Treating data set as 5''')
        
        Rate012=400;     %Rate per minute
        TimeStart012=3;
        TimeEnd012=7;
        
        Rate013=400;     %Rate per minute
        TimeStart013=5;
        TimeEnd013=10;
        
        Rate014=400;     %Rate per minute
        TimeStart014=5;
        TimeEnd014=1000;
    elseif exist('StemLoopEnd', 'var') && strcmp(StemLoopEnd, '3''')
        Delay=GeneLength3/ElongationRate;
        display('Treating data set as 3''')
        
        Rate012=400;     %Rate per minute
        TimeStart012=3;
        TimeEnd012=7;
        
        Rate013=400;     %Rate per minute
        TimeStart013=7.5;
        TimeEnd013=10;
        
        Rate014=400;     %Rate per minute
        TimeStart014=7.5;
        TimeEnd014=1000;
    else
        error('Could not recognize data type from the Prefix or from the value of StemLoopEnd in MovieDatabase.')
    end
end





%Set the first guess for the parameters for each AP bin and also
%disapproved the ones that did not have enough data points. Fit results
%is a structure with the fits corresponding to each AP position and nc13
%or nc14

if exist([DropboxFolder,filesep,Prefix,filesep,'FitResultsMultiple.mat'])
    load([DropboxFolder,filesep,Prefix,filesep,'FitResultsMultiple.mat']);
    if isempty(FitResultsMultiple)
        for i=1:length(FitResultsMultiple)
            FitResultsMultiple{i}.Rate0=[];
        end
    end
else
    FitResultsMultiple = cell(1,length(data));
    for i=1:length(FitResultsMultiple)
        FitResultsMultiple{i}.Rate0=[];
    end
    %Set default starting values for nc12, nc 13 and nc14
    for k = 1:3 %iterate over nuclear cycles
        for j = 1:length(data) %iterate over data sets
            for i=1:length(CP{j}) %iterate over particles
                    FitResultsMultiple{j}(i,k).Rate0=Rate012;    
                    FitResultsMultiple{j}(i,k).TimeStart0=TimeStart012;
                    FitResultsMultiple{j}(i,k).TimeEnd0=TimeEnd012;
                    FitResultsMultiple{j}(i,k).FrameFilter=[];
                    FitResultsMultiple{j}(i,k).FitFrameRange=[];
                    FitResultsMultiple{j}(i,k).Approved=0;
                    FitResultsMultiple{j}(i,k).APBin=0;
                    FitResultsMultiple{j}(i,k).TimeStart=0;
                    FitResultsMultiple{j}(i,k).TimeEnd=0;
                    FitResultsMultiple{j}(i,k).RateFit=0;
                    FitResultsMultiple{j}(i,k).CI=0;
                    FitResultsMultiple{j}(i,k).SDTimeStart=0;
                    FitResultsMultiple{j}(i,k).SDTimeEnd=0;
                    FitResultsMultiple{j}(i,k).SDRateFit=0;
                    for m=1:length(APbinID)
                        if abs(APbinID(m) - CP{j}(i).MeanAP) < APResolution
                            FitResultsMultiple{j}(i,k).APBin = m;
                        end
                    end
            end
        end
    end
end

%Figure out which nc we can use
for i=1:length(data)
    if data(i).nc14~=0
        CurrentNC=14;
        MinNC=14;       %We'll use this to keep track of the minimum nc
    end
end
for i=1:length(data)
    if data(i).nc13~=0
        CurrentNC=13;
        MinNC=13;
    end
end
for i=1:length(data)
    if data(i).nc12~=0
        CurrentNC=12;
        MinNC=12;
    end
end
if exist('CurrentNC','var')==0
    error('There is a problem. Are the ncs defined?')
end



%%
%Iterate over particles

FitFigure=figure;
i = 1; %iterates over particles
j = 1; %iterates over datasets
currentAPBin=1;
cc=1;

while (cc~=13)
    
    for m=1:length(APbinID)
        if abs(APbinID(m) - CP{j}(i).MeanAP) < APResolution
            currentAPBin = m;
        end
    end
    figure(FitFigure)
    clf
    
    if FitResultsMultiple{j}(i,CurrentNC-11).Approved==-1
        set(gcf,'Color','r')
    elseif FitResultsMultiple{j}(i,CurrentNC-11).Approved==1
        set(gcf,'Color','g')
    else
        set(gcf,'Color','default')
    end
    
    
    if data(j).APDivision(CurrentNC,currentAPBin)
        if CurrentNC~=14
            FrameWindow=data(j).APDivision(CurrentNC,currentAPBin):data(j).APDivision(CurrentNC+1,currentAPBin);
        else
            FrameWindow=data(j).APDivision(CurrentNC,currentAPBin):length(data(j).ElapsedTime);
        end

        %Check that we have the minimum number of particles for a minimum
        %amount of time
        
            VectorAP = NaN(length(data(j).ElapsedTime), length(APbinID));
%             ErrorAP = NaN(length(data(j).ElapsedTime), length(APbinID));
            for n = 1:length(Frames{j,i})
                VectorAP(Frames{j,i}(n),currentAPBin) = Fluos{j,i}(n);
%                 ErrorAP(Frames{j,i}(n),currentAPBin) = FluoErrors{j,i}(n);
            end


            %Extract the data for this range of frames
            FluoData=VectorAP(FrameWindow,currentAPBin);
%             SDFluoData=ErrorAP(FrameWindow,i);
            %NData=NParticlesAP(FrameWindow,i);
            TimeData=data(j).ElapsedTime(FrameWindow);
          
            %Now filter them according the number of particles
%             FrameFilter=NData>=MinParticles;


            %As an initial guess, use FrameFilter to determine the range of the
            %fit
            if isempty(FitResultsMultiple{j}(i,CurrentNC-11).FitFrameRange)
                FitFrameRange=FrameWindow;
                if CurrentNC==14
                    FitFrameRange=FitFrameRange((data(j).ElapsedTime(FitFrameRange)-data(j).ElapsedTime(data(j).APDivision(CurrentNC,currentAPBin)))<12);
                end
                FitResultsMultiple{j}(i,CurrentNC-11).FitFrameRange=FitFrameRange;
            else
                FitFrameRange=FitResultsMultiple{j}(i,CurrentNC-11).FitFrameRange;
            end

            %Filter the frames according to FitFrameRange
            FitFrameFilter=ismember(FrameWindow,FitFrameRange);
            
           

            FluoDataForFit=FluoData;
%             SDFluoDataForFit=SDFluoData;
%             NDataForFit=NData;
            TimeDataForFit=TimeData;
            


    

            if CurrentNC~=14
                %Do the fit
                x0=[FitResultsMultiple{j}(i,CurrentNC-11).TimeStart0,FitResultsMultiple{j}(i,CurrentNC-11).TimeEnd0,FitResultsMultiple{j}(i,CurrentNC-11).Rate0];


                
                %Get rid of any NaN in the data
                NanFilter=~isnan(FluoDataForFit);

                if ~isempty(TimeData(NanFilter))

                    [xFit,resnorm,residual,exitflag,output,lambda,jacobian]=...
                        lsqnonlin(@(x) lsqnonlinFitFluorescenceCurve(TimeDataForFit(NanFilter)-...
                        data(j).ElapsedTime(FrameWindow(1)),...
                        FluoDataForFit(NanFilter),Delay,...
                        data(j).ElapsedTime(FrameWindow(end))-data(j).ElapsedTime(FrameWindow(1)),x),x0);

                    FitResultsMultiple{j}(i,CurrentNC-11).TimeStart=xFit(1);
                    FitResultsMultiple{j}(i,CurrentNC-11).TimeEnd=xFit(2);
                    FitResultsMultiple{j}(i,CurrentNC-11).RateFit=xFit(3);

                    %Estimate an error bar out of the confidence intervals
                    FitResultsMultiple{j}(i,CurrentNC-11).CI=nlparci(xFit,residual,'jacobian',jacobian);

                    FitResultsMultiple{j}(i,CurrentNC-11).SDTimeStart=(FitResultsMultiple{j}(i,CurrentNC-11).CI(1,2)-FitResultsMultiple{j}(i,CurrentNC-11).CI(1,1))/2;
                    FitResultsMultiple{j}(i,CurrentNC-11).SDTimeEnd=(FitResultsMultiple{j}(i,CurrentNC-11).CI(2,2)-FitResultsMultiple{j}(i,CurrentNC-11).CI(2,1))/2;
                    FitResultsMultiple{j}(i,CurrentNC-11).SDRateFit=(FitResultsMultiple{j}(i,CurrentNC-11).CI(3,2)-FitResultsMultiple{j}(i,CurrentNC-11).CI(3,1))/2;




                    %Plot the results
                    %Get the corresponding fitted curve
                    [TimeFit,FluoFit]=FluorescenceCurve(data(j).ElapsedTime(FrameWindow(end))-...
                        data(j).ElapsedTime(FrameWindow(1)),...
                        xFit(1),xFit(2),xFit(3),Delay);
                    %Plot all the data
                    PlotHandle=plot(data(j).ElapsedTime(FrameWindow)-data(j).ElapsedTime(FrameWindow(1)),...
                        VectorAP(FrameWindow,currentAPBin),'.-k');
                    hold on
%                     Plot the data that could be used for the fit
                    PlotHandle(end+1)=plot(data(j).ElapsedTime(FrameWindow)-data(j).ElapsedTime(FrameWindow(1)),...
                        FluoData,'or');
%                     Plot the data that was actually used for the fit
                    PlotHandle(end+1)=plot(data(j).ElapsedTime(FitFrameRange)-data(j).ElapsedTime(FrameWindow(1)),...
                        FluoData(ismember(FrameWindow,FitFrameRange)),'or','MarkerFaceColor','r');
                    
                    %Plot the fit
                    PlotHandle(end+1)=plot(TimeFit,FluoFit,'-r');
                    hold off
                    %StandardFigure(PlotHandle,gca)
                    ylabel('Mean fluorescence nucleus')
                    xlabel('Time into nc (min)')
                    
                    try
                        ylim([0,max(VectorAP(FrameWindow,currentAPBin))]);
                    catch
                        display('Error in displaying the plot')
                    end

                    legend(['tON=',num2str(FitResultsMultiple{j}(i,CurrentNC-11).TimeStart),' \pm ',num2str(FitResultsMultiple{j}(i,CurrentNC-11).SDTimeStart)],...
                        ['tOFF=',num2str(FitResultsMultiple{j}(i,CurrentNC-11).TimeEnd),' \pm ',num2str(FitResultsMultiple{j}(i,CurrentNC-11).SDTimeEnd)],...
                        ['Rate=',num2str(FitResultsMultiple{j}(i,CurrentNC-11).RateFit),' \pm ',num2str(FitResultsMultiple{j}(i,CurrentNC-11).SDRateFit)],...
                        'Location','SouthOutside')
                end
            elseif CurrentNC==14
                
                %Do the fit
                x0=[FitResultsMultiple{j}(i,CurrentNC-11).TimeStart0,FitResultsMultiple{j}(i,CurrentNC-11).Rate0];


                
                %Get rid of any NaN in the data
                NanFilter=~isnan(FluoDataForFit);

                if ~isempty(TimeData(NanFilter))

                    [xFit,resnorm,residual,exitflag,output,lambda,jacobian]=...
                        lsqnonlin(@(x) lsqnonlinFitFluorescenceCurveNC14(TimeDataForFit(NanFilter)-...
                        data(j).ElapsedTime(FrameWindow(1)),...
                        FluoDataForFit(NanFilter),Delay,...
                        data(j).ElapsedTime(FrameWindow(end))-data(j).ElapsedTime(FrameWindow(1)),x),x0);

                    FitResultsMultiple{j}(i,CurrentNC-11).TimeStart=xFit(1);
                    FitResultsMultiple{j}(i,CurrentNC-11).RateFit=xFit(2);

                    %Estimate an error bar out of the confidence intervals
                    FitResultsMultiple{j}(i,CurrentNC-11).CI=nlparci(xFit,residual,'jacobian',jacobian);

                    FitResultsMultiple{j}(i,CurrentNC-11).SDTimeStart=(FitResultsMultiple{j}(i,CurrentNC-11).CI(1,2)-FitResultsMultiple{j}(i,CurrentNC-11).CI(1,1))/2;
                    FitResultsMultiple{j}(i,CurrentNC-11).SDRateFit=(FitResultsMultiple{j}(i,CurrentNC-11).CI(2,2)-FitResultsMultiple{j}(i,CurrentNC-11).CI(2,1))/2;




                    %Plot the results
                    %Get the corresponding fitted curve
                    [TimeFit,FluoFit]=FluorescenceCurve(data(j).ElapsedTime(FrameWindow(end))-...
                        data(j).ElapsedTime(FrameWindow(1)),...
                        xFit(1),1000,xFit(2),Delay);
                    %Plot all the data
                    PlotHandle=plot(data(j).ElapsedTime(FrameWindow)-data(j).ElapsedTime(FrameWindow(1)),...
                        VectorAP(FrameWindow,currentAPBin),'.-k');
                    hold on
                    %Plot the data that could be used for the fit
                    PlotHandle(end+1)=plot(data(j).ElapsedTime(FrameWindow)-data(j).ElapsedTime(FrameWindow(1)),...
                        FluoData,'or');
                    %Plot the data that was actually used for the fit
                    PlotHandle(end+1)=plot(data(j).ElapsedTime(FitFrameRange)-data(j).ElapsedTime(FrameWindow(1)),...
                        FluoData(ismember(FrameWindow,FitFrameRange)),'or','MarkerFaceColor','r');
                    
                    %Plot the fit
                    PlotHandle(end+1)=plot(TimeFit,FluoFit,'-r');
                    hold off
                    ylabel('Mean fluorescence nucleus')
                    xlabel('Time into nc (min)')
                    
                    
                    try
                        ylim([0,max(VectorAP(FrameWindow,currentAPBin))])
                    catch
                        display('Error in displaying the plot')
                    end

                    legend(['tON=',num2str(FitResultsMultiple{j}(i,CurrentNC-11).TimeStart),' \pm ',num2str(FitResultsMultiple{j}(i,CurrentNC-11).SDTimeStart)],...
                        ['Rate=',num2str(FitResultsMultiple{j}(i,CurrentNC-11).RateFit),' \pm ',num2str(FitResultsMultiple{j}(i,CurrentNC-11).SDRateFit)],...
                        'Location','SouthOutside')
                end
            end

        
    end
    
    title(['Data Set = ',num2str(j),', ','Particle = ',num2str(i),'/',num2str(length(CP{j})),', ',num2str(data(j).APbinID(currentAPBin)),' AP, TimeStart0=',num2str(FitResultsMultiple{j}(i,CurrentNC-11).TimeStart0),...
        ', TimeEnd0=',num2str(FitResultsMultiple{j}(i,CurrentNC-11).TimeEnd0),', Rate=',num2str(FitResultsMultiple{j}(i,CurrentNC-11).Rate0),...
        ', nc',num2str(CurrentNC)])
    
    
    %Set the limits on the x-axis
    if CurrentNC==14
        %xlim([0,data(j).ElapsedTime(end)])
        xlim([0,60])
%     else
%         xlim([0,data(j).ElapsedTime(eval(['nc',num2str(nc+1)]))-...
%             data(j).ElapsedTime(eval(['nc',num2str(nc)]))])
    end
    
    
    
    figure(FitFigure)
    try
        ct=waitforbuttonpress;
    catch
        error('Fits not saved.');
    end
    cc=get(FitFigure,'currentcharacter');
    cm=get(gca,'CurrentPoint');
    
    %Move between particles within a dataset
    if (ct~=0)&(cc=='.')&(i<length(CP{j}))
        i=i+1;
    elseif (ct~=0)&(cc==',')&(i>1)
        i=i-1;
    
    %Move between datasets
    elseif (ct~=0)&(cc=='>')&(j<length(data))
        j=j+1;
        i=1;
    elseif (ct~=0)&(cc=='<')&(j>1)
        j=j-1;
        i=1;
    
    
    
    %Approve, disapprove fit
    elseif (ct~=0)&(cc=='q')
        if FitResultsMultiple{j}(i,CurrentNC-11).Approved==0
            FitResultsMultiple{j}(i,CurrentNC-11).Approved=1;
        elseif FitResultsMultiple{j}(i,CurrentNC-11).Approved==1
            FitResultsMultiple{j}(i,CurrentNC-11).Approved=0;
        end

    
    %Disapprove, disapprove fit
    elseif (ct~=0)&(cc=='w')
        if FitResultsMultiple{j}(i,CurrentNC-11).Approved==0
            FitResultsMultiple{j}(i,CurrentNC-11).Approved=-1;
        elseif FitResultsMultiple{j}(i,CurrentNC-11).Approved==-1
            FitResultsMultiple{j}(i,CurrentNC-11).Approved=0;
        end
  
    
        
    %Move right range of fit
    elseif (ct~=0)&(cc=='k')&(length(FitResultsMultiple{j}(i,CurrentNC-11).FitFrameRange)>2)
        FitResultsMultiple{j}(i,CurrentNC-11).FitFrameRange=FitResultsMultiple{j}(i,CurrentNC-11).FitFrameRange(1:end-1);
    elseif (ct~=0)&(cc=='l')
        if ~isempty(find(~ismember(FrameWindow,FitResultsMultiple{j}(i,CurrentNC-11).FitFrameRange)))
            FilteredFramesTemp=FrameWindow;
            FitResultsMultiple{j}(i,CurrentNC-11).FitFrameRange(end+1)=...
                FilteredFramesTemp(min(find(~ismember(FilteredFramesTemp,FitResultsMultiple{j}(i,CurrentNC-11).FitFrameRange))));
        end
    %Move left range of fit
    elseif (ct~=0)&(cc=='j')&(length(FitResultsMultiple{j}(i,CurrentNC-11).FitFrameRange)>2)
        FitResultsMultiple{j}(i,CurrentNC-11).FitFrameRange=FitResultsMultiple{j}(i,CurrentNC-11).FitFrameRange(2:end);
    elseif (ct~=0)&(cc=='h')
        if ~isempty(find(~ismember(FrameWindow,FitResultsMultiple{j}(i,CurrentNC-11).FitFrameRange)))
            FilteredFramesTemp=FrameWindow;
            FitResultsMultiple{j}(i,CurrentNC-11).FitFrameRange=...
                [FilteredFramesTemp(max(find(~ismember(FilteredFramesTemp,FitResultsMultiple{j}(i,CurrentNC-11).FitFrameRange)))),...
                FitResultsMultiple{j}(i,CurrentNC-11).FitFrameRange];
        end
    %Reset frame fit range
     elseif (ct~=0)&(cc=='r')   
        FitResultsMultiple{j}(i,CurrentNC-11).FitFrameRange=FrameWindow;

        
        
    %Change the initial parameters
    %TimeStart
    elseif (ct~=0)&(cc=='a')&((CurrentNC==14)|(CurrentNC~=14&FitResultsMultiple{j}(i,CurrentNC-11).TimeStart0<FitResultsMultiple{j}(i,CurrentNC-11).TimeEnd0))
        FitResultsMultiple{j}(i,CurrentNC-11).TimeStart0=FitResultsMultiple{j}(i,CurrentNC-11).TimeStart0+1;
    elseif (ct~=0)&(cc=='z')&(FitResultsMultiple{j}(i,CurrentNC-11).TimeStart0>1)
        FitResultsMultiple{j}(i,CurrentNC-11).TimeStart0=FitResultsMultiple{j}(i,CurrentNC-11).TimeStart0-1;
    %TimeEnd
    elseif (ct~=0)&(cc=='s')&(FitResultsMultiple{j}(i,CurrentNC-11).TimeEnd0<data(j).ElapsedTime(FrameWindow(end))-data(j).ElapsedTime(FrameWindow(1)))
        FitResultsMultiple{j}(i,CurrentNC-11).TimeEnd0=FitResultsMultiple{j}(i,CurrentNC-11).TimeEnd0+1;
    elseif (ct~=0)&(cc=='x')&(FitResultsMultiple{j}(i,CurrentNC-11).TimeEnd0>FitResultsMultiple{j}(i,CurrentNC-11).TimeStart0)
        FitResultsMultiple{j}(i,CurrentNC-11).TimeEnd0=FitResultsMultiple{j}(i,CurrentNC-11).TimeEnd0-1;
    %Rate, fine
    elseif (ct~=0)&(cc=='c')&(FitResultsMultiple{j}(i,CurrentNC-11).Rate0>100)
        FitResultsMultiple{j}(i,CurrentNC-11).Rate0=FitResultsMultiple{j}(i,CurrentNC-11).Rate0-100;
    elseif (ct~=0)&(cc=='d')
        FitResultsMultiple{j}(i,CurrentNC-11).Rate0=FitResultsMultiple{j}(i,CurrentNC-11).Rate0+100;    
    %Rate, coarse
    elseif (ct~=0)&(cc=='C')&(FitResultsMultiple{j}(i,CurrentNC-11).Rate0>100)
        FitResultsMultiple{j}(i,CurrentNC-11).Rate0=FitResultsMultiple{j}(i,CurrentNC-11).Rate0-500;
    elseif (ct~=0)&(cc=='D')
        FitResultsMultiple{j}(i,CurrentNC-11).Rate0=FitResultsMultiple{j}(i,CurrentNC-11).Rate0+500; 
    
    %Switch NCs
    elseif (ct~=0)&(cc=='m')&CurrentNC<14
        CurrentNC=CurrentNC+1;
    elseif (ct~=0)&(cc=='n')&CurrentNC>MinNC
        CurrentNC=CurrentNC-1;
        
        
    %Save
    elseif (ct~=0)&(cc=='v')
        save([DropboxFolder,filesep,Prefix,filesep,'FitResultsMultiple.mat'],...
        'FitResultsMultiple')
    display('FitResultsMultiple.mat saved')
        
    %Debug mode
    elseif (ct~=0)&(cc=='9')
        keyboard
    
    end
    
end


%Save the information
save([DropboxFolder,filesep,Prefix,filesep,'FitResultsMultiple.mat'],...
    'FitResultsMultipleFitResultsMultiple')
display('FitResultsMultiple.mat saved')        
        
close(FitFigure)