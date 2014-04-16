function FitIndividualFluorescenceCurve(varargin)

%This function performs fits to the single particle fluorescence as a
%function of time of a particular dataset.


%Usage:
%, .: Move between particles
%n m: Move between points on a trace
%- =: Change the amount of points in a trace
%a z: Move the time point
%s x: Move the rate up and down
%SHIFT + a,z,s,x: Coarse move
%f  : Fit this trace
%b  : Switch between nc13 and nc14


close all

[SourcePath,FISHPath,DefaultDropboxFolder,MS2CodePath,SchnitzcellsFolder]=...
    DetermineLocalFolders;



if ~isempty(varargin)
    Prefix=varargin{1};
else
    FolderTemp=uigetdir(DefaultDropboxFolder,'Choose folder with files to analyze');
    Dashes=strfind(FolderTemp,'\');
    Prefix=FolderTemp((Dashes(end)+1):end);
end


%Get the actual folder using the Prefix
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,SchnitzcellsFolder]=...
    DetermineLocalFolders(Prefix);



%Load the compiled particles and the division information                                    
load([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat'])


%Determine what type of file format we're using. This is important in terms
%of the fluorescence units and steps
load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'])

if isfield(FrameInfo,'FileMode')
    FileMode=FrameInfo(1).FileMode;
else
    FileMode=[];
end

%Determine the step of rates for the manual fit
if strcmp(FileMode,'LSM')
    ScaleFactor=10;
else
    ScaleFactor=1;
end


%Figure out the gene length. 
[XLSNum,XLSTxt]=xlsread([DefaultDropboxFolder,filesep,'MovieDatabase.xlsx']);
DataFolderColumn=find(strcmp(XLSTxt(1,:),'DataFolder'));
StemLoopColumn=find(strcmp(XLSTxt(1,:),'StemLoop'));

% Convert the prefix into the string used in the XLS file
Dashes = strfind(Prefix, '-');
PrefixRow = find(strcmp(XLSTxt(:, DataFolderColumn),...
    [Prefix(1:Dashes(3)-1), filesep, Prefix(Dashes(3)+1:end)]));
StemLoop=XLSTxt(PrefixRow,StemLoopColumn);

%Distance from the first MS2 site to the end of the 3' UTR
if ~isempty(strfind(StemLoop,'Eve'))
    GeneLength=6.443;
elseif  ~isempty(strfind(StemLoop,'X1'))
    GeneLength=5.296;       %Distance from the first MS2 site to the end of the
                        %TUB3'UTR in kb.
else
    error('The gene length has not been defined for this construct')
end
    

%Parameters:
MinTimePoints=5;    %Minimum number of time points where we'll have per particle.
ElongationRate=1.54;    %In kb/minutes.
Delay=GeneLength/ElongationRate;    %Minutes for PolII to fall off after reaching
                                    %the first MS2 site.

              

ParticleNCFilter{1}=([CompiledParticles.nc]==13);
ParticlesNC{1}=find(ParticleNCFilter{1});

ParticleNCFilter{2}=([CompiledParticles.nc]==14);
ParticlesNC{2}=find(ParticleNCFilter{2});


%Initialize the fit results structure
for j=1:length(ParticlesNC)
    for i=1:length(ParticlesNC{j})
            FitResultsIndiv(i,j).nSteps=1;
            FitResultsIndiv(i,j).ManualTransitions=4;
            FitResultsIndiv(i,j).ManualRates=0.6E4;
            FitResultsIndiv(i,j).Approved=0;
            FitResultsIndiv(i,j).CompiledParticle=ParticlesNC{j}(i);
            FitResultsIndiv(i,j).FittedTransitions=[];
    end
end

%Load the fits if they exist
if exist([DropboxFolder,filesep,Prefix,filesep,'IndividualFits.mat'])
    load([DropboxFolder,filesep,Prefix,filesep,'IndividualFits.mat'])
end


%Go through each Particle
FitFigure=figure;
set(gcf,'Position',[35   226   560   420]);

RateFigure=figure;
set(gcf,'Position',[630   226   560   420])

i=1;
cc=1;
CurrentTransition=1;
nc=13;
while cc~=13
    figure(FitFigure)
    clf
    
    errorbar(ElapsedTime(CompiledParticles(ParticlesNC{nc-12}(i)).Frame)-...
        ElapsedTime(eval(['nc',num2str(nc)])),...
        CompiledParticles(ParticlesNC{nc-12}(i)).Fluo,...
        ones(size(CompiledParticles(ParticlesNC{nc-12}(i)).Fluo))*...
        CompiledParticles(ParticlesNC{nc-12}(i)).FluoError,'.-k')
    
    %Plot the fitting result
    [TimePrediction,FluoPrediction]=IndividualTrace(FitResultsIndiv(i,nc-12).ManualTransitions,...
        FitResultsIndiv(i,nc-12).ManualRates,Delay,80);
    hold on
    plot(TimePrediction,FluoPrediction,'-r')
    hold off

    %Locate the transition points
    hold on
    for j=1:FitResultsIndiv(i,nc-12).nSteps
        %Find the closest index to this time point
        TimeTransitionTemp=TimePrediction(max(find(TimePrediction<=FitResultsIndiv(i,nc-12).ManualTransitions(j))));
        FluoTransitionTemp=FluoPrediction(max(find(TimePrediction<=FitResultsIndiv(i,nc-12).ManualTransitions(j))));
        
        if CurrentTransition==j
            plot(TimeTransitionTemp,FluoTransitionTemp,'or','MarkerFaceColor','r')
        else
            plot(TimeTransitionTemp,FluoTransitionTemp,'ok','MarkerFaceColor','k')
        end
    end
    hold off
        
        
    %Plot the fitting results if they exist
    if FitResultsIndiv(i,nc-12).Approved==1
        if isfield(FitResultsIndiv,'FittedTransitions')
            if ~isempty(FitResultsIndiv(i,nc-12).FittedTransitions)
                hold on
                [TimeFit,FluoFit]=IndividualTrace(FitResultsIndiv(i,nc-12).FittedTransitions,...
                    FitResultsIndiv(i,nc-12).FittedRates,Delay,80);
                plot(TimeFit,FluoFit,'--b')
                hold off


                %Locate the transition points
                hold on
                for j=1:FitResultsIndiv(i,nc-12).nSteps
                    %Find the closest index to this time point
                    TimeTransitionTemp=TimeFit(max(find(TimePrediction<=FitResultsIndiv(i,nc-12).FittedTransitions(j))));
                    FluoTransitionTemp=FluoFit(max(find(TimePrediction<=FitResultsIndiv(i,nc-12).FittedTransitions(j))));

                    plot(TimeTransitionTemp,FluoTransitionTemp,'ob')
                end
                hold off

            end
        end
    end
        
        
    
    
    
    %Do the title
    if isfield(CompiledParticles,'MeanAP')
        title(['nc',num2str(nc),' (',num2str(i),'/',num2str(length(ParticlesNC{nc-12})),...
            '), Rate: ',num2str(FitResultsIndiv(i,nc-12).ManualRates(CurrentTransition)),...
            ', AP: ',num2str(CompiledParticles(ParticlesNC{nc-12}(i)).MeanAP)])
    else
        title(['nc',num2str(nc),' (',num2str(i),'/',num2str(length(ParticlesNC{nc-12})),...
            '), Rate: ',num2str(FitResultsIndiv(i,nc-12).ManualRates(CurrentTransition)),...
            ', x: ',num2str(mean(CompiledParticles(ParticlesNC{nc-12}(i)).xPos))])
    end
        
    %Show the flag status
    if FitResultsIndiv(i,nc-12).Approved==0
        set(gca,'Color','default')
    elseif FitResultsIndiv(i,nc-12).Approved==1
        set(gca,'Color','g')
    elseif FitResultsIndiv(i,nc-12).Approved==-1
        set(gca,'Color','r')
    end
    
    
    %Plot the rates we have
    figure(RateFigure)
    clf
    
    %These are the manual ones
    hold on
    for j=0:FitResultsIndiv(i,nc-12).nSteps
        if j<FitResultsIndiv(i,nc-12).nSteps
            if j==0
                xRange=linspace(0,FitResultsIndiv(i,nc-12).ManualTransitions(1));
                Rate=ones(size(xRange))*0;
            else
                xRangeTemp=linspace(FitResultsIndiv(i,nc-12).ManualTransitions(j),...
                    FitResultsIndiv(i,nc-12).ManualTransitions(j+1));
                RateTemp=ones(size(xRangeTemp))*FitResultsIndiv(i,nc-12).ManualRates(j);
                xRange=[xRange,xRangeTemp];
                Rate=[Rate,RateTemp];
            end
        else
            if nc==14
                xRangeTemp=linspace(FitResultsIndiv(i,nc-12).ManualTransitions(j),80);
            elseif nc==13
                xRangeTemp=linspace(FitResultsIndiv(i,nc-12).ManualTransitions(j),...
                    ElapsedTime(eval(['nc',num2str(nc+1)]))-ElapsedTime(eval(['nc',num2str(nc)])));
            end
            RateTemp=ones(size(xRangeTemp))*FitResultsIndiv(i,nc-12).ManualRates(j);
            xRange=[xRange,xRangeTemp];
            Rate=[Rate,RateTemp];
        end
    end
                
    plot(xRange,Rate,'-r','LineWidth',4)
    ylim([0,max(Rate)*1.1])
    hold off
    
    
    %These are the fitted ones
    hold on
    if ~isempty(FitResultsIndiv(i,nc-12).FittedTransitions)
        for j=0:length(FitResultsIndiv(i,nc-12).FittedTransitions)%nSteps
            if j<length(FitResultsIndiv(i,nc-12).FittedTransitions)%nSteps
                if j==0
                    xRange=linspace(0,FitResultsIndiv(i,nc-12).FittedTransitions(1));
                    Rate=ones(size(xRange))*0;
                else
                    xRangeTemp=linspace(FitResultsIndiv(i,nc-12).FittedTransitions(j),...
                        FitResultsIndiv(i,nc-12).FittedTransitions(j+1));
                    RateTemp=ones(size(xRangeTemp))*FitResultsIndiv(i,nc-12).FittedRates(j);
                    xRange=[xRange,xRangeTemp];
                    Rate=[Rate,RateTemp];
                end
            else
                if nc==14
                    xRangeTemp=linspace(FitResultsIndiv(i,nc-12).FittedTransitions(j),80);
                elseif nc==13
                    xRangeTemp=linspace(FitResultsIndiv(i,nc-12).FittedTransitions(j),...
                        ElapsedTime(eval(['nc',num2str(nc+1)]))-ElapsedTime(eval(['nc',num2str(nc)])));
                end
                RateTemp=ones(size(xRangeTemp))*FitResultsIndiv(i,nc-12).FittedRates(j);
                xRange=[xRange,xRangeTemp];
                Rate=[Rate,RateTemp];
            end
        end

        plot(xRange,Rate,'-b','LineWidth',4)
        ylim([0,max(Rate)*1.1])
        hold off
    end
    
    
    %Set the limits on the x-axis
    if nc==14
        %xlim([0,ElapsedTime(end)])
        xlim([0,80])
    else
        xlim([0,ElapsedTime(eval(['nc',num2str(nc+1)]))-...
            ElapsedTime(eval(['nc',num2str(nc)]))])
    end
    
    figure(FitFigure)
    ct=waitforbuttonpress;
    cc=get(FitFigure,'currentcharacter');
    cm=get(gca,'CurrentPoint');
    
    %Move between particles positions
    if (ct~=0)&(cc=='.')&(i<length(ParticlesNC{nc-12}))
        i=i+1;
        CurrentTransition=1;
    elseif (ct~=0)&(cc==',')&(i>1)
        i=i-1;
        CurrentTransition=1;
        
    %Add or remove transition points
    elseif (ct~=0)&(cc=='=')
        FitResultsIndiv(i,nc-12).nSteps=FitResultsIndiv(i,nc-12).nSteps+1;
        FitResultsIndiv(i,nc-12).ManualTransitions(end+1)=...
            FitResultsIndiv(i,nc-12).ManualTransitions(end)+3;
        FitResultsIndiv(i,nc-12).ManualRates(end+1)=0.6E4;
        CurrentTransition=FitResultsIndiv(i,nc-12).nSteps;
        FitResultsIndiv(i,nc-12).Approved=0;
    elseif (ct~=0)&(cc=='-')
        if FitResultsIndiv(i,nc-12).nSteps>1
            FitResultsIndiv(i,nc-12).nSteps=FitResultsIndiv(i,nc-12).nSteps-1;
            FitResultsIndiv(i,nc-12).ManualTransitions=FitResultsIndiv(i,nc-12).ManualTransitions(1:end-1);
            FitResultsIndiv(i,nc-12).ManualRates=FitResultsIndiv(i,nc-12).ManualRates(1:end-1);
            CurrentTransition=FitResultsIndiv(i,nc-12).nSteps;
            FitResultsIndiv(i,nc-12).Approved=0;
        end
        
    %Move between transition points
    elseif (ct~=0)&(cc=='m')&(CurrentTransition<FitResultsIndiv(i,nc-12).nSteps)
        CurrentTransition=CurrentTransition+1;
    elseif (ct~=0)&(cc=='n')&(CurrentTransition>1)
        CurrentTransition=CurrentTransition-1;
        
    %Change the parameters of the transition
    %Time, fine resolution
    elseif (ct~=0)&(cc=='a')
        FitResultsIndiv(i,nc-12).ManualTransitions(CurrentTransition)=...
            FitResultsIndiv(i,nc-12).ManualTransitions(CurrentTransition)+0.2;
        FitResultsIndiv(i,nc-12).Approved=0;
    elseif (ct~=0)&(cc=='z')
        FitResultsIndiv(i,nc-12).ManualTransitions(CurrentTransition)=...
            FitResultsIndiv(i,nc-12).ManualTransitions(CurrentTransition)-0.2;
        FitResultsIndiv(i,nc-12).Approved=0;
    %Time, coarse resolution
    elseif (ct~=0)&(cc=='A')
        FitResultsIndiv(i,nc-12).ManualTransitions(CurrentTransition)=...
            FitResultsIndiv(i,nc-12).ManualTransitions(CurrentTransition)+2;
        FitResultsIndiv(i,nc-12).Approved=0;
    elseif (ct~=0)&(cc=='Z')
        FitResultsIndiv(i,nc-12).ManualTransitions(CurrentTransition)=...
            FitResultsIndiv(i,nc-12).ManualTransitions(CurrentTransition)-2;
        FitResultsIndiv(i,nc-12).Approved=0;
        
        
    %Rate, fine resolution
    elseif (ct~=0)&(cc=='s')
        FitResultsIndiv(i,nc-12).ManualRates(CurrentTransition)=...
            FitResultsIndiv(i,nc-12).ManualRates(CurrentTransition)+200*ScaleFactor;
        FitResultsIndiv(i,nc-12).Approved=0;
    elseif (ct~=0)&(cc=='x')&(FitResultsIndiv(i,nc-12).ManualRates(CurrentTransition)>0)
        FitResultsIndiv(i,nc-12).ManualRates(CurrentTransition)=...
            FitResultsIndiv(i,nc-12).ManualRates(CurrentTransition)-200*ScaleFactor;
        FitResultsIndiv(i,nc-12).Approved=0;
        
    %Rate, coarse resolution
    elseif (ct~=0)&(cc=='S')
        FitResultsIndiv(i,nc-12).ManualRates(CurrentTransition)=...
            FitResultsIndiv(i,nc-12).ManualRates(CurrentTransition)+2000*ScaleFactor;
        FitResultsIndiv(i,nc-12).Approved=0;
    elseif (ct~=0)&(cc=='X')&(FitResultsIndiv(i,nc-12).ManualRates(CurrentTransition)>0)
        FitResultsIndiv(i,nc-12).ManualRates(CurrentTransition)=...
            FitResultsIndiv(i,nc-12).ManualRates(CurrentTransition)-2000*ScaleFactor;
        FitResultsIndiv(i,nc-12).Approved=0;
    
        
    %Do the fit
    elseif (ct~=0)&(cc=='f')
        x0=[FitResultsIndiv(i,nc-12).ManualTransitions,...
            FitResultsIndiv(i,nc-12).ManualRates];
        %Create lower and upper bounds for the rates
        x0Lower=zeros(size(x0));
        x0Upper=inf(size(x0));
        
        options = optimset('MaxFunEvals',2000*length(x0),'MaxIter',5000);
        
        [xFit,resnorm,residual,exitflag,output,lambda,jacobian]=...
            lsqnonlin(@(x) lsqnonlinFitIndividualFluorescenceCurve(ElapsedTime(CompiledParticles(ParticlesNC{nc-12}(i)).Frame)-...
            ElapsedTime(eval(['nc',num2str(nc)])),...
            CompiledParticles(ParticlesNC{nc-12}(i)).Fluo,Delay,FitResultsIndiv(i,nc-12).nSteps,x),x0,...
            x0Lower,x0Upper,options);

        FitResultsIndiv(i,nc-12).FittedTransitions=xFit(1:FitResultsIndiv(i,nc-12).nSteps);
        FitResultsIndiv(i,nc-12).FittedRates=xFit(FitResultsIndiv(i,nc-12).nSteps+1:end);
        
        FitResultsIndiv(i,nc-12).Approved=1;
        
    %Auto Fit
    elseif (ct~=0)&(cc=='g')
        [AutoTransitions,AutoRates] = AutoFitFluorescenceCurve(ElapsedTime(CompiledParticles(ParticlesNC{nc-12}(i)).Frame)-...
            ElapsedTime(eval(['nc',num2str(nc)])),...
            CompiledParticles(ParticlesNC{nc-12}(i)).Fluo,...
            CompiledParticles(ParticlesNC{nc-12}(i)).FluoError);

        FitResultsIndiv(i,nc-12).FittedTransitions=AutoTransitions;
        FitResultsIndiv(i,nc-12).FittedRates=AutoRates;
        
        FitResultsIndiv(i,nc-12).Approved=1;
        
        
%         
%                     FitResultsIndiv(i,CurrentNC-12).TimeStart=xFit(1);
%                     FitResultsIndiv(i,CurrentNC-12).TimeEnd=xFit(2);
%                     FitResultsIndiv(i,CurrentNC-12).RateFit=xFit(3);
% 
%                     %Estimate an error bar out of the confidence intervals
%                     FitResultsIndiv(i,CurrentNC-12).CI=nlparci(xFit,residual,'jacobian',jacobian);
% 
%                     FitResultsIndiv(i,CurrentNC-12).SDTimeStart=(FitResultsIndiv(i,CurrentNC-12).CI(1,2)-FitResultsIndiv(i,CurrentNC-12).CI(1,1))/2;
%                     FitResultsIndiv(i,CurrentNC-12).SDTimeEnd=(FitResultsIndiv(i,CurrentNC-12).CI(2,2)-FitResultsIndiv(i,CurrentNC-12).CI(2,1))/2;
%                     FitResultsIndiv(i,CurrentNC-12).SDRateFit=(FitResultsIndiv(i,CurrentNC-12).CI(3,2)-FitResultsIndiv(i,CurrentNC-12).CI(3,1))/2;
% 
%     
    %Approve/dissaprove
    elseif (ct~=0)&(cc=='q')
        if FitResultsIndiv(i,nc-12).Approved==0
            FitResultsIndiv(i,nc-12).Approved=1;
        elseif FitResultsIndiv(i,nc-12).Approved==1
            FitResultsIndiv(i,nc-12).Approved=-1;
        elseif FitResultsIndiv(i,nc-12).Approved==-1
            FitResultsIndiv(i,nc-12).Approved=0;
        end

    
    
    %Save results
    elseif (ct~=0)&(cc=='t')
        save([DropboxFolder,filesep,Prefix,filesep,'IndividualFits.mat'],...
            'FitResultsIndiv')
            display('IndividualFits.mat saved')
     
    %Switch ncs
    elseif (ct~=0)&(cc=='b')
        if nc==13
            if i>length(ParticlesNC{2})
                i=length(ParticlesNC{2});
            end
            nc=14;
            CurrentTransition=1;
        elseif nc==14
            if i>length(ParticlesNC{1})
                i=length(ParticlesNC{1});
            end
            nc=13;
            CurrentTransition=1;
        end
        
        
    %Debug mode
    elseif (ct~=0)&(cc=='9')
        keyboard
     
    
    %Jump to a particle
    elseif (ct~=0)&(cc=='j')
        try
            JumpParticle=input('Select particle to jump to: ');
            if (JumpParticle>0)&(JumpParticle<length(ParticlesNC{nc-12}))
                i=JumpParticle;
            end
        end
    end
   
end

            
save([DropboxFolder,filesep,Prefix,filesep,'IndividualFits.mat'],...
    'FitResultsIndiv')
    display('IndividualFits.mat saved')
