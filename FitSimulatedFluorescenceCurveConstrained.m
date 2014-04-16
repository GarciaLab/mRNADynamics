function FitSimulatedFluorescenceCurveConstrained(Time,Fluo,GeneLength,ElongationRate)

if nargin<3
    GeneLength=6.443;
end
if nargin<4
    ElongationRate=1.54;
end

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

ScaleFactor=1;     

FitResultsIndiv.nSteps=1;
FitResultsIndiv.ManualTransitions=4;
FitResultsIndiv.SingleRate=0.6E4;
FitResultsIndiv.NActive=2;
FitResultsIndiv.ManualRates=FitResultsIndiv.SingleRate...
    * FitResultsIndiv.NActive;
FitResultsIndiv.Approved=0;
FitResultsIndiv.FittedTransitions=[];

%Go through each Particle
FitFigure=figure;
set(gcf,'Position',[35   226   560   420]);

RateFigure=figure;
set(gcf,'Position',[630   226   560   420])

cc=1;
CurrentTransition=1;

while cc~=13
    figure(FitFigure)
    clf
    
    errorbar(Time, Fluo, ones(size(Fluo)) * 0, '.-k')
    
    %Plot the fitting result
    [TimePrediction,FluoPrediction]=IndividualTrace(FitResultsIndiv.ManualTransitions,...
        FitResultsIndiv.ManualRates,Delay,80);
    hold on
    plot(TimePrediction,FluoPrediction,'-r')
    hold off

    %Locate the transition points
    hold on
    for j=1:FitResultsIndiv.nSteps
        %Find the closest index to this time point
        TimeTransitionTemp=TimePrediction(max(find(TimePrediction<=FitResultsIndiv.ManualTransitions(j))));
        FluoTransitionTemp=FluoPrediction(max(find(TimePrediction<=FitResultsIndiv.ManualTransitions(j))));
        
        if CurrentTransition==j
            plot(TimeTransitionTemp,FluoTransitionTemp,'or','MarkerFaceColor','r')
        else
            plot(TimeTransitionTemp,FluoTransitionTemp,'ok','MarkerFaceColor','k')
        end
    end
    hold off
        
        
    %Plot the fitting results if they exist
    if FitResultsIndiv.Approved==1
        if isfield(FitResultsIndiv,'FittedTransitions')
            if ~isempty(FitResultsIndiv.FittedTransitions)
                hold on
                [TimeFit,FluoFit]=IndividualTrace(FitResultsIndiv.FittedTransitions,...
                    FitResultsIndiv.FittedRates,Delay,80);
                plot(TimeFit,FluoFit,'--b')
                hold off


                %Locate the transition points
                hold on
                for j=1:FitResultsIndiv.nSteps
                    %Find the closest index to this time point
                    TimeTransitionTemp=TimeFit(max(find(TimePrediction<=FitResultsIndiv.FittedTransitions(j))));
                    FluoTransitionTemp=FluoFit(max(find(TimePrediction<=FitResultsIndiv.FittedTransitions(j))));

                    plot(TimeTransitionTemp,FluoTransitionTemp,'ob')
                end
                hold off

            end
        end
    end
        
        
    
    
    
    %Do the title
    if isfield(CompiledParticles,'MeanAP')
        title(['nc',num2str(nc),' (',num2str(i),'/',num2str(length(ParticlesNC{nc-12})),...
            '), Rate: ',num2str(FitResultsIndiv.ManualRates(CurrentTransition)),...
            ', AP: ',num2str(CompiledParticles(ParticlesNC{nc-12}(i)).MeanAP)])
    else
        title(['nc',num2str(nc),' (',num2str(i),'/',num2str(length(ParticlesNC{nc-12})),...
            '), Rate: ',num2str(FitResultsIndiv.ManualRates(CurrentTransition)),...
            ', x: ',num2str(mean(CompiledParticles(ParticlesNC{nc-12}(i)).xPos))])
    end
        
    %Show the flag status
    if FitResultsIndiv.Approved==0
        set(gca,'Color','default')
    elseif FitResultsIndiv.Approved==1
        set(gca,'Color','g')
    elseif FitResultsIndiv.Approved==-1
        set(gca,'Color','r')
    end
    
    
    %Plot the rates we have
    figure(RateFigure)
    clf
    
    %These are the manual ones
    hold on
    for j=0:FitResultsIndiv.nSteps
        if j<FitResultsIndiv.nSteps
            if j==0
                xRange=linspace(0,FitResultsIndiv.ManualTransitions(1));
                Rate=ones(size(xRange))*0;
            else
                xRangeTemp=linspace(FitResultsIndiv.ManualTransitions(j),...
                    FitResultsIndiv.ManualTransitions(j+1));
                RateTemp=ones(size(xRangeTemp))*FitResultsIndiv.ManualRates(j);
                xRange=[xRange,xRangeTemp];
                Rate=[Rate,RateTemp];
            end
        else
            if nc==14
                xRangeTemp=linspace(FitResultsIndiv.ManualTransitions(j),80);
            elseif nc==13
                xRangeTemp=linspace(FitResultsIndiv.ManualTransitions(j),...
                    ElapsedTime(eval(['nc',num2str(nc+1)]))-ElapsedTime(eval(['nc',num2str(nc)])));
            end
            RateTemp=ones(size(xRangeTemp))*FitResultsIndiv.ManualRates(j);
            xRange=[xRange,xRangeTemp];
            Rate=[Rate,RateTemp];
        end
    end
                
    plot(xRange,Rate,'-r','LineWidth',4)
    if max(Rate) == 0
        ylim([0,10])
    else
        ylim([0,max(Rate)*1.1])
    end
    hold off
    
    
    %These are the fitted ones
    hold on
    if ~isempty(FitResultsIndiv.FittedTransitions)
        for j=0:FitResultsIndiv.nSteps
            if j<FitResultsIndiv.nSteps
                if j==0
                    xRange=linspace(0,FitResultsIndiv.FittedTransitions(1));
                    Rate=ones(size(xRange))*0;
                else
                    xRangeTemp=linspace(FitResultsIndiv.FittedTransitions(j),...
                        FitResultsIndiv.FittedTransitions(j+1));
                    RateTemp=ones(size(xRangeTemp))*FitResultsIndiv.FittedRates(j);
                    xRange=[xRange,xRangeTemp];
                    Rate=[Rate,RateTemp];
                end
            else
                if nc==14
                    xRangeTemp=linspace(FitResultsIndiv.FittedTransitions(j),80);
                elseif nc==13
                    xRangeTemp=linspace(FitResultsIndiv.FittedTransitions(j),...
                        ElapsedTime(eval(['nc',num2str(nc+1)]))-ElapsedTime(eval(['nc',num2str(nc)])));
                end
                RateTemp=ones(size(xRangeTemp))*FitResultsIndiv.FittedRates(j);
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
        FitResultsIndiv.nSteps=FitResultsIndiv.nSteps+1;
        FitResultsIndiv.ManualTransitions(end+1)=...
            FitResultsIndiv.ManualTransitions(end)+3;
        FitResultsIndiv.NActive(end+1)=1;
        FitResultsIndiv.ManualRates(end+1)=FitResultsIndiv.SingleRate...
            *FitResultsIndiv.NActive(end);
        CurrentTransition=FitResultsIndiv.nSteps;
        FitResultsIndiv.Approved=0;
    elseif (ct~=0)&(cc=='-')
        if FitResultsIndiv.nSteps>1
            FitResultsIndiv.nSteps=FitResultsIndiv.nSteps-1;
            FitResultsIndiv.ManualTransitions=FitResultsIndiv.ManualTransitions(1:end-1);
            FitResultsIndiv.NActive=FitResultsIndiv.NActive(1:end-1);
            FitResultsIndiv.ManualRates=FitResultsIndiv.ManualRates(1:end-1);
            CurrentTransition=FitResultsIndiv.nSteps;
            FitResultsIndiv.Approved=0;
        end
        
    %Move between transition points
    elseif (ct~=0)&(cc=='m')&(CurrentTransition<FitResultsIndiv.nSteps)
        CurrentTransition=CurrentTransition+1;
    elseif (ct~=0)&(cc=='n')&(CurrentTransition>1)
        CurrentTransition=CurrentTransition-1;
        
    %Change the parameters of the transition
    %Time, fine resolution
    elseif (ct~=0)&(cc=='a')
        FitResultsIndiv.ManualTransitions(CurrentTransition)=...
            FitResultsIndiv.ManualTransitions(CurrentTransition)+0.2;
        FitResultsIndiv.Approved=0;
    elseif (ct~=0)&(cc=='z')
        FitResultsIndiv.ManualTransitions(CurrentTransition)=...
            FitResultsIndiv.ManualTransitions(CurrentTransition)-0.2;
        FitResultsIndiv.Approved=0;
    %Time, coarse resolution
    elseif (ct~=0)&(cc=='A')
        FitResultsIndiv.ManualTransitions(CurrentTransition)=...
            FitResultsIndiv.ManualTransitions(CurrentTransition)+2;
        FitResultsIndiv.Approved=0;
    elseif (ct~=0)&(cc=='Z')
        FitResultsIndiv.ManualTransitions(CurrentTransition)=...
            FitResultsIndiv.ManualTransitions(CurrentTransition)-2;
        FitResultsIndiv.Approved=0;
        
        
    %Rate, fine resolution
    elseif (ct~=0)&(cc=='s')
        FitResultsIndiv.SingleRate=...
            FitResultsIndiv.SingleRate+200*ScaleFactor;
        FitResultsIndiv.ManualRates=FitResultsIndiv.SingleRate...
            * FitResultsIndiv.NActive;
        FitResultsIndiv.Approved=0;
    elseif (ct~=0)&(cc=='x')&(FitResultsIndiv.SingleRate>0)
        FitResultsIndiv.SingleRate=...
            FitResultsIndiv.SingleRate-200*ScaleFactor;
        FitResultsIndiv.ManualRates=FitResultsIndiv.SingleRate...
            * FitResultsIndiv.NActive;
        FitResultsIndiv.Approved=0;
        
    %Rate, coarse resolution
    elseif (ct~=0)&(cc=='S')
        FitResultsIndiv.SingleRate=...
            FitResultsIndiv.SingleRate+2000*ScaleFactor;
        FitResultsIndiv.ManualRates=FitResultsIndiv.SingleRate...
            * FitResultsIndiv.NActive;
        FitResultsIndiv.Approved=0;
    elseif (ct~=0)&(cc=='X')&(FitResultsIndiv.SingleRate>0)
        FitResultsIndiv.SingleRate=...
            FitResultsIndiv.SingleRate-2000*ScaleFactor;
        FitResultsIndiv.ManualRates=FitResultsIndiv.SingleRate...
            * FitResultsIndiv.NActive;
        FitResultsIndiv.Approved=0;
        
    %NActive
    elseif (ct~=0)&(cc=='d')&(FitResultsIndiv.NActive(CurrentTransition)<2)
        FitResultsIndiv.NActive(CurrentTransition)=...
            FitResultsIndiv.NActive(CurrentTransition)+1;
        FitResultsIndiv.ManualRates=FitResultsIndiv.SingleRate...
            * FitResultsIndiv.NActive;
        FitResultsIndiv.Approved=0;
    elseif (ct~=0)&(cc=='c')&(FitResultsIndiv.NActive(CurrentTransition)>0)
        FitResultsIndiv.NActive(CurrentTransition)=...
            FitResultsIndiv.NActive(CurrentTransition)-1;
        FitResultsIndiv.ManualRates=FitResultsIndiv.SingleRate...
            * FitResultsIndiv.NActive;
        FitResultsIndiv.Approved=0;
    
        
    %Do the fit
    elseif (ct~=0)&(cc=='f')
        x0=[FitResultsIndiv.ManualTransitions,...
            FitResultsIndiv.SingleRate];
        %Create lower and upper bounds for the rates
        x0Lower=zeros(size(x0));
        x0Upper=inf(size(x0));
        
        options = optimset('MaxFunEvals',2000*length(x0),'MaxIter',5000);
        
        [xFit,resnorm,residual,exitflag,output,lambda,jacobian]=...
            lsqnonlin(@(x) lsqnonlinFitIndividualFluorescenceCurveConstrained(ElapsedTime(CompiledParticles(ParticlesNC{nc-12}(i)).Frame)-...
            ElapsedTime(eval(['nc',num2str(nc)])),...
            CompiledParticles(ParticlesNC{nc-12}(i)).Fluo,Delay,...
            FitResultsIndiv.nSteps,FitResultsIndiv.NActive,x),...
            x0,x0Lower,x0Upper,options);

        FitResultsIndiv.FittedTransitions=xFit(1:FitResultsIndiv.nSteps);
        FitResultsIndiv.FittedRates=xFit(end)*FitResultsIndiv.NActive;
        
        FitResultsIndiv.Approved=1;
        
        
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
        if FitResultsIndiv.Approved==0
            FitResultsIndiv.Approved=1;
        elseif FitResultsIndiv.Approved==1
            FitResultsIndiv.Approved=-1;
        elseif FitResultsIndiv.Approved==-1
            FitResultsIndiv.Approved=0;
        end

    
    
    %Save results
    elseif (ct~=0)&(cc=='t')
        save([DropboxFolder,filesep,Prefix,filesep,'IndividualFitsConstrained.mat'],...
            'FitResultsIndiv')
            display('IndividualFitsConstrained.mat saved')
     
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

            
save([DropboxFolder,filesep,Prefix,filesep,'IndividualFitsConstrained.mat'],...
    'FitResultsIndiv')
    display('IndividualFitsConstrained.mat saved')
