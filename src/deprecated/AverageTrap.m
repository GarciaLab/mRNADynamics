function AverageTrap(DataType)

[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;
prompt='How many data sets are we analyzing today, Sir? \n';
NumberOfFiles=input(prompt);
disp('Very good, Sir')
pause(.5)
disp('Please select one of the files to analyze.')
FolderTemp=uigetdir(DropboxFolder,'Please select one of the files to analyze.');
disp('Excellent choice, Sir')
pause(.5)
Dashes=strfind(FolderTemp,'\');
Prefix=FolderTemp((Dashes(end)+1):end);
Dashes2 = strfind(Prefix, '-');
                                    
%Load the Fit Results                                  
load([DropboxFolder,filesep,Prefix,'\FitResultsMultiple.mat']);
TotalParticles=0;

%Combine the Fit Results into one concatenation
for dsn=1:length(FitResultsMultiple); %Data Set number = dsn
    if dsn >1
        FitResultsMultipleCombined=cat(1,FitResultsMultiple{1,1},FitResultsMultiple{1,dsn});
    end
end


FRC=FitResultsMultipleCombined; %Shortened the name of FitResultsMultipleCombined
FRCShape=size(FRC); 
NumberOfNCs=FRCShape(2);
MeanTraps=struct('Approved',1,'TimeStart',0,'TimeEnd',0,'RateFit',0,'SDTimeStart',0,'SDTimeEnd',0,'SDRateFit',0);
MeanTraps=repmat(MeanTraps,41,NumberOfNCs);
ParticlesPerBin=zeros(41,NumberOfNCs);
for nc=1:NumberOfNCs;
        for pn=1:length(FRC); %Particle Number = pn 
            particle=FRC(pn,nc); %particle is now referring to one specific particle
            if (particle.Approved==1);
                TotalParticles=TotalParticles+1;
                for binnum=1:41;
                    if binnum==particle.APBin;
                        if ((MeanTraps(binnum,nc).RateFit==0)... 
                                && (MeanTraps(binnum,nc).TimeStart==0) && (MeanTraps(binnum,nc).TimeEnd==0));
                            MeanTraps(binnum,nc).RateFit=particle.RateFit;
                            MeanTraps(binnum,nc).TimeStart=particle.TimeStart;
                            MeanTraps(binnum,nc).TimeEnd=particle.TimeEnd;
                            MeanTraps(binnum,nc).TimeEnd=particle.TimeEnd;
                            MeanTraps(binnum,nc).SDTimeStart=(particle.TimeStart)^2;
                            MeanTraps(binnum,nc).SDTimeEnd=(particle.TimeEnd)^2;
                            MeanTraps(binnum,nc).SDRateFit=(particle.RateFit)^2;
                            ParticlesPerBin(binnum,nc)=ParticlesPerBin(binnum,nc)+1;
                        else 
                            MeanTraps(binnum,nc).RateFit=particle.RateFit+MeanTraps(binnum,nc).RateFit;
                            MeanTraps(binnum,nc).TimeStart=particle.TimeStart+MeanTraps(binnum,nc).TimeStart;
                            MeanTraps(binnum,nc).TimeEnd=particle.TimeEnd+MeanTraps(binnum,nc).TimeEnd;
                            MeanTraps(binnum,nc).SDRateFit=(particle.RateFit)^2+MeanTraps(binnum,nc).SDRateFit;
                            MeanTraps(binnum,nc).SDTimeStart=(particle.TimeStart)^2+MeanTraps(binnum,nc).SDTimeStart;
                            MeanTraps(binnum,nc).SDTimeEnd=(particle.TimeEnd)^2+MeanTraps(binnum,nc).SDTimeEnd;
                            ParticlesPerBin(binnum,nc)=ParticlesPerBin(binnum,nc)+1;
                        end
                    end
                end
            end
        end
end
for binnum=1:41;
    for nc=1:NumberOfNCs;
        if ParticlesPerBin(binnum,nc)>0;
            MeanTraps(binnum,nc).RateFit=MeanTraps(binnum,nc).RateFit/ParticlesPerBin(binnum,nc);
            MeanTraps(binnum,nc).TimeStart=MeanTraps(binnum,nc).TimeStart/ParticlesPerBin(binnum,nc);
            MeanTraps(binnum,nc).TimeEnd=MeanTraps(binnum,nc).TimeEnd/ParticlesPerBin(binnum,nc);
            MeanTraps(binnum,nc).SDRateFit=sqrt(MeanTraps(binnum,nc).SDRateFit/ParticlesPerBin(binnum,nc)-(MeanTraps(binnum,nc).RateFit)^2);
            MeanTraps(binnum,nc).SDTimeStart=sqrt(MeanTraps(binnum,nc).SDTimeStart/ParticlesPerBin(binnum,nc)-(MeanTraps(binnum,nc).TimeStart)^2);
            MeanTraps(binnum,nc).SDTimeEnd=sqrt(MeanTraps(binnum,nc).SDTimeEnd/ParticlesPerBin(binnum,nc)-(MeanTraps(binnum,nc).TimeEnd)^2);
        end
    end
end

%Save the Mean Trapezoids
FitResults=MeanTraps;
save([DropboxFolder,filesep,Prefix,filesep,'MeanFits.mat'],...
    'FitResults')
save([DropboxFolder,filesep,Prefix,filesep,'FitResultsMultiple.mat'],...
    'FitResultsMultiple')
disp('I''ve taken the liberty of saving the results for you in the initial folder')
pause(1)
for i=1:(NumberOfFiles-1)
    disp('Please select another data set where I should save the results')
    Additional=uigetdir(DropboxFolder,'Please select another data set where I should save the results');
    save([Additional,filesep,'MeanFits.mat'],...
    'FitResults')
    save([Additional,filesep,'FitResultsMultiple.mat'],...
    'FitResultsMultiple')
    disp('Thank you, Sir')
end
TotalParticles

%s
    

