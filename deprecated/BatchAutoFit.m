function BatchAutoFit(varargin)
% BatchAutoFit: Auto fit all profiles in an embryo using the multi-rate
% fitting algorithm.

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
DropboxFolder
%Instead just use folder I configured
% DropboxFolder = 'C:\Users\Mark\Dropbox\LivemRNAData';

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
[Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
nc9, nc10, nc11, nc12, nc13, nc14, CF] = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder)

%Distance from the first MS2 site to the end of the 3' UTR
if ~isempty(strfind(StemLoop,'Eve'))
    GeneLength=6.443;
elseif  ~isempty(strfind(StemLoop,'X1'))
    GeneLength=5.296;       %Distance from the first MS2 site to the end of the
                        %TUB3'UTR in kb.
else
%     error('The gene length has not been defined for this construct')
    GeneLength=6.443;
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
            FitResultsIndiv(i,j).FittedRates=[];
            FitResultsIndiv(i,j).AutoTransitions=[];
            FitResultsIndiv(i,j).AutoRates=[];
    end
end

%Load the fits if they exist
if exist([DropboxFolder,filesep,Prefix,filesep,'IndividualFits.mat'])
    load([DropboxFolder,filesep,Prefix,filesep,'IndividualFits.mat'])
end

nc=14;
RedoAllPast=0; % If true, redo any autofits that have been performed already

for i=1:length(ParticlesNC{nc-12})
    i
    if (i>RedoAllPast || ~isfield(FitResultsIndiv(i,nc-12),'AutoTransitions') ||...
            isempty(FitResultsIndiv(i,nc-12).AutoTransitions))...
            && length(CompiledParticles(ParticlesNC{nc-12}(i)).Fluo) > MinTimePoints
        disp('############### NEW PARTICLE ################')
        disp(['### ',num2str(i),' of ',num2str(length(ParticlesNC{nc-12}))])
        [AutoTransitions,AutoRates] = AutoFitFluorescenceCurveLinear(ElapsedTime(CompiledParticles(ParticlesNC{nc-12}(i)).Frame)-...
            ElapsedTime(eval(['nc',num2str(nc)])),...
            CompiledParticles(ParticlesNC{nc-12}(i)).Fluo,...
            CompiledParticles(ParticlesNC{nc-12}(i)).FluoError,...
            GeneLength, ElongationRate);
        FitResultsIndiv(i,nc-12).AutoTransitions=AutoTransitions;
        FitResultsIndiv(i,nc-12).AutoRates=AutoRates;


        %Save results
        save([DropboxFolder,filesep,Prefix,filesep,'IndividualFits.mat'],...
            'FitResultsIndiv')
            display('IndividualFits.mat saved')
    end
end

            
save([DropboxFolder,filesep,Prefix,filesep,'IndividualFits.mat'],...
    'FitResultsIndiv')
    display('IndividualFits.mat saved')
