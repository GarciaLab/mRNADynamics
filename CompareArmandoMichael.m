function CompareArmandoMichael

%Compare the tracking of particles between Armando and Michael's code

%% Load both datasets

%Information about about folders
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;

ADataPrefix='2015-07-25-test_P2P_75uW_bi_test';
MDataPrefix='2015-07-25-P2P_75uW_bi';

%Load the raw particle traces
AData=load([DropboxFolder,filesep,ADataPrefix,filesep,'Particles.mat']);
MData=load([DropboxFolder,filesep,MDataPrefix,filesep,'Particles.mat']);

%Load the spots information
ASpots=load([DropboxFolder,filesep,ADataPrefix,filesep,'Spots.mat']);

%Load the CompiledParticles
ACompiled=load([DropboxFolder,filesep,ADataPrefix,filesep,'CompiledParticles.mat']);
MCompiled=load([DropboxFolder,filesep,MDataPrefix,filesep,'CompiledParticles.mat']);

%% Compare the means

figure(1)
plot(ACompiled.ElapsedTime,ACompiled.MeanVectorAP(:,11)','-k')
hold on
plot(MCompiled.ElapsedTime,MCompiled.MeanVectorAP(:,11)','-r')
hold off
legend('Armando','Michael')


%% Look at individual particles

%Compare a given particle in nc14

%First, find a good particle in Michael's data set
MParticle=50;

%Find the corresponding particle in Armando's code by looking at the
%associated nucleus
AParticle=find([ACompiled.CompiledParticles.Nucleus]==...
    MCompiled.CompiledParticles(MParticle).Nucleus);

%What original particle does this one correspond to?
AOriginalParticle=ACompiled.CompiledParticles(AParticle).OriginalParticle

%Find the missing frames between Armando and Michael's
figure(1)
plot(ACompiled.CompiledParticles(AParticle).Frame,...
    ACompiled.CompiledParticles(AParticle).Fluo,'.-k')
hold on
plot(MCompiled.CompiledParticles(MParticle).Frame,...
    MCompiled.CompiledParticles(MParticle).Fluo,'.-r')
hold off

%Is this a result of tracking or of compiling particles?
[Frame,AmpIntegral,AmpGaussian,Offset,...
    ErrorIntegral,ErrorGauss,optFit,FitType,noIntensityFlag]=...
    GetParticleTrace(AOriginalParticle,AData.Particles,ASpots.Spots);

figure(2)
plot(ACompiled.CompiledParticles(AParticle).Frame,...
    ACompiled.CompiledParticles(AParticle).Fluo,'.-k')
hold on
plot(Frame,AmpGaussian,'o-r')
hold off
legend('CompiledParticle','Particle')


