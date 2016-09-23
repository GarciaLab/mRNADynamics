function CompareArmandoMichael

%Compare the tracking of particles between Armando and Michael's code

%Information about about folders
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;

ADataPrefix='2015-07-25-test_P2P_75uW_bi_test';
MDataPrefix='2015-07-25-P2P_75uW_bi';

%Load the raw particle traces
AData=load([DropboxFolder,filesep,ADataPrefix,filesep,'Particles.mat']);
MData=load([DropboxFolder,filesep,MDataPrefix,filesep,'Particles.mat']);


%Load the CompiledParticles
ACompiled=load([DropboxFolder,filesep,ADataPrefix,filesep,'CompiledParticles.mat']);
MCompiled=load([DropboxFolder,filesep,MDataPrefix,filesep,'CompiledParticles.mat']);
