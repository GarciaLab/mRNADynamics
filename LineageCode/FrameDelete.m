function [ Particles] = FrameDelete(Prefix,Start, End)

% Marks all particles between Start and End as disapproved


[~,~,DropboxFolder,~,~]=DetermineLocalFolders(Prefix);
Particles=struct;
try
    dummy1=load([DropboxFolder,filesep,Prefix,filesep,'Particles.mat']);
    Particles=dummy1.Particles;
catch
    warning('No Particles.mat found. No changes will be made to the Particles structure.');
end

for i=1:size(Particles,2)
    Particles(i).FrameApproved(Particles(i).Frame<End && Particles(i).Frame>Start)=false;
end
save([DataFolder,filesep,'Particles.mat'],'Particles','fad','fad2','Threshold1','Threshold2')
end

