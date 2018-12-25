function [Particles, schnitzcells] = redoTracking(DataFolder, ...
    UseHistoneOverlay, FrameInfo, DropboxFolder, FilePrefix, schnitzcells, ...
    Particles, Threshold1, Threshold2, NChannels, CurrentChannel, numParticles)
%REDOTRACKING Summary of this function goes here
%   Detailed explanation goes here

Answer=lower(input('Are you sure you want to redo the tracking?  (y/n) ','s'));
if Answer=='y'
    warning('HG: Not clear that this feature will work with the multiple channels')

    %We need to save the data
    save([DataFolder,filesep,'FrameInfo.mat'],'FrameInfo')
    if UseHistoneOverlay
        save([DataFolder,filesep,'Particles.mat'],'Particles','Threshold1','Threshold2', '-v7.3')
        save([DropboxFolder,filesep,FilePrefix(1:end-1),filesep,FilePrefix(1:end-1),'_lin.mat'],'schnitzcells', '-v7.3')
    else
        save([DataFolder,filesep,'Particles.mat'],'Particles','Threshold1','Threshold2', '-v7.3')
    end
    disp('Particles saved.')
    if NChannels==1
        Particles=Particles{1};
    end

    [Particles,schnitzcells]=TrackmRNADynamics(FilePrefix(1:end-1),...
        Threshold1,Threshold2);
    if NChannels==1
        Particles={Particles};
    end
    %Check the FrameApproved field
    for i=1:numParticles
        if isempty(Particles{CurrentChannel}(i).FrameApproved)
            Particles{CurrentChannel}(i).FrameApproved=true(size(Particles{CurrentChannel}(i).Frame));
        end
    end
end
end

