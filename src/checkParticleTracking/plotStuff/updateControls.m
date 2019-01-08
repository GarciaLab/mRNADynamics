function updateControls(frame_num, z_num, particle_num, CurrentFrame, ...
    CurrentZ, CurrentParticle)
%UPDATECONTROLS Summary of this function goes here
%   Detailed explanation goes here

frame_num.Value = num2str(CurrentFrame);
z_num.Value = num2str(CurrentZ);
particle_num.Value = num2str(CurrentParticle);
end

