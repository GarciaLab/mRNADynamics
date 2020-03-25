function [EllipsePosAP, APAngle,EllipsePosDV]...
    = convertToFractionalEmbryoLength(Prefix)

thisExperiment = liveExperiment(Prefix);
resultsFolder = thisExperiment.resultsFolder;
load([resultsFolder,filesep,'APDetection.mat'])
Ellipses = getEllipses(thisExperiment);
correctDV = exist([resultsFolder, filesep,'DV',filesep,'DV_correction.mat'], 'file');
if correctDV
    load([DropboxFolder,filesep,Prefix,filesep,'DV',filesep,'DV_correction.mat'], 'DV_correction');
end

%Angle between the x-axis and the AP-axis
if exist('coordPZoom', 'var')
    APAngle=atan2((coordPZoom(2)-coordAZoom(2)),...
        (coordPZoom(1)-coordAZoom(1)));
else error('coordPZoom not defined. Was AddParticlePosition.m run?'); end

APLength=sqrt((coordPZoom(2)-coordAZoom(2))^2+(coordPZoom(1)-coordAZoom(1))^2);
DVLength = APLength/2;

%The information in Ellipses is
%(x, y, a, b, theta, maxcontourvalue, time, particle_id)
for i=1:length(Ellipses)
    for j=1:size(Ellipses{i},1)
        
        %Angle between the x-axis and the nucleus using the A position as a
        %zero
        
        nucleusAngles_AP=atan2((Ellipses{i}(j,2)-coordAZoom(2)),...
            (Ellipses{i}(j,1)-coordAZoom(1)));
        
        %Distance between the points and the A point
        Distances=sqrt((coordAZoom(2)-Ellipses{i}(j,2)).^2 ...
            + (coordAZoom(1)-Ellipses{i}(j,1)).^2);
        
        APPositions=Distances.*cos(nucleusAngles_AP-APAngle);
        EllipsePosAP{i}(j)=APPositions/APLength;
        
        if DVExperiment && correctDV
            DVPositions=Distances.*sin(nucleusAngles_AP-APAngle);
            EllipsePosDV{i}(j)=abs(DVPositions-DV_correction)/DVLength;
        else
            DVPositions=Distances.*sin(nucleusAngles_AP-APAngle);
            EllipsePosDV{i}(j)=DVPositions/DVLength;
        end
        
    end
end



end