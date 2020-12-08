function particle = add3DTrackingInfo(particle, SpotsCh, trackingOptions, kalmanOptions)

    particle.xPos3D = NaN(size(particle.xPos));
    particle.yPos3D = NaN(size(particle.xPos));
    particle.zPos3D = NaN(size(particle.xPos));
    particle.zPosDetrended3D = NaN(size(particle.xPos));            

    for f = 1:length(particle.Frame)
        spot = SpotsCh(particle.Frame(f)).Fits(particle.Index(f));
        posAdjustment = 0; 
%         if isfield(spot, 'updated3DOffset')
%           if spot.updated3DOffset == 1
%             posAdjustment = -.5;
%           end
%         end
        particle.xPos3D(f) = spot.GaussPos3D(2)+posAdjustment;
        particle.yPos3D(f) = spot.GaussPos3D(1)+posAdjustment;
        particle.zPos3D(f) = spot.GaussPos3D(3)+posAdjustment;
        particle.zPosDetrended3D(f) = particle.zPos3D(f)-trackingOptions.zPosStage(particle.Frame(f));
%                 FluoVec3D(f) = spot.FixedAreaIntensity(brightestZIndex);
    end

    % initialize dummy structure
    tempStruct = struct;
    tempStruct.Frame = particle.Frame;
    tempStruct.MeasurementVec = [particle.xPos3D particle.yPos3D...
            particle.zPosDetrended3D];

%             trackingStruct.MeasurementVec(:,end+1) = FluoVec3D./trackingOptions.kalmanOptions.fluoFactor;

    % make particle path predictions
    particle = pathPrediction(particle, tempStruct, trackingOptions, kalmanOptions, 1);