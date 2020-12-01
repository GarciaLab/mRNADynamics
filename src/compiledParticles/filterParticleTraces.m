function particleFiltered = filterParticleTraces(particle, trackingOptions, SpotsCh)

    % remove problematic frames
    FrameFilter = particle.FrameApproved;
    particleFiltered = particle;    
    particleFiltered.FirstFrame = find(FrameFilter,1);
    particleFiltered.LastFrame = find(FrameFilter,1,'last');
    particleFiltered.Frame = particle.Frame(FrameFilter);
    particleFiltered.Index = particle.Index(FrameFilter);
    particleFiltered.Index = particle.Index(FrameFilter);
    particleFiltered.xPos = particle.xPos(FrameFilter);
    particleFiltered.yPos = particle.yPos(FrameFilter);
    
    if isfield(particle, 'zPosDetrended')
        particleFiltered.zPosDetrended = particle.zPosDetrended(FrameFilter);
    end
    if isfield(particle, 'zPos')
        particleFiltered.zPos = particle.zPos(FrameFilter);
    end
    if isfield(particle, 'xPos3D')
        particleFiltered.xPos3D = particle.xPos3D(FrameFilter);
        particleFiltered.yPos3D = particle.yPos3D(FrameFilter);
        if isfield(particleFiltered, 'zPosDetrended')
            particleFiltered.zPosDetrended3D = particle.zPosDetrended3D(FrameFilter);
        end
        if isfield(particle, 'zPos')
            particleFiltered.zPos3D = particle.zPos3D(FrameFilter);
        end
    end
    
    particleFiltered.FrameApproved = particle.FrameApproved(FrameFilter);
    particleFiltered.ManuallyReviewed = particle.ManuallyReviewed(FrameFilter);    
    particleFiltered.nucleusProbability = particle.nucleusProbability(FrameFilter);
    
    if any(particleFiltered.FrameApproved)
        % update path predictions    
        trackingStruct.MeasurementVec = [particleFiltered.xPos particleFiltered.yPos...
                                                  particleFiltered.zPosDetrended];

        trackingStruct.Frame = particleFiltered.Frame;
    %     FluoVec = [];
    %     for f = 1:length(trackingStruct.Frame)
    %         spot = Spots(trackingStruct.Frame(f)).Fits(Particles(OriginalParticle).Index(f));
    %         brightestZ = spot.brightestZ;
    %         brightestZIndex = spot.z == brightestZ;
    %         FluoVec(f) = spot.FixedAreaIntensity(brightestZIndex);
    %     end
    %     trackingStruct.MeasurementVec(:,end+1) = FluoVec./trackingOptions.kalmanOptions.fluoFactor;

        % make particle path predictions
        particleFiltered = pathPrediction(particleFiltered, trackingStruct, trackingOptions, trackingOptions.kalmanOptions,0);
%         particleFiltered.framesFull = particle.FirstFrame:particleFiltered.LastFrame;
        % add 3D info
        if trackingOptions.has3DInfo
            particleFiltered = add3DTrackingInfo(particleFiltered, SpotsCh, trackingOptions, trackingOptions.kalmanOptions);
        end
    end