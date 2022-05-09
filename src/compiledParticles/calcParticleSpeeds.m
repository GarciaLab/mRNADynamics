function calcParticleSpeeds(NChannels, Particles, ...
    Spots, ElapsedTime, schnitzcells, Ellipses)
%
%CALCPARTICLESPEEDS Summary of this function goes here
%   TO DO: produce output and plots. 

    if NChannels > 1
%         disp('Speeds could not be calculated for your data. Please contact Emma.')
        % implement crude? initial slope and time on calculation
    else
        numberOfParticles = size(Particles{:},2);

        currentChannel = 1;
        for currentParticle = 1:numberOfParticles
            if ~isempty(Particles{currentChannel}(currentParticle).Nucleus)
                % getting the frames of the particle of interest --------------------------
                plotTraceSettings = PlotTraceSettings();
                [frames,~,~,~,~,~,~,~,~,~] =...
                    GetParticleTrace(currentParticle,...
                    Particles{currentChannel}, Spots{currentChannel}, plotTraceSettings, false);

                % getting tau (average time increment [minutes])--------------------------
                currentTimeArray = ElapsedTime(frame); % Units to seconds
                %AR 1/2/2019- what was this originally for?

                % getting x,y,z positions -------------------------------------------------
                xArray = zeros(1,length(frames));
                yArray = zeros(1,length(frames));
                zArray = zeros(1,length(frames));

                % coordinates with the Nucleus Center as the origin (NO).
                xArrayNO = zeros(1,length(frames));
                yArrayNO = zeros(1,length(frames));

                frameCounter = 0;
                for currentFrame = frames
                    frameCounter = frameCounter + 1;

                    % getting current particle index
                    currentParticleIndex=...
                        Particles{currentChannel}(currentParticle).Index(frameCounter);

                    % getting x and y positions of the brightest z slice
                    [x,y,z]=getSpotsXYZ(Spots{currentChannel}(currentFrame));

                    % storing x, y, and z positions in pixels relative to the image
                    % frame
                    xArray(frameCounter) = x(currentParticleIndex);
                    yArray(frameCounter) = y(currentParticleIndex);
                    zArray(frameCounter) = z(currentParticleIndex);

                    %getting current nucleus and schnitz
                    schnitzIndex=find(schnitzcells(Particles{currentChannel}(currentParticle).Nucleus).frames==currentFrame);
                    nucleusIndex=schnitzcells(Particles{currentChannel}(currentParticle).Nucleus).cellno(schnitzIndex);
                    nucleusCenterCoordinates = [Ellipses{currentFrame}(nucleusIndex,1:2)];

                    % storing x and y positions with the assigned nucleus as the origin.
                    xArrayNO(frameCounter) = x(currentParticleIndex) - nucleusCenterCoordinates(1);
                    yArrayNO(frameCounter) = y(currentParticleIndex) - nucleusCenterCoordinates(2);
                end

                % speed calculation and distance relative to the center of the nucleus
                dz = diff(zArray); %AR 1/2/2019- why isn't the z coordinate also NO?
                dxNO = diff(xArrayNO);
                dyNO = diff(yArrayNO);
                drNO = sqrt(dxNO.^2 + dyNO.^2 + dz.^2);
                rNO = sqrt(xArrayNO.^2 + yArrayNO.^2 + zArray.^2);

                %plot graph without showing and save it in a folder with the yy
            end
        end
    end
end