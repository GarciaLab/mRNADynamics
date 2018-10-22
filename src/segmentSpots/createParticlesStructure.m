% Creates a useful structure that can be fed into pipeline
function [Particles, particleFields] = createParticlesStructure(initialFrame, numFrames, allFrames, zSize, snippet_size, useIntegralCenter)
  Particles = struct('FixedAreaIntensity', [], 'xFit', [], 'yFit', [], 'Offset', [], ...
  'GaussianIntensity', [], 'CentralIntensity', [], 'xDoG', [], 'yDoG', [], ...
  'DOGIntensity', [], 'ConfidenceIntervals', {}, ...
  'gaussParams', [], 'z', [], 'discardThis', [], 'frame', [], 'r', []);
  particleFields = fieldnames(Particles);

  particlesIndex = 1;

  particleInformationWaitbar = waitbar(0, 'Saving particle information');

  for frameIndex = initialFrame:numFrames
    waitbar(frameIndex / (numFrames - initialFrame), particleInformationWaitbar);

    for zIndex = 1:zSize

      for spotIndex = 1:length(allFrames{frameIndex, zIndex})%Spots{channelIndex} within particular image
        
        spot = allFrames{frameIndex, zIndex}{spotIndex};

        if ~ isempty(spot)
          Particles(particlesIndex).FixedAreaIntensity(1) = cell2mat(spot(1));
          Particles(particlesIndex).xFit(1) = cell2mat(spot(2));
          Particles(particlesIndex).yFit(1) = cell2mat(spot(3));
          Particles(particlesIndex).Offset(1) = cell2mat(spot(4));
          Particles(particlesIndex).yDoG(1) = cell2mat(spot(9));
          Particles(particlesIndex).xDoG(1) = cell2mat(spot(10));
          Particles(particlesIndex).GaussianIntensity(1) = cell2mat(spot(11));
          Particles(particlesIndex).CentralIntensity(1) = cell2mat(spot(12));
          Particles(particlesIndex).DOGIntensity(1) = cell2mat(spot(13));
          Particles(particlesIndex).ConfidenceIntervals{1} = cell2mat(spot(19));
          Particles(particlesIndex).gaussParams = spot(22);
          Particles(particlesIndex).z(1) = zIndex;
          Particles(particlesIndex).discardThis = 0;
          Particles(particlesIndex).frame(1) = frameIndex;
          Particles(particlesIndex).r = 0;
          Particles(particlesIndex).intArea = cell2mat(spot(23));
          Particles(particlesIndex).cylIntensity = cell2mat(spot(24));
          Particles(particlesIndex).IntegralZ = useIntegralCenter;
          Particles(particlesIndex).snippet_size = snippet_size;
          particlesIndex = particlesIndex + 1;
        end

      end

    end

  end
  
  close(particleInformationWaitbar);

end
