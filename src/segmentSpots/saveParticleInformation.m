function [Particles, fields] = saveParticleInformation(numFrames, all_frames, zSize)
  % Create a useful structure that can be fed into pipeline
  Particles = struct;
  n = 1;

  waitbarFigure = waitbar(0, 'Saving particle information');

  for frameIndex = 1:numFrames
      
    waitbar(frameIndex / numFrames, waitbarFigure)

    for zIndex = 1:zSize

      for spot = 1:length(all_frames{frameIndex, zIndex}) % spots within particular image

        if ~isempty(all_frames{frameIndex, zIndex}{spot})
            
          Particles(n).FixedAreaIntensity(1) = cell2mat(all_frames{frameIndex, zIndex}{spot}(1));
          Particles(n).xFit(1) = cell2mat(all_frames{frameIndex, zIndex}{spot}(2));
          Particles(n).yFit(1) = cell2mat(all_frames{frameIndex, zIndex}{spot}(3));
          Particles(n).Offset(1) = cell2mat(all_frames{frameIndex, zIndex}{spot}(4));
          Particles(n).Area{1} = cell2mat(all_frames{frameIndex, zIndex}{spot}(6));
          Particles(n).xFitWidth{1} = cell2mat(all_frames{frameIndex, zIndex}{spot}(7));
          Particles(n).yFitWidth{1} = cell2mat(all_frames{frameIndex, zIndex}{spot}(8));
          Particles(n).yDoG(1) = cell2mat(all_frames{frameIndex, zIndex}{spot}(9));
          Particles(n).xDoG(1) = cell2mat(all_frames{frameIndex, zIndex}{spot}(10));
          Particles(n).GaussianIntensity(1) = cell2mat(all_frames{frameIndex, zIndex}{spot}(11));
          Particles(n).CentralIntensity(1) = cell2mat(all_frames{frameIndex, zIndex}{spot}(12));
          Particles(n).DOGIntensity(1) = cell2mat(all_frames{frameIndex, zIndex}{spot}(13));
          Particles(n).SisterDistance(1) = cell2mat(all_frames{frameIndex, zIndex}{spot}(17));
          Particles(n).ConfidenceIntervals{1} = cell2mat(all_frames{frameIndex, zIndex}{spot}(19));
          Particles(n).gaussParams = all_frames{frameIndex,zIndex}{spot}(22);
          Particles(n).intArea = cell2mat(all_frames{frameIndex, zIndex}{spot}(23));
          Particles(n).intArea = Particles(n).intArea(1); %this should be the same for all z slices.
          Particles(n).cylIntensity(1) = cell2mat(all_frames{frameIndex, zIndex}{spot}(24));
          Particles(n).z(1) = zIndex;
          Particles(n).discardThis = 0;
          Particles(n).frame(1) = frameIndex;
          Particles(n).r = 0;
          n = n + 1;
          
        end 

      end 

    end 

  end 

  close(waitbarFigure);

  fields = fieldnames(Particles);
  
end 
