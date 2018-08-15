function Spots = createSpotsStructure(Particles, numFrames, channelIndex)
  Spots{channelIndex} = [];
  fields = fieldnames(Particles);
  num_fields = length(fields);

  for framesIndex = 1:numFrames
    frames = find([Particles.frame] == framesIndex);

    if ~ isempty(frames)

      for j = frames(1):frames(end)

        if ~ Particles(j).discardThis
          Spots{channelIndex}(framesIndex).Fits(j - frames(1) + 1) = Particles(j);
        end 

        %Sometimes, all spots are discarded in a frame. In that
        %case, create an empty Spots entry in that frame.

        if length(Spots{channelIndex}) < framesIndex

          for l = 1:num_fields
            Spots{channelIndex}(framesIndex).Fits.(fields{l}) = [];
          end 

        end 

      end 

    else 

      for l = 1:num_fields
        Spots{channelIndex}(framesIndex).Fits.(fields{l}) = [];
      end 

    end 

  end 

  %Clean up Spots to remove empty rows
  Dots{channelIndex} = struct('Fits', []);

  for spotsIndex = 1:length(Spots{channelIndex})
    Dots{channelIndex}(spotsIndex).Fits = [];

    for j = 1:length(Spots{channelIndex}(spotsIndex).Fits)

      if j ~= 1

        if ~ isempty(Spots{channelIndex}(spotsIndex).Fits(j).z)...
          && ~ isequal(Spots{channelIndex}(spotsIndex).Fits(j).CentralIntensity, ...
          Spots{channelIndex}(spotsIndex).Fits(j - 1).CentralIntensity)
          
          Dots{channelIndex}(spotsIndex).Fits = [Dots{channelIndex}(spotsIndex).Fits, Spots{channelIndex}(spotsIndex).Fits(j)];
        end 

      else 

        if ~ isempty(Spots{channelIndex}(spotsIndex).Fits(j).z)
          Dots{channelIndex}(spotsIndex).Fits = [Dots{channelIndex}(spotsIndex).Fits, Spots{channelIndex}(spotsIndex).Fits(j)];
        end 

      end 

    end 

  end 

  for dotsIndex = 1:length(Dots{channelIndex})

    if isstruct(Dots{channelIndex}(dotsIndex).Fits)
      Spots{channelIndex}(dotsIndex).Fits = rmfield(Dots{channelIndex}(dotsIndex).Fits, 'r');
      Spots{channelIndex}(dotsIndex).Fits = rmfield(Spots{channelIndex}(dotsIndex).Fits, 'discardThis');
    else 
      Spots{channelIndex}(dotsIndex).Fits = [];
    end 

  end 

end 
