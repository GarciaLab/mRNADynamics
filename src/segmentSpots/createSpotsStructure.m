function Spots = createSpotsStructure(Particles, numFrames, initialFrame)
  Spots = [];
  fields = fieldnames(Particles);
  num_fields = length(fields);

  for framesIndex = initialFrame:numFrames
    frames = find([Particles.frame] == framesIndex);

    if ~ isempty(frames)

      for j = frames(1):frames(end)

        if ~ Particles(j).discardThis
          Spots(framesIndex).Fits(j - frames(1) + 1) = Particles(j);
        end 

        %Sometimes, all spots are discarded in a frame. In that
        %case, create an empty Spots entry in that frame.

        if length(Spots) < framesIndex

          for l = 1:num_fields
            Spots(framesIndex).Fits.(fields{l}) = [];
          end 

        end 

      end 

    else 

      for l = 1:num_fields
        Spots(framesIndex).Fits.(fields{l}) = [];
      end 

    end 

  end 

  %Clean up Spots to remove empty rows
  Dots = struct('Fits', []);

  for spotsIndex = 1:length(Spots)
    Dots(spotsIndex).Fits = [];

    for j = 1:length(Spots(spotsIndex).Fits)

      if j ~= 1

        if ~ isempty(Spots(spotsIndex).Fits(j).z)...
          && ~ isequal(Spots(spotsIndex).Fits(j).CentralIntensity, ...
          Spots(spotsIndex).Fits(j - 1).CentralIntensity)
          
          Dots(spotsIndex).Fits = [Dots(spotsIndex).Fits, Spots(spotsIndex).Fits(j)];
        end 

      else 

        if ~ isempty(Spots(spotsIndex).Fits(j).z)
          Dots(spotsIndex).Fits = [Dots(spotsIndex).Fits, Spots(spotsIndex).Fits(j)];
        end 

      end 

    end 

  end 

  for dotsIndex = 1:length(Dots)

    if isstruct(Dots(dotsIndex).Fits)
      Spots(dotsIndex).Fits = rmfield(Dots(dotsIndex).Fits, 'r');
      Spots(dotsIndex).Fits = rmfield(Spots(dotsIndex).Fits, 'discardThis');
    else 
      Spots(dotsIndex).Fits = [];
    end 

  end 

end 
