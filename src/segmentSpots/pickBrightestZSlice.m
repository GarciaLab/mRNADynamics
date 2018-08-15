function [Particles, falsePositives] = pickBrightestZSlice(Particles, fields, trackSpots, numShadows)
  falsePositives = 0;
  %pick the brightest z-slice

  for i = 1:length(Particles)
    [~, max_index] = max(Particles(i).CentralIntensity);

    if trackSpots

      for j = 1:numel(fields) - 2 % do not include fields 'r' or 'frame'
        Particles(i).(fields{j}) = Particles(i).(fields{j})(max_index);
      end 

    else 
      Particles(i).brightestZ = Particles(i).z(max_index);

      if numShadows == 1

        if length(Particles(i).z) <= 1
          Particles(i).discardThis = 1;
          Particles(i).noIntensityAnalysis = 1;
          falsePositives = falsePositives + 1;
        elseif Particles(i).brightestZ == Particles(i).z(end)

          if Particles(i).z(max_index - 1) ~= Particles(i).brightestZ - 1
            Particles(i).discardThis = 1;
            Particles(i).noIntensityAnalysis = 1;
            falsePositives = falsePositives + 1;
          end 

        elseif Particles(i).brightestZ == Particles(i).z(1)

          if Particles(i).z(max_index + 1) ~= Particles(i).brightestZ + 1
            Particles(i).discardThis = 1;
            Particles(i).noIntensityAnalysis = 1;
            falsePositives = falsePositives + 1;
          end 

        elseif Particles(i).z(max_index - 1) ~= Particles(i).brightestZ - 1 ...
          && Particles(i).z(max_index + 1) ~= Particles(i).brightestZ + 1
          Particles(i).discardThis = 1;
          Particles(i).noIntensityAnalysis = 1;
          falsePositives = falsePositives + 1;
        end 

      elseif numShadows == 2

        if Particles(i).brightestZ == Particles(i).z(end) || ...
          Particles(i).brightestZ == Particles(i).z(1)
          Particles(i).discardThis = 1;
          Particles(i).noIntensityAnalysis = 1;
          falsePositives = falsePositives + 1;
        elseif Particles(i).z(max_index - 1) ~= Particles(i).brightestZ - 1 ...
          || Particles(i).z(max_index + 1) ~= Particles(i).brightestZ + 1
          Particles(i).discardThis = 1;
          Particles(i).noIntensityAnalysis = 1;
          falsePositives = falsePositives + 1;
        end 

      end 

    end 

  end 

end 
