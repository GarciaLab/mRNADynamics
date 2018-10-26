function Particles = findZColumns(Particles, particleFields, initialFrame, numFrames, neighborhoodZ)
  changes = 1;

  while changes ~= 0
    changes = 0;
    i = 1;
    zColumnsWaitbar = waitbar(0, 'Finding z-columns');

    for frameIndex = initialFrame:numFrames
      waitbar(frameIndex / (numFrames - initialFrame), zColumnsWaitbar)
      l = length(Particles([Particles.frame] == frameIndex));
      i = i + length(Particles([Particles.frame] == (frameIndex - 1)));

      for j = i:i + l - 1

        for k = j + 1:i + l - 1
          dist = sqrt((Particles(j).xFit(end) - Particles(k).xFit(end))^2 + (Particles(j).yFit(end) - Particles(k).yFit(end))^2);

          if dist < neighborhoodZ && Particles(j).z(end) ~= Particles(k).z(end)

            for fieldIndex = 1:numel(particleFields) - 2 % do not include fields 'r' or 'frame'
              Particles(j).(particleFields{fieldIndex}) = [Particles(j).(particleFields{fieldIndex}),...
                Particles(k).(particleFields{fieldIndex})];
            end

            Particles(k).r = 1;
            changes = changes + 1;
          end

        end

      end

    end

    Particles = Particles([Particles.r] ~= 1);
    close(zColumnsWaitbar);

  end

end
