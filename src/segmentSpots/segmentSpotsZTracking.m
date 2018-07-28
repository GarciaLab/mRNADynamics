function [neighborhood, Particles] = segmentSpotsZTracking(pixelSize, numFrames, Particles, fields)
  changes = 1;

  while changes ~= 0
    changes = 0;
    i = 1;
    waitbarFigure = waitbar(0, 'Finding z-columns');

    neighborhood = round(1300 / pixelSize);
    for framesIndex = 1:numFrames
      waitbar(framesIndex / numFrames, waitbarFigure)
      l = length(Particles([Particles.frame] == framesIndex));
      i = i + length(Particles([Particles.frame] == (framesIndex - 1)));

      for j = i:i + l - 1

        for k = j + 1:i + l - 1
          dist = sqrt((Particles(j).xFit(end) - Particles(k).xFit(end))^2 + (Particles(j).yFit(end) - Particles(k).yFit(end))^2);

          if dist < neighborhood && Particles(j).z(end) ~= Particles(k).z(end)

            for m = 1:numel(fields) - 2 % do not include fields 'r' or 'frame'
              Particles(j).(fields{m}) = [Particles(j).(fields{m}), Particles(k).(fields{m})];
            end 

            Particles(k).r = 1;
            changes = changes + 1;
          end 

        end 

      end 

    end 

    Particles = Particles([Particles.r] ~= 1);
    close(waitbarFigure)
  end 

end 
