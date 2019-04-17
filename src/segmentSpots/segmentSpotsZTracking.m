function [neighborhood, Particles, Spots] = segmentSpotsZTracking(pixelSize, numFrames, Particles, Spots)


fields = fieldnames(Particles);

%%
changes = 1;

while changes ~= 0
    changes = 0;
    i = 1;
    waitbarFigure = waitbar(0, 'Finding z-columns');
    
    neighborhood = round(1300 / pixelSize); %empirically, the spot should be within ~1.3um of the same spot in other z-planes
    
    for framesIndex = 1:numFrames
        waitbar(framesIndex / numFrames, waitbarFigure)
        n = length(Particles([Particles.frame] == framesIndex));
        i = i + length(Particles([Particles.frame] == (framesIndex - 1)));
        
        for j = i:i + n - 1
            
            for k = j + 1:i + n - 1
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

%%
numFrames = length(Spots);
neighborhood = round(1300 / pixelSize); %empirically, the spot should be within ~1.3um of the same spot in other z-planes
waitbarFigure = waitbar(0, 'Finding z-columns');

if ~isempty([Spots.Fits])
    fields = fieldnames([Spots.Fits]);
    for currentFrame = 1:numFrames
        waitbar(currentFrame / numFrames, waitbarFigure)
        if ~isempty(Spots(currentFrame).Fits)
            changes = 1;
            while changes ~= 0
                changes = 0;

                lastSpot = length(Spots(currentFrame).Fits);

                for spotIndex = 1:lastSpot
                    for spotIndex2 = (spotIndex + 1):lastSpot
                        xDist = Spots(currentFrame).Fits(spotIndex).xFit(end) - Spots(currentFrame).Fits(spotIndex2).xFit(end);
                        yDist = Spots(currentFrame).Fits(spotIndex).yFit(end) - Spots(currentFrame).Fits(spotIndex2).yFit(end);
                        dist = sqrt(xDist^2 + yDist^2);

                        if dist < neighborhood && Spots(currentFrame).Fits(spotIndex).z(end) ~= Spots(currentFrame).Fits(spotIndex2).z(end)
                            for m = 1:numel(fields)
                                if strcmpi(fields{m}, 'r')
                                    Spots(currentFrame).Fits(spotIndex2).r = 1;
                                elseif strcmpi(fields{m}, 'frame')                                
                                    %do nothing
                                else
                                    Spots(currentFrame).Fits(spotIndex).(fields{m}) = [Spots(currentFrame).Fits(spotIndex).(fields{m}),Spots(currentFrame).Fits(spotIndex2).(fields{m}) ];
                                end
                                changes = changes + 1;
                            end
                        end
                    end
                end
                Spots(currentFrame).Fits = Spots(currentFrame).Fits([Spots(currentFrame).Fits.r] ~= 1);
            end
        end
    end
else
    fields = {};
end
    
    close(waitbarFigure)

end
