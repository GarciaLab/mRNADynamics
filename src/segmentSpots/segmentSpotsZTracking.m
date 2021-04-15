function Spots = segmentSpotsZTracking(pixelSize, Spots, varargin)

spotsTracking = Spots;

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end

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

                nSpots = length(Spots(currentFrame).Fits);

                for spotIndex = 1:nSpots
                    for spotIndex2 = (spotIndex + 1):nSpots
                           xDist = Spots(currentFrame).Fits(spotIndex).xFit(end) - Spots(currentFrame).Fits(spotIndex2).xFit(end);
                        yDist = Spots(currentFrame).Fits(spotIndex).yFit(end) - Spots(currentFrame).Fits(spotIndex2).yFit(end);
%                          xDist = Spots(currentFrame).Fits(spotIndex).xDoG(end) - Spots(currentFrame).Fits(spotIndex2).xDoG(end);
%                         yDist = Spots(currentFrame).Fits(spotIndex).yDoG(end) - Spots(currentFrame).Fits(spotIndex2).yDoG(end);
                        
%                          dist = sqrt(double(xDist)^2 + double(yDist)^2);
                        dist = sqrt(xDist^2 + yDist^2);


                        if dist < neighborhood &&...
                                Spots(currentFrame).Fits(spotIndex).z(end) ~= Spots(currentFrame).Fits(spotIndex2).z(end)
                            for m = 1:numel(fields)
                                if strcmpi(fields{m}, 'r')
                                    Spots(currentFrame).Fits(spotIndex2).r = 1;
                                elseif strcmpi(fields{m}, 'frame') ||...
                                        strcmpi(fields{m}, 'IntegralZ') ||...
                                        strcmpi(fields{m}, 'intArea') ||...
                                        strcmpi(fields{m}, 'snippet_size')                     
                                    %do nothing
                                else
                                    Spots(currentFrame).Fits(spotIndex).(fields{m}) =...
                                        [Spots(currentFrame).Fits(spotIndex).(fields{m}),Spots(currentFrame).Fits(spotIndex2).(fields{m}) ];
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
