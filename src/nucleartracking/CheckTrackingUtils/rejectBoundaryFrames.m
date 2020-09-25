function schnitzcells =...
    rejectBoundaryFrames(schnitzcells, padding)

if nargin == 1
    padding = 2;
end
Nslices = size(schnitzcells(1).Fluo, 2);
for i=1:length(schnitzcells)
    MaxFluos = max(schnitzcells(i).Fluo, [], 2);
    schnitzcells(i).FrameApproved(MaxFluos == 0) = 0;
    for j = 1:length(schnitzcells(i).frames)
        if isempty(find(schnitzcells(i).Fluo(j,(2+padding):(Nslices-1-padding)) == MaxFluos(j), 1))
            schnitzcells(i).FrameApproved(j) = 0;
        end
    end



end
end