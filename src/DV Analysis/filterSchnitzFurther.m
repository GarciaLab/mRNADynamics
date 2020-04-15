function schnitzcells = filterSchnitzFurther(schnitzcells)

% Disapprove nuclei without
%fluos.

for s = 1:length(schnitzcells)
    
    if isempty(schnitzcells(s).Fluo)
        schnitzcells(s).Approved = false;
    end
    
end

end