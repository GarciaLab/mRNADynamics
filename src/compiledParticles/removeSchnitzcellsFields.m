function schnitzcells = removeSchnitzcellsFields(schnitzcells)

    fieldsToRm = {'P', 'E', 'D'};
    for field = 1:length(fieldsToRm)
        if isfield(schnitzcells,fieldsToRm{field})
            schnitzcells = rmfield(schnitzcells, fieldsToRm{field});
        end
    end
    

end