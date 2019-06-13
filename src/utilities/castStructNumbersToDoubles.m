function s = castStructNumbersToDoubles(s)

nCh = 1;
if iscell(s)
    nCh = length(s);
end

for ch = 1:nCh
    
    if iscell(s)
        sCh = s{ch};
    else
        sCh = s;
    end
    
    for frame = 1:length(sCh)
        spots = sCh(frame).Fits;
        if ~isempty(spots)
            fieldnames = fields(spots);
            numericFields = {};
            for field = 1:numel(fieldnames)
                if isnumeric([spots.(fieldnames{field})])
                    numericFields = [numericFields, fieldnames{field}];
                end
            end
            for spotInd = 1:length(spots)
                %there's one special case. 
                for k = 1:length(spots(spotInd).gaussParams)
                    spots(spotInd).gaussParams{k} = double(spots(spotInd).gaussParams{k});
                end
                for j = 1:length(numericFields)
                    spots(spotInd).(numericFields{j}) = double(spots(spotInd).(numericFields{j}));
                end
            end
        end
        
        sCh(frame).Fits = spots;
    end
    
    
    
    if iscell(s)
        s{ch} = sCh;
    else
        s = sCh;
    end
    
    
end

end