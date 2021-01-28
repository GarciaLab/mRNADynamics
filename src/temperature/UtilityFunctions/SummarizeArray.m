function [BinCounts, BinEdges] = SummarizeArray(InputArray, BinNumber)
if ~exist('BinNumber', 'var')
    BinNumber = 50;
end
if isa(InputArray, 'integer')
    UniqueValues = unique(InputArray);
    BinCounts = [];
    if length(UniqueValues) <=  BinNumber
        for u=UniqueValues
            MatchingValues =  find(InputArray == u);
            BinCounts = [BinCounts, length(MatchingValues)];
            disp([num2str(u), ': ', num2str(length(MatchingValues))]);
            
        end
        BinEdges = UniqueValues;
    else
        BinWidth = ceil(length(UniqueValues)/BinNumber);
        BinEdges = min(UniqueValues):BinWidth:max(UniqueValues);
        
        if BinEdges(end) <= max(UniqueValues)
            BinEdges = [BinEdges, BinEdges(end)+BinWidth];
        end
        for bin_index = 1:length(BinEdges)-1
            IndexesMatchingBin = find((InputArray >= BinEdges(bin_index)) & (InputArray < BinEdges(bin_index+1)));
            NumberValuesMatchingBin = length(IndexesMatchingBin);
            disp([num2str(BinEdges(bin_index)),'-', num2str(BinEdges(bin_index+1)), ': ', num2str(NumberValuesMatchingBin)]);
            BinCounts = [BinCounts, NumberValuesMatchingBin];
        end
    end
    
    disp(['NaN: ', num2str(length(find(isnan(InputArray))))]);
    disp(['Total: ', num2str(length(InputArray))]);
elseif isa(InputArray, 'float')
    [BinnedValues, BinEdges] = discretize(InputArray, BinNumber);
    BinCounts = [];
    for bin_index = 1:BinNumber
        IndexesMatchingBin = find(BinnedValues == bin_index);
        NumberValuesMatchingBin = length(IndexesMatchingBin);
        disp([num2str(BinEdges(bin_index)),'-', num2str(BinEdges(bin_index+1)), ': ', num2str(NumberValuesMatchingBin)]);
        BinCounts = [BinCounts, NumberValuesMatchingBin];
    end
    disp(['NaN: ', num2str(length(find(isnan(InputArray))))]);
    disp(['Total: ', num2str(length(InputArray))]);
else
    error('Must be a float or integer array')
end
end