function SubplotPositionInfo = GetSubplotPositioningParameters(NumSets, includeLegend, forTheoryPlots, hasTitle)
if ~exist('includeLegend', 'var')
    includeLegend = true;
end
if ~exist('forTheoryPlots', 'var')
    forTheoryPlots = false;
end

if ~exist('hasTitle', 'var')
    hasTitle = true;
end


SubplotPositionInfo = {};
if NumSets == 2
    SubplotPositionInfo.SubplotDims = [1, 2];
elseif NumSets == 3
    SubplotPositionInfo.SubplotDims = [1, 3];
elseif NumSets == 4
    SubplotPositionInfo.SubplotDims = [2, 2];
elseif NumSets == 9
    if forTheoryPlots
        SubplotPositionInfo.SubplotDims = [3, 3];% GM 4/29
    else
        SubplotPositionInfo.SubplotDims = [2, 5];
    end
elseif NumSets == 12
    SubplotPositionInfo.SubplotDims = [2, 6];
    if ~includeLegend
        SubplotPositionInfo.SubplotDims = [3, 4];
    end
elseif NumSets == 16
    SubplotPositionInfo.SubplotDims = [3, 6];
elseif NumSets == 17
    SubplotPositionInfo.SubplotDims = [3, 6];
elseif NumSets == 18
    SubplotPositionInfo.SubplotDims = [3, 6];
elseif NumSets == 19
    SubplotPositionInfo.SubplotDims = [4, 5];
elseif NumSets == 20
    SubplotPositionInfo.SubplotDims = [4, 5];
elseif NumSets >= 21 &  NumSets <= 24
    SubplotPositionInfo.SubplotDims = [4, 6];
else
    SubplotPositionInfo.SubplotDims = numSubplots(NumSets);
end

SubplotPositionInfo.SubFigDims = [0, 0];

if includeLegend
    SubplotPositionInfo.SubplotDims(2) = SubplotPositionInfo.SubplotDims(2)+1;
    SubplotPositionInfo.LegendSubplots = (1:SubplotPositionInfo.SubplotDims(1))*SubplotPositionInfo.SubplotDims(2);
   
    % SubFigDims(1)
    if SubplotPositionInfo.SubplotDims(2) < 4
        SubplotPositionInfo.SubFigDims(1) = 0.7;
    elseif (SubplotPositionInfo.SubplotDims(2) >= 4) & (SubplotPositionInfo.SubplotDims(2) < 6)
        SubplotPositionInfo.SubFigDims(1) = 0.85;
    elseif SubplotPositionInfo.SubplotDims(2) >= 6
        SubplotPositionInfo.SubFigDims(1) = 0.9;
    end
    
    % SubFigDims(2)
    if SubplotPositionInfo.SubplotDims(1) == 1
        SubplotPositionInfo.SubFigDims(2) = 0.5;
    elseif (SubplotPositionInfo.SubplotDims(1) >= 2) & (SubplotPositionInfo.SubplotDims(1) < 4)
        SubplotPositionInfo.SubFigDims(2) = 0.8;
    elseif SubplotPositionInfo.SubplotDims(1) >= 4
        SubplotPositionInfo.SubFigDims(2) = 0.9;
    end
    
    % LegendXPosition
    if SubplotPositionInfo.SubplotDims(2) < 4
        SubplotPositionInfo.LegendXPosition = 0.88;
    elseif (SubplotPositionInfo.SubplotDims(2) >= 4) & (SubplotPositionInfo.SubplotDims(2) <  7)
        SubplotPositionInfo.LegendXPosition = 0.91;
    elseif (SubplotPositionInfo.SubplotDims(2) >= 7)
        SubplotPositionInfo.LegendXPosition = 0.93;
    end
    
    % SubplotXBuffer
    SubplotPositionInfo.SubplotXBuffer = 0.07;
    
    % SubplotYBuffer
    if SubplotPositionInfo.SubplotDims(1) < 3
        if ~forTheoryPlots
            SubplotPositionInfo.SubplotYBuffer = 0.175; 
        else
           SubplotPositionInfo.SubplotYBuffer = 0.125; % GM 4/29/21
        end
    elseif SubplotPositionInfo.SubplotDims(1) == 3
        SubplotPositionInfo.SubplotYBuffer = 0.125; 
    elseif SubplotPositionInfo.SubplotDims(1) > 3
        SubplotPositionInfo.SubplotYBuffer = 0.14;
    end
    
    % SubplotBottomEdge
    if SubplotPositionInfo.SubplotDims(1) < 4
        if hasTitle
            SubplotPositionInfo.SubplotBottomEdge = 0.15;
        else
            SubplotPositionInfo.SubplotBottomEdge = 0.1;
        end
    elseif SubplotPositionInfo.SubplotDims(1) >= 4
        SubplotPositionInfo.SubplotBottomEdge = 0.07;   
    end
    
    % SubplotTopEdge
    if SubplotPositionInfo.SubplotDims(1) == 1
        SubplotPositionInfo.SubplotTopEdge = 0.15;
    elseif (SubplotPositionInfo.SubplotDims(1) >= 2) & (SubplotPositionInfo.SubplotDims(1) < 4)
        if ~forTheoryPlots & hasTitle
           SubplotPositionInfo.SubplotTopEdge = 0.15; 
        elseif ~forTheoryPlots & ~hasTitle
            SubplotPositionInfo.SubplotTopEdge = 0.075; 
        else
           SubplotPositionInfo.SubplotTopEdge = 0.05; 
        end
    elseif SubplotPositionInfo.SubplotDims(1) >= 4
        SubplotPositionInfo.SubplotTopEdge = 0.07;
    end
    
    % SubplotLeftEdge
    SubplotPositionInfo.SubplotLeftEdge = 0.06;
    
    % SubplotRightEdge
    SubplotPositionInfo.SubplotRightEdge = 0.03;
    
    
     
    SubplotPositionInfo.SubplotWidth = (SubplotPositionInfo.LegendXPosition-SubplotPositionInfo.SubplotRightEdge-...
        SubplotPositionInfo.SubplotLeftEdge-SubplotPositionInfo.SubplotXBuffer*(SubplotPositionInfo.SubplotDims(2)-2))/(SubplotPositionInfo.SubplotDims(2)-1);
    SubplotPositionInfo.XPositions = SubplotPositionInfo.SubplotLeftEdge+(0:(SubplotPositionInfo.SubplotDims(2)-2))*(SubplotPositionInfo.SubplotWidth+SubplotPositionInfo.SubplotXBuffer);
else
    
    % SubFigDims(1)
    if SubplotPositionInfo.SubplotDims(2) < 4
        SubplotPositionInfo.SubFigDims(1) = 0.7;
    elseif (SubplotPositionInfo.SubplotDims(2) >= 4) & (SubplotPositionInfo.SubplotDims(2) < 6)
        SubplotPositionInfo.SubFigDims(1) = 0.85;
    elseif SubplotPositionInfo.SubplotDims(2) >= 6
        SubplotPositionInfo.SubFigDims(1) = 0.9;
    end
    
     % SubFigDims(2)
    if SubplotPositionInfo.SubplotDims(1) == 1
        SubplotPositionInfo.SubFigDims(2) = 0.5;
    elseif (SubplotPositionInfo.SubplotDims(1) >= 2) & (SubplotPositionInfo.SubplotDims(1) < 4)
        SubplotPositionInfo.SubFigDims(2) = 0.8;
    elseif SubplotPositionInfo.SubplotDims(1) >= 4
        SubplotPositionInfo.SubFigDims(2) = 0.9;
    end
    
    % SubplotXBuffer
    SubplotPositionInfo.SubplotXBuffer = 0.07;
    
    
    % SubplotYBuffer
    SubplotPositionInfo.SubplotYBuffer = 0.14;
    
    % SubplotBottomEdge
     SubplotPositionInfo.SubplotBottomEdge = 0.1;
     
     % SubplotTopEdge
    if SubplotPositionInfo.SubplotDims(2) <= 2
        SubplotPositionInfo.SubplotTopEdge = 0.1;
    elseif SubplotPositionInfo.SubplotDims(2) > 2
        SubplotPositionInfo.SubplotTopEdge = 0.125;
    end
    
    % SubplotLeftEdge
    SubplotPositionInfo.SubplotLeftEdge = 0.06;
    
    % SubplotRightEdge
    SubplotPositionInfo.SubplotRightEdge = 0.03;
    
    SubplotPositionInfo.SubplotWidth = (1-SubplotPositionInfo.SubplotRightEdge-...
        SubplotPositionInfo.SubplotLeftEdge-SubplotPositionInfo.SubplotXBuffer*(SubplotPositionInfo.SubplotDims(2)-1))/(SubplotPositionInfo.SubplotDims(2));
    SubplotPositionInfo.XPositions = SubplotPositionInfo.SubplotLeftEdge+(0:(SubplotPositionInfo.SubplotDims(2)-1))*(SubplotPositionInfo.SubplotWidth+SubplotPositionInfo.SubplotXBuffer);
end
    







SubplotPositionInfo.SubplotHeight = (1-SubplotPositionInfo.SubplotTopEdge-SubplotPositionInfo.SubplotBottomEdge-SubplotPositionInfo.SubplotYBuffer*(SubplotPositionInfo.SubplotDims(1)-1))/SubplotPositionInfo.SubplotDims(1);

SubplotPositionInfo.YPositions = SubplotPositionInfo.SubplotBottomEdge+(0:(SubplotPositionInfo.SubplotDims(1)-1))*(SubplotPositionInfo.SubplotHeight+SubplotPositionInfo.SubplotYBuffer);
SubplotPositionInfo.YPositions = fliplr(SubplotPositionInfo.YPositions);
SubplotPositionInfo.SubplotXPositions = NaN(1, NumSets);
SubplotPositionInfo.SubplotYPositions = NaN(1, NumSets);
SubplotPositionInfo.SubplotIndexList = NaN(1, NumSets);
SubplotIndex = 0;
for i = 1:NumSets
    SubplotIndex = SubplotIndex + 1;
    if includeLegend & (mod(SubplotIndex,  SubplotPositionInfo.SubplotDims(2)) == 0)
        SubplotIndex = SubplotIndex+1;
    end
    SubplotPositionInfo.SubplotIndexList(i) = SubplotIndex;
    ColumnIndex = mod(SubplotIndex, SubplotPositionInfo.SubplotDims(2));
    if ColumnIndex == 0
        ColumnIndex = SubplotPositionInfo.SubplotDims(2);
    end
    if includeLegend
        RowIndex = fix(SubplotIndex/SubplotPositionInfo.SubplotDims(2))+1;
    else
        RowIndex = ceil(SubplotIndex/SubplotPositionInfo.SubplotDims(2));
    end
    SubplotPositionInfo.SubplotYPositions(i) = SubplotPositionInfo.YPositions(RowIndex);
    if includeLegend
        if RowIndex ~= SubplotPositionInfo.SubplotDims(1) | (SubplotPositionInfo.SubplotDims(1)*(SubplotPositionInfo.SubplotDims(2)-1) == NumSets)
            SubplotPositionInfo.SubplotXPositions(i) = SubplotPositionInfo.XPositions(ColumnIndex);
            
        else
            XShift = (SubplotPositionInfo.SubplotDims(1)*(SubplotPositionInfo.SubplotDims(2)-1) -NumSets)*...
                (SubplotPositionInfo.SubplotXBuffer+SubplotPositionInfo.SubplotWidth)/2;
            SubplotPositionInfo.SubplotXPositions(i) = SubplotPositionInfo.XPositions(ColumnIndex)+XShift;
        end
    elseif RowIndex ~= SubplotPositionInfo.SubplotDims(1) | (SubplotPositionInfo.SubplotDims(1)*(SubplotPositionInfo.SubplotDims(2)) == NumSets)
        SubplotPositionInfo.SubplotXPositions(i) = SubplotPositionInfo.XPositions(ColumnIndex);
        
    else
        XShift = (SubplotPositionInfo.SubplotDims(1)*(SubplotPositionInfo.SubplotDims(2)) -NumSets)*...
            (SubplotPositionInfo.SubplotXBuffer+SubplotPositionInfo.SubplotWidth)/2;
        SubplotPositionInfo.SubplotXPositions(i) = SubplotPositionInfo.XPositions(ColumnIndex)+XShift;
    end
    
end

