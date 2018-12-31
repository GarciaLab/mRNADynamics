function pointGroupNumbers = groupingPoints(currentTimeArray,smoothedAmp)

dydx = diff([eps; smoothedAmp'])./diff([eps; currentTimeArray']);
pointsSlopeState = zeros(1,length(dydx));
% assigning points' slope state
pointsSlopeState(dydx<0) = -1;
pointsSlopeState(dydx>=0) = 1;

% Assigning points to groups based on states ----------------------
% done by looking for crude version of point of inflections
pointGroupNumbers = zeros(1,length(dydx));
% checking points up until the last one
pointRangeToSearch = 1:length(dydx)-1;
searchingPoint = pointRangeToSearch(1);
currentGroupNumber = 1;
% Grouping the points based on change in the points' slope
% state (ignores fluctions of duration of one point)
for currentPoint = pointRangeToSearch
    searchingPoint = searchingPoint + 1;
    pointGroupNumbers(currentPoint) = currentGroupNumber;
    currentState = pointsSlopeState(currentPoint);
    nextState = pointsSlopeState(searchingPoint);
    
    currentToNeighborEquality = isequal(currentState,nextState);
    
    if ~currentToNeighborEquality
        if currentPoint > 1
            previousState = pointsSlopeState(currentPoint-1);
            if ~isequal(previousState,nextState)
                currentGroupNumber = currentGroupNumber + 1;
            end
        end
    end
end
% assigning the last point to the last group made
pointGroupNumbers(end) = currentGroupNumber;

% d2ydx2 = diff([eps; dydx])./diff([eps; currentTimeArray']);
% plot(currentTimeArray,d2ydx2,'DisplayName','d2')
% hold on
% plot(currentTimeArray,dydx,'DisplayName','d1')
% plot(currentTimeArray,smoothedAmp,'DisplayName','Original')
% legend()

            
end
