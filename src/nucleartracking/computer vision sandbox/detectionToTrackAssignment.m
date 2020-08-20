function [assignments, unassignedTracks, unassignedDetections] = ...
    detectionToTrackAssignment(tracks, measurements)

nTracks = length(tracks);
nDetections = size(measurements, 1);

% Compute the cost of assigning each detection to each track.
cost = zeros(nTracks, nDetections);
for k = 1:nTracks
    cost(k, :) = distance(tracks(k).kalmanFilter, measurements);
end

% cost(cost > gatingThresh) = 1 + gatingCost;

% %visualize cost matrix
% figure(3) ; imagesc(cost); colorbar;
% title('cost matrix')
% xlabel('prediction')
% ylabel('detection')
% set(gca,'ColorScale','log')



% tracks.kalmanFilter.State

% Solve the assignment problem.

%AR- i don't know how to properly adjust
%this number.
costOfNonAssignment = 65;
[assignments, unassignedTracks, unassignedDetections] = ...
    assignDetectionsToTracks(cost, costOfNonAssignment);
% 
% try 
% predictions = tracks(1).kalmanFilter.State;
% %debugging figure
% figure(3);
% plot(predictions(:,1),predictions(:,2),'*',detections(:,1),...
%     detections(:,2),'ro');
% hold on;
% legend('predictions','detections');
% for k = 1:size(assignment,1)
%     text(predictions(assignment(k, 1),1)+0.1,...
%         predictions(assignment(k,1),2)-0.1,num2str(k));
%     text(detections(assignment(k, 2),1)+0.1,...
%         detections(assignment(k,2),2)-0.1,num2str(k));
% end
% for k = 1:length(unassignedDetections)
%     text(detections(unassignedDetections(k),1)+0.1,...
%         detections(unassignedDetections(k),2)+0.1,'unassigned');
% end
% xlim([0,4]);
% ylim([0,4]);
% end

% figure(4); imagesc(assignments); colorbar;

end