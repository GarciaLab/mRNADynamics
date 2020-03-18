function resegmentAllFrames(Prefix, varargin)

thisExperiment = liveExperiment(Prefix);

hisMat = getHisMat(thisExperiment);
Ellipses = getEllipses(thisExperiment);

for CurrentFrame = 1:thisExperiment.nFrames 
        
        HisImage = hisMat(:, :, CurrentFrame);
        [~, circles] = kSnakeCircles(HisImage, thisExperiment.pixelSize_nm/1000);    
        circles(:, 4) = circles(:, 3);
        circles(:, 5:9) = zeros(size(circles, 1), 5);
        Ellipses{CurrentFrame} = circles;
        
end

save([thisExperiment.DropboxFolder, 'Ellipses.mat'],'Ellipses', '-v6');

TrackNuclei(Prefix,'NoBulkShift','ExpandedSpaceTolerance', 1.5, 'retrack');


end