function resegmentAllFrames(Prefix, varargin)

thisExperiment = liveExperiment(Prefix);

hisMat = getHisMat(thisExperiment);
% if exist([thisExperiment.resultsFolder, filesep, 'Ellipses.mat'], 'file')
%     Ellipses = getEllipses(thisExperiment);
% else
Ellipses = cell(thisExperiment.nFrames, 1);
% end

parfor CurrentFrame = 1:thisExperiment.nFrames
        
        HisImage = hisMat(:, :, CurrentFrame);
        [~, circles] = kSnakeCircles(HisImage, thisExperiment.pixelSize_nm/1000);    
        circles(:, 4) = circles(:, 3);
        circles(:, 5:9) = zeros(size(circles, 1), 5);
        Ellipses{CurrentFrame} = circles;
                
        %report progress every tenth frame
        if ~mod(CurrentFrame, 10), disp(num2str(CurrentFrame)); end
        
end

save([thisExperiment.resultsFolder, 'Ellipses.mat'],'Ellipses', '-v6');

TrackNuclei(Prefix,'retrack');


end