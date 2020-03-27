function resegmentAllFrames(Prefix, varargin)

hisMat = [];
min_rad_um = 2; %default. good for fly embryos. 
max_rad_um = 6;

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end


thisExperiment = liveExperiment(Prefix);

if isempty(hisMat)    
    hisMat = getHisMat(thisExperiment);
end
% if exist([thisExperiment.resultsFolder, filesep, 'Ellipses.mat'], 'file')
%     Ellipses = getEllipses(thisExperiment);
% else
Ellipses = cell(thisExperiment.nFrames, 1);
% end

parfor CurrentFrame = 1:thisExperiment.nFrames
        
        HisImage = hisMat(:, :, CurrentFrame);
        [~, circles] = kSnakeCircles(HisImage,...
            thisExperiment.pixelSize_nm/1000, 'min_rad_um', min_rad_um,...
            'max_rad_um', max_rad_um);    
%         circles(:, 4) = circles(:, 3);
        circles(:, 6:9) = zeros(size(circles, 1), 4);
        Ellipses{CurrentFrame} = circles;
                
        %report progress every tenth frame
        if ~mod(CurrentFrame, 10), disp(['Segmenting frame: ' num2str(CurrentFrame), '...']); end
        
end

save([thisExperiment.resultsFolder, 'Ellipses.mat'],'Ellipses', '-v6');

TrackNuclei(Prefix,'retrack');


end