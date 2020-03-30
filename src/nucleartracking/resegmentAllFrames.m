function resegmentAllFrames(Prefix, varargin)

hisMat = [];
min_rad_um = 2; %default. good for fly embryos. 
max_rad_um = 6;

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for k = 1:2:(numel(varargin)-1)
    if k ~= numel(varargin)
        eval([varargin{k} '=varargin{k+1};']);
    end
end

thisExperiment = liveExperiment(Prefix);

if isempty(hisMat)    
    hisMat = getHisMat(thisExperiment);
end

Ellipses = cell(thisExperiment.nFrames, 1);

parfor CurrentFrame = 1:thisExperiment.nFrames
        
        HisImage = hisMat(:, :, CurrentFrame);
        [~, circles] = kSnakeCircles(HisImage,...
            thisExperiment.pixelSize_nm/1000, 'min_rad_um', min_rad_um,...
            'max_rad_um', max_rad_um);    
        circles(:, 6:9) = zeros(size(circles, 1), 4);
        Ellipses{CurrentFrame} = circles;
                
        %report progress every tenth frame
        if ~mod(CurrentFrame, 10), disp(['Segmenting frame: ' num2str(CurrentFrame), '...']); end
        
end

save([thisExperiment.resultsFolder, 'Ellipses.mat'],'Ellipses', '-v6');

end