function Ellipses = resegmentAllFrames(Prefix, varargin)

hisMat = [];
%defaults. good for fly embryos (Dorsal synthetics settings)
% to change use them as options in TrackNuclei
sigmaK_um = 0.2;
min_rad_um = 2; 
max_rad_um = 6;
mu = 0.1;
nIterSnakes = 100;

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for k = 1:2:(numel(varargin)-1)
    if k ~= numel(varargin)
        eval([varargin{k} '=varargin{k+1};']);
    end
end

thisExperiment = liveExperiment(Prefix);
pixelSize_um = thisExperiment.pixelSize_um;
hisMat = getHisMat(thisExperiment);


Ellipses = cell(thisExperiment.nFrames, 1);

parfor CurrentFrame = 1:thisExperiment.nFrames
        
        
        HisImage = hisMat(:, :, CurrentFrame);
        
        %segment current frame if nonempty
        %and add ellipse information to Ellipses
        if max(max(HisImage))~=0
            [~, ellipseFrame] = kSnakeCircles(HisImage,...
                pixelSize_um, 'min_rad_um', min_rad_um,...
                'max_rad_um', max_rad_um,'sigmaK_um',sigmaK_um,'mu', mu,...
            'nIterSnakes',nIterSnakes);    
            ellipseFrame(:, 6:9) = zeros(size(ellipseFrame, 1), 4);
            Ellipses{CurrentFrame} = ellipseFrame;
        else
            Ellipses{CurrentFrame} = [];
        end
                
        %report progress every tenth frame
        if ~mod(CurrentFrame, 10),...
                disp(['Segmenting frame: ' num2str(CurrentFrame), '...']); end
        
end

if whos(var2str(Ellipses)).bytes < 2E9
    save([thisExperiment.resultsFolder, 'Ellipses.mat'],var2str(Ellipses), '-v6');
else
    save([thisExperiment.resultsFolder, 'Ellipses.mat'],var2str(Ellipses), '-v7.3', '-nocompression');
end

end