function GenerateSegmentationLabels(varargin)
% 
if isempty(varargin)
    warning('No Prefix defined. Will output default Dropbox folder')
end

%%%Default options 
n_sigma=2;
WriteWeights = 1;

Prefix = varargin{1};
for i=1:length(varargin)
    if strcmp(varargin{i},'n_sigma')
        n_sigma=varargin{i+1};
    end
    if strcmp(varargin{i},'WriteWeights')
        WriteWeights = 1;
    end
end

%Get filepaths
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders(Prefix);
%Load Spots Data Set
load([DropboxFolder filesep Prefix filesep 'Spots.mat'])
%Load Frame Info (for image dim info)
load([DropboxFolder filesep Prefix filesep 'FrameInfo.mat'])
%Load first image in 'ProcessedData' folder to get image Dim info
files = dir([FISHPath, filesep, Prefix, '_', filesep, 'dogs',filesep, '*.tif']);
firstim = imread([FISHPath, filesep, Prefix, '_', filesep, 'dogs', filesep, files(1).name]);
%Set dimensions
zDim = FrameInfo(1).NumberSlices + 2;
xDim = size(firstim,2);
yDim = size(firstim,1);

%Set write path
BinaryPath = [FISHPath filesep Prefix '_' filesep 'binary_masks'];
if exist(BinaryPath) ~= 7
    mkdir(BinaryPath);
end

WeightPath = [FISHPath filesep Prefix '_' filesep 'weight_masks'];
if exist(WeightPath) ~= 7
    mkdir(WeightPath);
end
[X_mesh,Y_mesh] = meshgrid(1:xDim,1:yDim);
%Define position of each pixel as its centroid
X_mesh = X_mesh + .5;
Y_mesh = Y_mesh + .5;
f_max = Spots(end).Fits(1).frame;
digits = floor(log10(f_max));
zero_string = repelem('0',digits);
%Iterate through each time point
for t = 1:f_max
    Fits = Spots(t).Fits;
    zStack_binary = cell(1,zDim);
    zStack_potential = cell(1,zDim);
    %Intilize
    for i = 1:zDim
        zStack_binary{i} = zeros(yDim,xDim);  
        zStack_potential{i} = zeros(yDim,xDim);  
    end
    %Iterate through spots
    for s = 1:length(Fits)
        spot = Fits(s);
        %Iterate through z slices within a spot
        for z = 1:length(spot.z)
            zPos = spot.z(z);
            xPos = spot.xFit(z);
            yPos = spot.yFit(z);
            %set width of spot
            xTol = n_sigma*spot.xFitWidth{z};
            yTol= n_sigma*spot.yFitWidth{z};
            % Generate 2D array with radial distances from current spot            
            R_array = (X_mesh - floor(xPos)).^2 / xTol^2 + (Y_mesh - floor(yPos)).^2 / yTol^2;            
            % Add spot region to appropriate z slice for time point
            spot_slice = zStack_binary{zPos}; 
            % At a minumum add the pixel in which the spot is centered
            spot_slice(floor(yPos),floor(xPos)) = 1;
            % Now add proximate pixes
            spot_slice(R_array <=1) = 1; 
            zStack_binary{zPos} = spot_slice;
            %Addproximity score
            r_slice = zStack_potential{zPos}; 
            r_slice = r_slice + R_array.^.5; 
            r_slice(spot_slice == 1) = 0;
            zStack_potential{zPos} = r_slice;            
        end
    end
    %Write labeled stack to Tif 
    d = floor(log10(t));
    imwrite(uint16(zStack_binary{1}), [BinaryPath '/' 'binary_stack_' Prefix '_' zero_string(d+1:digits) num2str(t) '.tiff']);
    for z = 2:length(zStack_binary)
        imwrite(uint16(zStack_binary{z}), [BinaryPath '/' 'binary_stack_' Prefix '_' zero_string(d+1:digits) num2str(t) '.tiff'], 'WriteMode', 'append');
    end
    if WriteWeights
        %Write weights (Just testing for now)
        imwrite(uint16(zStack_potential{1}), [WeightPath '/' 'weight_stack_' Prefix '_' zero_string(d+1:digits) num2str(t) '.tiff']);
        for z = 2:length(zStack_potential)
            imwrite(uint16(zStack_potential{z}), [WeightPath '/' 'weight_stack_' Prefix '_' zero_string(d+1:digits) num2str(t) '.tiff'], 'WriteMode', 'append');
        end
    end
end