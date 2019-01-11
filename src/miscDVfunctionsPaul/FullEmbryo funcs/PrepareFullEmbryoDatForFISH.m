function PrepareFullEmbryoDatForFISH(Prefix)

%Get the folders, including the default Dropbox one
[SourcePath, FISHPath, DefaultDropboxFolder, DropboxFolder, MS2CodePath, PreProcPath,...
configValues, movieDatabasePath] = DetermineAllLocalFolders(Prefix);

%Determine division times
%Load the information about the nc from moviedatabase file
[Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
nc9, nc10, nc11, nc12, nc13, nc14, CF] = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder);

copyfile([SourcePath '/' Prefix(1:10) '/' Prefix(12:end) '/' 'FullEmbryo' '/' 'surf*'],...
    [SourcePath '/' Prefix(1:10) '/' Prefix(12:end)])

mkdir([SourcePath '/' Prefix(1:10) '/' Prefix(12:end) '/' 'MetaData'])
copyfile([SourcePath '/' Prefix(1:10) '/' Prefix(12:end) '/' 'FullEmbryo' '/' 'MetaData'],...
    [SourcePath '/' Prefix(1:10) '/' Prefix(12:end) '/' 'MetaData'])

list = dir([SourcePath '/' Prefix(1:10) '/' Prefix(12:end)]);

logical = zeros(1,length(list));

for q = 1:length(list)
    if length(list(q).name) > 3
        logical(q) = strcmp(list(q).name(end-3:end),'.lif');
    end
end

lif_ind = find(logical == 1);

movefile([SourcePath '/' Prefix(1:10) '/' Prefix(12:end) '/' list(lif_ind).name],...
    [SourcePath '/' Prefix(1:10) '/' Prefix(12:end) '/' Prefix '.lif']);

new_logical = zeros(1,length(list));

for r = 1:length(list)
    if length(list(r).name) > 3
        new_logical(r) = strcmp(list(r).name(end-3:end),'.tif');
    end
end

tif_ind = find(new_logical == 1);

for p = 1:length(tif_ind)

    movefile([SourcePath '/' Prefix(1:10) '/' Prefix(12:end) '/' list(tif_ind(p)).name],...
    [SourcePath '/' Prefix(1:10) '/' Prefix(12:end) '/' Prefix '_001_z0' num2str(p) '_ch01' '.tif']);

end

newlist = dir([SourcePath '/' Prefix(1:10) '/' Prefix(12:end)]);

newest_logical = zeros(1,length(newlist));

for s = 1:length(newlist)
    if length(newlist(s).name) > 3
        newest_logical(s) = strcmp(newlist(s).name(end-3:end),'.tif');
    end
end

newtif_ind = find(newest_logical == 1);

for i = 1:2
    for m = 1:length(newtif_ind)
        
        copyfile([SourcePath '/' Prefix(1:10) '/' Prefix(12:end) '/' newlist(newtif_ind(m)).name],...
            [SourcePath '/' Prefix(1:10) '/' Prefix(12:end) '/' Prefix '_00' num2str(i+1) '_z0' num2str(m) '_ch01' '.tif']);
        
    end
end

%Create the output folder
OutputFolder = [PreProcPath, filesep, Prefix];
mkdir(OutputFolder)

res_info = imfinfo([SourcePath '/' Prefix(1:10) '/' Prefix(12:end) '/' Prefix '_00' num2str(i) '_z0' num2str(m) '_ch01' '.tif']);

%Generate FrameInfo
FrameInfo = struct('LinesPerFrame', {}, 'PixelsPerLine', {}, ...
    'NumberSlices', {}, 'FileMode', {}, 'Time', {}, ...
    'PixelSize', {},'ZStep', {}, 'zPosition', {}, 'NChInput', {});

for x = 1:3
    
    FrameInfo(x).LinesPerFrame = res_info.Height;
    FrameInfo(x).PixelsPerLine = res_info.Width;
    FrameInfo(x).NumberSlices = length(newtif_ind) - 3;
    FrameInfo(x).FileMode = 'LIFExport';
    FrameInfo(x).PixelSize = 0.240600195694716;
    FrameInfo(x).ZStep = 0.499725800000000;
    FrameInfo(x).Time = 0;
    FrameInfo(x).zPosition = -9.00007152564194e-06;
    FrameInfo(x).NChInput = 1;
    
end

mkdir([DropboxFolder, filesep, Prefix]);
save([DropboxFolder, filesep, Prefix, filesep, 'FrameInfo.mat'], 'FrameInfo');

mkdir([PreProcPath, filesep, Prefix]);

newnewlist = dir([SourcePath '/' Prefix(1:10) '/' Prefix(12:end)]);

newnew_logical = zeros(1,length(newnewlist));

for y = 1:length(newnewlist)
    if length(newnewlist(y).name) > 3
        newnew_logical(y) = strcmp(newnewlist(y).name(end-3:end),'.tif');
    end
end

newnewtif_ind = find(newnew_logical == 1);

for v = 1:length(newnewtif_ind)
    
    movefile([SourcePath '/' Prefix(1:10) '/' Prefix(12:end) '/' newnewlist(newnewtif_ind(v)).name],...
        [PreProcPath, filesep, Prefix])
    
end

for c = 1:3
    
    movefile([PreProcPath, filesep, Prefix '/' Prefix '_00' num2str(c) '_z0' num2str(length(newnewtif_ind) / 3) '_ch01' '.tif'],...
        [PreProcPath, filesep, Prefix '/' Prefix '-His_0' num2str(c) '.tif']);
    
end

end