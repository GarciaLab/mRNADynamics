function ExtractAPprofileInfoNC14(Prefix, varargin)
% ExtractAPprofileInfo.m
% author: Gabriella Martini
% 12/06/19

[~,~,DefaultDropboxFolder,~,~]=...
    DetermineLocalFolders;


ConvertCompiledNucleiToTableArray(Prefix);
%Get the folders, including the default Dropbox one
[SourcePath, FISHPath, DefaultDropboxFolder, DropboxFolder, MS2CodePath, PreProcPath,...
configValues, movieDatabasePath] = DetermineAllLocalFolders(Prefix);



%Determine division times
%Load the information about the nc from moviedatabase file
[Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
nc9, nc10, nc11, nc12, nc13, nc14, CF] = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder);
DataFolder=[DropboxFolder,filesep,Prefix];
FilePrefix=[DataFolder(length(DropboxFolder)+2:end),'_'];
load([DropboxFolder,filesep, Prefix,'\CompiledNucleiTable.mat'])
cnt14  = CompiledNucleiTable(CompiledNucleiTable.nc == 14, :);
cnt14  = cnt14(cnt14.Fluo > 0,:);
%%
close all
numFrames = max(cnt14.FrameNC)+1;
numBins = 1/APResolution+1;
MeanFluoAP = zeros(numFrames, numBins);
StdFluoAP = zeros(numFrames, numBins);
NumNucAP = zeros(numFrames, numBins);
profFolder = 'APprofiles';
timeFolder = 'Timeprofiles';
% First, get the name of the folder you're using.
% For example if your folder is 'D:\photos\Info', parentFolder  would = 'D:\photos, and deepestFolder would = 'Info'.
[parentFolder deepestFolder] = fileparts([DropboxFolder,filesep,Prefix]);
% Next, create a name for a subfolder within that.
% For example 'D:\photos\Info\DATA-Info'
newSubFolder = sprintf('%s/%s/%s', parentFolder, deepestFolder, profFolder);
% Finally, create the folder if it doesn't exist already.
if ~exist(newSubFolder, 'dir')
  mkdir(newSubFolder);
end
newSubFolder = sprintf('%s/%s/%s', parentFolder, deepestFolder, timeFolder);
% Finally, create the folder if it doesn't exist already.
if ~exist(newSubFolder, 'dir')
  mkdir(newSubFolder);
end
APbins = 0:0.025:1;
fc = 1;
for f = 0:max(cnt14.FrameNC)
    for b = 0:(numBins-1)
        APbin = b*APResolution;
        sub14 = cnt14((cnt14.FrameNC == f) & (abs(cnt14.APbin- APbin) < 0.0001), :);
        cutoff = quantile(sub14.Fluo, .1);
        fluoVector = sub14.Fluo(sub14.Fluo >= cutoff);
        MeanFluo = mean(fluoVector);
        StdFluo = std(fluoVector);
        MeanFluoAP(f+1, b+1) = MeanFluo;
        StdFluoAP(f+1, b+1) = StdFluo;
        NumNucAP(f+1, b+1) = size(fluoVector, 1);
    end
    fig = figure(fc);
    ax = axes('Parent', fig);
    plot(APbins, MeanFluoAP(f+1,:), '.');
    
    %histogram(sub14.Fluo)
    xlabel('Fraction Embryo Length')
    ylabel('Fluo (AU)')
    ylim([0, 1600])
    xlim([0, 0.95])
    %CurrentFrame = sub14.Frame(1);

%     ImageHis=imread([PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
%             FilePrefix(1:end-1),'-His_',iIndex(CurrentFrame,3),'.tif']);
%     imshow(ImageHis,'DisplayRange',[]);
    title(ax, ['Frame: ',num2str(f)])
%     hold on 

    
%     scatter(sub14.xPos, sub14.yPos,40, sub14.Fluo, 'filled')
%     plot(sub14.MedianDV, sub14.Fluo,'o')


    %ylabel('Fluo (AU)')
    %xlabel('Median DV (pixels)')
    saveas(fig,[DropboxFolder,filesep,Prefix,filesep, profFolder, filesep, 'FluoProfile_frame', num2str(f),'.png']);
    fc = fc+1;
end

close all
frames = 0:max(cnt14.FrameNC);
fc = 1;
for b = 0:(numBins-1)
    fig = figure(fc);
    ax = axes('Parent', fig);
    plot(frames, MeanFluoAP(:,b+1), '.');
    
    %histogram(sub14.Fluo)
    xlabel('Frame')
    ylabel('Fluo (AU)')
    ylim([0, 1600])
    %xlim([0.05, 0.95])
    %CurrentFrame = sub14.Frame(1);

%     ImageHis=imread([PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
%             FilePrefix(1:end-1),'-His_',iIndex(CurrentFrame,3),'.tif']);
%     imshow(ImageHis,'DisplayRange',[]);
    title(ax, ['Fraction Embryo Length: ',num2str(b*APResolution)])
%     hold on 

    saveas(fig,[DropboxFolder,filesep,Prefix,filesep, timeFolder, filesep,'TimeProfile_AP', num2str(b),'.png']);
    fc = fc+1;
end

close all

save([DropboxFolder,filesep,Prefix,filesep,'MeanFluoAPNC14.mat'],...
        'MeanFluoAP');
    
save([DropboxFolder,filesep,Prefix,filesep,'StdFluoAPNC14.mat'],...
        'StdFluoAP');
    
save([DropboxFolder,filesep,Prefix,filesep,'NumNucAPNC14.mat'],...
        'NumNucAP');
