function ExtractNuclearFluoProfilesAP(Prefix, varargin)
% ExtractNuclearFluoProfilesAP.m
% author: Gabriella Martini
% date created: 8/27/20
% date last modified: 9/26/20
%% 

nuclear_cycles = [12, 13, 14];
quantile_cutoff = 0;
min_nuc_count = 5;
min_fluo_pos_range = 0.55:0.025:0.7;

if ~isempty(varargin)
    x = 1;
    while x <= length(varargin)
        switch varargin{x}
            case{'NuclearCycles'}
                 nuclear_cycles = varargin{x+1};
                 x = x+1;
            case{'QuantileCutoff'}
                quantile_cutoff = varargin{x+1};
                x = x+1;
            case{'MinNucBinCount'}
                min_nuc_count = varargin{x+1};
                x = x+1;
            case{'MinFluoPosRange'}
                min_fluo_pos_range = varargin{x+1};
                x = x+1;
            otherwise
                error(['Flag "', varargin{x},'" not valid'])
        end
        x = x +1;
    end
end
%% 

[~,~,DefaultDropboxFolder,~,~]=...
    DetermineLocalFolders;

[SourcePath, FISHPath, DefaultDropboxFolder, DropboxFolder, MS2CodePath, PreProcPath,...
configValues, movieDatabasePath] = DetermineAllLocalFolders(Prefix);


if exist([DropboxFolder,filesep,Prefix,filesep,'CompiledNucleiTable.mat'], 'file') 
    load([DropboxFolder,filesep,Prefix,filesep,'CompiledNucleiTable.mat']);
else
    CompiledNucleiTable = ConvertCompiledNucleiToTableArray(Prefix);
end
%Get the folders, including the default Dropbox one
load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat']);
load([DropboxFolder,filesep,Prefix,filesep,'APDivision.mat']);

CompiledNucleiTable((CompiledNucleiTable.nc >= min(nuclear_cycles)) & (CompiledNucleiTable.nc <= max(nuclear_cycles)), :);
%Determine division times
%Load the information about the nc from moviedatabase file
[Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
nc9, nc10, nc11, nc12, nc13, nc14, CF] = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder);
DataFolder=[DropboxFolder,filesep,Prefix];
FilePrefix=[DataFolder(length(DropboxFolder)+2:end),'_'];
Prefix_label =  strrep(Prefix, '_5C', '.5C');
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
%% 
nc_frames = [nc9, nc10, nc11, nc12, nc13, nc14];
min_fluo_pos_idx = uint16(min_fluo_pos_range/APResolution+1);
numBins = uint16(1/APResolution+1);
for nc=nuclear_cycles
    table_lengths = [];
    disp(['Nuclear cycle: ', num2str(nc)])
    cnt  = CompiledNucleiTable(CompiledNucleiTable.nc == nc, :);
    
    cnt  = cnt(cnt.Fluo > 0,:);
    table_lengths(length(table_lengths) + 1) = height(cnt);
    disp(['Compiled Nuclei Table row count: ', num2str(height(cnt))])
    cnt = cnt(cnt.FrameApproved == 1,:);
    table_lengths(length(table_lengths) + 1) = height(cnt);
    disp(['Compiled Nuclei Table new row count after Frame Approved: ', num2str(height(cnt)),...
        ' (',num2str(height(cnt)/table_lengths(1)), ')'])
    cnt  = cnt(cnt.Flag1 == 0,:);
    table_lengths(length(table_lengths) + 1) = height(cnt);
    disp(['Compiled Nuclei Table new row count after Flag 1 filter: ', num2str(height(cnt)),...
        ' (',num2str(height(cnt)/table_lengths(1)), ')'])
    cnt  = cnt(cnt.Flag2 <  0.5,:);
    table_lengths(length(table_lengths) + 1) = height(cnt);
    disp(['Compiled Nuclei Table new row count after Flag 2 filter: ', num2str(height(cnt)),...
        ' (',num2str(height(cnt)/table_lengths(1) ), ')'])
    cnt  = cnt(cnt.Flag3 <  0.5,:);
    table_lengths(length(table_lengths) + 1) = height(cnt);
    disp(['Compiled Nuclei Table new row count after Flag 3 filter: ', num2str(height(cnt)),...
        ' (',num2str(height(cnt)/table_lengths(1)), ')'])
    cnt  = cnt(cnt.Flag4 <  0.5,:);
    table_lengths(length(table_lengths) + 1) = height(cnt);
    disp(['Compiled Nuclei Table new row count after Flag 4 filter: ', num2str(height(cnt)),...
        ' (',num2str(height(cnt)/table_lengths(1)), ')'])
    cnt = cnt(cnt.Flag5 < 0.75,:);
    table_lengths(length(table_lengths) + 1) = height(cnt);
    disp(['Compiled Nuclei Table new row count after Flag 5 filter: ', num2str(height(cnt)),...
        ' (',num2str(height(cnt)/table_lengths(1)), ')'])
    if nc == 14
        cnt = cnt(cnt.Flag6 == 0,:);
        table_lengths(length(table_lengths) + 1) = height(cnt);
        disp(['Compiled Nuclei Table new row count after Flag 6 filter: ', num2str(height(cnt)),...
            ' (',num2str(height(cnt)/table_lengths(1)), ')'])
    end
    cnt  = cnt(cnt.Flag7 < 0.2,:);
    table_lengths(length(table_lengths) + 1) = height(cnt);
    disp(['Compiled Nuclei Table new row count after Flag 7 filter: ', num2str(height(cnt)),...
        ' (',num2str(height(cnt)/table_lengths(1)), ')'])
    if isempty(cnt)
        continue
    end
    %%
    numFrames = max(cnt.FrameNC)+1;

    MeanFluoAP = NaN(numFrames, numBins);
    StdFluoAP = NaN(numFrames, numBins);
    NumNucAP = zeros(numFrames, numBins);
    DiffMeanFluoAP = NaN(numFrames, numBins);
    DiffStdErrorAP = NaN(numFrames, numBins);
    APbins = 0:APResolution:1;
    NCFrames = unique(cnt.FrameNC);
    NCTimes = [];
    for i=1:length(NCFrames)
        NCTimes(i) =  mean(cnt.TimeNC(cnt.FrameNC == NCFrames(i)));
    end
    for i = 1:length(NCFrames)
        f = NCFrames(i);
        t = NCTimes(i);
        for APbin = min(cnt.APbinIdx):max(cnt.APbinIdx)
            
            sub_table = cnt((cnt.FrameNC == f) & (cnt.APbinIdx == APbin), :);
            if (height(sub_table) < 5)
                continue 
            end
            if (height(sub_table) >= 10) & (quantile_cutoff > 0)
                cutoff = quantile(sub_table.Fluo, quantile_cutoff);
                fluoVector = sub_table.Fluo(sub_table.Fluo >= cutoff);
            else
                fluoVector = sub_table.Fluo;
            end
            MeanFluo = mean(fluoVector);
            StdFluo = std(fluoVector);
            MeanFluoAP(i, APbin) = MeanFluo;
            StdFluoAP(i, APbin) = StdFluo;
            NumNucAP(i, APbin) = size(fluoVector, 1);
        end
        if sum(~isnan(MeanFluoAP(i,  min_fluo_pos_idx ))) > 0
            MinFluoAPidx = find(MeanFluoAP(i,  min_fluo_pos_idx ) == min(MeanFluoAP(i,  min_fluo_pos_idx )));
            if length(MinFluoAPidx) > 1
                MinFluoAPidx = MinFluoAPidx(end);
            end
            idx = min_fluo_pos_idx(MinFluoAPidx);
            DiffMeanFluoAP(i,:) = MeanFluoAP(i,:) - MeanFluoAP(i,idx); 
            stdErrorAPframe = StdFluoAP(i,:)./sqrt(NumNucAP(i,:));
            DiffStdErrorAP(i,:) = sqrt(stdErrorAPframe.^2 + stdErrorAPframe(idx).^2);
        end      
    end
    
    
% 

 

    save([DropboxFolder,filesep,Prefix,filesep,'MeanFluoAPNC', num2str(nc),'_FramesNC.mat'],...
            'MeanFluoAP');

    save([DropboxFolder,filesep,Prefix,filesep,'StdFluoAPNC',num2str(nc),'_FramesNC.mat'],...
            'StdFluoAP');

    save([DropboxFolder,filesep,Prefix,filesep,'NumNucAPNC',num2str(nc),'_FramesNC.mat'],...
            'NumNucAP');
    save([DropboxFolder,filesep,Prefix,filesep,'DiffMeanFluoAPNC',num2str(nc),'_FramesNC.mat'],...
            'DiffMeanFluoAP');
    save([DropboxFolder,filesep,Prefix,filesep,'DiffStdErrorAPNC',num2str(nc),'_FramesNC.mat'],...
            'DiffStdErrorAP');
    save([DropboxFolder,filesep,Prefix,filesep,'TimeNC',num2str(nc),'_FramesNC.mat'],...
            'NCTimes');

        
        
        
%% 
        

        


    if max(unique(cnt.timeSinceAnaphase)) == 0
        disp('Cannot create profiles using time Since Anaphase')
        continue
    end
    AllNCtimesSinceAnaphase = unique(cnt.timeSinceAnaphase);
    BinnedNCtimesSinceAnaphase = zeros(1, length(AllNCtimesSinceAnaphase));
    for i = 2:length(AllNCtimesSinceAnaphase)
        if AllNCtimesSinceAnaphase(i) - BinnedNCtimesSinceAnaphase(i-1) < 0.2 
            BinnedNCtimesSinceAnaphase(i) = BinnedNCtimesSinceAnaphase(i-1);
        else
            BinnedNCtimesSinceAnaphase(i) = AllNCtimesSinceAnaphase(i);
        end
    end
        
    BinnedTimeSinceAnaphase = zeros(1, length(cnt.timeSinceAnaphase));
    timesSinceAnaphase = zeros(1, length(cnt.timeSinceAnaphase));
    timesSinceAnaphase(:) = cnt.timeSinceAnaphase;
    for idx=1:height(cnt)
       BinnedTimeSinceAnaphase(idx) =  BinnedNCtimesSinceAnaphase(AllNCtimesSinceAnaphase == timesSinceAnaphase(idx));
    end
    cnt.BinnedTimeSinceAnaphase = BinnedTimeSinceAnaphase.';
    NCTimes = unique(BinnedNCtimesSinceAnaphase);
    numTimes = length(NCTimes);
    MeanFluoAP = NaN(numTimes, numBins);
    StdFluoAP = NaN(numTimes, numBins);
    NumNucAP = zeros(numTimes, numBins);
    DiffMeanFluoAP = NaN(numTimes, numBins);
    DiffStdErrorAP = NaN(numTimes, numBins);
   
    
    %NCTimes = 
    for i = 1:numTimes
        t = NCTimes(i);
        for APbin = min(cnt.APbinIdx):max(cnt.APbinIdx)
            sub_table = cnt((cnt.BinnedTimeSinceAnaphase == t) & (cnt.APbinIdx == APbin), :);
            if height(sub_table) < 5
                continue
            end
            if (height(sub_table) >= 10) & (quantile_cutoff > 0)
                cutoff = quantile(sub_table.Fluo, quantile_cutoff);
                fluoVector = sub_table.Fluo(sub_table.Fluo >= cutoff);
            else
                fluoVector = sub_table.Fluo;
            end
            MeanFluo = mean(fluoVector);
            StdFluo = std(fluoVector);
            MeanFluoAP(i, APbin) = MeanFluo;
            StdFluoAP(i,APbin) = StdFluo;
            NumNucAP(i, APbin) = size(fluoVector, 1);
        end
        if sum(~isnan(MeanFluoAP(i,  min_fluo_pos_idx )))
            MinFluoAPidx = find(MeanFluoAP(i,  min_fluo_pos_idx ) == min(MeanFluoAP(i,  min_fluo_pos_idx )));
            if length(MinFluoAPidx) > 1
                MinFluoAPidx = MinFluoAPidx(end);
            end
            idx = min_fluo_pos_idx(MinFluoAPidx);
            DiffMeanFluoAP(i,:) = MeanFluoAP(i,:) - MeanFluoAP(i,idx); 
            stdErrorAPframe = StdFluoAP(i,:)./sqrt(NumNucAP(i,:));
            DiffStdErrorAP(i,:) = sqrt(stdErrorAPframe.^2 + stdErrorAPframe(idx).^2);
        end   
       
    end



    save([DropboxFolder,filesep,Prefix,filesep,'MeanFluoAPNC', num2str(nc),'_timesSinceAnaphase.mat'],...
            'MeanFluoAP');

    save([DropboxFolder,filesep,Prefix,filesep,'StdFluoAPNC',num2str(nc),'__timesSinceAnaphase.mat'],...
            'StdFluoAP');

    save([DropboxFolder,filesep,Prefix,filesep,'NumNucAPNC',num2str(nc),'__timesSinceAnaphase.mat'],...
            'NumNucAP');
    save([DropboxFolder,filesep,Prefix,filesep,'DiffMeanFluoAPNC',num2str(nc),'__timesSinceAnaphase.mat'],...
            'DiffMeanFluoAP');
    save([DropboxFolder,filesep,Prefix,filesep,'DiffStdErrorAPNC',num2str(nc),'__timesSinceAnaphase.mat'],...
            'DiffStdErrorAP');
    save([DropboxFolder,filesep,Prefix,filesep,'BinnedTimesSinceAnaphaseNC',num2str(nc),'.mat'],...
        'NCTimes');

end


%%
table_lengths = [];
disp(['All Nuclear Cycles'])

cnt  = CompiledNucleiTable(CompiledNucleiTable.Fluo > 0,:);
table_lengths(length(table_lengths) + 1) = height(cnt);
disp(['Compiled Nuclei Table row count: ', num2str(height(cnt))])
cnt = cnt(cnt.FrameApproved == 1,:);
table_lengths(length(table_lengths) + 1) = height(cnt);
disp(['Compiled Nuclei Table new row count after Frame Approved: ', num2str(height(cnt)),...
    ' (',num2str(height(cnt)/table_lengths(1)), ')'])
cnt  = cnt(cnt.Flag1 == 0,:);
table_lengths(length(table_lengths) + 1) = height(cnt);
disp(['Compiled Nuclei Table new row count after Flag 1 filter: ', num2str(height(cnt)),...
    ' (',num2str(height(cnt)/table_lengths(1)), ')'])
cnt  = cnt(cnt.Flag2 <  0.5,:);
table_lengths(length(table_lengths) + 1) = height(cnt);
disp(['Compiled Nuclei Table new row count after Flag 2 filter: ', num2str(height(cnt)),...
    ' (',num2str(height(cnt)/table_lengths(1) ), ')'])
cnt  = cnt(cnt.Flag3 <  0.5,:);
table_lengths(length(table_lengths) + 1) = height(cnt);
disp(['Compiled Nuclei Table new row count after Flag 3 filter: ', num2str(height(cnt)),...
    ' (',num2str(height(cnt)/table_lengths(1)), ')'])
cnt  = cnt(cnt.Flag4 <  0.5,:);
table_lengths(length(table_lengths) + 1) = height(cnt);
disp(['Compiled Nuclei Table new row count after Flag 4 filter: ', num2str(height(cnt)),...
    ' (',num2str(height(cnt)/table_lengths(1)), ')'])
cnt = cnt(cnt.Flag5 < 0.75,:);
table_lengths(length(table_lengths) + 1) = height(cnt);
disp(['Compiled Nuclei Table new row count after Flag 5 filter: ', num2str(height(cnt)),...
    ' (',num2str(height(cnt)/table_lengths(1)), ')'])
cnt  = cnt(cnt.Flag7 < 0.2,:);
table_lengths(length(table_lengths) + 1) = height(cnt);
disp(['Compiled Nuclei Table new row count after Flag 7 filter: ', num2str(height(cnt)),...
    ' (',num2str(height(cnt)/table_lengths(1)), ')'])

%%


MeanFluoAP = NaN(numFrames, numBins);
StdFluoAP = NaN(numFrames, numBins);
NumNucAP = zeros(numFrames, numBins);
DiffMeanFluoAP = NaN(numFrames, numBins);
DiffStdErrorAP = NaN(numFrames, numBins);
Frames = unique(cnt.Frame);
numFrames = length(Frames);
Times = [];
for i=1:length(Frames)
    Times(i) =  mean(cnt.Time(cnt.Frame == Frames(i)));
end
for i = 1:length(Frames)
    f = Frames(i);

    for APbin = min(cnt.APbinIdx):max(cnt.APbinIdx)

        sub_table = cnt((cnt.Frame == f) & (cnt.APbinIdx == APbin), :);
        if (height(sub_table) < 5)
            continue 
        end
        if (height(sub_table) >= 10) & (quantile_cutoff > 0)
            cutoff = quantile(sub_table.Fluo, quantile_cutoff);
            fluoVector = sub_table.Fluo(sub_table.Fluo >= cutoff);
        else
            fluoVector = sub_table.Fluo;
        end
        MeanFluo = mean(fluoVector);
        StdFluo = std(fluoVector);
        MeanFluoAP(i, APbin) = MeanFluo;
        StdFluoAP(i, APbin) = StdFluo;
        NumNucAP(i, APbin) = size(fluoVector, 1);
    end
    if sum(~isnan(MeanFluoAP(i,  min_fluo_pos_idx ))) > 0
        MinFluoAPidx = find(MeanFluoAP(i,  min_fluo_pos_idx ) == min(MeanFluoAP(i,  min_fluo_pos_idx )));
        if length(MinFluoAPidx) > 1
            MinFluoAPidx = MinFluoAPidx(end);
        end
        idx = min_fluo_pos_idx(MinFluoAPidx);
        DiffMeanFluoAP(i,:) = MeanFluoAP(i,:) - MeanFluoAP(i,idx); 
        stdErrorAPframe = StdFluoAP(i,:)./sqrt(NumNucAP(i,:));
        DiffStdErrorAP(i,:) = sqrt(stdErrorAPframe.^2 + stdErrorAPframe(idx).^2);
    end      
end


% 



save([DropboxFolder,filesep,Prefix,filesep,'MeanFluoAP.mat'],...
        'MeanFluoAP');

save([DropboxFolder,filesep,Prefix,filesep,'StdFluoAP.mat'],...
        'StdFluoAP');

save([DropboxFolder,filesep,Prefix,filesep,'NumNucAP.mat'],...
        'NumNucAP');
save([DropboxFolder,filesep,Prefix,filesep,'DiffMeanFluoAP.mat'],...
        'DiffMeanFluoAP');
save([DropboxFolder,filesep,Prefix,filesep,'DiffStdErrorAP.mat'],...
        'DiffStdErrorAP');
save([DropboxFolder,filesep,Prefix,filesep,'AllTimes.mat'],...
        'Times');

