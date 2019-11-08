% Script to sample MCP levels from image stacks using nucleus center
% locations

function MCPLevelsSummary = MCPLevelAnalysis(Prefix)

% get path info, spot channel, and movie database
[rawDataPath,ProcPath,DropboxFolder,MS2CodePath, PreProcPath,...
    rawDataFolder, Prefix, ExperimentType,Channel1,Channel2,OutputFolder,...
    Channel3, spotChannels, MovieDataBaseFolder, movieDatabase]= readMovieDatabase(Prefix);

% get nc timing info
[~, ~, ~, ~, ~,~,~, ~, ~, ~, ~, ~, ~, nc9, nc10, nc11, nc12, nc13, nc14] = ...
    getExperimentDataFromMovieDatabase(Prefix, movieDatabase);

% load FrameInfo
load([DropboxFolder, filesep, Prefix, filesep, 'FrameInfo.mat'], 'FrameInfo');

% load Ellipse info
load([DropboxFolder, filesep, Prefix, filesep, 'Ellipses.mat'], 'Ellipses');

% load particles set if it exists
if exist([DropboxFolder, filesep, Prefix, filesep, 'Particles.mat'])
    load([DropboxFolder, filesep, Prefix, filesep, 'Particles.mat'],'Particles')
end

% generate time vector normalized to nc14
TimeVec = [FrameInfo.Time];
TimeVec = TimeVec - TimeVec(nc14);

%%
% shape parameters
FrameVec = 1:numel(FrameInfo);
PixelSize = FrameInfo(1).PixelSize;
MinRad = round(1/PixelSize);
IntRad = ceil(.5/PixelSize);
zDim = FrameInfo(1).NumberSlices;
xDim = FrameInfo(1).PixelsPerLine;
yDim = FrameInfo(1).LinesPerFrame;
% initialize
MCPLevelsData = struct;
% check to for "ch" suffix
chFiles = dir([PreProcPath filesep Prefix '\*ch01.tif']);
multiChannel = ~isempty(chFiles); 
tic
parfor f = FrameVec    
    % load data
    Frame = FrameVec(f);    
    ImStack = zeros(yDim,xDim,zDim-2);
    for z = 2:zDim-1
        if multiChannel
            fname = [Prefix filesep Prefix '_' sprintf('%03d',Frame) '_z' sprintf('%02d',z) '_ch' sprintf('%02d',spotChannels(1)) '.tif'];
        else
            fname = [Prefix filesep Prefix '_' sprintf('%03d',Frame) '_z' sprintf('%02d',z) '.tif'];
        end
        im = imread([PreProcPath filesep fname]);
        ImStack(:,:,z-1) = im;
    end
    % process stack
    med_frame = median(ImStack,3);
    med_frame_sm = imgaussfilt(med_frame,IntRad);
    min_frame = min(ImStack,[],3);
    min_frame_sm = imgaussfilt(min_frame,IntRad);
    % get list of nucleus coordinates and linearize
    nc_loc_array = round(Ellipses{f}(:,1:2));
    keep_flags_x = nc_loc_array(:,1) > MinRad & nc_loc_array(:,1) < xDim-MinRad;
    keep_flags_y = nc_loc_array(:,2) > MinRad & nc_loc_array(:,2) < yDim-MinRad;
    nc_loc_ind = sub2ind([yDim,xDim],nc_loc_array(keep_flags_x&keep_flags_y,2),nc_loc_array(keep_flags_x&keep_flags_y,1));    
    % pull samples     
    med_mcp_vec = med_frame_sm(nc_loc_ind);
    min_mcp_vec = min_frame_sm(nc_loc_ind);
    % save results
    MCPLevelsData(f).med_mcp_vec = med_mcp_vec;
    MCPLevelsData(f).min_mcp_vec = min_mcp_vec;
    MCPLevelsData(f).nc_indices = nc_loc_ind;
    MCPLevelsData(f).yxz = [yDim, xDim, zDim];
    MCPLevelsData(f).time = TimeVec(f);    
end
toc
%% analyze results
% generate flags indicating middle 68%. Removing outliers to gurad against
% influence of active spots and other random shit
for f= FrameVec
    [~, med_rank] = sort(MCPLevelsData(f).med_mcp_vec);
    med_prctile = med_rank / numel(med_rank);
    MCPLevelsData(f).med_analysis_flags = med_prctile >=  .16 & med_prctile <= .84;
    [~, min_rank] = sort(MCPLevelsData(f).min_mcp_vec);
    min_prctile = min_rank / numel(min_rank);
    MCPLevelsData(f).min_analysis_flags = min_prctile >=  .16 & min_prctile <= .84;
end

% Bootstrap to estimate mean and standard error
n_boots = 100;
min_mcp_array = NaN(n_boots,numel(FrameVec));
med_mcp_array = NaN(n_boots,numel(FrameVec));
for f = FrameVec
    med_analysis_indices = find(MCPLevelsData(f).med_analysis_flags);
    med_mcp_vec = MCPLevelsData(f).med_mcp_vec;
    min_analysis_indices = find(MCPLevelsData(f).min_analysis_flags);
    min_mcp_vec = MCPLevelsData(f).min_mcp_vec;
    for n = 1:n_boots
        med_ids = randsample(med_analysis_indices,numel(med_analysis_indices),true);
        min_ids = randsample(min_analysis_indices,numel(min_analysis_indices),true);
        med_mcp_array(n,f) = nanmean(med_mcp_vec(med_ids));
        min_mcp_array(n,f) = nanmean(min_mcp_vec(min_ids));
    end
end

% calculate mean and se
MCPLevelsSummary.med_mcp_mean = nanmean(med_mcp_array);
MCPLevelsSummary.med_mcp_se = nanstd(med_mcp_array);
MCPLevelsSummary.min_mcp_mean = nanmean(min_mcp_array);
MCPLevelsSummary.min_mcp_se = nanstd(min_mcp_array);

% save
save([DropboxFolder, filesep, Prefix, filesep, 'MCPLevelsData.mat'], 'MCPLevelsData','MCPLevelsSummary');
