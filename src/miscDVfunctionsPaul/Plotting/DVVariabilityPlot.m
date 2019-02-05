% Plots variability in DV Gradient. Using data from full embryo is still a
% work in progress, but for now set 'im_type' to be any string that is not
% 'FullEmbryo'. Also, this function uses the downloaded function
% 'notBoxPlot,' make sure to cite it appropriately if you use this function
% (info should be in the same folder as this function).
% Paul Marchando
% 1-28-19

function [plot_data, stats] = DVVariabilityPlot(Prefix, para_values, fake_dv_comb, im_type)

%Get the folders, including the default Dropbox one
[SourcePath, FISHPath, DefaultDropboxFolder, DropboxFolder, MS2CodePath, PreProcPath,...
configValues, movieDatabasePath] = DetermineAllLocalFolders(Prefix);

%Determine division times
%Load the information about the nc from moviedatabase file
[Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
nc9, nc10, nc11, nc12, nc13, nc14, CF] = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder);

load([DropboxFolder '/' Prefix '/' 'Ellipses.mat'])

if strcmp(im_type, 'FullEmbryo') == 0

true_frames = para_values(:,1);
plot_data = cell(length(true_frames),4);
plot_data(:,1) = num2cell(true_frames);

for i = 1:length(true_frames)
    
    bin_width = Ellipses{true_frames(i)}(1,3) * 2;
    
    DV_position_vec = fake_dv_comb(fake_dv_comb(:,3) == true_frames(i),1);
    new_DV_position = para_values(i,3) - DV_position_vec;
    fluo_vec = fake_dv_comb(fake_dv_comb(:,3) == true_frames(i),2);
    
    
    bin_edges = min(new_DV_position):bin_width:max(new_DV_position);
    plot_data{i,3} = discretize(new_DV_position, bin_edges);
    plot_data{i,2} = bin_edges + (((bin_edges(1) + bin_edges(2)) / 2) - bin_edges(1));
    plot_data{i,4} = fluo_vec;
    
end

stats = cell(length(plot_data),1);

for m = 1:length(plot_data)
    
    temp_cell = cell(2,length(plot_data{m,2}));
    temp_cell(1,:) = num2cell(plot_data{m,2});
    
    for n = 1:max(plot_data{m,3})
        
        temp_cell{2,n} = plot_data{m,4}(plot_data{m,3} == n);
        
    end
    
    temp_stats = cell(max(plot_data{m,3}),1);
    
    for o = 1:max(plot_data{m,3})
        
        [~, temp_stats{o}] = notBoxPlot(temp_cell{2,o},temp_cell{1,o},'style','patch','jitter',5);
        
    end
    
    push_to_stats = zeros(length(temp_stats),1);
    
    for p = 1:length(temp_stats)
        
        push_to_stats(p) = ((2 * temp_stats{p}.sd) / temp_stats{p}.mu) * 100;
        
    end
    
    stats{m} = push_to_stats;
    xlabel('DV Position')
    ylabel('Fluorescence')
    title(['Frame: ' num2str(plot_data{m,1})])
    xlim([-300 300])
    ylim([0 2400])
    waitforbuttonpress
    clf
    
end

stats = [plot_data(:,1) stats];

else

load([DropboxFolder '/' Prefix '/' Prefix '_lin.mat'])
[working_dat] = CleanDVDataFullEmbryo(schnitzcells);
temp_ellipse = Ellipses{1,1};
radii_dat = [temp_ellipse(:,3); temp_ellipse(:,4)];
bin_width = mode(radii_dat) * 2;
DV_pos_vec = working_dat(:,2);
fluo_vec = working_dat(:,1);
bin_edges = min(DV_pos_vec):bin_width:max(DV_pos_vec);
bin_centers = bin_edges + (((bin_edges(1) + bin_edges(2)) / 2) - bin_edges(1));
bin_indices = discretize(DV_pos_vec,bin_edges);


end



end