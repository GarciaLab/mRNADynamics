% script to explore spot movemnt characteristics
clear
close all
dataPath = 'D:\Data\Nick\projects\spot_tracking\dat\';
figPath = '../fig/';
mkdir(figPath)
% load data
load([dataPath 'Dl-Ven_sna-mCherry_201909026.mat'])
% find time periods with 30+ contiguous non-missing frames
min_len = 30;
iter = 1;
% iterate through structure
for i = 1:numel(nucleus_struct)
    fluo = nucleus_struct(i).fluo;
    time = nucleus_struct(i).time;
    x = nucleus_struct(i).xPosParticle';
    y = nucleus_struct(i).yPosParticle';
    z = nucleus_struct(i).zPosParticle';
    % find continuous obs periods    
    nan_frames = unique([1 find(isnan(fluo)) numel(fluo)]);
    gap_vec = diff(nan_frames);
    select_indices = find(gap_vec > min_len);
    % iterate through long active periods
    for j = 1:numel(select_indices)
        start = nan_frames(select_indices(j));
        indices = start+1:start+min_len;
        % concatenate
        xyz = [x(indices) y(indices) z(indices)];
        if iter == 1
            xyz_snip_array = xyz;
            fluo_snip_array = fluo(indices)';
            time_snip_array = time(indices)';
        else
            xyz_snip_array = cat(3,xyz,xyz_snip_array);
            fluo_snip_array = [fluo_snip_array fluo(indices)'];
            time_snip_array = [time_snip_array time(indices)'];
        end
        iter = iter + 1;
    end
end


%%
test_pos = [xyz_snip_array(:,1,500),xyz_snip_array(:,2,500)];
test_pos_sm = [imgaussfilt(xyz_snip_array(:,1,500),1),imgaussfilt(xyz_snip_array(:,2,500),1)];

close all
plot(test_pos(:,1),test_pos(:,2))
hold on
plot(test_pos_sm(:,1),test_pos_sm(:,2))
% kalmanFilter = configureKalmanFilter('ConstantVelocity', test_pos(1,:),...
%          20*[1,1], [1 1], 30);
% 
% kal_pos_array = NaN(size(test_pos));
% 
% kal_pos_array(1,:) = correct(kalmanFilter, test_pos(1,:));
% for i = 2:size(test_pos,1)
%     curr_pos = test_pos(i,:);
%     predict(kalmanFilter);
%     kal_pos_array(i,:) = correct(kalmanFilter, curr_pos);
% end
% 
% 
% close all
% plot(kal_pos_array(:,1),kal_pos_array(:,2))
% hold on
% plot(test_pos(:,1),test_pos(:,2))