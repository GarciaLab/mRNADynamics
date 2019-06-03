% Cleans up data out of the pipeline by splitting nuclei from the full
% embryo and removing NaNs for use with Gaussian fitting
% Paul Marchando
% 1-04-2019

function [working_dat_final] = CleanDVDataFullEmbryo(schnitzcells)

tic
schnitz = struct2cell(schnitzcells);
schnitz = transpose(schnitz(:,:));

Fluo = schnitz(:,9);
DVpos = schnitz(:,12);
frame = schnitz(:,1);
APpos = schnitz(:,11);

[fluo_size, ~] = cellfun(@size,Fluo);
filter_size = fluo_size(:) == 1;
mult_vec = fluo_size(fluo_size ~= 1);
additional_rows = sum(mult_vec) - length(mult_vec);
new_vec = cell(length(Fluo) + additional_rows,4);

for i = 1:length(Fluo)
    
    temp_ind = find(cellfun('isempty',new_vec(:,1)));
    
    if filter_size(i) == 1
        
        new_vec{temp_ind(1),1} = Fluo{i};
        new_vec{temp_ind(1),2} = DVpos{i};
        new_vec{temp_ind(1),3} = frame{i};
        new_vec{temp_ind(1),4} = APpos{i};
        
    else
        
        temp_dat_fluo = Fluo{i};
        temp_dat_DVpos = DVpos{i};
        temp_dat_frame = frame{i};
        temp_dat_APpos = APpos{i};
        temp_size = size(temp_dat_fluo);
        row_dist = ones(1,temp_size(1));
        temp_dat_fluo = mat2cell(temp_dat_fluo,row_dist);
        temp_dat_DVpos = mat2cell(transpose(temp_dat_DVpos),row_dist);
        temp_dat_frame = mat2cell(transpose(temp_dat_frame),row_dist);
        temp_dat_APpos = mat2cell(transpose(temp_dat_APpos),row_dist);
        
        new_vec(temp_ind(1):temp_ind(1) + temp_size(1) - 1,1) = temp_dat_fluo;
        new_vec(temp_ind(1):temp_ind(1) + temp_size(1) - 1,2) = temp_dat_DVpos;
        new_vec(temp_ind(1):temp_ind(1) + temp_size(1) - 1,3) = temp_dat_frame;
        new_vec(temp_ind(1):temp_ind(1) + temp_size(1) - 1,4) = temp_dat_APpos;
        
    end
    
end

fluo_mean = cellfun(@mean,new_vec(:,1));
filter_NAN = isnan(fluo_mean);
working_dat_final = new_vec(~filter_NAN,:);
temp_final = working_dat_final(:,1);
temp_final = cellfun(@(x) max(x),temp_final(:,1),'un',0);
working_dat_final(:,1) = temp_final;
working_dat_final = cell2mat(working_dat_final);
working_dat_final = working_dat_final(working_dat_final(:,3) == 1, :);
working_dat_final(:,3) = [];

toc

end
