% Cleans up data out of the pipeline by splitting nuclei and removing NaNs
% for use with DVpolyfit
% Paul Marchando
% 11-28-2018

function [working_dat_final] = CleanDVData(schnitzcells,CompiledNuclei)

% Selects relevant info from CompiledNuclei and schnitzcells, converts each
% struct into a cell array, splits up linked traces so each row corresponds
% to a single nucleus in a single frame, and then removes NaNs and removes
% the dummy frames from the beginning and end of each z-stack.

tic
schnitz = struct2cell(schnitzcells);
schnitz = transpose(schnitz(:,:));

comp_nuc = struct2cell(CompiledNuclei);
comp_nuc = transpose(comp_nuc(:,:));

Fluo = schnitz(:,9);
DVpos = schnitz(:,12);
frame = schnitz(:,4);
ypos = comp_nuc(:,6);

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
        new_vec{temp_ind(1),4} = ypos{i};
        
    else
        
        temp_dat_fluo = Fluo{i};
        temp_dat_DVpos = DVpos{i};
        temp_dat_frame = frame{i};
        temp_dat_ypos = ypos{i};
        temp_size = size(temp_dat_fluo);
        row_dist = ones(1,temp_size(1));
        temp_dat_fluo = mat2cell(temp_dat_fluo,row_dist);
        temp_dat_DVpos = mat2cell(transpose(temp_dat_DVpos),row_dist);
        temp_dat_frame = mat2cell(temp_dat_frame,row_dist);
        temp_dat_ypos = mat2cell(transpose(temp_dat_ypos),row_dist);
        
        new_vec(temp_ind(1):temp_ind(1) + temp_size(1) - 1,1) = temp_dat_fluo;
        new_vec(temp_ind(1):temp_ind(1) + temp_size(1) - 1,2) = temp_dat_DVpos;
        new_vec(temp_ind(1):temp_ind(1) + temp_size(1) - 1,3) = temp_dat_frame;
        new_vec(temp_ind(1):temp_ind(1) + temp_size(1) - 1,4) = temp_dat_ypos;
        
    end
    
end

fluo_mean = cellfun(@mean,new_vec(:,1));
filter_NAN = isnan(fluo_mean);
working_dat_final = new_vec(~filter_NAN,:);
temp_final = working_dat_final(:,1);
temp_final = cellfun(@(x) x(2:end-1),temp_final(:,1),'un',0);
working_dat_final(:,1) = temp_final;
toc

end
