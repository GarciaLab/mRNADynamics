function [nucleus_struct] = interpolateTraces(nucleus_struct, minDP, ...
    TresInterp, num_outputs)
%CLEANTRACE Summary of this function goes here
%   Detailed explanation goes here

    
    if ~exist('TresInterp')
        TresInterp = round(med_time);
    end
    interpGrid = 0:TresInterp:50*60;
    %%% Cleaning Parameters
    jump_thresholds = zeros(1, num_outputs);
    for cidx = 1:num_outputs
        big_jump = prctile([nucleus_struct.(['fluo' num2str(cidx)])],99);
        jump_thresholds(cidx) = big_jump/1.5; % this should be a conservative threshold for single time step increase
    end
    for i = 1:length(nucleus_struct) 
        temp = nucleus_struct(i);
        for cidx = 1:num_outputs
            trace1 = temp.(['fluo' num2str(cidx)]); %Load full trace, including intervening NaN's    
            pt_time = temp.time;      
            quality_flag = 1; % indicates whether trace suitable for inference

            if sum(~isnan(trace1)) == 0 
                t_start = 0;
                t_stop = -1;
            else
                t_start = interpGrid(find(interpGrid>=min(pt_time(~isnan(trace1))),1));
                t_stop = interpGrid(find(interpGrid<=max(pt_time(~isnan(trace1))),1,'last'));
            end
            time_interp = t_start:TresInterp:t_stop;

            if sum(~isnan(trace1)) < minDP
                quality_flag = 0;
                trace1_interp = NaN(size(time_interp));
                time_interp = NaN(size(time_interp));
            else
                %Null assumption is that all clusters of 6 or more NaNs are 0s. Smaller
                %clusters are assumed to have been missed nonzero dps
                trace1_nans = isnan(trace1);      
                %Look for clusters of 6 or more NaNs
                kernel = [1,1,1,1,1];
                tn_conv = conv(kernel,trace1_nans);
                tn_conv = tn_conv(3:end-2);
                z_ids = find(tn_conv==5);
                z_ids = unique([z_ids-1 z_ids z_ids+1]); % get set of z_ids    
                trace1(z_ids) = 0; % set clusters to zeros    
                trace1(trace1<0) = 0; % deal with negative values    
                % find single dp "blips". These will be replaced via interpolation
                tr_dd1 = abs([0 diff(diff(trace1)) 0]);
                trace1(tr_dd1>2*jump_thresholds(cidx)) = NaN;    

                % interpolate remaining NaNs    
                query_points1 = pt_time(isnan(trace1));
                interp_t1 = pt_time(~isnan(trace1));
                interp_f1 = trace1(~isnan(trace1));

                new_f1 = interp1(interp_t1,interp_f1,query_points1);  

                trace1(ismember(pt_time,query_points1)) = new_f1;        

                %%% flag traces with unreasonably large rises or falls    
                tr_d1 = diff(trace1);    
                if max(abs(tr_d1)) >= jump_thresholds(cidx)         
                    quality_flag = 0;        
                end

                % Interpolate to standardize spacing    
                trace1_interp = interp1(pt_time,trace1,time_interp);    
            end
            nucleus_struct(i).(['fluo_interp' num2str(cidx)]) = trace1_interp;    
            nucleus_struct(i).(['time_interp' num2str(cidx)]) = time_interp;
            nucleus_struct(i).(['inference_flag' num2str(cidx)]) = quality_flag;
        end
        nucleus_struct(i).TresInterp = TresInterp;
        
        interp_fields = {'xPos','yPos','ap_vector'};
        % interpolate other vector fields
        for j = 1:length(interp_fields)
            field_string = interp_fields{j};
            if isfield(nucleus_struct,field_string)
                init_vec = temp.(field_string);
                init_time = temp.time;   
                if numel(init_time) < 2 || numel(time_interp) < 2
                    nucleus_struct(i).([interp_fields{j} '_interp']) = init_vec;                
                else
                    nucleus_struct(i).([interp_fields{j} '_interp']) = interp1(init_time,init_vec,time_interp);
                end
            end
        end     
        
    end
end

