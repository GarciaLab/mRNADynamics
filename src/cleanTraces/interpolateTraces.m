function [nucleus_struct] = interpolateTraces(nucleus_struct, minDP, ...
    TresInterp, num_outputs, nan_buffer)
%CLEANTRACE Summary of this function goes here
%   Detailed explanation goes here    
    interpGrid = 0:TresInterp:100*60;
    %%% Cleaning Parameters
    jump_thresholds = zeros(1, num_outputs);
    jump_thresholds3D = zeros(1, num_outputs);
    for cidx = 1:num_outputs
        fluo_together = [nucleus_struct.fluo];
        big_jump = prctile(fluo_together(cidx,:),99);
        jump_thresholds(cidx) = big_jump/1.5; % this should be a conservative threshold for single time step increase
        fluo3D_together = [nucleus_struct.fluo3D];
        big_jump3D = prctile(fluo3D_together(cidx,:),99);
        jump_thresholds3D(cidx) = big_jump3D/1.5;
    end
    for i = 1:length(nucleus_struct) 
        temp = nucleus_struct(i);
        for cidx = 1:num_outputs
            trace = temp.fluo(cidx,:);
            trace3D = temp.fluo3D(cidx,:);
            pt_time = temp.time;      

%             if sum(~isnan(trace)) == 0 
%                 t_start = 0;
%                 t_stop = -1;
%             else
%                 t_start = interpGrid(find(interpGrid>=min(pt_time(~isnan(trace))),1));
%                 t_stop = interpGrid(find(interpGrid<=max(pt_time(~isnan(trace))),1,'last'));
%             end
            t_start = interpGrid(find(interpGrid>=min(pt_time),1));
            t_stop = interpGrid(find(interpGrid<=max(pt_time),1,'last'));            
            time_interp = t_start:TresInterp:t_stop;
            
            if sum(~isnan(trace)) > 1
                [trace_interp, quality_flag] = processTrace(trace,nan_buffer,pt_time,...
                    time_interp,jump_thresholds(cidx),minDP);
            else
                trace_interp = nan(1, length(time_interp));
                quality_flag = 0;
            end            
            if sum(~isnan(trace3D)) > 1
                [trace3D_interp,~] = processTrace(trace3D,nan_buffer,pt_time,...
                    time_interp,jump_thresholds3D(cidx),minDP);
            else
                trace3D_interp = nan(1, length(time_interp));
            end
            
            nucleus_struct(i).fluo_interp(cidx,:) = trace_interp; 
            nucleus_struct(i).fluo3D_interp(cidx,:) = trace3D_interp;
            nucleus_struct(i).time_interp(cidx,:) = time_interp;
            % indicates whether trace suitable for inference
            nucleus_struct(i).inference_flag(cidx,:) = quality_flag;
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
                    nucleus_struct(i).([interp_fields{j} '_interp'])(cidx,:) = init_vec;                
                else
                    try 
                        nucleus_struct(i).([interp_fields{j} '_interp'])(cidx,:) = interp1(init_time,init_vec,time_interp);
                    catch
                        error('asfa')
                    end
                end
            end
        end     
    end
end

