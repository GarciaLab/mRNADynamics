function [trace1_interp, approved] = processTrace(trace1,nan_buffer,pt_time,...
    time_interp,jump_thresh,minDP)
%PROCESSTRACE Summary of this function goes here
%   Detailed explanation goes here
    
    % edited version of trace1
    trace1_clean = trace1;
    tr_start = find(~isnan(trace1_clean),1);
    tr_end = find(~isnan(trace1_clean),1,'last');
    trace1_clean = trace1_clean(tr_start:tr_end);
    pt_time = pt_time(tr_start:tr_end);
    %Null assumption is that all clusters of 6 or more NaNs are 0s. Smaller
    %clusters are assumed to have been missed nonzero dps
    trace1_nans = isnan(trace1_clean);      
    %Look for clusters of 6 or more NaNs
    kernel = ones(1, nan_buffer * 2 + 1);
    tn_conv = conv(trace1_nans,kernel,'same');
    % finds nan points surrounded by at least nan_buffer nan points
    z_ids = find(tn_conv==length(kernel)); 
    trace1_clean(z_ids) = 0; % set clusters to zeros    
    trace1_clean(trace1_clean<0) = 0; % deal with negative values    
    % find single dp "blips". These will be replaced via interpolation
    tr_dd1 = abs([0 diff(diff(trace1_clean)) 0]);
    trace1_clean(tr_dd1>2*jump_thresh) = NaN;    

    % interpolate remaining NaNs    
    query_points1 = pt_time(isnan(trace1_clean));
    interp_t1 = pt_time(~isnan(trace1_clean));
    interp_f1 = trace1_clean(~isnan(trace1_clean));   
    new_f1 = interp1(interp_t1,interp_f1,query_points1);    
    trace1_clean(ismember(pt_time,query_points1)) = new_f1;        

    % Interpolate to standardize spacing   
    trace1_interp = interp1(pt_time,trace1_clean,time_interp); 

    % says if trace is good
    approved =  standardQualityControl(trace1,trace1_clean,minDP,jump_thresh);
end

