function approved = standardQualityControl(trace1,trace1_clean,minDP,jump_thresh)
%STANDARDQUALITYCONTROL Summary of this function goes here
%   Detailed explanation goes here

    approved = true;
    
    %%% flag traces that are too small
    if sum(~isnan(trace1)) < minDP
        approved = false;
    end
    
    %%% flag traces with unreasonably large rises or falls    
    tr_d1 = diff(trace1_clean);    
    if max(abs(tr_d1)) >= jump_thresh        
        approved = false;       
    end
end

