function [FluoTimeTrace, FrameApproved, FluoZValue] = CalculateNuclearFluorescence(FluoMatrix, method)
% author: Gabriella Martini
% date created: 9/25/20
% date last modified: 9/26/20


FluoTimeTrace = [];
FrameApproved = [];
FluoZValue = [];
padding = 0;
if isempty(FluoMatrix)
    return;
end

Frames = size(FluoMatrix,1);

for f = 1:Frames
    FluoZTrace = FluoMatrix(f,2:end-1);
    if strcmp(method, 'max')
        FluoTimeTrace(f) = max(FluoZTrace);
        if ~isnan(FluoTimeTrace(f))
            FluoZValue(f) = find(FluoZTrace == FluoTimeTrace(f), 1);
        else
            FluoZValue(f) = nan;
        end
        if isnan(FluoTimeTrace(f))
            FrameApproved(f) = 0;
        elseif sum(FluoZTrace(padding+1:end-2) == FluoTimeTrace(f)) == 0
            FrameApproved(f) = 0;
        else
            FrameApproved(f) =0;
        end
        
    elseif strcmp(method, 'median')
        FluoTimeTrace(f) = median(FluoZTrace);
        if ~isnan(FluoTimeTrace(f))
            FrameApproved(f) = 1;
            FluoZ = find(FluoZTrace == FluoTimeTrace(f), 1);
            if ~isempty(FluoZ)
                FluoZValue(f) = FluoZ;
            else
                FluoZValue(f) = nan;
            end
        else
            FrameApproved(f) = 0;
             FluoZValue(f) = nan;
        end
        
    elseif strcmp(method, 'midmedian')
        MaxFluo = max(FluoZTrace);
        MaxFluoPos = find(FluoZTrace == MaxFluo, 1);
        FluoTimeTrace(f) = median(FluoZTrace(max([1, MaxFluoPos-5]):min([length(FluoZTrace), MaxFluoPos+5])));
        if ~isnan(FluoTimeTrace(f))
            FrameApproved(f) = 1;
            FluoZ = find(FluoZTrace == FluoTimeTrace(f), 1);
            if ~isempty(FluoZ)
                FluoZValue(f) = FluoZ;
            else
                FluoZValue(f) = nan;
            end
        else
            FrameApproved(f) = 0;
             FluoZValue(f) = nan;
        end
    elseif strcmp(method, 'meansmooth')
        smoothedFluoZTrace = smoothdata(FluoZTrace, 'movmean', 5);
        FluoTimeTrace(f) = max(smoothedFluoZTrace);
        if ~isnan(FluoTimeTrace(f))
            FrameApproved(f) = 1;
            FluoZ = find(smoothedFluoZTrace == FluoTimeTrace(f), 1);
            if ~isempty(FluoZ)
                FluoZValue(f) = FluoZ;
            else
                FluoZValue(f) = nan;
            end
        else
            FrameApproved(f) = 0;
            FluoZValue(f) = nan;
        end
    elseif strcmp(method, 'gaussiansmooth')
        smoothedFluoZTrace = smoothdata(FluoZTrace, 'gaussian', 5);
        FluoTimeTrace(f) = max(smoothedFluoZTrace);
        if ~isnan(FluoTimeTrace(f))
            FrameApproved(f) = 1;
            FluoZ = find(smoothedFluoZTrace == FluoTimeTrace(f), 1);
            if ~isempty(FluoZ)
                FluoZValue(f) = FluoZ;
            else
                FluoZValue(f) = nan;
            end
        else
            FrameApproved(f) = 0;
            FluoZValue(f) = nan;
        end
    else
        error('invalid choice of fluo calculation method')
    end
    
    
end