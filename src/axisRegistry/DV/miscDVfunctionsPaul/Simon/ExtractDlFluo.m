function FluoTimeTrace = extractDorsalFluo(FluoMatrix, Thresh)

% Arguments:
% FLUOMATRIX is a txz matrix of fluorescence values where t is frames and z
% is Z slice. It corresponds to the 'schnitzcells(i).Fluo' field in the schnitzcells struct
% THRESH corresponds to the quadratic coefficient in the parabola we fit to
% fluo(z). .5 works. 

% Output 
% FLUOTIMETRACE is a tx1 vector of nuclear fluorescence values.

% it classifies each frame in 'schnitzcells(i).Fluo' in three different categores,
% which represent the three different ways the Dl nuclear fluorescence compares to the
% surrounding cytoplasm:

% Case 1: the nuclear fluorescence is higher than the cytoplasm. In this
% case a line of y = fluo(z) looks kind of like a parabola pointing up. In this
% case we calculate max(fluo(z)). This is basically what we do for any
% other nuclear protein.

%Case 2: the nuclear fluorescence is lower than the cytoplasm. In this
% case a line of y = fluo(z) looks like a parabola pointing down. In this
% case we calculate min(fluo(z)). 

%Case 3: the nuclear fluorescence is comparable to the cytoplasm. In this
% case a line of y = fluo(z) looks roughly flat. In this
% case we calculate the median(fluo(z)). 

% The way a fluo trace is assigned to one of these categories is by
% fitting a parabola such that fluo(z) = az^2 + bz + c.
% the first coefficient 'a' indicates the orientation of the parabola and
% is the one that we apply 'THRESH' to.

% Contact: Simon Alamos simon.alamos@berkeley.edu
%          Jiaxi Zhao jiaxi.zhao@berkeley.edu

%% Fit fluo(z) to a parabola, plot the data and the fit

%FluoMatrix = schnitzcells(1086).Fluo;

FluoTimeTrace = [];

if isempty(FluoMatrix)
    return;
end

Frames = size(FluoMatrix,1);

for f = 1:Frames
    FluoZTrace = FluoMatrix(f,2:end-1);
    Xvals = [1:length(FluoZTrace)];
    coefficients = polyfit(Xvals,FluoZTrace,2);    
    % This section here is for visualization/debugging
    %\/\/\/\/\/\/\/
%     figure(1)
%     plot(Xvals,FluoZTrace,'o','MarkerEdgeColor','none','MarkerFaceColor',[1 0.7 0.7])
%     hold on
%     fit = polyval(coefficients,Xvals);
%     plot(Xvals,fit,'LineWidth',1.5,'Color',[0.7 0.7 1])
%     title(['frame ' num2str(f)])
%     xlabel('z slice')
%     ylabel('nuclear fluorescence (a.u)')
%     legend('nuclear fluorescence data','fit')
%     hold off
%     waitforbuttonpress
    %\/\/\/\/\/\/\/


    % Classify the trace and calculate the fluo
    %if the a coefficient in the model fluo(z) = az^2 + bz + c is negative, the
    %parabola opens up. If its positive, it opens down. If it's slightly
    %negative or slightly positive it's a very shallow parabola.
    if ~isnan(coefficients(1))
        if coefficients(1) < -Thresh
            FluoFrame = nanmax(FluoZTrace);
        elseif coefficients(1) > Thresh
            FluoFrame = nanmin(FluoZTrace);
        elseif -Thresh <= coefficients(1) || coefficients(1) <= Thresh
            FluoFrame = nanmedian(FluoZTrace);
        end
        FluoTimeTrace(f) = FluoFrame;
    else
        FluoTimeTrace(f) = nan;
    % visualization/debugging
    % \/\/\/\/\/\/
%     figure(2)
%     plot(FluoTimeTrace)
    % \/\/\/\/\/\/ 
    end

end
