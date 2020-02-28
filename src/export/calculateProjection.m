function Projection = calculateProjection(ProjectionType, NSlices, HisSlices, varargin)
% Calculate the projection (either Maximum or Median)

if strcmpi(ProjectionType, 'medianprojection')
    Projection = median(HisSlices, 3);
elseif strcmpi(ProjectionType, 'middleprojection')
    Projection = max(HisSlices(:, :, round(NSlices * .50):round(NSlices * .75)), [], 3);
elseif strcmpi(ProjectionType, 'maxprojection')
    Projection = max(HisSlices, [], 3);
else
    SortedHisSlices = sort(HisSlices, 3, 'descend');
    if ~isempty(varargin)
        max_custom = varargin{1};
        min_custom = varargin{2};
        if length(size(SortedHisSlices)) > 2
            Projection = mean(SortedHisSlices(:, :, max_custom:min_custom), 3);
        else
            Projection = SortedHisSlices;
        end
    else
        ProjectionBounds = strsplit(ProjectionType, ':');
        max_custom = str2double(ProjectionBounds{2});
        min_custom = str2double(ProjectionBounds{3});
        Projection = mean(SortedHisSlices(:, :, max_custom:min_custom), 3);
    end
end

end