function Projection = calculateProjection(ProjectionType, NSlices, HisSlices, varargin)

% Calculate projection for a nuclear channel
lowerSlice = 1;
upperSlice = NSlices;

if strcmpi(ProjectionType, 'middleprojection')
    lowerSlice = round(NSlices * .50);
    upperSlice =round(NSlices * .75);
elseif strcmpi(ProjectionType, 'midsumprojection')
    lowerSlice = round(NSlices * .33);
    upperSlice =round(NSlices * .66);
end
% 
% if ~isempty(varargin)
%     upperSlice = varargin{1};
%     lowerSlice = varargin{2};
% end

HisSlices = HisSlices(:, :, lowerSlice:upperSlice);
    

if strcmpi(ProjectionType, 'medianprojection')
    Projection = median(HisSlices, 3);
elseif strcmpi(ProjectionType, 'maxprojection') || strcmpi(ProjectionType, 'middleprojection')
    Projection = max(HisSlices, [], 3);
elseif strcmpi(ProjectionType, 'midsumprojection')
    Projection = sum(HisSlices, 3);
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