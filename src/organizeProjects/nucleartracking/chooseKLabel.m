function k = chooseKLabel(kLabel, varargin)

subt = false; %i forgot what this flag is supposed to represent, but i think it's important
areaFilter = [];

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for k = 1:2:(numel(varargin)-1)
    if k ~= numel(varargin)
        eval([varargin{k} '=varargin{k+1};']);
    end
end

if islogical(kLabel)
    kLabel = uint8(kLabel) + 1;
    subt = true;
end

maxLabel = uint8(max(kLabel(:)));
stats = cell(1, maxLabel);
descriptorCell = cell(1, maxLabel);
stds = zeros(1, maxLabel);

for label = 1:maxLabel
    
    descriptor= 'Area';
    stats{label} = regionprops(kLabel==label, descriptor);
    descriptorCell{label} = [stats{label}.(descriptor)];
    stds(label) = std(descriptorCell{label});
    
end

%pick the label that minimizes the standard deviation 
%of our chosen descriptor. We'll decide this is the
%least weird label. Unless the standard deviation is zero. Then
% it's likely to be the background

stds(stds==0) = [];

if subt
    stds(stds==0) = nan;
end


[~, k] = min(stds);

if subt
    k = k - 1;
end



end