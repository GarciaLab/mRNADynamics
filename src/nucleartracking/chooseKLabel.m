function k = chooseKLabel(kLabel)

subt = false; %i forgot what this flag is supposed to represent, but i think it's important


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
%least weird label. 

if subt
    stds(stds==0) = nan;
end


[~, k] = min(stds);

if subt
    k = k - 1;
end



end