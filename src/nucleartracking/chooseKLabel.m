function k = chooseKLabel(kLabel)

subt = false;
if islogical(kLabel)
    kLabel = uint8(kLabel) + 1;
    subt = true;
end
kMax = uint8(max(kLabel(:)));
stats = cell(1, kMax);
descriptorCell = cell(1, kMax);
stds = zeros(1, kMax);

for i = 1:kMax
    
    descriptor= 'Area';
    stats{i} = regionprops(kLabel==i, descriptor);
    descriptorCell{i} = [stats{i}.(descriptor)];
    stds(i) = std(descriptorCell{i});
    
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