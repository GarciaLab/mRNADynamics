function [cleanedAttributes, rejectedAttributes,...
    keepIndices, removeIndices] = validateAttributes(attributes, dim)

%pass dummy images through "filterImage" to see if they can be successfully
%filtered. 

cleanedAttributes = {};
rejectedAttributes = {};
keepIndices = [];
removeIndices = [];

if dim == 2
    sampleImage = zeros(256, 256);
elseif dim == 3
    sampleImage = zeros(256, 256, 8);
end

numAttributes = numel(attributes);

for i = 1:numAttributes
    
    att = attributes{i};   
    
    if ~strcmpi(att, 'original') && ~strcmpi(att, 'class')
        [~, successFlag]  = filterAttribute(att, sampleImage);
    else
        successFlag = true;
    end

    if successFlag
        cleanedAttributes = [cleanedAttributes, att];
        keepIndices = [keepIndices, i];
    else
        rejectedAttributes = [rejectedAttributes, att];
        removeIndices = [removeIndices, i];
    end
    
end

keepIndices = keepIndices - 1;
removeIndices = removeIndices - 1;


end