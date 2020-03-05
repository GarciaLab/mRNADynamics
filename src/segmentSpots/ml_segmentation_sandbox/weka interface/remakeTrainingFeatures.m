function trainingMat = remakeTrainingFeatures(trainingMat, cleanedAttributes)

attributes = cleanedAttributes;


numAttributes = numel(classifier.PredictorNames);
testMatrix = zeros(numInstances, numAttributes);

usedFeatures = {};
ignoredFeatures = {};

for i = 1:numAttributes
    
    att = attributes{i};
    %     disp(['Generating feature ', num2str((i)), '/', num2str(numAttributes) , ': ', att]);
    
    
    if i > 1
        [filteredIm, successFlag]  = filterAttribute(att, im);
    else
        filteredIm = im;
        successFlag = true;
    end
    
    if successFlag
        testMatrix(:,i) = filteredIm(:);
        usedFeatures = [usedFeatures, att];
    else
        testMatrix(:,i) = nan(numel(im),1);
        ignoredFeatures = [ignoredFeatures, att];
    end
    
end

end