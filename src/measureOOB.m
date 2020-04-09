function measureOOB(classifier, trainingData)

figure; 
oob = [];
nTrees = 0:1:30;
k = 0;
for t = nTrees
    k = k + 1;
    classifier = weka.classifiers.trees.RandomForest;
     options = {'-I', '64', '-K', '2', '-S', '-1650757608', '-depth', num2str(t),'-O'};
    classifier.setOptions(options);
    classifier.buildClassifier(trainingData);
    oob = [oob, classifier.measureOutOfBagError];
    plot(gca, nTrees(1:k), oob);
    drawnow;
end

end