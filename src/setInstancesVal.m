function dataSet = setInstancesVal(dataSet, instInd, attInd, val)

%in the arff file, instInd<->rows and attInd<->cols
% if dataSet.numInstances < instInd
%     
%     
%     attributes = {};
%     nAtt = trainingData.numAttributes()-2;
%     for i = 0:nAtt
%         attributes{i+1} = char(trainingData.attribute(i).name());
%     end
%     
%     inst =  weka.core.Instance(numel(attributes));
% end
inst = dataSet.instance(instInd);

inst.setValue(attInd, val);

dataSet.set(instInd, inst);

end