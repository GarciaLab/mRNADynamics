function [data, arffLoader] = loadArff(file, varargin)

balance = false;

for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end

arffLoader = weka.core.converters.ArffLoader; %this constructs an object of the arffloader class
arffLoader.setFile(javaObject('java.io.File',file)); %construct an arff file object
data= arffLoader.getDataSet;
data.setClassIndex(data.numAttributes - 1);

if balance
    data = balanceClasses(data);
end

end