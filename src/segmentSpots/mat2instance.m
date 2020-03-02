function [dataSet, attList] = mat2instances(attributes, name, data)

%attributes- cell array. include the class label last. 
%val - linear array of numeric values. one value for each feature

attList = cell2AttributeList(attributes);

nAtt = attList.size;
dataSet = weka.core.Instances(name,attList, nAtt);
dataSet.setClassIndex(end);

nRows = size(data, 1);
nCols = size(data, 2);

for row = 1:nRows
    
    dataRow = data(row, :);
    
    inst = mat2instance(attList, dataRow);

    dataSet = weka.core.Instances(name,attList, nCols);

end

function list = cell2AttributeList(cell)
    list = java.util.ArrayList;
    for i = 1:numel(cell)
        list.add(weka.core.Attribute(cell{i}));
    end
end

function mat2instance(attList, val)

    nAtt = attList.size;

    inst=weka.core.DenseInstance(nAtt+1);

    for i = 0:nAtt
        for j = 1:numel(val)
            inst.setValue(attList.get(i), val);
        end
    end
    
end


end
