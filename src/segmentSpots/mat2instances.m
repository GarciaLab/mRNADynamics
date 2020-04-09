function [dataSet, attList] = mat2instances(attributes, name, data)

%attributes- cell array. include the class label last. 
%val - linear array of numeric values. one value for each feature

attList = cell2AttributeList(attributes);
nAtt = attList.size;

nRows = size(data, 1);
nCols = size(data, 2);

dataSet = weka.core.Instances(name,attList, nCols);
dataSet.setClassIndex(nAtt);

for row = 1:nRows
    
    dataRow = data(row, :);
    
    inst = mat2instance(attList, dataRow);

    dataSet.add(inst);

end

function list = cell2AttributeList(cell)
    list = java.util.ArrayList;
    for i = 1:numel(cell)
        list.add(weka.core.Attribute(cell{i}));
    end
end



end
