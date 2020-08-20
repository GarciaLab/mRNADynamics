function attributeCell = getAttributes(dataSet)

%get the attributes of a weka data set in a matlab readable cell
    attributeCell = {};
    nAtt = dataSet.numAttributes - 2;
    for i = 0:nAtt+1
        attributeCell{i+1} = char(dataSet.attribute(i).name());
    end

end