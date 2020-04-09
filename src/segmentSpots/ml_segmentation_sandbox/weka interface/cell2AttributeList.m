function list = cell2AttributeList(cell)
    list = java.util.ArrayList;
    for i = 1:numel(cell)
        list.add(weka.core.Attribute(cell{i}));
    end
end