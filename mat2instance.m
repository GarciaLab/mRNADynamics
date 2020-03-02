function mat2instance(attList, val)

    nAtt = attList.size;

    inst=weka.core.DenseInstance(nAtt+1);

    for i = 0:nAtt
        inst.setValue(attList.get(i), val(i));
    end
    
end

row = 1;
inst = isTrainingSet.instance(row);
col = 4;
val = 6;
inst.setValue(col, val)
isTrainingSet.set(row, inst)