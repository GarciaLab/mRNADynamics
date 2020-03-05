function mat2instance(attList, val)

    nAtt = attList.size;

    inst=weka.core.DenseInstance(nAtt+1);

    for i = 0:nAtt
        inst.setValue(attList.get(i), val(i));
    end
    
end
