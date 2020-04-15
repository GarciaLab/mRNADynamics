function dataSet = setInstancesVal(dataSet, instInd, attInd, val)

    inst = dataSet.instance(instInd);

    inst.setValue(attInd, val);

    dataSet.set(instInd, inst);

end