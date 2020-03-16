function cleanedData= cleanArff(data, keepIndices)

    removeFilter = weka.filters.unsupervised.attribute.Remove;
    removeFilter.setAttributeIndicesArray(keepIndices);
    removeFilter.setInvertSelection(true);
    removeFilter.setInputFormat(data);
    cleanedData= weka.filters.Filter.useFilter(data, removeFilter);
    
end