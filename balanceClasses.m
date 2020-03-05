function data = balanceClasses(data)

    subSampleFilter = weka.filters.supervised.instance.SpreadSubsample;
    subSampleFilter.setOptions({'-M', '1'});
    subSampleFilter.setInputFormat(data);
    data = weka.filters.Filter.useFilter(data, subSampleFilter);
    
end