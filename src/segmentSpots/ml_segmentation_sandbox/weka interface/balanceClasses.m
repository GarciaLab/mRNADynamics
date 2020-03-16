function data = balanceClasses(data, varargin)

saveAsArff = false;
filename = '';

subSampleFilter = weka.filters.supervised.instance.SpreadSubsample;
subSampleFilter.setOptions({'-M', '1'});
subSampleFilter.setInputFormat(data);
data = weka.filters.Filter.useFilter(data, subSampleFilter);

if ~isempty(varargin)
    saveAsArff = true;
    filename = varargin{1};
end


if saveAsArff
    
    saver = weka.core.converters.ArffSaver();
    saver.setInstances(data);
    saver.setFile(java.io.File(filename));
    saver.writeBatch();
    
end


end