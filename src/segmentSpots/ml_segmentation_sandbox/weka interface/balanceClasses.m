function data = balanceClasses(data, varargin)

saveAsArff = false;
filename = '';
sampler = 'Resample';

if strcmpi(sampler, 'SpreadSubSample')
    
    subSampleFilter = weka.filters.supervised.instance.SpreadSubsample;
    subSampleFilter.setOptions({'-M', '1'});
    subSampleFilter.setInputFormat(data);
    data = weka.filters.Filter.useFilter(data, subSampleFilter);
    
elseif strcmpi(sampler, 'Resample')
    
    resampleFilter = weka.filters.supervised.instance.Resample;
    resampleFilter.setOptions({'-B', '1.0','-S', '1', '-Z','100.0'});
    resampleFilter.setInputFormat(data);
    data = weka.filters.Filter.useFilter(data, resampleFilter);
    
end

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