function Value = getValueFromOMEMetadata(MetaData, index)
    Value = split(MetaData(index), '=');
    Value = split(Value(2), '/');
    
    Value = Value(1);
    Value = Value{1};
end
