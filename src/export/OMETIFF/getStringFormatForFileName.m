function Format = getStringFormatForFileName(Value)
  Format = ['%0', num2str(length(num2str(Value)) + 1), 'd'];
end