function index = findColumnIndex(row, value)
  indexMatches = strfind(row, value);
  index = find(not(cellfun('isempty', indexMatches)));
end