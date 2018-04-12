function value = getValueFromMovieDatabase(movieDatabase, row, fieldName)
try
  movieDatabaseHeaderRow = movieDatabase(1, :);
  fieldColumn = findColumnIndex(movieDatabaseHeaderRow, fieldName);
  
  value = movieDatabase(row, fieldColumn);
  value = value{1};
catch
    value = [];
end
end