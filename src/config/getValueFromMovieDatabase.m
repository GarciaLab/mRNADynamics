function value = getValueFromMovieDatabase(movieDatabase, row, fieldName)

  movieDatabaseHeaderRow = movieDatabase(1, :);
  fieldColumn = findColumnIndex(movieDatabaseHeaderRow, fieldName);
  value = movieDatabase(row, fieldColumn);
  value = value{1};
  
end