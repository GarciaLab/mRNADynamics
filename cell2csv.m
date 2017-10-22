function cell2csv(filename, cellArray, delimiter, types)
% Writes cell array content into a *.csv file.
% 
% CELL2CSV(filename,cellArray,delimiter)
%
% filename      = Name of the file to save. [ i.e. 'text.csv' ]
% cellarray    = Name of the Cell Array where the data is in
% delimiter = seperating sign, normally:',' (it's default)
%
% by Sylvain Fiedler, KA, 2004
% modified by Rob Kohr, Rutgers, 2005 - changed to english and fixed delimiter
% UC Berkeley, 2017, added support for column types, specifically excel dates.
  if nargin < 3
    delimiter = ',';
  end

  file = fopen(filename, 'w');
  for z = 1:size(cellArray, 1)
    for s = 1:size(cellArray, 2)
      var = eval(['cellArray{z,s}']);

      if size(var, 1) == 0
        var = '';
      end

      if isnumeric(var) == 1
        if nargin < 4 || s > length(types)
          if(isnan(var))
            var = '';
          else
            var = num2str(var);
          end
        else
          format = types{s};
          if (format == 'date')
            % dates in excel are integer numbers from 30/dec/1899
            var = datestr(var+datenum('30-Dec-1899'));
          else
            var = num2str(var);
          end
        end
      end

      var = strrep(var, '"', '""'); %escape double quotes
      var = regexprep(var, '[\n\r]+' , ''); %eliminate newlines
      var = ['"', var, '"']; %add text delimiters (double quotes)
      fprintf(file, '%s', var); %print to file as literal string %s, to prevent interpreting var as a formatSpec for fprintf

      if s ~= size(cellArray,2)
        fprintf(file, [delimiter]);
      end
    end
    fprintf(file, '\n');
  end

  fclose(file);
end