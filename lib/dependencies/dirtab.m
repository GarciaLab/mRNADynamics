function files = dirtab(in)
%DIRTAB  Tabulate directory list (dir results as table also with extension).
%   files = dirtab(in) returns the directory listing of in as a table. MATLABs builtin DIR command generates directory
%   listing but the struct it would return is converted into a table. The variable names of the resulting table are the
%   field names as returned by DIR, supplemented by 'ext' which identifies the file-extension, and folder (returned by
%   DIR only since R2016b). The differences between dirtab and dir are as follows:
%      DIRTAB returns a table instead of a struct
%      DIRTAB adds the variable folder, even before R2016b. A trailing filesep is added (not done by DIR in R2016b)
%      The name column does not include the file extension, which appears in a separate column
%      The date column now contains datetime representation of file dates instead of cellstr text.
%      the pseudo-dirs '.' and '..' don't appear in table rows.
%   Tabulating the directory with dirtab makes it easier to sort by file size or file extension, say, using sortrows.
%   The first 3 columns are folder, file, and ext which can be catenated to obtain full file paths.
%
%   Example    
%      dirtab(which('dirtab')) % produces a single match
%      dirtab('') % likely produces no match
%      files = dirtab;  % likely produces multiple matches
%      files = sortrows(files, 'bytes')
%      list = strcat(files.folder, files.name, files.ext)
%      assert(isequal( ~files.isdir, cellfun(@exist, list)==2))
%
%   See also DIR, WHAT, LS, TABLE.

%% Obtain directory
if nargin==0, in='.'; end
files = dir(in);
f = ~ismember( {files.name}, {'.' '..'}); % logical filter for entries except pseudo-files '.' and '..'

%% Convert into table
if numel(files)>1 
  files = struct2table(files); % we had at least 2 examples, we can go in 1 step
else % files is empty or scalar
  % in this case we create shims, based on a single example of a match, which we remove below
  shim = dir([mfilename('fullpath') '.m']); % guaranteed scalar struct of correct forms
  files = struct2table([files shim([1 1])]); % need 2 shims to get cellstr instead of char array
end
files.Properties.Description = sprintf('dir(''%s'') obtained on %s for %s at %s', ...
                      in, getenv('COMPUTERNAME'), getenv('USERNAME'), datestr(now)); % metadata for PC
files = files(f, :);  % this will drop shims if used, because f is short

%% Split file-name and extension
files.ext = cell(height(files), 1); % Add a column for the file extension
f = ~files.isdir; % logical filter to exclude directory names (files only)
in2 = regexp([{} files.name(f)], '^(.*?)(\.[^\\\./]*)$', 'tokens', 'once'); % split name/ext (i.e. in two)
ok = ~cellfun(@isempty, in2); % where regex matched; we should have {name ext} in each cell
f(f) = ok; % f, overall files filter, now excludes non-match (ie. files without ext)
in2 = vertcat(cell(0,2), in2{ok}); % includes sizer shim (expect 2 columns, name & ext)
files{f, {'name' 'ext'}} = in2; % replace name and ext cellstr arrays
files.ext(~f) = {''}; % any files without ext

%% put datenum as datetime 
dt = datetime(files.datenum, 'ConvertFrom', 'datenum');
% assert( isequal(files.date, cellstr(dt)) ) % check that representation matches
files.date = dt; % write back date with datetime representation
arrange = {'folder' 'name' 'ext'}; % arrange in consecutive columns in order required
if ~any(strcmp('folder', files.Properties.VariableNames))
  in = {in};
  files.folder = in(ones(size(ext)));
end
% Add trailing filesep to folder column to facilitate full path names
files.folder = regexprep(files.folder, '(\\|/)?$', filesep, 'emptymatch', 'once');
files = files(:, [arrange setdiff(files.Properties.VariableNames, arrange, 'stable' )]);
