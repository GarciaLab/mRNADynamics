function prefix_cell = compilePrefixNames(project,userDynResPath)

% The following code is creditted to Nick Lammers and Meghan Turner and
% modified for this purpose by Emma Luu. 

% find sheet
sheet_path = [userDynResPath,filesep, 'DataStatus.xlsx'];
[~,sheet_names]=xlsfinfo(sheet_path);
sheet_index = find(ismember(sheet_names,project));
if isempty(sheet_index)
    error('no tab matching "project" string found in DataStatus')
end
[~,~,sheet_cell] = xlsread(sheet_path,sheet_index);
name_col = sheet_cell(1:33,1); % hard coded for now
ready_ft = contains(name_col,'Ran ExportDataForFISH');
ready_cols = 1 + find([sheet_cell{ready_ft,2:end}]==1);
sheet_cell = sheet_cell(:,[1 ready_cols]);
% get list of project names
prefix_ft = contains(name_col,'Prefix');
prefix_cell_raw = sheet_cell(prefix_ft,2:end);
prefix_cell = {};
for i = 1:numel(prefix_cell_raw)
    if ~isempty(prefix_cell_raw{i})
        eval([prefix_cell_raw{i} ';'])
        prefix_cell = [prefix_cell{:} {Prefix}];
    end
end

end
