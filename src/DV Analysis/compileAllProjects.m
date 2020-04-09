function compileAllProjects(DataType)

[~, ~, prefixes] = getDorsalPrefixes(DataType);

compiledProjects = cell(1, length(prefixes));
for i = 1:length(prefixes)
    compiledProjects{i} = makeCompiledProject(prefixes{i});
end


% binDorsal(dataType)
binDorsal(DataType, true)


end