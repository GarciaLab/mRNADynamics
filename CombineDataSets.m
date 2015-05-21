function CombineDataSets(DataType)

[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;

all_data = LoadMS2Sets(DataType);
Combined_cp = [];
Combined_apdiv = {};

for i=1:length(all_data)
    Combined_cp = [Combined_cp, all_data(i).CompiledParticles];
    Combined_apdiv = [Combined_apdiv, all_data(i).APDivision];
end

save([DropboxFolder,filesep,[DataType,'_Combined_cp.mat']],'Combined_cp');
display('Combined_cp saved');
save([DropboxFolder,filesep,[DataType,'_Combined_apdiv.mat']],'Combined_apdiv');
display('Combined_apdiv saved'); 
