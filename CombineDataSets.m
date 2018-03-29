function CombineDataSets(DataType)

[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;

all_data = LoadMS2Sets(DataType);
all_data2 = struct();
%Loop over all the fields in the data structure
for str = fieldnames(all_data)'
    if all(size(all_data(1).(char(str))) == size(all_data(2).(char(str))))
        all_data2 = setfield(all_data2, char(str), 0);
        all_data2.(char(str)) = [all_data(1).(char(str)), all_data(2).(char(str))]
    end
end

MeanVectorAPTot = vertcat(all_data(1).MeanVectorAP, all_data(2).MeanVectorAP);
ElapsedTimeTot = horzcat(all_data(1).ElapsedTime,all_data(2).ElapsedTime);
NParticlesAPTot = vertcat(all_data(1).NParticlesAP,all_data(2).NParticlesAP);
SDVectorAPTot = vertcat(all_data(1).SDVectorAP,all_data(2).SDVectorAP);
OnRatioAPTot = vertcat(all_data(1).OnRatioAP,all_data(2).OnRatioAP);
APDivTot = all_data(1).APDivision+all_data(2).APDivision;

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
