function ImportFromAWS(Prefix)

%Get folders
    [~,UserProcPath,UserDynResPath,~,UserPreProcPath] = ...
        DetermineLocalFolders(Prefix);
    UserPreProcPath_Prefix = [UserPreProcPath,filesep,Prefix];
    UserProcPath_Prefix = [UserProcPath,filesep,Prefix,'_'];
    UserDynResPath_Prefix =[UserDynResPath,filesep,Prefix];
    [~, username] = system('echo %USERNAME%');
    username = strrep(username, sprintf('\n'),''); %removes new line
    HGlabDataFolder = ['S:\HGlab\Dropbox\',username,'\LivemRNA\Data'];
    HGlabPreProcPath_Prefix = [HGlabDataFolder,filesep,'PreProcessedData\',Prefix];
    HGlabProcPath_Prefix = [HGlabDataFolder,filesep,'ProcessedData\',Prefix,'_'];
    HGlabDynResPath_Prefix = [HGlabDataFolder,filesep,'DynamicsResults\',Prefix];
    
%Move files from HGlab Dropbox to users Data folder
    movefile(HGlabPreProcPath_Prefix, UserPreProcPath)
    movefile(HGlabProcPath_Prefix, UserProcPath)
    movefile(HGlabDynResPath_Prefix, UserDynResPath)
    
    disp(['Import from AWS for ' Prefix ' complete.']);

end