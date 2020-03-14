function ExportToAWS(Prefix)

%Get folders
    [~,UserProcPath,UserDynResPath,~,UserPreProcPath] = ...
        DetermineLocalFolders(Prefix);
    UserPreProcPath_Prefix = [UserPreProcPath,filesep,Prefix];
    UserPreProcPath_Prefix_stacks = [UserPreProcPath_Prefix,filesep,'stacks'];
    UserProcPath_Prefix = [UserProcPath,filesep,Prefix,'_'];
    UserDynResPath_Prefix =[UserDynResPath,filesep,Prefix];
    CONFIG_CSV_PATH = 'ComputerFolders.csv';
    configValues = csv2cell(CONFIG_CSV_PATH, 'fromfile');
    DefaultDropboxFolder = getConfigValue(configValues, 'DropboxFolder');
    
%Generate AWS folders on the HGlab user
    [~, username] = system('echo %USERNAME%');
    username = strrep(username, sprintf('\n'),''); %removes new line
    HGlabLivemRNAFolder = ['S:\HGlab\Dropbox\',username,'\LivemRNA'];
    HGlabDataFolder = [HGlabLivemRNAFolder,'\Data'];
    warning('off','MATLAB:MKDIR:DirectoryExists');
    warning('off','MATLAB:legend:IgnoringExtraEntries');
    mkdir(HGlabDataFolder)
    mkdir(HGlabDataFolder, ['\PreProcessedData\',Prefix,'\stacks'])
    mkdir(HGlabDataFolder, ['\ProcessedData\',Prefix,'_'])
    mkdir(HGlabDataFolder, ['\DynamicsResults\',Prefix])
    
%Copy Relevant Folders 
    %Should I change this to movefile() instead?
    copyfile(UserPreProcPath_Prefix_stacks, [HGlabDataFolder,'\PreProcessedData\',Prefix,'\stacks'])
    if exist(UserProcPath_Prefix)
        copyfile(UserProcPath_Prefix, [HGlabDataFolder,'\ProcessedData\',Prefix,'_'])
    end
    copyfile(UserDynResPath_Prefix, [HGlabDataFolder,'\DynamicsResults\',Prefix])
    copyfile([DefaultDropboxFolder,'\MovieDatabase.csv'],...
        [HGlabDataFolder,'\DynamicsResults\'])
    
%Generate ComputerFolders.csv
    if exist([HGlabLivemRNAFolder,filesep,'ComputerFolders.csv'])
        warning([HGlabLivemRNAFolder,filesep,'ComputerFolders.csv already exists. Not overwriting.'])
    else
        AWSLivemRNAFolder = ['C:\Users\',username,'\Dropbox\',username,'\LivemRNA']; %CHANGE ME%
        AWSDataFolder = [AWSLivemRNAFolder,'\Data'];
        ComputerFolder{1,1}= 'Computer Name';
        ComputerFolder{2,1}= 'User Name';
        ComputerFolder{3,1}= 'SourcePath';
        ComputerFolder{4,1}= 'PreProcPath';
        ComputerFolder{5,1}= 'FISHPath';
        ComputerFolder{6,1}= 'DropboxFolder';
        ComputerFolder{7,1}= 'MS2CodePath';
        ComputerFolder{1,2}= 'AWS Server';
        ComputerFolder{2,2}= ['C:\Users\',username];
        ComputerFolder{3,2}= [AWSDataFolder,filesep,'RawDynamicsData'];
        ComputerFolder{4,2}= [AWSDataFolder,filesep,'PreProcessedData'];
        ComputerFolder{5,2}= [AWSDataFolder,filesep,'ProcessedData'];
        ComputerFolder{6,2}= [AWSDataFolder,filesep,'DynamicsResults'];
        ComputerFolder{7,2}= [AWSLivemRNAFolder,filesep,'mRNADynamics'];

        cell2csv([HGlabLivemRNAFolder,filesep,'ComputerFolders.csv'],ComputerFolder);
    end
    
    disp(['Export for AWS for ' Prefix ' complete.']);
    
end
