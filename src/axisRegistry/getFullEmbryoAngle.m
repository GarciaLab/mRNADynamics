function full_embryo_angle = getFullEmbryoAngle(...
    fullEmbryoPath, surfFile, Prefix)

    %deprecated
%     xml_file_dir_surf = dir([fullEmbryoPath,...
%         'MetaData', filesep,'*', surfOrMidStr, '*.xml']);
%     xml_file_surf = xml_file_dir_surf(1).name;
%     
%     
%     try
        
        xml_file_surf = [fullEmbryoPath, filesep, 'lifMeta.xml'];
        
        generateLIFMetaDataXML(surfFile, xml_file_surf);
        
        [~, settingStruct] =...
            readSettingsFromLIFMetaDataXML(Prefix, xml_file_surf);
        
        full_embryo_angle = settingStruct.rotatorAngle;
%     catch
%         warning('leaving this here. ignore for now.');
%     end
%     
%     
%     
%     xDoc2 = searchXML([fullEmbryoPath,'MetaData', filesep, xml_file_surf]);
%     full_embryo_angle = str2double(evalin('base','rot'));
% else
%     warning('No full embryo metadata found.')
% end
% 
% evalin('base','clear rot')

end
