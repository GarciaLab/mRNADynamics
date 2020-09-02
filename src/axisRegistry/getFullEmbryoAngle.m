function full_embryo_angle = getFullEmbryoAngle(...
    fullEmbryoPath, surfFile, Prefix)


%If the full embryo has been exported in LASX, use that.

if isfolder([fullEmbryoPath, filesep, 'MetaData']) &&...
        ~isempty(dir([fullEmbryoPath, filesep, 'MetaData', filesep, '*.xml']))
    surfStr = 'Surf';
    xml_file_dir_surf = dir([fullEmbryoPath,...
        'MetaData', filesep,'*', surfStr, '*.xml']);
    xml_file_surf = xml_file_dir_surf(1).name;
    xDoc2 = searchXML([fullEmbryoPath,'MetaData', filesep, xml_file_surf]);
    full_embryo_angle = str2double(evalin('base','rot'));
    evalin('base','clear rot')
    
else
    
    xml_file_surf = [fullEmbryoPath, filesep, 'lifMeta.xml'];
    
    generateLIFMetaDataXML(surfFile, xml_file_surf);
    
    [~, settingStruct] =...
        readSettingsFromLIFMetaDataXML(Prefix, xml_file_surf);
    
    full_embryo_angle = settingStruct.rotatorAngle;
    
end

end
