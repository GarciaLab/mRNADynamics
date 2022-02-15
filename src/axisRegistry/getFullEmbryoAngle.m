function full_embryo_angle = getFullEmbryoAngle(...
    fullEmbryoPath, surfFile, Prefix)


%If the full embryo has been exported in LASX, use that.

xml_file_surf = [fullEmbryoPath, filesep, 'lifMeta.xml'];

generateLIFMetaDataXML(surfFile, xml_file_surf);

[~, settingStruct] =...
    readSettingsFromLIFMetaDataXML(Prefix, xml_file_surf);

full_embryo_angle = settingStruct.rotatorAngle;
