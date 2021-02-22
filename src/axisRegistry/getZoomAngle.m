function zoom_angle = getZoomAngle(Prefix, rawPrefixPath)


%If the full embryo has been exported in LASX, use that. 

if isfolder([rawPrefixPath, 'MetaData']) && ~isempty(dir([rawPrefixPath, 'MetaData', filesep, '*.xml']))
    
    xml_file_dir = dir([rawPrefixPath,'MetaData', filesep, '*.xml']);
    xml_file = xml_file_dir(1).name;
    xDoc = searchXML([rawPrefixPath,'MetaData', filesep, xml_file]);
    zoom_angle = str2double(evalin('base','rot'));
    evalin('base','clear rot')
    
else
     xml_file = [rawPrefixPath, filesep, 'lifMeta.xml'];

    generateLIFMetaDataXML(Prefix, xml_file);

    [~, settingStruct] =...
        readSettingsFromLIFMetaDataXML(Prefix, xml_file);

    zoom_angle = settingStruct.rotatorAngle;
end

end