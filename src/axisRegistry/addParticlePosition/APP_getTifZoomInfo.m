function [ZoomImage, SurfImage, ZoomRatio, SurfInfo, SurfColumns, SurfRows, Rows, Columns] =...
	APP_getTifZoomInfo(DropboxFolder, Prefix, rawPrefixPath, fullEmbryoPath, ChannelToLoad)

        D=dir([rawPrefixPath,'*.tif']);
        ImageInfo = imfinfo([rawPrefixPath,D(1).name]);
        
        %Figure out the zoom factor
        MovieZoom=ExtractInformationField(ImageInfo(1),'state.acq.zoomFactor=');
        MovieZoom=str2double(MovieZoom);
        
        
        %Get the zoomed out surface image and its dimensions from the FullEmbryo folder
        D=dir([fullEmbryoPath,filesep,'*.tif']);
        SurfName=D(find(~cellfun('isempty',strfind(lower({D.name}),'surf')))).name;
        SurfImage=imread([rawPrefixPath,...
            'FullEmbryo',filesep,SurfName],ChannelToLoad);
        
        %Get the size of the zoom image
        Rows = str2double(ExtractInformationField(ImageInfo(1), 'state.acq.linesPerFrame='));
        Columns = str2double(ExtractInformationField(ImageInfo(1), 'state.acq.pixelsPerLine='));
        
        % TO-DO: use rawPrefixPath and fullEmbryoPath to concat this
        SurfInfo = imfinfo([RawDynamicsPath, filesep, projectDate, filesep, EmbryoName, filesep, 'FullEmbryo', filesep, SurfName]);
        SurfZoom = ExtractInformationField(SurfInfo(1), 'state.acq.zoomFactor=');
        SurfZoom = str2double(SurfZoom);
        
        SurfRows = str2double(ExtractInformationField(SurfInfo(1), 'state.acq.linesPerFrame='));
        SurfColumns = str2double(ExtractInformationField(SurfInfo(1), 'state.acq.pixelsPerLine='));
        
        %HG: I had to add this for some surface images that were edited
        %with ImageJ and lost their metadata. Note that I'm hardcoding the
        %zoom of the low magnification images.
        if isnan(SurfRows)
            SurfRows=SurfInfo(1).Height;
            SurfColumns=SurfInfo(1).Width;
            SurfZoom=1;
        end
        
        
        %Get the full embryo image
        FullEmbryo=imread([DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'FullEmbryo.tif']);
        
        %Ratio between the two zoom levels
        ZoomRatio = MovieZoom / SurfZoom;
        ResizeFactor = max([Rows/SurfRows*ZoomRatio, Columns/SurfColumns*ZoomRatio]);
        % ES 2013-10-30: the reason I have to define ResizeFactor differently
        % from ZoomRatio is because you can't necessarily infer the
        % microns-per-pixel resolution from the zoom alone: it also depends on
        % the dimensions of the image. This may not work for all possible
        % resolutions, though...
        ZoomRatio = ResizeFactor;
end