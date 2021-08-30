function [ZoomImage, SurfImage, ZoomRatio, SurfInfo, SurfColumns, SurfRows, Rows, Columns] =...
	APP_getTifZoomInfo(DropboxFolder, Prefix, rawPrefixPath, fullEmbryoPath, ChannelToLoad)

        FileModeObj = TIFFileMode(rawPrefixPath, DropboxFolder, Prefix);
        
        ImageInfo = FileModeObj.getImageInfo();
        
        MovieZoom = FileModeObj.getMovieZoom(ImageInfo);
        
        %Get the zoomed out surface image and its dimensions from the FullEmbryo folder
        [Rows, Columns, SurfInfo, SurfZoom, SurfRows, SurfColumns] = FileModeObj.getSurfValues(ChannelToLoad);
        
        % HG: I had to add this for some surface images that were edited
        % with ImageJ and lost their metadata. Note that I'm hardcoding the
        % zoom of the low magnification images.
        if isnan(SurfRows)
            SurfRows = SurfInfo(1).Height;
            SurfColumns = SurfInfo(1).Width;
            SurfZoom = 1;
        end
        
        FullEmbryo = FileModeObj.getFullEmbryoImage();
        
        % Ratio between the two zoom levels
        ZoomRatio = MovieZoom / SurfZoom;
        ResizeFactor = max([Rows / SurfRows * ZoomRatio, Columns / SurfColumns * ZoomRatio]);

        % ES 2013-10-30: the reason I have to define ResizeFactor differently
        % from ZoomRatio is because you can't necessarily infer the
        % microns-per-pixel resolution from the zoom alone: it also depends on
        % the dimensions of the image. This may not work for all possible
        % resolutions, though...
        ZoomRatio = ResizeFactor;
end
