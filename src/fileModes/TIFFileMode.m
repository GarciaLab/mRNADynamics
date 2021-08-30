classdef TIFFileMode < FileMode
   methods
   	  function obj = TIFFileMode(rawPrefixPathm, DropboxFolder, Prefix)
         obj = obj@FileMode('TIF', rawPrefixPath, DropboxFolder, Prefix);
      end 

      function D = readMovieDir(this)
      	fprintf('Reading dir for filemode %s\n', this.identifier);
      	D = dir([this.rawPrefixPath, '*.tif']);
      end

      function ImageInfo = getImageInfo(this)
      	D = this.readMovieDir();
        ImageInfo = imfinfo([this.rawPrefixPath, D(1).name]);
      end

      function MovieZoom = getMovieZoom(this, ImageInfo)
      	%Figure out the zoom factor
        MovieZoom = ExtractInformationField(ImageInfo(1),'state.acq.zoomFactor=');
        MovieZoom = str2double(MovieZoom);
      end

      function [SurfName, SurfImage] = getSurfs(this, ChannelToLoad)
      	D = dir([this.fullEmbryoPath,filesep,'*.tif']);
        SurfName = D(find(~cellfun('isempty', strfind(lower({D.name}), 'surf')))).name;
        SurfImage = imread([this.fullEmbryoPath,filesep,SurfName],ChannelToLoad);
      end

      function [Rows, Columns, SurfInfo, SurfZoom, SurfRows, SurfColumns, SurfImage] = getSurfValues(this, ChannelToLoad)
      	[SurfName, SurfImage] = this.getSurfs(ChannelToLoad);

     	% Get the size of the zoom image
        Rows = str2double(ExtractInformationField(ImageInfo(1), 'state.acq.linesPerFrame='));
        Columns = str2double(ExtractInformationField(ImageInfo(1), 'state.acq.pixelsPerLine='));
        
        SurfInfo = imfinfo([this.fullEmbryoPath, filesep, SurfName]);
        SurfZoom = ExtractInformationField(SurfInfo(1), 'state.acq.zoomFactor=');
        SurfZoom = str2double(SurfZoom);
        
        SurfRows = str2double(ExtractInformationField(SurfInfo(1), 'state.acq.linesPerFrame='));
        SurfColumns = str2double(ExtractInformationField(SurfInfo(1), 'state.acq.pixelsPerLine='));
      end

      function FullEmbryo = getFullEmbryoImage(this)
      	FullEmbryo = imread([this.DropboxFolder, filesep, this.Prefix, filesep, 'APDetection', filesep, 'FullEmbryo.tif']);
      end
   end
end
