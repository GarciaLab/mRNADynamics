classdef (Abstract) FileMode

	properties
		identifier = '';
	end

	methods
		function this = FileMode(identifier, rawPrefixPath, DropboxFolder, Prefix)
			this.identifier = identifier;
			this.rawPrefixPath = rawPrefixPath;
			this.fullEmbryoPath = [rawPrefixPath, 'FullEmbryo', filesep];
			this.DropboxFolder = DropboxFolder;
			this.Prefix = Prefix;
		end
	end
	
	methods (Abstract)
		D = readMovieDir()
	end

end
