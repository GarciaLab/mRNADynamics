classdef DSPINFileMode < FileMode
   methods
      function obj = DSPINFileMode()
         obj = obj@FileMode('DSPIN');
      end 

      function D = readMovieDir(this)
      	fprintf('Reading dir for filemode %s\n', this.identifier);
      end
   end
end
