function [NIndices, MaxShift, MaxHistone, ProjectionType, PrefixOverrideFlag, SkipFrames] = ExportDataForFISH_processArguments(varargin);
  varargin = varargin{1};
  NIndices = 3;  %Number of indices ScanImage used to save the files
  MaxShift = 9;  %Maximum shift in pixels corresponding to image shift and alignment
  MaxHistone = 1000; %Maximum intensity for the histone channel. Anything above this will be capped.
  ProjectionType = 'maxprojection'; %Default setting for z-projection is maximum-based. This may fail when high intensity reflections are present.

  PrefixOverrideFlag = 0;
  SkipFrames = [];

  k = 1;
  while k <= length(varargin)
      if strcmpi(varargin{k}, 'skipframes')
          SkipFrames = varargin{k+1};
          k = k + 1;
          warning('SkipFrame mode.')
      elseif strcmpi(varargin{k}, 'medianprojection')
          ProjectionType = 'medianprojection';
      else
          Prefix = varargin{k};
          PrefixOverrideFlag = 1;
      end
      k = k + 1;
  end
end
