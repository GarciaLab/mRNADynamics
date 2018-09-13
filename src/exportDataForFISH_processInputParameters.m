function [Prefix, SkipFrames, ProjectionType, PreferredFileNameForTest, keepTifs] = exportDataForFISH_processInputParameters(varargin)
  %Look at parameters
  SkipFrames = [];
  Prefix = '';
  %Default setting for z-projection is maximum-based. 
  %This may fail when high intensity reflections are present
  ProjectionType = 'maxprojection'; 
  %Added new argument to specify a preferred file name and enable automatic testing
  PreferredFileNameForTest = '';
  keepTifs = false;
  
  k=1;
  while k<=length(varargin)
    if strcmpi(varargin{k},'skipframes')
      SkipFrames=varargin{k+1};
      k=k+1;
      warning('SkipFrame mode.')
    elseif strcmpi(varargin{k},'medianprojection')
      ProjectionType = 'medianprojection';
    elseif strcmpi(varargin{k},'middleprojection')
      ProjectionType = 'middleprojection';
    elseif isobject(varargin{k}) && isa(varargin{k}, 'PreferredFileForTest')
      PreferredFileForTest = varargin{k};
      PreferredFileNameForTest = PreferredFileForTest.fileName;
    elseif strcmpi(varargin{k}, 'keepTifs')
      keepTifs = true;
    else
      Prefix = varargin{k};
    end
    k=k+1;
  end
end