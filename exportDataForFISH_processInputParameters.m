function [Prefix, SkipFrames, ProjectionType] = exportDataForFISH_processInputParameters(varargin)
  %Look at parameters
  SkipFrames=[];
  Prefix='';

  %Default setting for z-projection is maximum-based...
  %This may fail when high intensity reflections are present
  ProjectionType = 'maxprojection'; 

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
    else
      Prefix = varargin{k};
    end
    k=k+1;
  end
end