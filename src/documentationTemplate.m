%[This is a template/example of the standard documentation. Anything that
%appears in [] is a comment about the template or something that should be
%replaced with the appropriate information about your code.]
%[Created by Meghan & Armando on 6/1/17]

function segmentSpotsTemp(Prefix,Threshold,varargin)
% segmentSpotsTemp(Prefix, Threshold, [Options])
%
% DESCRIPTION
% Identify and segregate individual transcription spots. [Make the
% description as short or long as you feel is necessary to communicate the
% purpose of the script/function.]
%
% ARGUMENTS
% Prefix: Prefix of the data set to analyze
% Threshold: Threshold to be used. Should be kept at ~90-200 for lattice
%           light-sheet data, and at ~5-10 for confocal data (Leica SP8).
%           If left empty, then the code just generates the DoG files.
% [Any other parameters you have]
%
% OPTIONS
% [In this example, and many of our other scripts, the options are the
% different parameters you can pass to the function (typically via varargin)
% to adjust how the code runs.]
% 'displayFigures':   If you want to display plots and images.
%                
% 'TrackSpots':   Do you want to use this code to track the particles instead
%                of using TrackmRNADynamics? 
% 'Frames',N:     Run the code from frame 1 to frame N. Defaults to all
%                frames. It's suggested to run 5-20 frames for debugging.
% 'NoShadows':    Spots without valid (peaked) z-profiles are normally
%                discarded. This option overrides that.
% [Any other options you may have]
%
% CONTROLS
% [For interactive scripts with keyboard controls]
%
% OUTPUT
% [Any output that your script provides with a description of what it is. If 
% applicable, please note where this output is saved.]
%
% Author (contact): Armando Reimer (areimer@berkeley.edu)
% Created: 01/01/2016
% Last Updated: 12/31/2016
%
% Documented by: Meghan Turner (meghan_turner@berkeley.edu)

%[Leave at least one blank line (without a comment symbol) before beginning
%any additonal commenting or beginning the code. Everything contained in
%the inital comment block will be included; everything after the first 
%blank, uncommented line will be excluded. Thus, do not add any blank,
%uncommented lines within the documentation header.]

%[CODE]

