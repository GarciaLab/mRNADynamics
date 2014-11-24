%These parameters override those in setFishDefaultParams.m

%Trying this to be able to track better
params.shadowN=2;

%params.matlabWorkersToUse=2;        %Otherwise it tries with as many as the computer supports.
params.stopAfter='fits';
params.saveCompactFadInPreanalysisMode=true;

params.fit_prefitMode = FITMODE_ELLIPTICAL;

params.useGUIprogressbar=false;
params.paramID='preanalysis';
params.ap_fullyAutomatic=true;
params.usePreanalysis=false;
params.DoG_center=1.5;
params.DoG_surround=2.5;
params.storeShadowsInCompactFAD=true;
params.processColumnsPeakingAtStackEdge=true;

%Change the separation in xy between slices that defines the maximum
%separation between slices that make a shadow
params.shadow_dist=4;
params.align_mode='none';       %Don't do any alignment

%Get rid of the chopping of the image at the margins. This was there because
%of the alignment
params.align_maxShiftOverStack=0;


%Output file:
params.outputFileID=fopen('AnalysisLog.txt','a');


%params.IAI_mode=1;      %Load images one at a time instead of having everything in memory
%params.ap_zoomFactorRange=5:0.05:5.1;       %This is for the 40x objective

% 
% 'DoG_center',1.5,...         % the center size in pixels of the mexican hat
% 'DoG_surround',2.5,...       % the surround size in pixels of the mexican hat
% 'DoG_filterSize',15,...      % the size of the DoG filter in pixels, 
% ...                          % filterSize x filterSize.
% ...                          % Must be odd (this avoids a half-pixel bias) and preferably substantially bigger than DoG_surround.
% 'DoG_neighborhood',1,...     % when searching for local maxima in the DoG-filtered image, 
% ...                          % this parameter defines what we understand by "local" 
% ...%'DoG_threshold',30,...   % NOW PART OF stackDescription
% 'shadowN',3,...              % The required number of shadows for a spot to make it from brightSpots to candidateSpots
% 
