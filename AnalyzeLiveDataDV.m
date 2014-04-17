function AnalyzeLiveDataDV

%This function does all the required data exporting and runs the FISH
%analysis code


%TO-DO:
%1) Code to do automated threshold detection. This might only have to
%be done once per construct and setup. 
%2) Integrate FindAPAxis.m and FindAPAxisFullEmbryo.m.
%3) Put the segmentation and tracking into an independent function and not
%in this script.


%2013-08-08: Modified to support Laurent's segmentation and tracking code






%% Copy the images and do the initial FISH analysis


%Process the raw images from the microscopes so that they can be analyzed
% %by the FISH code.
% Prefix=ExportDataForFISH;
%Optional parameters:
%TAGOnly: Generate the TAG file only

Prefix = '2014-04-08-P2P_B';


%Get the relevant folders now:
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,SchnitzcellsFolder]=...
    DetermineLocalFolders(Prefix);


% %Start the matlab workers for the FISH analysis code
% try
%     matlabpool
% catch
%     display('matlabpool already running')
% end
% 
% %First do an analysis without a threshold to generate the DoG images.
% cd([FISHPath,filesep,'Analysis'])
% analyzeShawnLibrary('fad',@(x)tagged(x,'id',[Prefix,'_']),'params_mRNADynamics',inf)
% cd([MS2CodePath])


%% Look at the dog-filtered images and decide on a threshold to use
% 
% %We will keep the threshold low and then increase it after the fact.
% 
% %For LSM settings
 Threshold=250;   
% 
% %For Princeton
% Threshold=90;



%Now, do an analysis with an actual threshold
cd([FISHPath,filesep,'Analysis'])
analyzeShawnLibrary('fad',@(x)tagged(x,'id',[Prefix,'_']),'params_mRNADynamics',Threshold)
cd([MS2CodePath])



%% Track the nuclei and check their segmentation
TrackNuclei(Prefix)

%Check the segmentation
CheckNucleiSegmentation(Prefix)
%If segmentation was modified then we need to rerun the tracking
TrackNuclei(Prefix)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Track the particles and check the tracking. We also have the option to
%% check the nuclei tracking in this round.

%Track the particles, the two numbers are Threshold1 and Threshold2
TrackmRNADynamics(Prefix,300,300)


CheckParticleTracking(Prefix)


%% Done


CompileParticles(Prefix,'ApproveAll')




    