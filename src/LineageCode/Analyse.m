function Prefix=AnalyseDV()

%This function does all the required data exporting and runs the FISH and
%nuclear analysis code for the DV axis


%TO-DO:
%1) Code to locate DV position in embryo, just as is done in AP

%% Copy the images and do the initial FISH analysis
%Process the raw images from the microscopes so that they can be analyzed
%by the FISH code.
Prefix=ExportDataForFISH();
%Optional parameters:
%TAGOnly: Generate the TAG file only

%First do an analysis without a threshold to generate the DoG images.
RunFISHToolbox(Prefix)


%% Look at the dog-filtered images and decide on a threshold to use

%It is advised to keep the threshold low and then increase it after the fact.

%For power of 10mW
Threshold=input('Please input threshold (it is advised to keep it low):       ')

%Now, do an analysis with an actual threshold
RunFISHToolbox(Prefix,Threshold)

%% Track the nuclei and check their segmentation
TrackNuclei(Prefix)

%Check the segmentation
CheckNucleiSegmentation(Prefix)

%If segmentation was modified then we need to rerun the tracking
TrackNuclei(Prefix)

%% Track the particles and check the tracking. 
% We also have the option to check the nuclei tracking in this round.

%Track the particles, the two numbers are Threshold1 and Threshold2
TrackmRNADynamics(Prefix,Threshold,Threshold*1.5)

%% Lineage - Optional Lineage Step
FixLineage(Prefix)
CheckNuclei(Prefix)

%% Final manual curation
CheckParticleTracking(Prefix)

%% And Analyse!
CompileParticles(Prefix);
PostAnalysis(Prefix);
