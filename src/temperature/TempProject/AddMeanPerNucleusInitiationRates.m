 function this = AddMeanPerNucleusInitiationRates(this)
     alpha = this.alpha;
     Temperatures = fliplr(unique(this.Temp_sps));
     NumTemperatures = length(Temperatures);
     NumExperiments = length(this.ExperimentPrefixes);
     APResolution = this.Experiments{1}.APResolution;
     NumAPbins = uint16(1/APResolution)+1;
     this.PerNucleusProfileParameters.MeanInitiationRates = {};
     this.PerNucleusProfileParameters.TimeOns = {};
     this.PerNucleusProfileParameters.TimeOffs = {};
     this.PerNucleusProfileParameters.ElongationTimes = {};
     this.PerNucleusProfileParameters.UnloadingRates = {};
     this.PerNucleusProfileParameters.MeanFitR2s = {};

     this.PerNucleusProfileParameters.MeanFitR2s.AnaphaseAligned = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.MeanFitR2s.AnaphaseAligned3D = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.MeanFitR2s.Tbinned = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.MeanFitR2s.Tbinned3D = NaN(NumExperiments, NumAPbins, 6);

     this.PerNucleusProfileParameters.MeanInitiationRates.AnaphaseAligned = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.MeanInitiationRates.AnaphaseAlignedStdError = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.MeanInitiationRates.AnaphaseAlignedCIhigh = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.MeanInitiationRates.AnaphaseAlignedCIlow= NaN(NumExperiments, NumAPbins, 6);

     this.PerNucleusProfileParameters.TimeOns.AnaphaseAligned = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.TimeOns.AnaphaseAlignedStdError = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.TimeOns.AnaphaseAlignedCIhigh = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.TimeOns.AnaphaseAlignedCIlow= NaN(NumExperiments, NumAPbins, 6);

     this.PerNucleusProfileParameters.TimeOffs.AnaphaseAligned = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.TimeOffs.AnaphaseAlignedStdError = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.ElongationTimes.AnaphaseAligned = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.ElongationTimes.AnaphaseAlignedStdError = NaN(NumExperiments, NumAPbins, 6);

     this.PerNucleusProfileParameters.UnloadingRates.AnaphaseAligned = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.UnloadingRates.AnaphaseAlignedStdError = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.UnloadingRates.AnaphaseAlignedCIhigh = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.UnloadingRates.AnaphaseAlignedCIlow= NaN(NumExperiments, NumAPbins, 6);
     
     this.PerNucleusProfileParameters.MeanFitR2s.Unaligned = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.MeanFitR2s.Unaligned3D = NaN(NumExperiments, NumAPbins, 6);
     
     
     this.PerNucleusProfileParameters.MeanInitiationRates.Unaligned = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.MeanInitiationRates.UnalignedStdError = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.MeanInitiationRates.UnalignedCIhigh = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.MeanInitiationRates.UnalignedCIlow= NaN(NumExperiments, NumAPbins, 6);

     this.PerNucleusProfileParameters.TimeOns.Unaligned = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.TimeOns.UnalignedStdError = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.TimeOns.UnalignedCIhigh = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.TimeOns.UnalignedCIlow= NaN(NumExperiments, NumAPbins, 6);

     this.PerNucleusProfileParameters.TimeOffs.Unaligned = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.TimeOffs.UnalignedStdError = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.ElongationTimes.Unaligned = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.ElongationTimes.UnalignedStdError = NaN(NumExperiments, NumAPbins, 6);

     this.PerNucleusProfileParameters.UnloadingRates.Unaligned = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.UnloadingRates.UnalignedStdError = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.UnloadingRates.UnalignedCIhigh = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.UnloadingRates.UnalignedCIlow= NaN(NumExperiments, NumAPbins, 6);


     this.PerNucleusProfileParameters.MeanInitiationRates.Tbinned = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.MeanInitiationRates.TbinnedStdError = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.MeanInitiationRates.TbinnedCIhigh = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.MeanInitiationRates.TbinnedCIlow= NaN(NumExperiments, NumAPbins, 6);

     this.PerNucleusProfileParameters.TimeOns.Tbinned = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.TimeOns.TbinnedStdError = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.TimeOns.TbinnedCIhigh = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.TimeOns.TbinnedCIlow= NaN(NumExperiments, NumAPbins, 6);

     this.PerNucleusProfileParameters.TimeOffs.Tbinned = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.TimeOffs.TbinnedStdError = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.ElongationTimes.Tbinned = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.ElongationTimes.TbinnedStdError = NaN(NumExperiments, NumAPbins, 6);

     this.PerNucleusProfileParameters.UnloadingRates.Tbinned = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.UnloadingRates.TbinnedStdError = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.UnloadingRates.TbinnedCIhigh = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.UnloadingRates.TbinnedCIlow= NaN(NumExperiments, NumAPbins, 6);

     this.PerNucleusProfileParameters.MeanInitiationRates.AnaphaseAligned3D = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.MeanInitiationRates.AnaphaseAligned3DStdError = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.MeanInitiationRates.AnaphaseAligned3DCIhigh = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.MeanInitiationRates.AnaphaseAligned3DCIlow= NaN(NumExperiments, NumAPbins, 6);

     this.PerNucleusProfileParameters.TimeOns.AnaphaseAligned3D = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.TimeOns.AnaphaseAligned3DStdError = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.TimeOns.AnaphaseAligned3DCIhigh = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.TimeOns.AnaphaseAligned3DCIlow= NaN(NumExperiments, NumAPbins, 6);

     this.PerNucleusProfileParameters.TimeOffs.AnaphaseAligned3D = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.ElongationTimes.AnaphaseAligned3D = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.TimeOffs.AnaphaseAligned3DStdError = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.ElongationTimes.AnaphaseAligned3DStdError = NaN(NumExperiments, NumAPbins, 6);

     this.PerNucleusProfileParameters.UnloadingRates.AnaphaseAligned3D = NaN(NumTemperatures, NumAPbins, 6);
     this.PerNucleusProfileParameters.UnloadingRates.AnaphaseAligned3DStdError = NaN(NumTemperatures, NumAPbins, 6);
     this.PerNucleusProfileParameters.UnloadingRates.AnaphaseAligned3DCIhigh = NaN(NumTemperatures, NumAPbins, 6);
     this.PerNucleusProfileParameters.UnloadingRates.AnaphaseAligned3DCIlow= NaN(NumTemperatures, NumAPbins, 6);
     
     this.PerNucleusProfileParameters.MeanInitiationRates.Unaligned3D = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.MeanInitiationRates.Unaligned3DStdError = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.MeanInitiationRates.Unaligned3DCIhigh = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.MeanInitiationRates.Unaligned3DCIlow= NaN(NumExperiments, NumAPbins, 6);

     this.PerNucleusProfileParameters.TimeOns.Unaligned3D = NaN(NumTemperatures, NumAPbins, 6);
     this.PerNucleusProfileParameters.TimeOns.Unaligned3DStdError = NaN(NumTemperatures, NumAPbins, 6);
     this.PerNucleusProfileParameters.TimeOns.Unaligned3DCIhigh = NaN(NumTemperatures, NumAPbins, 6);
     this.PerNucleusProfileParameters.TimeOns.Unaligned3DCIlow= NaN(NumTemperatures, NumAPbins, 6);

     this.PerNucleusProfileParameters.TimeOffs.Unaligned3D = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.TimeOffs.Unaligned3DStdError = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.ElongationTimes.Unaligned3D = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.ElongationTimes.Unaligned3DStdError = NaN(NumExperiments, NumAPbins, 6);

     this.PerNucleusProfileParameters.UnloadingRates.Unaligned3D = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.UnloadingRates.Unaligned3DStdError = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.UnloadingRates.Unaligned3DCIhigh = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.UnloadingRates.Unaligned3DCIlow= NaN(NumExperiments, NumAPbins, 6);



     this.PerNucleusProfileParameters.MeanInitiationRates.Tbinned3D= NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.MeanInitiationRates.Tbinned3DStdError = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.MeanInitiationRates.Tbinned3DCIhigh = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.MeanInitiationRates.Tbinned3DCIlow= NaN(NumExperiments, NumAPbins, 6);

     this.PerNucleusProfileParameters.TimeOns.Tbinned3D= NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.TimeOns.Tbinned3DStdError = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.TimeOns.Tbinned3DCIhigh = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.TimeOns.Tbinned3DCIlow= NaN(NumExperiments, NumAPbins, 6);

     this.PerNucleusProfileParameters.TimeOffs.Tbinned3D= NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.ElongationTimes.Tbinned3D= NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.TimeOffs.Tbinned3DStdError= NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.ElongationTimes.Tbinned3DStdError= NaN(NumExperiments, NumAPbins, 6);

     this.PerNucleusProfileParameters.UnloadingRates.Tbinned3D= NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.UnloadingRates.Tbinned3DStdError = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.UnloadingRates.Tbinned3DCIhigh = NaN(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.UnloadingRates.Tbinned3DCIlow= NaN(NumExperiments, NumAPbins, 6);

     this.PerNucleusProfileParameters.Fits.AnaphaseAligned = cell(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.Fits.AnaphaseAligned3D = cell(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.Fits.Tbinned = cell(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.Fits.Tbinned3D = cell(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.Fits.Unaligned = cell(NumExperiments, NumAPbins, 6);
     this.PerNucleusProfileParameters.Fits.Unaligned3D = cell(NumExperiments, NumAPbins, 6);

     TraceTypes = {'fluo', 'fluo3d', 'anaphaseAligned', 'anaphaseAligned3d',  'tbinned', 'tbinned3d'};
     disp('Fitting trapezoids to mean traces.')
     for SetIndex = 1:length(this.ExperimentPrefixes)
         disp([num2str(SetIndex),'/',num2str(length(this.ExperimentPrefixes))])
         if ~ismember(SetIndex, this.ProcessedExperiments)
             continue
         end
         for NC = 9:14
             for APindex = 1:NumAPbins
                 %disp([num2str(SetIndex),', ', num2str(NC),', ', num2str(APindex)])
                 for tindex = 1:length(TraceTypes)
                     %disp(['NC: ', num2str(NC), ', APbin: ', num2str(APindex), ', TraceType: ', TraceTypes{tindex}])
                     this = FitTrapezoidsForPerNucleusMS2Traces(this, SetIndex, NC,...
                         APindex, TraceTypes{tindex});
                 end
             end
         end
     end
 end
        