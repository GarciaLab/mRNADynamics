 function this = AddMeanInitiationRates(this)
     alpha = this.alpha;
     NumExperiments = length(this.ExperimentPrefixes);
     APResolution = this.Experiments{1}.APResolution;
     NumAPbins = uint16(1/APResolution)+1;
     this.MeanInitiationRates = {};
     this.TimeOns = {};
     this.TimeOffs = {};
     this.ElongationTimes = {};
     this.UnloadingRates = {};
     this.MeanFitR2s = {};

     this.MeanFitR2s.AnaphaseAligned = NaN(NumExperiments, NumAPbins, 6);
     this.MeanFitR2s.AnaphaseAligned3D = NaN(NumExperiments, NumAPbins, 6);
     this.MeanFitR2s.Unaligned = NaN(NumExperiments, NumAPbins, 6);
     this.MeanFitR2s.Unaligned3D = NaN(NumExperiments, NumAPbins, 6);
     this.MeanFitR2s.Tbinned = NaN(NumExperiments, NumAPbins, 6);
     this.MeanFitR2s.Tbinned3D = NaN(NumExperiments, NumAPbins, 6);

     this.MeanInitiationRates.AnaphaseAligned = NaN(NumExperiments, NumAPbins, 6);
     this.MeanInitiationRates.AnaphaseAlignedStdError = NaN(NumExperiments, NumAPbins, 6);
     this.MeanInitiationRates.AnaphaseAlignedCIhigh = NaN(NumExperiments, NumAPbins, 6);
     this.MeanInitiationRates.AnaphaseAlignedCIlow= NaN(NumExperiments, NumAPbins, 6);

     this.TimeOns.AnaphaseAligned = NaN(NumExperiments, NumAPbins, 6);
     this.TimeOns.AnaphaseAlignedStdError = NaN(NumExperiments, NumAPbins, 6);
     this.TimeOns.AnaphaseAlignedCIhigh = NaN(NumExperiments, NumAPbins, 6);
     this.TimeOns.AnaphaseAlignedCIlow= NaN(NumExperiments, NumAPbins, 6);

     this.TimeOffs.AnaphaseAligned = NaN(NumExperiments, NumAPbins, 6);
     this.TimeOffs.AnaphaseAlignedStdError = NaN(NumExperiments, NumAPbins, 6);
     this.ElongationTimes.AnaphaseAligned = NaN(NumExperiments, NumAPbins, 6);
     this.ElongationTimes.AnaphaseAlignedStdError = NaN(NumExperiments, NumAPbins, 6);

     this.UnloadingRates.AnaphaseAligned = NaN(NumExperiments, NumAPbins, 6);
     this.UnloadingRates.AnaphaseAlignedStdError = NaN(NumExperiments, NumAPbins, 6);
     this.UnloadingRates.AnaphaseAlignedCIhigh = NaN(NumExperiments, NumAPbins, 6);
     this.UnloadingRates.AnaphaseAlignedCIlow= NaN(NumExperiments, NumAPbins, 6);

     this.MeanInitiationRates.Tbinned = NaN(NumExperiments, NumAPbins, 6);
     this.MeanInitiationRates.TbinnedStdError = NaN(NumExperiments, NumAPbins, 6);
     this.MeanInitiationRates.TbinnedCIhigh = NaN(NumExperiments, NumAPbins, 6);
     this.MeanInitiationRates.TbinnedCIlow= NaN(NumExperiments, NumAPbins, 6);

     this.TimeOns.Tbinned = NaN(NumExperiments, NumAPbins, 6);
     this.TimeOns.TbinnedStdError = NaN(NumExperiments, NumAPbins, 6);
     this.TimeOns.TbinnedCIhigh = NaN(NumExperiments, NumAPbins, 6);
     this.TimeOns.TbinnedCIlow= NaN(NumExperiments, NumAPbins, 6);

     this.TimeOffs.Tbinned = NaN(NumExperiments, NumAPbins, 6);
     this.TimeOffs.TbinnedStdError = NaN(NumExperiments, NumAPbins, 6);
     this.ElongationTimes.Tbinned = NaN(NumExperiments, NumAPbins, 6);
     this.ElongationTimes.TbinnedStdError = NaN(NumExperiments, NumAPbins, 6);

     this.UnloadingRates.Tbinned = NaN(NumExperiments, NumAPbins, 6);
     this.UnloadingRates.TbinnedStdError = NaN(NumExperiments, NumAPbins, 6);
     this.UnloadingRates.TbinnedCIhigh = NaN(NumExperiments, NumAPbins, 6);
     this.UnloadingRates.TbinnedCIlow= NaN(NumExperiments, NumAPbins, 6);

     this.MeanInitiationRates.Unaligned = NaN(NumExperiments, NumAPbins, 6);
     this.MeanInitiationRates.UnalignedStdError = NaN(NumExperiments, NumAPbins, 6);
     this.MeanInitiationRates.UnalignedCIhigh = NaN(NumExperiments, NumAPbins, 6);
     this.MeanInitiationRates.UnalignedCIlow= NaN(NumExperiments, NumAPbins, 6);

     this.TimeOns.Unaligned = NaN(NumExperiments, NumAPbins, 6);
     this.TimeOns.UnalignedStdError = NaN(NumExperiments, NumAPbins, 6);
     this.TimeOns.UnalignedCIhigh = NaN(NumExperiments, NumAPbins, 6);
     this.TimeOns.UnalignedCIlow= NaN(NumExperiments, NumAPbins, 6);

     this.TimeOffs.Unaligned = NaN(NumExperiments, NumAPbins, 6);
     this.ElongationTimes.Unaligned = NaN(NumExperiments, NumAPbins, 6);
     this.TimeOffs.UnalignedStdError = NaN(NumExperiments, NumAPbins, 6);
     this.ElongationTimes.UnalignedStdError = NaN(NumExperiments, NumAPbins, 6);

     this.UnloadingRates.Unaligned = NaN(NumExperiments, NumAPbins, 6);
     this.UnloadingRates.UnalignedStdError = NaN(NumExperiments, NumAPbins, 6);
     this.UnloadingRates.UnalignedCIhigh = NaN(NumExperiments, NumAPbins, 6);
     this.UnloadingRates.UnalignedCIlow= NaN(NumExperiments, NumAPbins, 6);

     this.MeanInitiationRates.AnaphaseAligned3D = NaN(NumExperiments, NumAPbins, 6);
     this.MeanInitiationRates.AnaphaseAligned3DStdError = NaN(NumExperiments, NumAPbins, 6);
     this.MeanInitiationRates.AnaphaseAligned3DCIhigh = NaN(NumExperiments, NumAPbins, 6);
     this.MeanInitiationRates.AnaphaseAligned3DCIlow= NaN(NumExperiments, NumAPbins, 6);

     this.TimeOns.AnaphaseAligned3D = NaN(NumExperiments, NumAPbins, 6);
     this.TimeOns.AnaphaseAligned3DStdError = NaN(NumExperiments, NumAPbins, 6);
     this.TimeOns.AnaphaseAligned3DCIhigh = NaN(NumExperiments, NumAPbins, 6);
     this.TimeOns.AnaphaseAligned3DCIlow= NaN(NumExperiments, NumAPbins, 6);

     this.TimeOffs.AnaphaseAligned3D = NaN(NumExperiments, NumAPbins, 6);
     this.ElongationTimes.AnaphaseAligned3D = NaN(NumExperiments, NumAPbins, 6);
     this.TimeOffs.AnaphaseAligned3DStdError = NaN(NumExperiments, NumAPbins, 6);
     this.ElongationTimes.AnaphaseAligned3DStdError = NaN(NumExperiments, NumAPbins, 6);

     this.UnloadingRates.AnaphaseAligned3D = NaN(NumExperiments, NumAPbins, 6);
     this.UnloadingRates.AnaphaseAligned3DStdError = NaN(NumExperiments, NumAPbins, 6);
     this.UnloadingRates.AnaphaseAligned3DCIhigh = NaN(NumExperiments, NumAPbins, 6);
     this.UnloadingRates.AnaphaseAligned3DCIlow= NaN(NumExperiments, NumAPbins, 6);

     this.MeanInitiationRates.Tbinned3D= NaN(NumExperiments, NumAPbins, 6);
     this.MeanInitiationRates.Tbinned3DStdError = NaN(NumExperiments, NumAPbins, 6);
     this.MeanInitiationRates.Tbinned3DCIhigh = NaN(NumExperiments, NumAPbins, 6);
     this.MeanInitiationRates.Tbinned3DCIlow= NaN(NumExperiments, NumAPbins, 6);

     this.TimeOns.Tbinned3D= NaN(NumExperiments, NumAPbins, 6);
     this.TimeOns.Tbinned3DStdError = NaN(NumExperiments, NumAPbins, 6);
     this.TimeOns.Tbinned3DCIhigh = NaN(NumExperiments, NumAPbins, 6);
     this.TimeOns.Tbinned3DCIlow= NaN(NumExperiments, NumAPbins, 6);

     this.TimeOffs.Tbinned3D= NaN(NumExperiments, NumAPbins, 6);
     this.ElongationTimes.Tbinned3D= NaN(NumExperiments, NumAPbins, 6);
     this.TimeOffs.Tbinned3DStdError= NaN(NumExperiments, NumAPbins, 6);
     this.ElongationTimes.Tbinned3DStdError= NaN(NumExperiments, NumAPbins, 6);

     this.UnloadingRates.Tbinned3D= NaN(NumExperiments, NumAPbins, 6);
     this.UnloadingRates.Tbinned3DStdError = NaN(NumExperiments, NumAPbins, 6);
     this.UnloadingRates.Tbinned3DCIhigh = NaN(NumExperiments, NumAPbins, 6);
     this.UnloadingRates.Tbinned3DCIlow= NaN(NumExperiments, NumAPbins, 6);

     this.MeanInitiationRates.Unaligned3D = NaN(NumExperiments, NumAPbins, 6);
     this.MeanInitiationRates.Unaligned3DStdError = NaN(NumExperiments, NumAPbins, 6);
     this.MeanInitiationRates.Unaligned3DCIhigh = NaN(NumExperiments, NumAPbins, 6);
     this.MeanInitiationRates.Unaligned3DCIlow= NaN(NumExperiments, NumAPbins, 6);

     this.TimeOns.Unaligned3D = NaN(NumExperiments, NumAPbins, 6);
     this.TimeOns.Unaligned3DStdError = NaN(NumExperiments, NumAPbins, 6);
     this.TimeOns.Unaligned3DCIhigh = NaN(NumExperiments, NumAPbins, 6);
     this.TimeOns.Unaligned3DCIlow= NaN(NumExperiments, NumAPbins, 6);

     this.TimeOffs.Unaligned3D = NaN(NumExperiments, NumAPbins, 6);
     this.ElongationTimes.Unaligned3D = NaN(NumExperiments, NumAPbins, 6);
     this.TimeOffs.Unaligned3DStdError = NaN(NumExperiments, NumAPbins, 6);
     this.ElongationTimes.Unaligned3DStdError = NaN(NumExperiments, NumAPbins, 6);

     this.UnloadingRates.Unaligned3D = NaN(NumExperiments, NumAPbins, 6);
     this.UnloadingRates.Unaligned3DStdError = NaN(NumExperiments, NumAPbins, 6);
     this.UnloadingRates.Unaligned3DCIhigh = NaN(NumExperiments, NumAPbins, 6);
     this.UnloadingRates.Unaligned3DCIlow= NaN(NumExperiments, NumAPbins, 6);

     this.Fits.Unaligned = cell(NumExperiments, NumAPbins, 6);
     this.Fits.Unaligned3D = cell(NumExperiments, NumAPbins, 6);
     this.Fits.AnaphaseAligned = cell(NumExperiments, NumAPbins, 6);
     this.Fits.AnaphaseAligned3D = cell(NumExperiments, NumAPbins, 6);
     this.Fits.Tbinned = cell(NumExperiments, NumAPbins, 6);
     this.Fits.Tbinned3D = cell(NumExperiments, NumAPbins, 6);

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
                     this = FitTrapezoidsForMeanMS2Traces(this, SetIndex, NC,...
                         APindex, TraceTypes{tindex}, UseWeights, alpha);
                 end
             end
         end
     end
     disp('Done')
 end
        