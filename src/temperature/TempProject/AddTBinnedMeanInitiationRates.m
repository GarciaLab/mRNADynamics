 function this = AddTBinnedMeanInitiationRates(this)
     alpha = this.alpha;
     Temperatures = fliplr(unique(this.Temp_sps));
     NumTemperatures = length(Temperatures);
     NumExperiments = length(this.ExperimentPrefixes);
     APResolution = this.Experiments{1}.APResolution;
     NumAPbins = uint16(1/APResolution)+1;
     this.BinnedProfileParameters.MeanInitiationRates = {};
     this.BinnedProfileParameters.TimeOns = {};
     this.BinnedProfileParameters.TimeOffs = {};
     this.BinnedProfileParameters.ElongationTimes = {};
     this.BinnedProfileParameters.UnloadingRates = {};
     this.BinnedProfileParameters.MeanFitR2s = {};

     this.BinnedProfileParameters.MeanFitR2s.AnaphaseAligned = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.MeanFitR2s.AnaphaseAligned3D = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.MeanFitR2s.Tbinned = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.MeanFitR2s.Tbinned3D = NaN(NumTemperatures, NumAPbins, 6);

     this.BinnedProfileParameters.MeanInitiationRates.AnaphaseAligned = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.MeanInitiationRates.AnaphaseAlignedStdError = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.MeanInitiationRates.AnaphaseAlignedCIhigh = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.MeanInitiationRates.AnaphaseAlignedCIlow= NaN(NumTemperatures, NumAPbins, 6);

     this.BinnedProfileParameters.TimeOns.AnaphaseAligned = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.TimeOns.AnaphaseAlignedStdError = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.TimeOns.AnaphaseAlignedCIhigh = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.TimeOns.AnaphaseAlignedCIlow= NaN(NumTemperatures, NumAPbins, 6);

     this.BinnedProfileParameters.TimeOffs.AnaphaseAligned = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.TimeOffs.AnaphaseAlignedStdError = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.ElongationTimes.AnaphaseAligned = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.ElongationTimes.AnaphaseAlignedStdError = NaN(NumTemperatures, NumAPbins, 6);

     this.BinnedProfileParameters.UnloadingRates.AnaphaseAligned = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.UnloadingRates.AnaphaseAlignedStdError = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.UnloadingRates.AnaphaseAlignedCIhigh = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.UnloadingRates.AnaphaseAlignedCIlow= NaN(NumTemperatures, NumAPbins, 6);

     this.BinnedProfileParameters.MeanInitiationRates.Tbinned = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.MeanInitiationRates.TbinnedStdError = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.MeanInitiationRates.TbinnedCIhigh = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.MeanInitiationRates.TbinnedCIlow= NaN(NumTemperatures, NumAPbins, 6);

     this.BinnedProfileParameters.TimeOns.Tbinned = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.TimeOns.TbinnedStdError = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.TimeOns.TbinnedCIhigh = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.TimeOns.TbinnedCIlow= NaN(NumTemperatures, NumAPbins, 6);

     this.BinnedProfileParameters.TimeOffs.Tbinned = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.TimeOffs.TbinnedStdError = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.ElongationTimes.Tbinned = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.ElongationTimes.TbinnedStdError = NaN(NumTemperatures, NumAPbins, 6);

     this.BinnedProfileParameters.UnloadingRates.Tbinned = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.UnloadingRates.TbinnedStdError = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.UnloadingRates.TbinnedCIhigh = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.UnloadingRates.TbinnedCIlow= NaN(NumTemperatures, NumAPbins, 6);

     this.BinnedProfileParameters.MeanInitiationRates.AnaphaseAligned3D = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.MeanInitiationRates.AnaphaseAligned3DStdError = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.MeanInitiationRates.AnaphaseAligned3DCIhigh = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.MeanInitiationRates.AnaphaseAligned3DCIlow= NaN(NumTemperatures, NumAPbins, 6);

     this.BinnedProfileParameters.TimeOns.AnaphaseAligned3D = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.TimeOns.AnaphaseAligned3DStdError = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.TimeOns.AnaphaseAligned3DCIhigh = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.TimeOns.AnaphaseAligned3DCIlow= NaN(NumTemperatures, NumAPbins, 6);

     this.BinnedProfileParameters.TimeOffs.AnaphaseAligned3D = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.ElongationTimes.AnaphaseAligned3D = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.TimeOffs.AnaphaseAligned3DStdError = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.ElongationTimes.AnaphaseAligned3DStdError = NaN(NumTemperatures, NumAPbins, 6);

     this.BinnedProfileParameters.UnloadingRates.AnaphaseAligned3D = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.UnloadingRates.AnaphaseAligned3DStdError = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.UnloadingRates.AnaphaseAligned3DCIhigh = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.UnloadingRates.AnaphaseAligned3DCIlow= NaN(NumTemperatures, NumAPbins, 6);

     this.BinnedProfileParameters.MeanInitiationRates.Tbinned3D= NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.MeanInitiationRates.Tbinned3DStdError = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.MeanInitiationRates.Tbinned3DCIhigh = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.MeanInitiationRates.Tbinned3DCIlow= NaN(NumTemperatures, NumAPbins, 6);

     this.BinnedProfileParameters.TimeOns.Tbinned3D= NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.TimeOns.Tbinned3DStdError = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.TimeOns.Tbinned3DCIhigh = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.TimeOns.Tbinned3DCIlow= NaN(NumTemperatures, NumAPbins, 6);

     this.BinnedProfileParameters.TimeOffs.Tbinned3D= NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.ElongationTimes.Tbinned3D= NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.TimeOffs.Tbinned3DStdError= NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.ElongationTimes.Tbinned3DStdError= NaN(NumTemperatures, NumAPbins, 6);

     this.BinnedProfileParameters.UnloadingRates.Tbinned3D= NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.UnloadingRates.Tbinned3DStdError = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.UnloadingRates.Tbinned3DCIhigh = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.UnloadingRates.Tbinned3DCIlow= NaN(NumTemperatures, NumAPbins, 6);

     this.BinnedProfileParameters.Fits.AnaphaseAligned = cell(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.Fits.AnaphaseAligned3D = cell(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.Fits.Tbinned = cell(NumTemperatures, NumAPbins, 6);
     this.BinnedProfileParameters.Fits.Tbinned3D = cell(NumTemperatures, NumAPbins, 6);

     TraceTypes = {'anaphaseAligned', 'anaphaseAligned3d',  'tbinned', 'tbinned3d'};
     disp('Fitting trapezoids to mean traces.')
     for TempIndex = 1:NumTemperatures
         disp([num2str(TempIndex),'/',num2str(NumTemperatures)])
         for NC = 9:14
             for APindex = 1:NumAPbins
                 %disp([num2str(SetIndex),', ', num2str(NC),', ', num2str(APindex)])
                 for tindex = 1:length(TraceTypes)
                     this = FitTrapezoidsForTbinnedMS2Traces(this, TempIndex, NC,...
                         APindex, TraceTypes{tindex});
                 end
             end
         end
     end
     disp('Done')
 end
        