 function this = AddTBinnedMeanPerNucleusInitiationRates(this)
     alpha = this.alpha;
     Temperatures = fliplr(unique(this.Temp_sps));
     NumTemperatures = length(Temperatures);
     NumExperiments = length(this.ExperimentPrefixes);
     APResolution = this.Experiments{1}.APResolution;
     NumAPbins = uint16(1/APResolution)+1;
     this.BinnedPerNucleusProfileParameters.MeanInitiationRates = {};
     this.BinnedPerNucleusProfileParameters.TimeOns = {};
     this.BinnedPerNucleusProfileParameters.TimeOffs = {};
     this.BinnedPerNucleusProfileParameters.ElongationTimes = {};
     this.BinnedPerNucleusProfileParameters.UnloadingRates = {};
     this.BinnedPerNucleusProfileParameters.MeanFitR2s = {};

     this.BinnedPerNucleusProfileParameters.MeanFitR2s.AnaphaseAligned = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.MeanFitR2s.AnaphaseAligned3D = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.MeanFitR2s.Tbinned = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.MeanFitR2s.Tbinned3D = NaN(NumTemperatures, NumAPbins, 6);

     this.BinnedPerNucleusProfileParameters.MeanInitiationRates.AnaphaseAligned = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.MeanInitiationRates.AnaphaseAlignedStdError = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.MeanInitiationRates.AnaphaseAlignedCIhigh = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.MeanInitiationRates.AnaphaseAlignedCIlow= NaN(NumTemperatures, NumAPbins, 6);

     this.BinnedPerNucleusProfileParameters.TimeOns.AnaphaseAligned = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.TimeOns.AnaphaseAlignedStdError = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.TimeOns.AnaphaseAlignedCIhigh = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.TimeOns.AnaphaseAlignedCIlow= NaN(NumTemperatures, NumAPbins, 6);

     this.BinnedPerNucleusProfileParameters.TimeOffs.AnaphaseAligned = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.TimeOffs.AnaphaseAlignedStdError = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.ElongationTimes.AnaphaseAligned = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.ElongationTimes.AnaphaseAlignedStdError = NaN(NumTemperatures, NumAPbins, 6);

     this.BinnedPerNucleusProfileParameters.UnloadingRates.AnaphaseAligned = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.UnloadingRates.AnaphaseAlignedStdError = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.UnloadingRates.AnaphaseAlignedCIhigh = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.UnloadingRates.AnaphaseAlignedCIlow= NaN(NumTemperatures, NumAPbins, 6);

     this.BinnedPerNucleusProfileParameters.MeanInitiationRates.Tbinned = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.MeanInitiationRates.TbinnedStdError = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.MeanInitiationRates.TbinnedCIhigh = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.MeanInitiationRates.TbinnedCIlow= NaN(NumTemperatures, NumAPbins, 6);

     this.BinnedPerNucleusProfileParameters.TimeOns.Tbinned = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.TimeOns.TbinnedStdError = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.TimeOns.TbinnedCIhigh = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.TimeOns.TbinnedCIlow= NaN(NumTemperatures, NumAPbins, 6);

     this.BinnedPerNucleusProfileParameters.TimeOffs.Tbinned = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.TimeOffs.TbinnedStdError = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.ElongationTimes.Tbinned = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.ElongationTimes.TbinnedStdError = NaN(NumTemperatures, NumAPbins, 6);

     this.BinnedPerNucleusProfileParameters.UnloadingRates.Tbinned = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.UnloadingRates.TbinnedStdError = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.UnloadingRates.TbinnedCIhigh = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.UnloadingRates.TbinnedCIlow= NaN(NumTemperatures, NumAPbins, 6);

     this.BinnedPerNucleusProfileParameters.MeanInitiationRates.AnaphaseAligned3D = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.MeanInitiationRates.AnaphaseAligned3DStdError = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.MeanInitiationRates.AnaphaseAligned3DCIhigh = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.MeanInitiationRates.AnaphaseAligned3DCIlow= NaN(NumTemperatures, NumAPbins, 6);

     this.BinnedPerNucleusProfileParameters.TimeOns.AnaphaseAligned3D = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.TimeOns.AnaphaseAligned3DStdError = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.TimeOns.AnaphaseAligned3DCIhigh = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.TimeOns.AnaphaseAligned3DCIlow= NaN(NumTemperatures, NumAPbins, 6);

     this.BinnedPerNucleusProfileParameters.TimeOffs.AnaphaseAligned3D = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.ElongationTimes.AnaphaseAligned3D = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.TimeOffs.AnaphaseAligned3DStdError = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.ElongationTimes.AnaphaseAligned3DStdError = NaN(NumTemperatures, NumAPbins, 6);

     this.BinnedPerNucleusProfileParameters.UnloadingRates.AnaphaseAligned3D = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.UnloadingRates.AnaphaseAligned3DStdError = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.UnloadingRates.AnaphaseAligned3DCIhigh = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.UnloadingRates.AnaphaseAligned3DCIlow= NaN(NumTemperatures, NumAPbins, 6);

     this.BinnedPerNucleusProfileParameters.MeanInitiationRates.Tbinned3D= NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.MeanInitiationRates.Tbinned3DStdError = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.MeanInitiationRates.Tbinned3DCIhigh = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.MeanInitiationRates.Tbinned3DCIlow= NaN(NumTemperatures, NumAPbins, 6);

     this.BinnedPerNucleusProfileParameters.TimeOns.Tbinned3D= NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.TimeOns.Tbinned3DStdError = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.TimeOns.Tbinned3DCIhigh = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.TimeOns.Tbinned3DCIlow= NaN(NumTemperatures, NumAPbins, 6);

     this.BinnedPerNucleusProfileParameters.TimeOffs.Tbinned3D= NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.ElongationTimes.Tbinned3D= NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.TimeOffs.Tbinned3DStdError= NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.ElongationTimes.Tbinned3DStdError= NaN(NumTemperatures, NumAPbins, 6);

     this.BinnedPerNucleusProfileParameters.UnloadingRates.Tbinned3D= NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.UnloadingRates.Tbinned3DStdError = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.UnloadingRates.Tbinned3DCIhigh = NaN(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.UnloadingRates.Tbinned3DCIlow= NaN(NumTemperatures, NumAPbins, 6);

     this.BinnedPerNucleusProfileParameters.Fits.AnaphaseAligned = cell(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.Fits.AnaphaseAligned3D = cell(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.Fits.Tbinned = cell(NumTemperatures, NumAPbins, 6);
     this.BinnedPerNucleusProfileParameters.Fits.Tbinned3D = cell(NumTemperatures, NumAPbins, 6);

     TraceTypes = {'anaphaseAligned', 'anaphaseAligned3d',  'tbinned', 'tbinned3d'};
     disp('Fitting trapezoids to mean traces.')
     for TempIndex = 1:NumTemperatures
         disp([num2str(TempIndex),'/',num2str(NumTemperatures)])
         for NC = 9:14
             for APindex = 1:NumAPbins
                 %disp([num2str(SetIndex),', ', num2str(NC),', ', num2str(APindex)])
                 for tindex = 1:length(TraceTypes)
                     this = FitTrapezoidsForTbinnedPerNucleusMS2Traces(this, TempIndex, NC,...
                         APindex, TraceTypes{tindex});
                 end
             end
         end
     end
     disp('Done')
 end
        