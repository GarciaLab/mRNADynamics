function TimingInfo = ExtractTimingCoeffs(ltm)
%%
EmbryoDivisionTimes = ltm.EmbryoStats.NCDivisionInfo;
NumEmbryos = size(EmbryoDivisionTimes, 1);
EmbryoTemperatures = ltm.Temp_sps;
IncludeEmbryos = ismember(1:NumEmbryos, ltm.IncludedExperiments);
Temperatures = fliplr(unique(ltm.Temp_sps));
NumTemperatures = length(Temperatures);
PropNCDivisionInfo = NaN(NumTemperatures, 6);
MeanNCDivisionInfo = NaN(NumTemperatures, 6);
StdNCDivisionInfo = NaN(NumTemperatures, 6);
CountNCDivisionInfo = NaN(NumTemperatures, 6);
TimingCoeffs = NaN(1, NumTemperatures);
for t_index = 1:NumTemperatures
    TemperatureDivisionTimes = EmbryoDivisionTimes(EmbryoTemperatures == Temperatures(t_index) & IncludeEmbryos,:);
    CountNCDivisionInfo(t_index,:) = sum(~isnan(TemperatureDivisionTimes), 1);
    MeanNCDivisionInfo(t_index,:) = mean(TemperatureDivisionTimes, 1, 'omitnan');
    StdNCDivisionInfo(t_index,:) = std(TemperatureDivisionTimes, 1, 'omitnan'); 
end
PropNCDivisionInfo = MeanNCDivisionInfo./MeanNCDivisionInfo(find(Temperatures == 25, 1),:);
TimeCoeffs = mean(1./PropNCDivisionInfo(:,4:5), 2).';
TimingInfo = {};
TimingInfo.Temperatures = Temperatures;
TimingInfo.PropNCDivisionInfo = PropNCDivisionInfo;
TimingInfo.MeanNCDivisionInfo = MeanNCDivisionInfo;
TimingInfo.StdNCDivisionInfo = StdNCDivisionInfo;
TimingInfo.CountNCDivisionInfo = CountNCDivisionInfo;
TimingInfo.TimeCoeffs = TimeCoeffs;

