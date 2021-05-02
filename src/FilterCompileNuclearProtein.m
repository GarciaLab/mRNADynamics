function CompiledNuclei = FilterCompileNuclearProtein(Prefix, FluoLabel, varargin)
%liveExperiment = LiveExperiment(Prefix);

liveExperiment = LiveExperiment(Prefix);
SummarizeEmbryoHealth(Prefix, false);
CompiledNuclei = CompileNuclearProteinFlexFluo(Prefix);
if ~exist('FluoLabel', 'var')
    CompiledNuclei = SetDefaultFluo(Prefix);
else
    CompiledNuclei = SetDefaultFluo(Prefix, FluoLabel);
end

CompiledNuclei = AddQCInfoToCompiledNuclei(Prefix, varargin{:});

try
    save([liveExperiment.resultsFolder,filesep,'CompiledNuclei.mat'],...
        'CompiledNuclei','-v6');
catch
    save([liveExperiment.resultsFolder,filesep,'CompiledNuclei.mat'],...
        'CompiledNuclei','-v7.3', '-nocompression');
end
