function [BinnedParams, BinnedSEParams, Counts, ParamTemperatures, ParamSETemperatures] = ...
    getBinnedPlottingVariables(this, PlottedParams, PlottedParamSEs,R2s, R2bound)
%%
if ~exist('R2s', 'var')
    R2s = ones(size(PlottedParams));
end
if ~exist('R2bound', 'var')
    R2bound = 0;
end

%%
if isfield(this, 'UniqueTemperatures')
    temperatures = this.UniqueTemperatures;
else
    temperatures = flip(unique(this.Temp_sps(this.ProcessedExperiments)));
end
temp_obs = this.Temp_obs;
se_temp_obs = 0.5*ones(1, length(temp_obs));
APResolution = this.Experiments{1}.APResolution;
Nbins = uint8(1/APResolution+1);

UseSetVector = ismember(1:size(PlottedParams, 1), this.ProcessedExperiments);

if ndims(PlottedParams) == 3
    BinnedParams = NaN(length(temperatures), Nbins, 6);
    BinnedSEParams = NaN(length(temperatures), Nbins, 6);
    ParamTemperatures = NaN(length(temperatures), Nbins, 6);
    ParamSETemperatures =NaN(length(temperatures), Nbins, 6);
    Counts = NaN(length(temperatures), Nbins, 6);
    for temp_idx = 1:length(temperatures)
        for NC = 9:14
            for APidx = 1:Nbins
                d1_matches = (this.Temp_sps == temperatures(temp_idx)) & ...
                    (UseSetVector == 1) & (R2s(:,APidx, NC-8).' >= R2bound);
                if sum(d1_matches) > 0
                    Counts(temp_idx, APidx, NC-8) = sum(d1_matches);
                    means = PlottedParams(d1_matches, APidx, NC-8).';
                    ses = PlottedParamSEs(d1_matches, APidx, NC-8).';
                    if ~all(isnan(ses))
                        BinnedParams(temp_idx, APidx, NC-8) = ...
                            (1/sum(1./(ses.^2)))*(sum(means./(ses.^2)));
                        BinnedSEParams(temp_idx, APidx, NC-8) = ...
                            sqrt((1/sum(1./(ses.^2))));
                    else
                        BinnedParams(temp_idx, APidx, NC-8) = mean(means);
                        BinnedSEParams(temp_idx, APidx, NC-8) = std(means)/sum(d1_matches);
                    end
                    
                    ParamTemperatures(temp_idx, APidx, NC-8) = mean(temp_obs(d1_matches));
                    ParamSETemperatures(temp_idx, APidx, NC-8) = std(temp_obs(d1_matches))/sum(d1_matches);
                    
                    
                end
            end
        end
        
    end
elseif ndims(PlottedParams) == 2
    BinnedParams = NaN(length(temperatures),  6);
    BinnedSEParams = NaN(length(temperatures), 6);
    ParamTemperatures = NaN(length(temperatures), 6);
    ParamSETemperatures =NaN(length(temperatures), 6);
    Counts = NaN(length(temperatures), 6);
    for temp_idx = 1:length(temperatures)
        for NC = 9:14
            
            d1_matches = (this.Temp_sps == temperatures(temp_idx)) & ...
                (UseSetVector == 1) & (~isnan(PlottedParams(:,NC-8).'));
            if sum(d1_matches) > 0
                Counts(temp_idx, NC-8) = sum(d1_matches);
                means = PlottedParams(d1_matches, NC-8).';
                ses = PlottedParamSEs(d1_matches,NC-8).';
                if ~all(isnan(ses))
                    BinnedParams(temp_idx, NC-8) = ...
                        (1/sum(1./(ses.^2)))*(sum(means./(ses.^2)));
                    BinnedSEParams(temp_idx, NC-8) = ...
                        sqrt((1/sum(1./(ses.^2))));
                else
                    BinnedParams(temp_idx, NC-8) = mean(means);
                    BinnedSEParams(temp_idx, NC-8) = std(means)/sum(d1_matches);
                end
                
                ParamTemperatures(temp_idx, NC-8) = mean(temp_obs(d1_matches));
                ParamSETemperatures(temp_idx, NC-8) = std(temp_obs(d1_matches))/sum(d1_matches);
                
                
                
            end
        end
        
    end
end
