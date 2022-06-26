function CompiledEmbryos = AddBinnedProfiles(CompiledEmbryos, deltafc_binwidth, dubuistime_binwidth)
%%
if ~exist('deltafc_binwidth', 'var')
    deltafc_binwidth = 2.5;
end
if ~exist('dubuistime_binwidth', 'var')
    dubuistime_binwidth = 3;
end
if ~exist('yw25Ctime_binwidth', 'var')
    yw25Ctime_binwidth = 5;
end
if ~exist('hisrfp25Ctime_binwidth', 'var')
    hisrfp25Ctime_binwidth = 5;
end
APbins = 0:0.025:1;
NumAPbins = length(APbins);
NarrowAPbins = 0:0.0125:1;
NumNarrowAPbins = length(NarrowAPbins);
DorsalProfiles = CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles;
NChannels = size(DorsalProfiles, 3);
CompiledEmbryos.BinnedProfiles = {};
CompiledEmbryos.WindowedProfiles = {};
CompiledEmbryos.NarrowBinnedProfiles = {};
CompiledEmbryos.NarrowWindowedProfiles = {};
%% Add DeltaFC-based Binning - standard AP
AllDeltaValidProfilesTestTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.TestSetEmbryos & ~isnan([CompiledEmbryos.FixCorrectedDeltaFC_um.mean(:)].');% &...
AllDeltaValidProfilesControlTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.ControlSetEmbryos & ~isnan([CompiledEmbryos.FixCorrectedDeltaFC_um.mean(:)].');% &...


x = 0:deltafc_binwidth:45;
NumTbins = length(x);

DorsalProfiles = CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles;
CompiledEmbryos.BinnedProfiles.DeltaFC = {};
CompiledEmbryos.BinnedProfiles.DeltaFC.BinWidth = deltafc_binwidth;
CompiledEmbryos.BinnedProfiles.DeltaFC.x = x;
CompiledEmbryos.BinnedProfiles.DeltaFC.Test = {};
CompiledEmbryos.BinnedProfiles.DeltaFC.Test.mean = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.BinnedProfiles.DeltaFC.Test.std = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.BinnedProfiles.DeltaFC.Test.se = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.BinnedProfiles.DeltaFC.Test.count = zeros(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.BinnedProfiles.DeltaFC.Control = {}; 
CompiledEmbryos.BinnedProfiles.DeltaFC.Control.mean = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.BinnedProfiles.DeltaFC.Control.std = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.BinnedProfiles.DeltaFC.Control.se = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.BinnedProfiles.DeltaFC.Control.count = zeros(NumTbins, NumAPbins, NChannels);


for i = 1:NumTbins
    TFtimeBin = [CompiledEmbryos.FixCorrectedDeltaFC_um.mean(:)].' >= x(i)-deltafc_binwidth/2 & [CompiledEmbryos.FixCorrectedDeltaFC_um.mean(:)].'< x(i)+deltafc_binwidth/2;
    TFtestBin = TFtimeBin & AllDeltaValidProfilesTestTF;
    TFcontrolBin = TFtimeBin & AllDeltaValidProfilesControlTF;
    for ch_index = 2:NChannels
        if sum(TFtestBin) >= 1
            for ap_index = 1:NumAPbins
                TFGoodAP = ~isnan(DorsalProfiles(TFtestBin,ap_index,ch_index));
                if sum(TFGoodAP) > 1
                    CompiledEmbryos.BinnedProfiles.DeltaFC.Test.mean(i,ap_index,ch_index) = mean(DorsalProfiles(TFtestBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.BinnedProfiles.DeltaFC.Test.std(i,ap_index,ch_index) = std(DorsalProfiles(TFtestBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.BinnedProfiles.DeltaFC.Test.count(i,ap_index,ch_index) = sum(~isnan(DorsalProfiles(TFtestBin,ap_index,ch_index)));
                    CompiledEmbryos.BinnedProfiles.DeltaFC.Test.se(i,ap_index,ch_index) = CompiledEmbryos.BinnedProfiles.DeltaFC.Test.std(i,ap_index,ch_index)/...
                        sqrt(CompiledEmbryos.BinnedProfiles.DeltaFC.Test.count(i,ap_index,ch_index));
                elseif sum(TFGoodAP)  == 1
                    TestProfs = DorsalProfiles(TFtestBin,ap_index,ch_index);
                    CompiledEmbryos.BinnedProfiles.DeltaFC.Test.mean(i,ap_index,ch_index) = TestProfs(TFGoodAP);
                    CompiledEmbryos.BinnedProfiles.DeltaFC.Test.std(i,ap_index,ch_index) = NaN;
                    CompiledEmbryos.BinnedProfiles.DeltaFC.Test.count(i,ap_index,ch_index) = 1;
                    CompiledEmbryos.BinnedProfiles.DeltaFC.Test.se(i,ap_index,ch_index) = NaN;
                end
            end
        end
        
        if sum(TFcontrolBin) >= 1
            for ap_index = 1:NumAPbins
                TFGoodAP = ~isnan(DorsalProfiles(TFcontrolBin,ap_index,ch_index));
                if sum(TFGoodAP) > 1
                    CompiledEmbryos.BinnedProfiles.DeltaFC.Control.mean(i,ap_index,ch_index) = mean(DorsalProfiles(TFcontrolBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.BinnedProfiles.DeltaFC.Control.std(i,ap_index,ch_index) = std(DorsalProfiles(TFcontrolBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.BinnedProfiles.DeltaFC.Control.count(i,ap_index,ch_index) = sum(~isnan(DorsalProfiles(TFcontrolBin,ap_index,ch_index)));
                    CompiledEmbryos.BinnedProfiles.DeltaFC.Control.se(i,ap_index,ch_index) = CompiledEmbryos.BinnedProfiles.DeltaFC.Control.std(i,ap_index,ch_index)/...
                        sqrt(CompiledEmbryos.BinnedProfiles.DeltaFC.Control.count(i,ap_index,ch_index));
                elseif sum(TFGoodAP)  == 1
                    ControlProfs = DorsalProfiles(TFcontrolBin,ap_index,ch_index);
                    CompiledEmbryos.BinnedProfiles.DeltaFC.Control.mean(i,ap_index,ch_index) = ControlProfs(TFGoodAP);
                    CompiledEmbryos.BinnedProfiles.DeltaFC.Control.std(i,ap_index,ch_index) = NaN;
                    CompiledEmbryos.BinnedProfiles.DeltaFC.Control.count(i,ap_index,ch_index) = 1;
                    CompiledEmbryos.BinnedProfiles.DeltaFC.Control.se(i,ap_index,ch_index) = NaN;
                end
            end
        end
    end
end

%% Add DeltaFC-based Binning - narrow AP
NarrowDorsalProfiles = CompiledEmbryos.SlideRescaledDorsalAvgNarrowAPProfiles;
CompiledEmbryos.NarrowBinnedProfiles.DeltaFC = {};
CompiledEmbryos.NarrowBinnedProfiles.DeltaFC.x = x;
CompiledEmbryos.NarrowBinnedProfiles.DeltaFC.BinWidth = deltafc_binwidth;
CompiledEmbryos.NarrowBinnedProfiles.DeltaFC.Test = {};
CompiledEmbryos.NarrowBinnedProfiles.DeltaFC.Test.mean = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowBinnedProfiles.DeltaFC.Test.std = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowBinnedProfiles.DeltaFC.Test.se = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowBinnedProfiles.DeltaFC.Test.count = zeros(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowBinnedProfiles.DeltaFC.Control = {}; 
CompiledEmbryos.NarrowBinnedProfiles.DeltaFC.Control.mean = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowBinnedProfiles.DeltaFC.Control.std = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowBinnedProfiles.DeltaFC.Control.se = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowBinnedProfiles.DeltaFC.Control.count = zeros(NumTbins, NumNarrowAPbins, NChannels);


for i = 1:NumTbins
    TFtimeBin = [CompiledEmbryos.FixCorrectedDeltaFC_um.mean(:)].' >= x(i)-deltafc_binwidth/2  & [CompiledEmbryos.FixCorrectedDeltaFC_um.mean(:)].'< x(i)+deltafc_binwidth/2 ;
    TFtestBin = TFtimeBin & AllDeltaValidProfilesTestTF;
    TFcontrolBin = TFtimeBin & AllDeltaValidProfilesControlTF;
    for ch_index = 2:NChannels
        if sum(TFtestBin) >= 1
            for ap_index = 1:NumNarrowAPbins
                TFGoodAP = ~isnan(NarrowDorsalProfiles(TFtestBin,ap_index,ch_index));
                if sum(TFGoodAP) > 1
                    CompiledEmbryos.NarrowBinnedProfiles.DeltaFC.Test.mean(i,ap_index,ch_index) = mean(NarrowDorsalProfiles(TFtestBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.NarrowBinnedProfiles.DeltaFC.Test.std(i,ap_index,ch_index) = std(NarrowDorsalProfiles(TFtestBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.NarrowBinnedProfiles.DeltaFC.Test.count(i,ap_index,ch_index) = sum(~isnan(NarrowDorsalProfiles(TFtestBin,ap_index,ch_index)));
                    CompiledEmbryos.NarrowBinnedProfiles.DeltaFC.Test.se(i,ap_index,ch_index) = CompiledEmbryos.NarrowBinnedProfiles.DeltaFC.Test.std(i,ap_index,ch_index)/...
                        sqrt(CompiledEmbryos.NarrowBinnedProfiles.DeltaFC.Test.count(i,ap_index,ch_index));
                elseif sum(TFGoodAP)  == 1
                    NarrowTestProfs = NarrowDorsalProfiles(TFtestBin,ap_index,ch_index);
                    CompiledEmbryos.NarrowBinnedProfiles.DeltaFC.Test.mean(i,ap_index,ch_index) = NarrowTestProfs(TFGoodAP);
                    CompiledEmbryos.NarrowBinnedProfiles.DeltaFC.Test.std(i,ap_index,ch_index) = NaN;
                    CompiledEmbryos.NarrowBinnedProfiles.DeltaFC.Test.count(i,ap_index,ch_index) = 1;
                    CompiledEmbryos.NarrowBinnedProfiles.DeltaFC.Test.se(i,ap_index,ch_index) = NaN;
                end
            end
        end
        
        if sum(TFcontrolBin) >= 1
            for ap_index = 1:NumNarrowAPbins
                TFGoodAP = ~isnan(NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index));
                if sum(TFGoodAP) > 1
                    CompiledEmbryos.NarrowBinnedProfiles.DeltaFC.Control.mean(i,ap_index,ch_index) = mean(NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.NarrowBinnedProfiles.DeltaFC.Control.std(i,ap_index,ch_index) = std(NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.NarrowBinnedProfiles.DeltaFC.Control.count(i,ap_index,ch_index) = sum(~isnan(NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index)));
                    CompiledEmbryos.NarrowBinnedProfiles.DeltaFC.Control.se(i,ap_index,ch_index) = CompiledEmbryos.NarrowBinnedProfiles.DeltaFC.Control.std(i,ap_index,ch_index)/...
                        sqrt(CompiledEmbryos.NarrowBinnedProfiles.DeltaFC.Control.count(i,ap_index,ch_index));
                elseif sum(TFGoodAP)  == 1
                    NarrowControlProfs = NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index);
                    CompiledEmbryos.NarrowBinnedProfiles.DeltaFC.Control.mean(i,ap_index,ch_index) = NarrowControlProfs(TFGoodAP);
                    CompiledEmbryos.NarrowBinnedProfiles.DeltaFC.Control.std(i,ap_index,ch_index) = NaN;
                    CompiledEmbryos.NarrowBinnedProfiles.DeltaFC.Control.count(i,ap_index,ch_index) = 1;
                    CompiledEmbryos.NarrowBinnedProfiles.DeltaFC.Control.se(i,ap_index,ch_index) = NaN;
                end
            end
        end
    end
end

%% Windowed Profiles - Delta FC

x = 0:1:45;
NumTbins = length(x);

DorsalProfiles = CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles;
CompiledEmbryos.WindowedProfiles.DeltaFC = {};
CompiledEmbryos.WindowedProfiles.DeltaFC.BinWidth = deltafc_binwidth;
CompiledEmbryos.WindowedProfiles.DeltaFC.x = x;
CompiledEmbryos.WindowedProfiles.DeltaFC.Test = {};
CompiledEmbryos.WindowedProfiles.DeltaFC.Test.mean = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.WindowedProfiles.DeltaFC.Test.std = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.WindowedProfiles.DeltaFC.Test.se = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.WindowedProfiles.DeltaFC.Test.count = zeros(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.WindowedProfiles.DeltaFC.Control = {}; 
CompiledEmbryos.WindowedProfiles.DeltaFC.Control.mean = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.WindowedProfiles.DeltaFC.Control.std = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.WindowedProfiles.DeltaFC.Control.se = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.WindowedProfiles.DeltaFC.Control.count = zeros(NumTbins, NumAPbins, NChannels);


for i = 1:NumTbins
    TFtimeBin = [CompiledEmbryos.FixCorrectedDeltaFC_um.mean(:)].' >= x(i)-deltafc_binwidth/2 & [CompiledEmbryos.FixCorrectedDeltaFC_um.mean(:)].'< x(i)+deltafc_binwidth/2;
    TFtestBin = TFtimeBin & AllDeltaValidProfilesTestTF;
    TFcontrolBin = TFtimeBin & AllDeltaValidProfilesControlTF;
    for ch_index = 2:NChannels
        if sum(TFtestBin) >= 1
            for ap_index = 1:NumAPbins
                TFGoodAP = ~isnan(DorsalProfiles(TFtestBin,ap_index,ch_index));
                if sum(TFGoodAP) > 1
                    CompiledEmbryos.WindowedProfiles.DeltaFC.Test.mean(i,ap_index,ch_index) = mean(DorsalProfiles(TFtestBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.WindowedProfiles.DeltaFC.Test.std(i,ap_index,ch_index) = std(DorsalProfiles(TFtestBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.WindowedProfiles.DeltaFC.Test.count(i,ap_index,ch_index) = sum(~isnan(DorsalProfiles(TFtestBin,ap_index,ch_index)));
                    CompiledEmbryos.WindowedProfiles.DeltaFC.Test.se(i,ap_index,ch_index) = CompiledEmbryos.WindowedProfiles.DeltaFC.Test.std(i,ap_index,ch_index)/...
                        sqrt(CompiledEmbryos.WindowedProfiles.DeltaFC.Test.count(i,ap_index,ch_index));
                elseif sum(TFGoodAP)  == 1
                    TestProfs = DorsalProfiles(TFtestBin,ap_index,ch_index);
                    CompiledEmbryos.WindowedProfiles.DeltaFC.Test.mean(i,ap_index,ch_index) = TestProfs(TFGoodAP);
                    CompiledEmbryos.WindowedProfiles.DeltaFC.Test.std(i,ap_index,ch_index) = NaN;
                    CompiledEmbryos.WindowedProfiles.DeltaFC.Test.count(i,ap_index,ch_index) = 1;
                    CompiledEmbryos.WindowedProfiles.DeltaFC.Test.se(i,ap_index,ch_index) = NaN;
                end
            end
        end
        
        if sum(TFcontrolBin) >= 1
            for ap_index = 1:NumAPbins
                TFGoodAP = ~isnan(DorsalProfiles(TFcontrolBin,ap_index,ch_index));
                if sum(TFGoodAP) > 1
                    CompiledEmbryos.WindowedProfiles.DeltaFC.Control.mean(i,ap_index,ch_index) = mean(DorsalProfiles(TFcontrolBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.WindowedProfiles.DeltaFC.Control.std(i,ap_index,ch_index) = std(DorsalProfiles(TFcontrolBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.WindowedProfiles.DeltaFC.Control.count(i,ap_index,ch_index) = sum(~isnan(DorsalProfiles(TFcontrolBin,ap_index,ch_index)));
                    CompiledEmbryos.WindowedProfiles.DeltaFC.Control.se(i,ap_index,ch_index) = CompiledEmbryos.WindowedProfiles.DeltaFC.Control.std(i,ap_index,ch_index)/...
                        sqrt(CompiledEmbryos.WindowedProfiles.DeltaFC.Control.count(i,ap_index,ch_index));
                elseif sum(TFGoodAP)  == 1
                    ControlProfs = DorsalProfiles(TFcontrolBin,ap_index,ch_index);
                    CompiledEmbryos.WindowedProfiles.DeltaFC.Control.mean(i,ap_index,ch_index) = ControlProfs(TFGoodAP);
                    CompiledEmbryos.WindowedProfiles.DeltaFC.Control.std(i,ap_index,ch_index) = NaN;
                    CompiledEmbryos.WindowedProfiles.DeltaFC.Control.count(i,ap_index,ch_index) = 1;
                    CompiledEmbryos.WindowedProfiles.DeltaFC.Control.se(i,ap_index,ch_index) = NaN;
                end
            end
        end
    end
end

%% Narrow windowed Profiles - Delta FC
NarrowDorsalProfiles = CompiledEmbryos.SlideRescaledDorsalAvgNarrowAPProfiles;
CompiledEmbryos.NarrowWindowedProfiles.DeltaFC = {};
CompiledEmbryos.NarrowWindowedProfiles.DeltaFC.x = x;
CompiledEmbryos.NarrowWindowedProfiles.DeltaFC.BinWidth = deltafc_binwidth;
CompiledEmbryos.NarrowWindowedProfiles.DeltaFC.Test = {};
CompiledEmbryos.NarrowWindowedProfiles.DeltaFC.Test.mean = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowWindowedProfiles.DeltaFC.Test.std = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowWindowedProfiles.DeltaFC.Test.se = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowWindowedProfiles.DeltaFC.Test.count = zeros(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowWindowedProfiles.DeltaFC.Control = {}; 
CompiledEmbryos.NarrowWindowedProfiles.DeltaFC.Control.mean = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowWindowedProfiles.DeltaFC.Control.std = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowWindowedProfiles.DeltaFC.Control.se = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowWindowedProfiles.DeltaFC.Control.count = zeros(NumTbins, NumNarrowAPbins, NChannels);


for i = 1:NumTbins
    TFtimeBin = [CompiledEmbryos.FixCorrectedDeltaFC_um.mean(:)].' >= x(i)-deltafc_binwidth/2  & [CompiledEmbryos.FixCorrectedDeltaFC_um.mean(:)].'< x(i)+deltafc_binwidth/2 ;
    TFtestBin = TFtimeBin & AllDeltaValidProfilesTestTF;
    TFcontrolBin = TFtimeBin & AllDeltaValidProfilesControlTF;
    for ch_index = 2:NChannels
        if sum(TFtestBin) >= 1
            for ap_index = 1:NumNarrowAPbins
                TFGoodAP = ~isnan(NarrowDorsalProfiles(TFtestBin,ap_index,ch_index));
                if sum(TFGoodAP) > 1
                    CompiledEmbryos.NarrowWindowedProfiles.DeltaFC.Test.mean(i,ap_index,ch_index) = mean(NarrowDorsalProfiles(TFtestBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.NarrowWindowedProfiles.DeltaFC.Test.std(i,ap_index,ch_index) = std(NarrowDorsalProfiles(TFtestBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.NarrowWindowedProfiles.DeltaFC.Test.count(i,ap_index,ch_index) = sum(~isnan(NarrowDorsalProfiles(TFtestBin,ap_index,ch_index)));
                    CompiledEmbryos.NarrowWindowedProfiles.DeltaFC.Test.se(i,ap_index,ch_index) = CompiledEmbryos.NarrowWindowedProfiles.DeltaFC.Test.std(i,ap_index,ch_index)/...
                        sqrt(CompiledEmbryos.NarrowWindowedProfiles.DeltaFC.Test.count(i,ap_index,ch_index));
                elseif sum(TFGoodAP)  == 1
                    NarrowTestProfs = NarrowDorsalProfiles(TFtestBin,ap_index,ch_index);
                    CompiledEmbryos.NarrowWindowedProfiles.DeltaFC.Test.mean(i,ap_index,ch_index) = NarrowTestProfs(TFGoodAP);
                    CompiledEmbryos.NarrowWindowedProfiles.DeltaFC.Test.std(i,ap_index,ch_index) = NaN;
                    CompiledEmbryos.NarrowWindowedProfiles.DeltaFC.Test.count(i,ap_index,ch_index) = 1;
                    CompiledEmbryos.NarrowWindowedProfiles.DeltaFC.Test.se(i,ap_index,ch_index) = NaN;
                end
            end
        end
        
        if sum(TFcontrolBin) >= 1
            for ap_index = 1:NumNarrowAPbins
                TFGoodAP = ~isnan(NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index));
                if sum(TFGoodAP) > 1
                    CompiledEmbryos.NarrowWindowedProfiles.DeltaFC.Control.mean(i,ap_index,ch_index) = mean(NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.NarrowWindowedProfiles.DeltaFC.Control.std(i,ap_index,ch_index) = std(NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.NarrowWindowedProfiles.DeltaFC.Control.count(i,ap_index,ch_index) = sum(~isnan(NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index)));
                    CompiledEmbryos.NarrowWindowedProfiles.DeltaFC.Control.se(i,ap_index,ch_index) = CompiledEmbryos.NarrowWindowedProfiles.DeltaFC.Control.std(i,ap_index,ch_index)/...
                        sqrt(CompiledEmbryos.NarrowWindowedProfiles.DeltaFC.Control.count(i,ap_index,ch_index));
                elseif sum(TFGoodAP)  == 1
                    NarrowControlProfs = NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index);
                    CompiledEmbryos.NarrowWindowedProfiles.DeltaFC.Control.mean(i,ap_index,ch_index) = NarrowControlProfs(TFGoodAP);
                    CompiledEmbryos.NarrowWindowedProfiles.DeltaFC.Control.std(i,ap_index,ch_index) = NaN;
                    CompiledEmbryos.NarrowWindowedProfiles.DeltaFC.Control.count(i,ap_index,ch_index) = 1;
                    CompiledEmbryos.NarrowWindowedProfiles.DeltaFC.Control.se(i,ap_index,ch_index) = NaN;
                end
            end
        end
    end
end

%% Add Dubuis Time-based Binning - standard AP
AllDubuisTimeValidProfilesTestTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.TestSetEmbryos & ~isnan(CompiledEmbryos.DubuisEmbryoTimes);% &...
AllDubuisTimeProfilesControlTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.ControlSetEmbryos & ~isnan(CompiledEmbryos.DubuisEmbryoTimes);% &...
        %CompiledEmbryos.FixCorrectedDeltaFC_um(:,1).' > 3;


x = 0:dubuistime_binwidth:65;

NumTbins = length(x);

DorsalProfiles = CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles;
CompiledEmbryos.BinnedProfiles.DubuisTime = {};
CompiledEmbryos.BinnedProfiles.DubuisTime.x = x;
CompiledEmbryos.BinnedProfiles.DubuisTime.BinWidth = dubuistime_binwidth;
CompiledEmbryos.BinnedProfiles.DubuisTime.Test = {};
CompiledEmbryos.BinnedProfiles.DubuisTime.Test.mean = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.BinnedProfiles.DubuisTime.Test.std = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.BinnedProfiles.DubuisTime.Test.se = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.BinnedProfiles.DubuisTime.Test.count = zeros(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.BinnedProfiles.DubuisTime.Test.count_abovecenter = zeros(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.BinnedProfiles.DubuisTime.Test.count_belowcenter = zeros(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.BinnedProfiles.DubuisTime.Control = {}; 
CompiledEmbryos.BinnedProfiles.DubuisTime.Control.mean = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.BinnedProfiles.DubuisTime.Control.std = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.BinnedProfiles.DubuisTime.Control.se = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.BinnedProfiles.DubuisTime.Control.count = zeros(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.BinnedProfiles.DubuisTime.Control.count_abovecenter = zeros(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.BinnedProfiles.DubuisTime.Control.count_belowcenter = zeros(NumTbins, NumAPbins, NChannels);


for i = 1:NumTbins
    TFtimeBin = CompiledEmbryos.DubuisEmbryoTimes >= x(i)-dubuistime_binwidth/2 & CompiledEmbryos.DubuisEmbryoTimes < x(i)+dubuistime_binwidth/2;
    TFtestBin = TFtimeBin & AllDubuisTimeValidProfilesTestTF;
    TFcontrolBin = TFtimeBin & AllDubuisTimeProfilesControlTF;
    for ch_index = 2:NChannels
        if sum(TFtestBin) >= 1
            
            for ap_index = 1:NumAPbins
                TFGoodAP = ~isnan(DorsalProfiles(TFtestBin,ap_index,ch_index));
                TestCountsAboveCenter = sum(~isnan(DorsalProfiles(CompiledEmbryos.DubuisEmbryoTimes > x(i) &  TFtestBin,ap_index,ch_index)));
                TestCountsBelowCenter = sum(~isnan(DorsalProfiles(CompiledEmbryos.DubuisEmbryoTimes < x(i) &  TFtestBin,ap_index,ch_index)));
                if sum(TFGoodAP) > 1
                    CompiledEmbryos.BinnedProfiles.DubuisTime.Test.mean(i,ap_index,ch_index) = mean(DorsalProfiles(TFtestBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.BinnedProfiles.DubuisTime.Test.std(i,ap_index,ch_index) = std(DorsalProfiles(TFtestBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.BinnedProfiles.DubuisTime.Test.count(i,ap_index,ch_index) = sum(~isnan(DorsalProfiles(TFtestBin,ap_index,ch_index)));
                    CompiledEmbryos.BinnedProfiles.DubuisTime.Test.se(i,ap_index,ch_index) = CompiledEmbryos.BinnedProfiles.DubuisTime.Test.std(i,ap_index,ch_index)/...
                        sqrt(CompiledEmbryos.BinnedProfiles.DubuisTime.Test.count(i,ap_index,ch_index));
                elseif sum(TFGoodAP)  == 1
                    TestProfs = DorsalProfiles(TFtestBin,ap_index,ch_index);
                    CompiledEmbryos.BinnedProfiles.DubuisTime.Test.mean(i,ap_index,ch_index) = TestProfs(TFGoodAP);
                    CompiledEmbryos.BinnedProfiles.DubuisTime.Test.std(i,ap_index,ch_index) = NaN;
                    CompiledEmbryos.BinnedProfiles.DubuisTime.Test.count(i,ap_index,ch_index) = 1;
                    CompiledEmbryos.BinnedProfiles.DubuisTime.Test.se(i,ap_index,ch_index) = NaN;
                end
                
                CompiledEmbryos.BinnedProfiles.DubuisTime.Test.count_abovecenter(i,ap_index,ch_index) = TestCountsAboveCenter;
                CompiledEmbryos.BinnedProfiles.DubuisTime.Test.count_belowcenter(i,ap_index,ch_index) = TestCountsBelowCenter;
            end
        end
        
        if sum(TFcontrolBin) >= 1
            ControlCountsAboveCenter = sum(~isnan(DorsalProfiles(CompiledEmbryos.DubuisEmbryoTimes > x(i) &  TFcontrolBin,ap_index,ch_index)));
            ControlCountsBelowCenter = sum(~isnan(DorsalProfiles(CompiledEmbryos.DubuisEmbryoTimes < x(i) &  TFcontrolBin,ap_index,ch_index)));
            for ap_index = 1:NumAPbins
                TFGoodAP = ~isnan(DorsalProfiles(TFcontrolBin,ap_index,ch_index));
                if sum(TFGoodAP) > 1
                    CompiledEmbryos.BinnedProfiles.DubuisTime.Control.mean(i,ap_index,ch_index) = mean(DorsalProfiles(TFcontrolBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.BinnedProfiles.DubuisTime.Control.std(i,ap_index,ch_index) = std(DorsalProfiles(TFcontrolBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.BinnedProfiles.DubuisTime.Control.count(i,ap_index,ch_index) = sum(~isnan(DorsalProfiles(TFcontrolBin,ap_index,ch_index)));
                    CompiledEmbryos.BinnedProfiles.DubuisTime.Control.se(i,ap_index,ch_index) = CompiledEmbryos.BinnedProfiles.DubuisTime.Control.std(i,ap_index,ch_index)/...
                        sqrt(CompiledEmbryos.BinnedProfiles.DubuisTime.Control.count(i,ap_index,ch_index));
                elseif sum(TFGoodAP)  == 1
                    ControlProfs = DorsalProfiles(TFcontrolBin,ap_index,ch_index);
                    CompiledEmbryos.BinnedProfiles.DubuisTime.Control.mean(i,ap_index,ch_index) = ControlProfs(TFGoodAP);
                    CompiledEmbryos.BinnedProfiles.DubuisTime.Control.std(i,ap_index,ch_index) = NaN;
                    CompiledEmbryos.BinnedProfiles.DubuisTime.Control.count(i,ap_index,ch_index) = 1;
                    CompiledEmbryos.BinnedProfiles.DubuisTime.Control.se(i,ap_index,ch_index) = NaN;
                end
                
                CompiledEmbryos.BinnedProfiles.DubuisTime.Control.count_abovecenter(i,ap_index,ch_index) = ControlCountsAboveCenter;
                CompiledEmbryos.BinnedProfiles.DubuisTime.Control.count_belowcenter(i,ap_index,ch_index) = ControlCountsBelowCenter;
            end
        end
    end
end

CompiledEmbryos.BinnedProfiles.DubuisTime.Test.count_balance = (CompiledEmbryos.BinnedProfiles.DubuisTime.Test.count_abovecenter)./(CompiledEmbryos.BinnedProfiles.DubuisTime.Test.count_abovecenter+CompiledEmbryos.BinnedProfiles.DubuisTime.Test.count_belowcenter);
CompiledEmbryos.BinnedProfiles.DubuisTime.Control.count_balance = (CompiledEmbryos.BinnedProfiles.DubuisTime.Control.count_abovecenter)./(CompiledEmbryos.BinnedProfiles.DubuisTime.Control.count_abovecenter+CompiledEmbryos.BinnedProfiles.DubuisTime.Control.count_belowcenter);

%% Add Dubuis Time-based Binning - narrow AP
NarrowDorsalProfiles = CompiledEmbryos.SlideRescaledDorsalAvgNarrowAPProfiles;
CompiledEmbryos.NarrowBinnedProfiles.DubuisTime = {};
CompiledEmbryos.NarrowBinnedProfiles.DubuisTime.BinWidth = deltafc_binwidth;
CompiledEmbryos.NarrowBinnedProfiles.DubuisTime.x = x;
CompiledEmbryos.NarrowBinnedProfiles.DubuisTime.Test = {};
CompiledEmbryos.NarrowBinnedProfiles.DubuisTime.Test.mean = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowBinnedProfiles.DubuisTime.Test.std = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowBinnedProfiles.DubuisTime.Test.se = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowBinnedProfiles.DubuisTime.Test.count = zeros(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowBinnedProfiles.DubuisTime.Control = {}; 
CompiledEmbryos.NarrowBinnedProfiles.DubuisTime.Control.mean = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowBinnedProfiles.DubuisTime.Control.std = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowBinnedProfiles.DubuisTime.Control.se = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowBinnedProfiles.DubuisTime.Control.count = zeros(NumTbins, NumNarrowAPbins, NChannels);


for i = 1:NumTbins
    TFtimeBin = CompiledEmbryos.DubuisEmbryoTimes >= x(i)-dubuistime_binwidth/2 & CompiledEmbryos.DubuisEmbryoTimes < x(i)+dubuistime_binwidth/2;
    TFtestBin = TFtimeBin & AllDubuisTimeValidProfilesTestTF;
    TFcontrolBin = TFtimeBin & AllDubuisTimeProfilesControlTF;
    for ch_index = 2:NChannels
        if sum(TFtestBin) >= 1
            for ap_index = 1:NumNarrowAPbins
                TFGoodAP = ~isnan(NarrowDorsalProfiles(TFtestBin,ap_index,ch_index));
                if sum(TFGoodAP) > 1
                    CompiledEmbryos.NarrowBinnedProfiles.DubuisTime.Test.mean(i,ap_index,ch_index) = mean(NarrowDorsalProfiles(TFtestBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.NarrowBinnedProfiles.DubuisTime.Test.std(i,ap_index,ch_index) = std(NarrowDorsalProfiles(TFtestBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.NarrowBinnedProfiles.DubuisTime.Test.count(i,ap_index,ch_index) = sum(~isnan(NarrowDorsalProfiles(TFtestBin,ap_index,ch_index)));
                    CompiledEmbryos.NarrowBinnedProfiles.DubuisTime.Test.se(i,ap_index,ch_index) = CompiledEmbryos.NarrowBinnedProfiles.DubuisTime.Test.std(i,ap_index,ch_index)/...
                        sqrt(CompiledEmbryos.NarrowBinnedProfiles.DubuisTime.Test.count(i,ap_index,ch_index));
                elseif sum(TFGoodAP)  == 1
                    NarrowTestProfs = NarrowDorsalProfiles(TFtestBin,ap_index,ch_index);
                    CompiledEmbryos.NarrowBinnedProfiles.DubuisTime.Test.mean(i,ap_index,ch_index) = NarrowTestProfs(TFGoodAP);
                    CompiledEmbryos.NarrowBinnedProfiles.DubuisTime.Test.std(i,ap_index,ch_index) = NaN;
                    CompiledEmbryos.NarrowBinnedProfiles.DubuisTime.Test.count(i,ap_index,ch_index) = 1;
                    CompiledEmbryos.NarrowBinnedProfiles.DubuisTime.Test.se(i,ap_index,ch_index) = NaN;
                end
            end
        end
        
        if sum(TFcontrolBin) >= 1
            for ap_index = 1:NumNarrowAPbins
                TFGoodAP = ~isnan(NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index));
                if sum(TFGoodAP) > 1
                    CompiledEmbryos.NarrowBinnedProfiles.DubuisTime.Control.mean(i,ap_index,ch_index) = mean(NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.NarrowBinnedProfiles.DubuisTime.Control.std(i,ap_index,ch_index) = std(NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.NarrowBinnedProfiles.DubuisTime.Control.count(i,ap_index,ch_index) = sum(~isnan(NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index)));
                    CompiledEmbryos.NarrowBinnedProfiles.DubuisTime.Control.se(i,ap_index,ch_index) = CompiledEmbryos.NarrowBinnedProfiles.DubuisTime.Control.std(i,ap_index,ch_index)/...
                        sqrt(CompiledEmbryos.NarrowBinnedProfiles.DubuisTime.Control.count(i,ap_index,ch_index));
                elseif sum(TFGoodAP)  == 1
                    NarrowControlProfs = NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index);
                    CompiledEmbryos.NarrowBinnedProfiles.DubuisTime.Control.mean(i,ap_index,ch_index) = NarrowControlProfs(TFGoodAP);
                    CompiledEmbryos.NarrowBinnedProfiles.DubuisTime.Control.std(i,ap_index,ch_index) = NaN;
                    CompiledEmbryos.NarrowBinnedProfiles.DubuisTime.Control.count(i,ap_index,ch_index) = 1;
                    CompiledEmbryos.NarrowBinnedProfiles.DubuisTime.Control.se(i,ap_index,ch_index) = NaN;
                end
            end
        end
    end
end

%% Windowed Profiles - Dubuis Times

x = 0:1:65;

NumTbins = length(x);


DorsalProfiles = CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles;
CompiledEmbryos.WindowedProfiles.DubuisTime = {};
CompiledEmbryos.WindowedProfiles.DubuisTime.x = x;
CompiledEmbryos.WindowedProfiles.DubuisTime.BinWidth = dubuistime_binwidth;
CompiledEmbryos.WindowedProfiles.DubuisTime.Test = {};
CompiledEmbryos.WindowedProfiles.DubuisTime.Test.mean = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.WindowedProfiles.DubuisTime.Test.std = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.WindowedProfiles.DubuisTime.Test.se = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.WindowedProfiles.DubuisTime.Test.count = zeros(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.WindowedProfiles.DubuisTime.Test.count_abovecenter = zeros(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.WindowedProfiles.DubuisTime.Test.count_belowcenter = zeros(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.WindowedProfiles.DubuisTime.Control = {}; 
CompiledEmbryos.WindowedProfiles.DubuisTime.Control.mean = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.WindowedProfiles.DubuisTime.Control.std = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.WindowedProfiles.DubuisTime.Control.se = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.WindowedProfiles.DubuisTime.Control.count = zeros(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.WindowedProfiles.DubuisTime.Control.count_abovecenter = zeros(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.WindowedProfiles.DubuisTime.Control.count_belowcenter = zeros(NumTbins, NumAPbins, NChannels);


for i = 1:NumTbins
    TFtimeBin = CompiledEmbryos.DubuisEmbryoTimes >= x(i)-dubuistime_binwidth/2 & CompiledEmbryos.DubuisEmbryoTimes < x(i)+dubuistime_binwidth/2;
    TFtestBin = TFtimeBin & AllDubuisTimeValidProfilesTestTF;
    TFcontrolBin = TFtimeBin & AllDubuisTimeProfilesControlTF;
    for ch_index = 2:NChannels
        if sum(TFtestBin) >= 1
            
            for ap_index = 1:NumAPbins
                TFGoodAP = ~isnan(DorsalProfiles(TFtestBin,ap_index,ch_index));
                TestCountsAboveCenter = sum(~isnan(DorsalProfiles(CompiledEmbryos.DubuisEmbryoTimes > x(i) &  TFtestBin,ap_index,ch_index)));
                TestCountsBelowCenter = sum(~isnan(DorsalProfiles(CompiledEmbryos.DubuisEmbryoTimes < x(i) &  TFtestBin,ap_index,ch_index)));
                if sum(TFGoodAP) > 1
                    CompiledEmbryos.WindowedProfiles.DubuisTime.Test.mean(i,ap_index,ch_index) = mean(DorsalProfiles(TFtestBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.WindowedProfiles.DubuisTime.Test.std(i,ap_index,ch_index) = std(DorsalProfiles(TFtestBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.WindowedProfiles.DubuisTime.Test.count(i,ap_index,ch_index) = sum(~isnan(DorsalProfiles(TFtestBin,ap_index,ch_index)));
                    CompiledEmbryos.WindowedProfiles.DubuisTime.Test.se(i,ap_index,ch_index) = CompiledEmbryos.WindowedProfiles.DubuisTime.Test.std(i,ap_index,ch_index)/...
                        sqrt(CompiledEmbryos.WindowedProfiles.DubuisTime.Test.count(i,ap_index,ch_index));
                elseif sum(TFGoodAP)  == 1
                    TestProfs = DorsalProfiles(TFtestBin,ap_index,ch_index);
                    CompiledEmbryos.WindowedProfiles.DubuisTime.Test.mean(i,ap_index,ch_index) = TestProfs(TFGoodAP);
                    CompiledEmbryos.WindowedProfiles.DubuisTime.Test.std(i,ap_index,ch_index) = NaN;
                    CompiledEmbryos.WindowedProfiles.DubuisTime.Test.count(i,ap_index,ch_index) = 1;
                    CompiledEmbryos.WindowedProfiles.DubuisTime.Test.se(i,ap_index,ch_index) = NaN;
                end
                
                CompiledEmbryos.WindowedProfiles.DubuisTime.Test.count_abovecenter(i,ap_index,ch_index) = TestCountsAboveCenter;
                CompiledEmbryos.WindowedProfiles.DubuisTime.Test.count_belowcenter(i,ap_index,ch_index) = TestCountsBelowCenter;
            end
        end
        
        if sum(TFcontrolBin) >= 1
            
            for ap_index = 1:NumAPbins
                TFGoodAP = ~isnan(DorsalProfiles(TFcontrolBin,ap_index,ch_index));
                ControlCountsAboveCenter = sum(~isnan(DorsalProfiles(CompiledEmbryos.DubuisEmbryoTimes > x(i) &  TFcontrolBin,ap_index,ch_index)));
                ControlCountsBelowCenter = sum(~isnan(DorsalProfiles(CompiledEmbryos.DubuisEmbryoTimes < x(i) &  TFcontrolBin,ap_index,ch_index)));
                if sum(TFGoodAP) > 1
                    CompiledEmbryos.WindowedProfiles.DubuisTime.Control.mean(i,ap_index,ch_index) = mean(DorsalProfiles(TFcontrolBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.WindowedProfiles.DubuisTime.Control.std(i,ap_index,ch_index) = std(DorsalProfiles(TFcontrolBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.WindowedProfiles.DubuisTime.Control.count(i,ap_index,ch_index) = sum(~isnan(DorsalProfiles(TFcontrolBin,ap_index,ch_index)));
                    CompiledEmbryos.WindowedProfiles.DubuisTime.Control.se(i,ap_index,ch_index) = CompiledEmbryos.WindowedProfiles.DubuisTime.Control.std(i,ap_index,ch_index)/...
                        sqrt(CompiledEmbryos.WindowedProfiles.DubuisTime.Control.count(i,ap_index,ch_index));
                elseif sum(TFGoodAP)  == 1
                    ControlProfs = DorsalProfiles(TFcontrolBin,ap_index,ch_index);
                    CompiledEmbryos.WindowedProfiles.DubuisTime.Control.mean(i,ap_index,ch_index) = ControlProfs(TFGoodAP);
                    CompiledEmbryos.WindowedProfiles.DubuisTime.Control.std(i,ap_index,ch_index) = NaN;
                    CompiledEmbryos.WindowedProfiles.DubuisTime.Control.count(i,ap_index,ch_index) = 1;
                    CompiledEmbryos.WindowedProfiles.DubuisTime.Control.se(i,ap_index,ch_index) = NaN;
                end
                
                CompiledEmbryos.WindowedProfiles.DubuisTime.Control.count_abovecenter(i,ap_index,ch_index) = ControlCountsAboveCenter;
                CompiledEmbryos.WindowedProfiles.DubuisTime.Control.count_belowcenter(i,ap_index,ch_index) = ControlCountsBelowCenter;
            end
        end
    end
end

CompiledEmbryos.WindowedProfiles.DubuisTime.Test.count_balance = (CompiledEmbryos.WindowedProfiles.DubuisTime.Test.count_abovecenter)./(CompiledEmbryos.WindowedProfiles.DubuisTime.Test.count_abovecenter+CompiledEmbryos.WindowedProfiles.DubuisTime.Test.count_belowcenter);
CompiledEmbryos.WindowedProfiles.DubuisTime.Control.count_balance = (CompiledEmbryos.WindowedProfiles.DubuisTime.Control.count_abovecenter)./(CompiledEmbryos.WindowedProfiles.DubuisTime.Control.count_abovecenter+CompiledEmbryos.WindowedProfiles.DubuisTime.Control.count_belowcenter);

%% Narrow windowed Profiles - Dubuis Times
NarrowDorsalProfiles = CompiledEmbryos.SlideRescaledDorsalAvgNarrowAPProfiles;
CompiledEmbryos.NarrowWindowedProfiles.DubuisTime = {};
CompiledEmbryos.NarrowWindowedProfiles.DubuisTime.BinWidth = deltafc_binwidth;
CompiledEmbryos.NarrowWindowedProfiles.DubuisTime.x = x;
CompiledEmbryos.NarrowWindowedProfiles.DubuisTime.Test = {};
CompiledEmbryos.NarrowWindowedProfiles.DubuisTime.Test.mean = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowWindowedProfiles.DubuisTime.Test.std = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowWindowedProfiles.DubuisTime.Test.se = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowWindowedProfiles.DubuisTime.Test.count = zeros(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowWindowedProfiles.DubuisTime.Control = {}; 
CompiledEmbryos.NarrowWindowedProfiles.DubuisTime.Control.mean = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowWindowedProfiles.DubuisTime.Control.std = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowWindowedProfiles.DubuisTime.Control.se = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowWindowedProfiles.DubuisTime.Control.count = zeros(NumTbins, NumNarrowAPbins, NChannels);


for i = 1:NumTbins
    TFtimeBin = CompiledEmbryos.DubuisEmbryoTimes >= x(i)-dubuistime_binwidth/2 & CompiledEmbryos.DubuisEmbryoTimes < x(i)+dubuistime_binwidth/2;
    TFtestBin = TFtimeBin & AllDubuisTimeValidProfilesTestTF;
    TFcontrolBin = TFtimeBin & AllDubuisTimeProfilesControlTF;
    for ch_index = 2:NChannels
        if sum(TFtestBin) >= 1
            for ap_index = 1:NumNarrowAPbins
                TFGoodAP = ~isnan(NarrowDorsalProfiles(TFtestBin,ap_index,ch_index));
                if sum(TFGoodAP) > 1
                    CompiledEmbryos.NarrowWindowedProfiles.DubuisTime.Test.mean(i,ap_index,ch_index) = mean(NarrowDorsalProfiles(TFtestBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.NarrowWindowedProfiles.DubuisTime.Test.std(i,ap_index,ch_index) = std(NarrowDorsalProfiles(TFtestBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.NarrowWindowedProfiles.DubuisTime.Test.count(i,ap_index,ch_index) = sum(~isnan(NarrowDorsalProfiles(TFtestBin,ap_index,ch_index)));
                    CompiledEmbryos.NarrowWindowedProfiles.DubuisTime.Test.se(i,ap_index,ch_index) = CompiledEmbryos.NarrowWindowedProfiles.DubuisTime.Test.std(i,ap_index,ch_index)/...
                        sqrt(CompiledEmbryos.NarrowWindowedProfiles.DubuisTime.Test.count(i,ap_index,ch_index));
                elseif sum(TFGoodAP)  == 1
                    NarrowTestProfs = NarrowDorsalProfiles(TFtestBin,ap_index,ch_index);
                    CompiledEmbryos.NarrowWindowedProfiles.DubuisTime.Test.mean(i,ap_index,ch_index) = NarrowTestProfs(TFGoodAP);
                    CompiledEmbryos.NarrowWindowedProfiles.DubuisTime.Test.std(i,ap_index,ch_index) = NaN;
                    CompiledEmbryos.NarrowWindowedProfiles.DubuisTime.Test.count(i,ap_index,ch_index) = 1;
                    CompiledEmbryos.NarrowWindowedProfiles.DubuisTime.Test.se(i,ap_index,ch_index) = NaN;
                end
            end
        end
        
        if sum(TFcontrolBin) >= 1
            for ap_index = 1:NumNarrowAPbins
                TFGoodAP = ~isnan(NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index));
                if sum(TFGoodAP) > 1
                    CompiledEmbryos.NarrowWindowedProfiles.DubuisTime.Control.mean(i,ap_index,ch_index) = mean(NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.NarrowWindowedProfiles.DubuisTime.Control.std(i,ap_index,ch_index) = std(NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.NarrowWindowedProfiles.DubuisTime.Control.count(i,ap_index,ch_index) = sum(~isnan(NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index)));
                    CompiledEmbryos.NarrowWindowedProfiles.DubuisTime.Control.se(i,ap_index,ch_index) = CompiledEmbryos.NarrowWindowedProfiles.DubuisTime.Control.std(i,ap_index,ch_index)/...
                        sqrt(CompiledEmbryos.NarrowWindowedProfiles.DubuisTime.Control.count(i,ap_index,ch_index));
                elseif sum(TFGoodAP)  == 1
                    NarrowControlProfs = NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index);
                    CompiledEmbryos.NarrowWindowedProfiles.DubuisTime.Control.mean(i,ap_index,ch_index) = NarrowControlProfs(TFGoodAP);
                    CompiledEmbryos.NarrowWindowedProfiles.DubuisTime.Control.std(i,ap_index,ch_index) = NaN;
                    CompiledEmbryos.NarrowWindowedProfiles.DubuisTime.Control.count(i,ap_index,ch_index) = 1;
                    CompiledEmbryos.NarrowWindowedProfiles.DubuisTime.Control.se(i,ap_index,ch_index) = NaN;
                end
            end
        end
    end
end


%% Add yw25C Time-based Binning - standard AP
Allyw25CTimeValidProfilesTestTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.TestSetEmbryos & ~isnan(CompiledEmbryos.yw25CEmbryoTimes);% &...
Allyw25CTimeProfilesControlTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.ControlSetEmbryos & ~isnan(CompiledEmbryos.yw25CEmbryoTimes);% &...
        %CompiledEmbryos.FixCorrectedDeltaFC_um(:,1).' > 3;


x = 0:yw25Ctime_binwidth:65;
NumTbins = length(x);

DorsalProfiles = CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles;
CompiledEmbryos.BinnedProfiles.yw25CTime = {};
CompiledEmbryos.BinnedProfiles.yw25CTime.BinWidth = yw25Ctime_binwidth;
CompiledEmbryos.BinnedProfiles.yw25CTime.x = x;
CompiledEmbryos.BinnedProfiles.yw25CTime.Test = {};
CompiledEmbryos.BinnedProfiles.yw25CTime.Test.mean = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.BinnedProfiles.yw25CTime.Test.std = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.BinnedProfiles.yw25CTime.Test.se = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.BinnedProfiles.yw25CTime.Test.count = zeros(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.BinnedProfiles.yw25CTime.Control = {}; 
CompiledEmbryos.BinnedProfiles.yw25CTime.Control.mean = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.BinnedProfiles.yw25CTime.Control.std = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.BinnedProfiles.yw25CTime.Control.se = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.BinnedProfiles.yw25CTime.Control.count = zeros(NumTbins, NumAPbins, NChannels);


for i = 1:NumTbins
    TFtimeBin = CompiledEmbryos.yw25CEmbryoTimes >= x(i)-yw25Ctime_binwidth/2 & CompiledEmbryos.yw25CEmbryoTimes < x(i)+yw25Ctime_binwidth/2;
    TFtestBin = TFtimeBin & Allyw25CTimeValidProfilesTestTF;
    TFcontrolBin = TFtimeBin & Allyw25CTimeProfilesControlTF;
    for ch_index = 2:NChannels
        if sum(TFtestBin) >= 1
            for ap_index = 1:NumAPbins
                TFGoodAP = ~isnan(DorsalProfiles(TFtestBin,ap_index,ch_index));
                if sum(TFGoodAP) > 1
                    CompiledEmbryos.BinnedProfiles.yw25CTime.Test.mean(i,ap_index,ch_index) = mean(DorsalProfiles(TFtestBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.BinnedProfiles.yw25CTime.Test.std(i,ap_index,ch_index) = std(DorsalProfiles(TFtestBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.BinnedProfiles.yw25CTime.Test.count(i,ap_index,ch_index) = sum(~isnan(DorsalProfiles(TFtestBin,ap_index,ch_index)));
                    CompiledEmbryos.BinnedProfiles.yw25CTime.Test.se(i,ap_index,ch_index) = CompiledEmbryos.BinnedProfiles.yw25CTime.Test.std(i,ap_index,ch_index)/...
                        sqrt(CompiledEmbryos.BinnedProfiles.yw25CTime.Test.count(i,ap_index,ch_index));
                elseif sum(TFGoodAP)  == 1
                    TestProfs = DorsalProfiles(TFtestBin,ap_index,ch_index);
                    CompiledEmbryos.BinnedProfiles.yw25CTime.Test.mean(i,ap_index,ch_index) = TestProfs(TFGoodAP);
                    CompiledEmbryos.BinnedProfiles.yw25CTime.Test.std(i,ap_index,ch_index) = NaN;
                    CompiledEmbryos.BinnedProfiles.yw25CTime.Test.count(i,ap_index,ch_index) = 1;
                    CompiledEmbryos.BinnedProfiles.yw25CTime.Test.se(i,ap_index,ch_index) = NaN;
                end
            end
        end
        
        if sum(TFcontrolBin) >= 1
            for ap_index = 1:NumAPbins
                TFGoodAP = ~isnan(DorsalProfiles(TFcontrolBin,ap_index,ch_index));
                if sum(TFGoodAP) > 1
                    CompiledEmbryos.BinnedProfiles.yw25CTime.Control.mean(i,ap_index,ch_index) = mean(DorsalProfiles(TFcontrolBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.BinnedProfiles.yw25CTime.Control.std(i,ap_index,ch_index) = std(DorsalProfiles(TFcontrolBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.BinnedProfiles.yw25CTime.Control.count(i,ap_index,ch_index) = sum(~isnan(DorsalProfiles(TFcontrolBin,ap_index,ch_index)));
                    CompiledEmbryos.BinnedProfiles.yw25CTime.Control.se(i,ap_index,ch_index) = CompiledEmbryos.BinnedProfiles.yw25CTime.Control.std(i,ap_index,ch_index)/...
                        sqrt(CompiledEmbryos.BinnedProfiles.yw25CTime.Control.count(i,ap_index,ch_index));
                elseif sum(TFGoodAP)  == 1
                    ControlProfs = DorsalProfiles(TFcontrolBin,ap_index,ch_index);
                    CompiledEmbryos.BinnedProfiles.yw25CTime.Control.mean(i,ap_index,ch_index) = ControlProfs(TFGoodAP);
                    CompiledEmbryos.BinnedProfiles.yw25CTime.Control.std(i,ap_index,ch_index) = NaN;
                    CompiledEmbryos.BinnedProfiles.yw25CTime.Control.count(i,ap_index,ch_index) = 1;
                    CompiledEmbryos.BinnedProfiles.yw25CTime.Control.se(i,ap_index,ch_index) = NaN;
                end
            end
        end
    end
end

%% Add yw25C Time-based Binning - narrow AP
NarrowDorsalProfiles = CompiledEmbryos.SlideRescaledDorsalAvgNarrowAPProfiles;
CompiledEmbryos.NarrowBinnedProfiles.yw25CTime = {};
CompiledEmbryos.NarrowBinnedProfiles.yw25CTime.BinWidth = deltafc_binwidth;
CompiledEmbryos.NarrowBinnedProfiles.yw25CTime.x = x;
CompiledEmbryos.NarrowBinnedProfiles.yw25CTime.Test = {};
CompiledEmbryos.NarrowBinnedProfiles.yw25CTime.Test.mean = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowBinnedProfiles.yw25CTime.Test.std = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowBinnedProfiles.yw25CTime.Test.se = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowBinnedProfiles.yw25CTime.Test.count = zeros(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowBinnedProfiles.yw25CTime.Control = {}; 
CompiledEmbryos.NarrowBinnedProfiles.yw25CTime.Control.mean = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowBinnedProfiles.yw25CTime.Control.std = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowBinnedProfiles.yw25CTime.Control.se = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowBinnedProfiles.yw25CTime.Control.count = zeros(NumTbins, NumNarrowAPbins, NChannels);


for i = 1:NumTbins
    TFtimeBin = CompiledEmbryos.yw25CEmbryoTimes >= x(i)-yw25Ctime_binwidth/2 & CompiledEmbryos.yw25CEmbryoTimes < x(i)+yw25Ctime_binwidth/2;
    TFtestBin = TFtimeBin & Allyw25CTimeValidProfilesTestTF;
    TFcontrolBin = TFtimeBin & Allyw25CTimeProfilesControlTF;
    for ch_index = 2:NChannels
        if sum(TFtestBin) >= 1
            for ap_index = 1:NumNarrowAPbins
                TFGoodAP = ~isnan(NarrowDorsalProfiles(TFtestBin,ap_index,ch_index));
                if sum(TFGoodAP) > 1
                    CompiledEmbryos.NarrowBinnedProfiles.yw25CTime.Test.mean(i,ap_index,ch_index) = mean(NarrowDorsalProfiles(TFtestBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.NarrowBinnedProfiles.yw25CTime.Test.std(i,ap_index,ch_index) = std(NarrowDorsalProfiles(TFtestBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.NarrowBinnedProfiles.yw25CTime.Test.count(i,ap_index,ch_index) = sum(~isnan(NarrowDorsalProfiles(TFtestBin,ap_index,ch_index)));
                    CompiledEmbryos.NarrowBinnedProfiles.yw25CTime.Test.se(i,ap_index,ch_index) = CompiledEmbryos.NarrowBinnedProfiles.yw25CTime.Test.std(i,ap_index,ch_index)/...
                        sqrt(CompiledEmbryos.NarrowBinnedProfiles.yw25CTime.Test.count(i,ap_index,ch_index));
                elseif sum(TFGoodAP)  == 1
                    NarrowTestProfs = NarrowDorsalProfiles(TFtestBin,ap_index,ch_index);
                    CompiledEmbryos.NarrowBinnedProfiles.yw25CTime.Test.mean(i,ap_index,ch_index) = NarrowTestProfs(TFGoodAP);
                    CompiledEmbryos.NarrowBinnedProfiles.yw25CTime.Test.std(i,ap_index,ch_index) = NaN;
                    CompiledEmbryos.NarrowBinnedProfiles.yw25CTime.Test.count(i,ap_index,ch_index) = 1;
                    CompiledEmbryos.NarrowBinnedProfiles.yw25CTime.Test.se(i,ap_index,ch_index) = NaN;
                end
            end
        end
        
        if sum(TFcontrolBin) >= 1
            for ap_index = 1:NumNarrowAPbins
                TFGoodAP = ~isnan(NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index));
                if sum(TFGoodAP) > 1
                    CompiledEmbryos.NarrowBinnedProfiles.yw25CTime.Control.mean(i,ap_index,ch_index) = mean(NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.NarrowBinnedProfiles.yw25CTime.Control.std(i,ap_index,ch_index) = std(NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.NarrowBinnedProfiles.yw25CTime.Control.count(i,ap_index,ch_index) = sum(~isnan(NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index)));
                    CompiledEmbryos.NarrowBinnedProfiles.yw25CTime.Control.se(i,ap_index,ch_index) = CompiledEmbryos.NarrowBinnedProfiles.yw25CTime.Control.std(i,ap_index,ch_index)/...
                        sqrt(CompiledEmbryos.NarrowBinnedProfiles.yw25CTime.Control.count(i,ap_index,ch_index));
                elseif sum(TFGoodAP)  == 1
                    NarrowControlProfs = NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index);
                    CompiledEmbryos.NarrowBinnedProfiles.yw25CTime.Control.mean(i,ap_index,ch_index) = NarrowControlProfs(TFGoodAP);
                    CompiledEmbryos.NarrowBinnedProfiles.yw25CTime.Control.std(i,ap_index,ch_index) = NaN;
                    CompiledEmbryos.NarrowBinnedProfiles.yw25CTime.Control.count(i,ap_index,ch_index) = 1;
                    CompiledEmbryos.NarrowBinnedProfiles.yw25CTime.Control.se(i,ap_index,ch_index) = NaN;
                end
            end
        end
    end
end

%% Windowed Profiles - yw25C Times

x = 0:1:65;
NumTbins = length(x);

DorsalProfiles = CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles;
CompiledEmbryos.WindowedProfiles.yw25CTime = {};
CompiledEmbryos.WindowedProfiles.yw25CTime.BinWidth = yw25Ctime_binwidth;
CompiledEmbryos.WindowedProfiles.yw25CTime.x = x;
CompiledEmbryos.WindowedProfiles.yw25CTime.Test = {};
CompiledEmbryos.WindowedProfiles.yw25CTime.Test.mean = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.WindowedProfiles.yw25CTime.Test.std = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.WindowedProfiles.yw25CTime.Test.se = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.WindowedProfiles.yw25CTime.Test.count = zeros(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.WindowedProfiles.yw25CTime.Control = {}; 
CompiledEmbryos.WindowedProfiles.yw25CTime.Control.mean = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.WindowedProfiles.yw25CTime.Control.std = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.WindowedProfiles.yw25CTime.Control.se = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.WindowedProfiles.yw25CTime.Control.count = zeros(NumTbins, NumAPbins, NChannels);


for i = 1:NumTbins
    TFtimeBin = CompiledEmbryos.yw25CEmbryoTimes >= x(i)-yw25Ctime_binwidth/2 & CompiledEmbryos.yw25CEmbryoTimes < x(i)+yw25Ctime_binwidth/2;
    TFtestBin = TFtimeBin & Allyw25CTimeValidProfilesTestTF;
    TFcontrolBin = TFtimeBin & Allyw25CTimeProfilesControlTF;
    for ch_index = 2:NChannels
        if sum(TFtestBin) >= 1
            for ap_index = 1:NumAPbins
                TFGoodAP = ~isnan(DorsalProfiles(TFtestBin,ap_index,ch_index));
                if sum(TFGoodAP) > 1
                    CompiledEmbryos.WindowedProfiles.yw25CTime.Test.mean(i,ap_index,ch_index) = mean(DorsalProfiles(TFtestBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.WindowedProfiles.yw25CTime.Test.std(i,ap_index,ch_index) = std(DorsalProfiles(TFtestBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.WindowedProfiles.yw25CTime.Test.count(i,ap_index,ch_index) = sum(~isnan(DorsalProfiles(TFtestBin,ap_index,ch_index)));
                    CompiledEmbryos.WindowedProfiles.yw25CTime.Test.se(i,ap_index,ch_index) = CompiledEmbryos.WindowedProfiles.yw25CTime.Test.std(i,ap_index,ch_index)/...
                        sqrt(CompiledEmbryos.WindowedProfiles.yw25CTime.Test.count(i,ap_index,ch_index));
                elseif sum(TFGoodAP)  == 1
                    TestProfs = DorsalProfiles(TFtestBin,ap_index,ch_index);
                    CompiledEmbryos.WindowedProfiles.yw25CTime.Test.mean(i,ap_index,ch_index) = TestProfs(TFGoodAP);
                    CompiledEmbryos.WindowedProfiles.yw25CTime.Test.std(i,ap_index,ch_index) = NaN;
                    CompiledEmbryos.WindowedProfiles.yw25CTime.Test.count(i,ap_index,ch_index) = 1;
                    CompiledEmbryos.WindowedProfiles.yw25CTime.Test.se(i,ap_index,ch_index) = NaN;
                end
            end
        end
        
        if sum(TFcontrolBin) >= 1
            for ap_index = 1:NumAPbins
                TFGoodAP = ~isnan(DorsalProfiles(TFcontrolBin,ap_index,ch_index));
                if sum(TFGoodAP) > 1
                    CompiledEmbryos.WindowedProfiles.yw25CTime.Control.mean(i,ap_index,ch_index) = mean(DorsalProfiles(TFcontrolBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.WindowedProfiles.yw25CTime.Control.std(i,ap_index,ch_index) = std(DorsalProfiles(TFcontrolBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.WindowedProfiles.yw25CTime.Control.count(i,ap_index,ch_index) = sum(~isnan(DorsalProfiles(TFcontrolBin,ap_index,ch_index)));
                    CompiledEmbryos.WindowedProfiles.yw25CTime.Control.se(i,ap_index,ch_index) = CompiledEmbryos.WindowedProfiles.yw25CTime.Control.std(i,ap_index,ch_index)/...
                        sqrt(CompiledEmbryos.WindowedProfiles.yw25CTime.Control.count(i,ap_index,ch_index));
                elseif sum(TFGoodAP)  == 1
                    ControlProfs = DorsalProfiles(TFcontrolBin,ap_index,ch_index);
                    CompiledEmbryos.WindowedProfiles.yw25CTime.Control.mean(i,ap_index,ch_index) = ControlProfs(TFGoodAP);
                    CompiledEmbryos.WindowedProfiles.yw25CTime.Control.std(i,ap_index,ch_index) = NaN;
                    CompiledEmbryos.WindowedProfiles.yw25CTime.Control.count(i,ap_index,ch_index) = 1;
                    CompiledEmbryos.WindowedProfiles.yw25CTime.Control.se(i,ap_index,ch_index) = NaN;
                end
            end
        end
    end
end

%% Narrow windowed Profiles - yw25C Times
NarrowDorsalProfiles = CompiledEmbryos.SlideRescaledDorsalAvgNarrowAPProfiles;
CompiledEmbryos.NarrowWindowedProfiles.yw25CTime = {};
CompiledEmbryos.NarrowWindowedProfiles.yw25CTime.BinWidth = deltafc_binwidth;
CompiledEmbryos.NarrowWindowedProfiles.yw25CTime.x = x;
CompiledEmbryos.NarrowWindowedProfiles.yw25CTime.Test = {};
CompiledEmbryos.NarrowWindowedProfiles.yw25CTime.Test.mean = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowWindowedProfiles.yw25CTime.Test.std = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowWindowedProfiles.yw25CTime.Test.se = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowWindowedProfiles.yw25CTime.Test.count = zeros(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowWindowedProfiles.yw25CTime.Control = {}; 
CompiledEmbryos.NarrowWindowedProfiles.yw25CTime.Control.mean = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowWindowedProfiles.yw25CTime.Control.std = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowWindowedProfiles.yw25CTime.Control.se = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowWindowedProfiles.yw25CTime.Control.count = zeros(NumTbins, NumNarrowAPbins, NChannels);


for i = 1:NumTbins
    TFtimeBin = CompiledEmbryos.yw25CEmbryoTimes >= x(i)-yw25Ctime_binwidth/2 & CompiledEmbryos.yw25CEmbryoTimes < x(i)+yw25Ctime_binwidth/2;
    TFtestBin = TFtimeBin & Allyw25CTimeValidProfilesTestTF;
    TFcontrolBin = TFtimeBin & Allyw25CTimeProfilesControlTF;
    for ch_index = 2:NChannels
        if sum(TFtestBin) >= 1
            for ap_index = 1:NumNarrowAPbins
                TFGoodAP = ~isnan(NarrowDorsalProfiles(TFtestBin,ap_index,ch_index));
                if sum(TFGoodAP) > 1
                    CompiledEmbryos.NarrowWindowedProfiles.yw25CTime.Test.mean(i,ap_index,ch_index) = mean(NarrowDorsalProfiles(TFtestBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.NarrowWindowedProfiles.yw25CTime.Test.std(i,ap_index,ch_index) = std(NarrowDorsalProfiles(TFtestBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.NarrowWindowedProfiles.yw25CTime.Test.count(i,ap_index,ch_index) = sum(~isnan(NarrowDorsalProfiles(TFtestBin,ap_index,ch_index)));
                    CompiledEmbryos.NarrowWindowedProfiles.yw25CTime.Test.se(i,ap_index,ch_index) = CompiledEmbryos.NarrowWindowedProfiles.yw25CTime.Test.std(i,ap_index,ch_index)/...
                        sqrt(CompiledEmbryos.NarrowWindowedProfiles.yw25CTime.Test.count(i,ap_index,ch_index));
                elseif sum(TFGoodAP)  == 1
                    NarrowTestProfs = NarrowDorsalProfiles(TFtestBin,ap_index,ch_index);
                    CompiledEmbryos.NarrowWindowedProfiles.yw25CTime.Test.mean(i,ap_index,ch_index) = NarrowTestProfs(TFGoodAP);
                    CompiledEmbryos.NarrowWindowedProfiles.yw25CTime.Test.std(i,ap_index,ch_index) = NaN;
                    CompiledEmbryos.NarrowWindowedProfiles.yw25CTime.Test.count(i,ap_index,ch_index) = 1;
                    CompiledEmbryos.NarrowWindowedProfiles.yw25CTime.Test.se(i,ap_index,ch_index) = NaN;
                end
            end
        end
        
        if sum(TFcontrolBin) >= 1
            for ap_index = 1:NumNarrowAPbins
                TFGoodAP = ~isnan(NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index));
                if sum(TFGoodAP) > 1
                    CompiledEmbryos.NarrowWindowedProfiles.yw25CTime.Control.mean(i,ap_index,ch_index) = mean(NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.NarrowWindowedProfiles.yw25CTime.Control.std(i,ap_index,ch_index) = std(NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.NarrowWindowedProfiles.yw25CTime.Control.count(i,ap_index,ch_index) = sum(~isnan(NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index)));
                    CompiledEmbryos.NarrowWindowedProfiles.yw25CTime.Control.se(i,ap_index,ch_index) = CompiledEmbryos.NarrowWindowedProfiles.yw25CTime.Control.std(i,ap_index,ch_index)/...
                        sqrt(CompiledEmbryos.NarrowWindowedProfiles.yw25CTime.Control.count(i,ap_index,ch_index));
                elseif sum(TFGoodAP)  == 1
                    NarrowControlProfs = NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index);
                    CompiledEmbryos.NarrowWindowedProfiles.yw25CTime.Control.mean(i,ap_index,ch_index) = NarrowControlProfs(TFGoodAP);
                    CompiledEmbryos.NarrowWindowedProfiles.yw25CTime.Control.std(i,ap_index,ch_index) = NaN;
                    CompiledEmbryos.NarrowWindowedProfiles.yw25CTime.Control.count(i,ap_index,ch_index) = 1;
                    CompiledEmbryos.NarrowWindowedProfiles.yw25CTime.Control.se(i,ap_index,ch_index) = NaN;
                end
            end
        end
    end
end


%% Add HisRFP 25C Time-based Binning - standard AP
AllHisRFP25CTimeValidProfilesTestTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.TestSetEmbryos & ~isnan(CompiledEmbryos.HisRFP25CEmbryoTimes);% &...
AllHisRFP25CTimeProfilesControlTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.ControlSetEmbryos & ~isnan(CompiledEmbryos.HisRFP25CEmbryoTimes);% &...
        %CompiledEmbryos.FixCorrectedDeltaFC_um(:,1).' > 3;


x = 0:hisrfp25Ctime_binwidth:65;
NumTbins = length(x);

DorsalProfiles = CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles;
CompiledEmbryos.BinnedProfiles.HisRFP25CTime = {};
CompiledEmbryos.BinnedProfiles.HisRFP25CTime.BinWidth = hisrfp25Ctime_binwidth;
CompiledEmbryos.BinnedProfiles.HisRFP25CTime.x = x;
CompiledEmbryos.BinnedProfiles.HisRFP25CTime.Test = {};
CompiledEmbryos.BinnedProfiles.HisRFP25CTime.Test.mean = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.BinnedProfiles.HisRFP25CTime.Test.std = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.BinnedProfiles.HisRFP25CTime.Test.se = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.BinnedProfiles.HisRFP25CTime.Test.count = zeros(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.BinnedProfiles.HisRFP25CTime.Control = {}; 
CompiledEmbryos.BinnedProfiles.HisRFP25CTime.Control.mean = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.BinnedProfiles.HisRFP25CTime.Control.std = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.BinnedProfiles.HisRFP25CTime.Control.se = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.BinnedProfiles.HisRFP25CTime.Control.count = zeros(NumTbins, NumAPbins, NChannels);


for i = 1:NumTbins
    TFtimeBin = CompiledEmbryos.HisRFP25CEmbryoTimes >= x(i)-hisrfp25Ctime_binwidth/2 & CompiledEmbryos.HisRFP25CEmbryoTimes < x(i)+hisrfp25Ctime_binwidth/2;
    TFtestBin = TFtimeBin & AllHisRFP25CTimeValidProfilesTestTF;
    TFcontrolBin = TFtimeBin & AllHisRFP25CTimeProfilesControlTF;
    for ch_index = 2:NChannels
        if sum(TFtestBin) >= 1
            for ap_index = 1:NumAPbins
                TFGoodAP = ~isnan(DorsalProfiles(TFtestBin,ap_index,ch_index));
                if sum(TFGoodAP) > 1
                    CompiledEmbryos.BinnedProfiles.HisRFP25CTime.Test.mean(i,ap_index,ch_index) = mean(DorsalProfiles(TFtestBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.BinnedProfiles.HisRFP25CTime.Test.std(i,ap_index,ch_index) = std(DorsalProfiles(TFtestBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.BinnedProfiles.HisRFP25CTime.Test.count(i,ap_index,ch_index) = sum(~isnan(DorsalProfiles(TFtestBin,ap_index,ch_index)));
                    CompiledEmbryos.BinnedProfiles.HisRFP25CTime.Test.se(i,ap_index,ch_index) = CompiledEmbryos.BinnedProfiles.HisRFP25CTime.Test.std(i,ap_index,ch_index)/...
                        sqrt(CompiledEmbryos.BinnedProfiles.HisRFP25CTime.Test.count(i,ap_index,ch_index));
                elseif sum(TFGoodAP)  == 1
                    TestProfs = DorsalProfiles(TFtestBin,ap_index,ch_index);
                    CompiledEmbryos.BinnedProfiles.HisRFP25CTime.Test.mean(i,ap_index,ch_index) = TestProfs(TFGoodAP);
                    CompiledEmbryos.BinnedProfiles.HisRFP25CTime.Test.std(i,ap_index,ch_index) = NaN;
                    CompiledEmbryos.BinnedProfiles.HisRFP25CTime.Test.count(i,ap_index,ch_index) = 1;
                    CompiledEmbryos.BinnedProfiles.HisRFP25CTime.Test.se(i,ap_index,ch_index) = NaN;
                end
            end
        end
        
        if sum(TFcontrolBin) >= 1
            for ap_index = 1:NumAPbins
                TFGoodAP = ~isnan(DorsalProfiles(TFcontrolBin,ap_index,ch_index));
                if sum(TFGoodAP) > 1
                    CompiledEmbryos.BinnedProfiles.HisRFP25CTime.Control.mean(i,ap_index,ch_index) = mean(DorsalProfiles(TFcontrolBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.BinnedProfiles.HisRFP25CTime.Control.std(i,ap_index,ch_index) = std(DorsalProfiles(TFcontrolBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.BinnedProfiles.HisRFP25CTime.Control.count(i,ap_index,ch_index) = sum(~isnan(DorsalProfiles(TFcontrolBin,ap_index,ch_index)));
                    CompiledEmbryos.BinnedProfiles.HisRFP25CTime.Control.se(i,ap_index,ch_index) = CompiledEmbryos.BinnedProfiles.HisRFP25CTime.Control.std(i,ap_index,ch_index)/...
                        sqrt(CompiledEmbryos.BinnedProfiles.HisRFP25CTime.Control.count(i,ap_index,ch_index));
                elseif sum(TFGoodAP)  == 1
                    ControlProfs = DorsalProfiles(TFcontrolBin,ap_index,ch_index);
                    CompiledEmbryos.BinnedProfiles.HisRFP25CTime.Control.mean(i,ap_index,ch_index) = ControlProfs(TFGoodAP);
                    CompiledEmbryos.BinnedProfiles.HisRFP25CTime.Control.std(i,ap_index,ch_index) = NaN;
                    CompiledEmbryos.BinnedProfiles.HisRFP25CTime.Control.count(i,ap_index,ch_index) = 1;
                    CompiledEmbryos.BinnedProfiles.HisRFP25CTime.Control.se(i,ap_index,ch_index) = NaN;
                end
            end
        end
    end
end

%% Add HisRFP Time-based Binning - narrow AP
NarrowDorsalProfiles = CompiledEmbryos.SlideRescaledDorsalAvgNarrowAPProfiles;
CompiledEmbryos.NarrowBinnedProfiles.HisRFP25CTime = {};
CompiledEmbryos.NarrowBinnedProfiles.HisRFP25CTime.BinWidth = deltafc_binwidth;
CompiledEmbryos.NarrowBinnedProfiles.HisRFP25CTime.x = x;
CompiledEmbryos.NarrowBinnedProfiles.HisRFP25CTime.Test = {};
CompiledEmbryos.NarrowBinnedProfiles.HisRFP25CTime.Test.mean = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowBinnedProfiles.HisRFP25CTime.Test.std = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowBinnedProfiles.HisRFP25CTime.Test.se = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowBinnedProfiles.HisRFP25CTime.Test.count = zeros(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowBinnedProfiles.HisRFP25CTime.Control = {}; 
CompiledEmbryos.NarrowBinnedProfiles.HisRFP25CTime.Control.mean = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowBinnedProfiles.HisRFP25CTime.Control.std = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowBinnedProfiles.HisRFP25CTime.Control.se = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowBinnedProfiles.HisRFP25CTime.Control.count = zeros(NumTbins, NumNarrowAPbins, NChannels);


for i = 1:NumTbins
    TFtimeBin = CompiledEmbryos.HisRFP25CEmbryoTimes >= x(i)-hisrfp25Ctime_binwidth/2 & CompiledEmbryos.HisRFP25CEmbryoTimes < x(i)+hisrfp25Ctime_binwidth/2;
    TFtestBin = TFtimeBin & AllHisRFP25CTimeValidProfilesTestTF;
    TFcontrolBin = TFtimeBin & AllHisRFP25CTimeProfilesControlTF;
    for ch_index = 2:NChannels
        if sum(TFtestBin) >= 1
            for ap_index = 1:NumNarrowAPbins
                TFGoodAP = ~isnan(NarrowDorsalProfiles(TFtestBin,ap_index,ch_index));
                if sum(TFGoodAP) > 1
                    CompiledEmbryos.NarrowBinnedProfiles.HisRFP25CTime.Test.mean(i,ap_index,ch_index) = mean(NarrowDorsalProfiles(TFtestBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.NarrowBinnedProfiles.HisRFP25CTime.Test.std(i,ap_index,ch_index) = std(NarrowDorsalProfiles(TFtestBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.NarrowBinnedProfiles.HisRFP25CTime.Test.count(i,ap_index,ch_index) = sum(~isnan(NarrowDorsalProfiles(TFtestBin,ap_index,ch_index)));
                    CompiledEmbryos.NarrowBinnedProfiles.HisRFP25CTime.Test.se(i,ap_index,ch_index) = CompiledEmbryos.NarrowBinnedProfiles.HisRFP25CTime.Test.std(i,ap_index,ch_index)/...
                        sqrt(CompiledEmbryos.NarrowBinnedProfiles.HisRFP25CTime.Test.count(i,ap_index,ch_index));
                elseif sum(TFGoodAP)  == 1
                    NarrowTestProfs = NarrowDorsalProfiles(TFtestBin,ap_index,ch_index);
                    CompiledEmbryos.NarrowBinnedProfiles.HisRFP25CTime.Test.mean(i,ap_index,ch_index) = NarrowTestProfs(TFGoodAP);
                    CompiledEmbryos.NarrowBinnedProfiles.HisRFP25CTime.Test.std(i,ap_index,ch_index) = NaN;
                    CompiledEmbryos.NarrowBinnedProfiles.HisRFP25CTime.Test.count(i,ap_index,ch_index) = 1;
                    CompiledEmbryos.NarrowBinnedProfiles.HisRFP25CTime.Test.se(i,ap_index,ch_index) = NaN;
                end
            end
        end
        
        if sum(TFcontrolBin) >= 1
            for ap_index = 1:NumNarrowAPbins
                TFGoodAP = ~isnan(NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index));
                if sum(TFGoodAP) > 1
                    CompiledEmbryos.NarrowBinnedProfiles.HisRFP25CTime.Control.mean(i,ap_index,ch_index) = mean(NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.NarrowBinnedProfiles.HisRFP25CTime.Control.std(i,ap_index,ch_index) = std(NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.NarrowBinnedProfiles.HisRFP25CTime.Control.count(i,ap_index,ch_index) = sum(~isnan(NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index)));
                    CompiledEmbryos.NarrowBinnedProfiles.HisRFP25CTime.Control.se(i,ap_index,ch_index) = CompiledEmbryos.NarrowBinnedProfiles.HisRFP25CTime.Control.std(i,ap_index,ch_index)/...
                        sqrt(CompiledEmbryos.NarrowBinnedProfiles.HisRFP25CTime.Control.count(i,ap_index,ch_index));
                elseif sum(TFGoodAP)  == 1
                    NarrowControlProfs = NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index);
                    CompiledEmbryos.NarrowBinnedProfiles.HisRFP25CTime.Control.mean(i,ap_index,ch_index) = NarrowControlProfs(TFGoodAP);
                    CompiledEmbryos.NarrowBinnedProfiles.HisRFP25CTime.Control.std(i,ap_index,ch_index) = NaN;
                    CompiledEmbryos.NarrowBinnedProfiles.HisRFP25CTime.Control.count(i,ap_index,ch_index) = 1;
                    CompiledEmbryos.NarrowBinnedProfiles.HisRFP25CTime.Control.se(i,ap_index,ch_index) = NaN;
                end
            end
        end
    end
end


%% Windowed Profiles - HisRFP 25C Times

x = 0:1:65;
NumTbins = length(x);

DorsalProfiles = CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles;
CompiledEmbryos.WindowedProfiles.HisRFP25CTime = {};
CompiledEmbryos.WindowedProfiles.HisRFP25CTime.BinWidth = hisrfp25Ctime_binwidth;
CompiledEmbryos.WindowedProfiles.HisRFP25CTime.x = x;
CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Test = {};
CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Test.mean = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Test.std = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Test.se = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Test.count = zeros(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Control = {}; 
CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Control.mean = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Control.std = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Control.se = NaN(NumTbins, NumAPbins, NChannels);
CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Control.count = zeros(NumTbins, NumAPbins, NChannels);


for i = 1:NumTbins
    TFtimeBin = CompiledEmbryos.HisRFP25CEmbryoTimes >= x(i)-hisrfp25Ctime_binwidth/2 & CompiledEmbryos.HisRFP25CEmbryoTimes < x(i)+hisrfp25Ctime_binwidth/2;
    TFtestBin = TFtimeBin & AllHisRFP25CTimeValidProfilesTestTF;
    TFcontrolBin = TFtimeBin & AllHisRFP25CTimeProfilesControlTF;
    for ch_index = 2:NChannels
        if sum(TFtestBin) >= 1
            for ap_index = 1:NumAPbins
                TFGoodAP = ~isnan(DorsalProfiles(TFtestBin,ap_index,ch_index));
                if sum(TFGoodAP) > 1
                    CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Test.mean(i,ap_index,ch_index) = mean(DorsalProfiles(TFtestBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Test.std(i,ap_index,ch_index) = std(DorsalProfiles(TFtestBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Test.count(i,ap_index,ch_index) = sum(~isnan(DorsalProfiles(TFtestBin,ap_index,ch_index)));
                    CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Test.se(i,ap_index,ch_index) = CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Test.std(i,ap_index,ch_index)/...
                        sqrt(CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Test.count(i,ap_index,ch_index));
                elseif sum(TFGoodAP)  == 1
                    TestProfs = DorsalProfiles(TFtestBin,ap_index,ch_index);
                    CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Test.mean(i,ap_index,ch_index) = TestProfs(TFGoodAP);
                    CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Test.std(i,ap_index,ch_index) = NaN;
                    CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Test.count(i,ap_index,ch_index) = 1;
                    CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Test.se(i,ap_index,ch_index) = NaN;
                end
            end
        end
        
        if sum(TFcontrolBin) >= 1
            for ap_index = 1:NumAPbins
                TFGoodAP = ~isnan(DorsalProfiles(TFcontrolBin,ap_index,ch_index));
                if sum(TFGoodAP) > 1
                    CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Control.mean(i,ap_index,ch_index) = mean(DorsalProfiles(TFcontrolBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Control.std(i,ap_index,ch_index) = std(DorsalProfiles(TFcontrolBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Control.count(i,ap_index,ch_index) = sum(~isnan(DorsalProfiles(TFcontrolBin,ap_index,ch_index)));
                    CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Control.se(i,ap_index,ch_index) = CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Control.std(i,ap_index,ch_index)/...
                        sqrt(CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Control.count(i,ap_index,ch_index));
                elseif sum(TFGoodAP)  == 1
                    ControlProfs = DorsalProfiles(TFcontrolBin,ap_index,ch_index);
                    CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Control.mean(i,ap_index,ch_index) = ControlProfs(TFGoodAP);
                    CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Control.std(i,ap_index,ch_index) = NaN;
                    CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Control.count(i,ap_index,ch_index) = 1;
                    CompiledEmbryos.WindowedProfiles.HisRFP25CTime.Control.se(i,ap_index,ch_index) = NaN;
                end
            end
        end
    end
end

%% Narrow windowed Profiles - HisRFP 25C Times
NarrowDorsalProfiles = CompiledEmbryos.SlideRescaledDorsalAvgNarrowAPProfiles;
CompiledEmbryos.NarrowWindowedProfiles.HisRFP25CTime = {};
CompiledEmbryos.NarrowWindowedProfiles.HisRFP25CTime.BinWidth = deltafc_binwidth;
CompiledEmbryos.NarrowWindowedProfiles.HisRFP25CTime.x = x;
CompiledEmbryos.NarrowWindowedProfiles.HisRFP25CTime.Test = {};
CompiledEmbryos.NarrowWindowedProfiles.HisRFP25CTime.Test.mean = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowWindowedProfiles.HisRFP25CTime.Test.std = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowWindowedProfiles.HisRFP25CTime.Test.se = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowWindowedProfiles.HisRFP25CTime.Test.count = zeros(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowWindowedProfiles.HisRFP25CTime.Control = {}; 
CompiledEmbryos.NarrowWindowedProfiles.HisRFP25CTime.Control.mean = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowWindowedProfiles.HisRFP25CTime.Control.std = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowWindowedProfiles.HisRFP25CTime.Control.se = NaN(NumTbins, NumNarrowAPbins, NChannels);
CompiledEmbryos.NarrowWindowedProfiles.HisRFP25CTime.Control.count = zeros(NumTbins, NumNarrowAPbins, NChannels);


for i = 1:NumTbins
    TFtimeBin = CompiledEmbryos.HisRFP25CEmbryoTimes >= x(i)-hisrfp25Ctime_binwidth/2 & CompiledEmbryos.HisRFP25CEmbryoTimes < x(i)+hisrfp25Ctime_binwidth/2;
    TFtestBin = TFtimeBin & AllHisRFP25CTimeValidProfilesTestTF;
    TFcontrolBin = TFtimeBin & AllHisRFP25CTimeProfilesControlTF;
    for ch_index = 2:NChannels
        if sum(TFtestBin) >= 1
            for ap_index = 1:NumNarrowAPbins
                TFGoodAP = ~isnan(NarrowDorsalProfiles(TFtestBin,ap_index,ch_index));
                if sum(TFGoodAP) > 1
                    CompiledEmbryos.NarrowWindowedProfiles.HisRFP25CTime.Test.mean(i,ap_index,ch_index) = mean(NarrowDorsalProfiles(TFtestBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.NarrowWindowedProfiles.HisRFP25CTime.Test.std(i,ap_index,ch_index) = std(NarrowDorsalProfiles(TFtestBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.NarrowWindowedProfiles.HisRFP25CTime.Test.count(i,ap_index,ch_index) = sum(~isnan(NarrowDorsalProfiles(TFtestBin,ap_index,ch_index)));
                    CompiledEmbryos.NarrowWindowedProfiles.HisRFP25CTime.Test.se(i,ap_index,ch_index) = CompiledEmbryos.NarrowWindowedProfiles.HisRFP25CTime.Test.std(i,ap_index,ch_index)/...
                        sqrt(CompiledEmbryos.NarrowWindowedProfiles.HisRFP25CTime.Test.count(i,ap_index,ch_index));
                elseif sum(TFGoodAP)  == 1
                    NarrowTestProfs = NarrowDorsalProfiles(TFtestBin,ap_index,ch_index);
                    CompiledEmbryos.NarrowWindowedProfiles.HisRFP25CTime.Test.mean(i,ap_index,ch_index) = NarrowTestProfs(TFGoodAP);
                    CompiledEmbryos.NarrowWindowedProfiles.HisRFP25CTime.Test.std(i,ap_index,ch_index) = NaN;
                    CompiledEmbryos.NarrowWindowedProfiles.HisRFP25CTime.Test.count(i,ap_index,ch_index) = 1;
                    CompiledEmbryos.NarrowWindowedProfiles.HisRFP25CTime.Test.se(i,ap_index,ch_index) = NaN;
                end
            end
        end
        
        if sum(TFcontrolBin) >= 1
            for ap_index = 1:NumNarrowAPbins
                TFGoodAP = ~isnan(NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index));
                if sum(TFGoodAP) > 1
                    CompiledEmbryos.NarrowWindowedProfiles.HisRFP25CTime.Control.mean(i,ap_index,ch_index) = mean(NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.NarrowWindowedProfiles.HisRFP25CTime.Control.std(i,ap_index,ch_index) = std(NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index), 'omitnan');
                    CompiledEmbryos.NarrowWindowedProfiles.HisRFP25CTime.Control.count(i,ap_index,ch_index) = sum(~isnan(NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index)));
                    CompiledEmbryos.NarrowWindowedProfiles.HisRFP25CTime.Control.se(i,ap_index,ch_index) = CompiledEmbryos.NarrowWindowedProfiles.HisRFP25CTime.Control.std(i,ap_index,ch_index)/...
                        sqrt(CompiledEmbryos.NarrowWindowedProfiles.HisRFP25CTime.Control.count(i,ap_index,ch_index));
                elseif sum(TFGoodAP)  == 1
                    NarrowControlProfs = NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index);
                    CompiledEmbryos.NarrowWindowedProfiles.HisRFP25CTime.Control.mean(i,ap_index,ch_index) = NarrowControlProfs(TFGoodAP);
                    CompiledEmbryos.NarrowWindowedProfiles.HisRFP25CTime.Control.std(i,ap_index,ch_index) = NaN;
                    CompiledEmbryos.NarrowWindowedProfiles.HisRFP25CTime.Control.count(i,ap_index,ch_index) = 1;
                    CompiledEmbryos.NarrowWindowedProfiles.HisRFP25CTime.Control.se(i,ap_index,ch_index) = NaN;
                end
            end
        end
    end
end










% disp('test button')




