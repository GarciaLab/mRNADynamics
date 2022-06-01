% AnalyzeSpotFluos.m
% author: Gabriella Martini
% date last modified: 1/26/22
%
% Compares 60mer 50uW 20C measurements to 120mer 50uW 20C measurements (all
% taken on 1/19/22)
clear all
plotdir =  ['S:/Gabriella/Dropbox/FixedTissueExperiments/StandardCandles/202202023_GoodSpotsNoFilters'];
mkdir(plotdir);

Prefixes_60mer_50uW = {'2022-01-19-eVasax5-60mer-GFP-lamin-T20C_856x856_LA1_Zoom2_zstep025_50uW-E1',...
    '2022-01-19-eVasax5-60mer-GFP-lamin-T20C_856x856_LA1_Zoom2_zstep025_50uW-E2',...
    '2022-01-19-eVasax5-60mer-GFP-lamin-T20C_856x856_LA1_Zoom2_zstep025_50uW-E3',...
    '2022-01-19-eVasax5-60mer-GFP-lamin-T20C_856x856_LA1_Zoom2_zstep025_50uW-E4',...
    '2022-01-19-eVasax5-60mer-GFP-lamin-T20C_856x856_LA1_Zoom2_zstep025_50uW-E5',...
    };

outfolder_60mer_50uW = '2022-01-19_20C_60mer_50uW';
outdir_60mer_50uW  = ['S:/Gabriella/Dropbox\StandardCandles\' outfolder_60mer_50uW];
disp('60mer 50uW Good Spots')
% for i = 1:5
%     disp(['i = ', num2str(i)]);
%     AddNewGaussianIntensityEstimates(Prefixes_60mer_50uW{i});
% end
% 
    CombineSpotData(Prefixes_60mer_50uW,outfolder_60mer_50uW);
    SpotInfo_60mer_50uW  = load([outdir_60mer_50uW '\SpotInfo.mat']);
% 
% 
Prefixes_120mer_50uW = {'2022-01-19-eVasax5-120mer-GFP-lamin-T20C_856x856_LA1_Zoom2_zstep025_50uW-E1',...
    '2022-01-19-eVasax5-120mer-GFP-lamin-T20C_856x856_LA1_Zoom2_zstep025_50uW-E2',...
    '2022-01-19-eVasax5-120mer-GFP-lamin-T20C_856x856_LA1_Zoom2_zstep025_50uW-E3',...
    '2022-01-19-eVasax5-120mer-GFP-lamin-T20C_856x856_LA1_Zoom2_zstep025_50uW-E4',...
    '2022-01-19-eVasax5-120mer-GFP-lamin-T20C_856x856_LA1_Zoom2_zstep025_50uW-E5',...
    };
outfolder_120mer_50uW = '2022-01-19_20C_120mer_50uW';
outdir_120mer_50uW  = ['S:/Gabriella/Dropbox\StandardCandles\' outfolder_120mer_50uW];


    CombineSpotData(Prefixes_120mer_50uW,outfolder_120mer_50uW);
    SpotInfo_120mer_50uW = load([outdir_120mer_50uW '\SpotInfo.mat']);

% disp('120mer 50uW Good Spots')
% for i = 1:5
%     disp(['i = ', num2str(i)]);
%     AddNewGaussianIntensityEstimates(Prefixes_120mer_50uW{i});
% end
% 
% disp('60mer 50uW All Spots')
% for i = 1:5
%     disp(['i = ', num2str(i)]);
%     AddNewGaussianIntensityEstimates(Prefixes_60mer_50uW{i},'useallspots');
% end
% 
% 
% disp('120mer 50uW All Spots')
% for i = 1:5
%     disp(['i = ', num2str(i)]);
%     AddNewGaussianIntensityEstimates(Prefixes_120mer_50uW{i},'useallspots');
% end
%%
sigx_limits=[0.45, 0.9];
sigy_limits = [0.5, 1.0];
sigsum_limits = [1.05, 1.75];
sigdiff_limits =[0 0.25];
Flags_120mer = ones(1, length(SpotInfo_120mer_50uW.GaussianIntensitySmallSnip), 'logical');
% Flags_120mer = Flags_120mer & SpotInfo_120mer_50uW.GaussianWidthXSmallSnip >= sigx_limits(1);
% Flags_120mer = Flags_120mer & SpotInfo_120mer_50uW.GaussianWidthXSmallSnip <= sigx_limits(2);
% Flags_120mer = Flags_120mer & SpotInfo_120mer_50uW.GaussianWidthYSmallSnip >= sigy_limits(1);
% Flags_120mer = Flags_120mer & SpotInfo_120mer_50uW.GaussianWidthYSmallSnip <= sigy_limits(2);
% Flags_120mer = Flags_120mer & SpotInfo_120mer_50uW.GaussianWidthXSmallSnip+SpotInfo_120mer_50uW.GaussianWidthYSmallSnip >= sigsum_limits(1);
% Flags_120mer = Flags_120mer & SpotInfo_120mer_50uW.GaussianWidthXSmallSnip+SpotInfo_120mer_50uW.GaussianWidthYSmallSnip <= sigsum_limits(2);
% Flags_120mer = Flags_120mer & abs(SpotInfo_120mer_50uW.GaussianWidthXSmallSnip-SpotInfo_120mer_50uW.GaussianWidthYSmallSnip) >= sigdiff_limits(1);
% Flags_120mer = Flags_120mer & abs(SpotInfo_120mer_50uW.GaussianWidthXSmallSnip-SpotInfo_120mer_50uW.GaussianWidthYSmallSnip) <= sigdiff_limits(2);
Flags_120mer = Flags_120mer & SpotInfo_120mer_50uW.nearest_neighbors_Corrected > 10;
Flags_120mer = Flags_120mer & SpotInfo_120mer_50uW.nearest_neighbors_SS > 10;
Flags_120mer = Flags_120mer & SpotInfo_120mer_50uW.nearest_neighbors_kernel > 10;

Flags_60mer = ones(1, length(SpotInfo_60mer_50uW.GaussianIntensitySmallSnip), 'logical');
% Flags_60mer = Flags_60mer & SpotInfo_60mer_50uW.GaussianWidthXSmallSnip >= sigx_limits(1);
% Flags_60mer = Flags_60mer & SpotInfo_60mer_50uW.GaussianWidthXSmallSnip <= sigx_limits(2);
% Flags_60mer = Flags_60mer & SpotInfo_60mer_50uW.GaussianWidthYSmallSnip >= sigy_limits(1);
% Flags_60mer = Flags_60mer & SpotInfo_60mer_50uW.GaussianWidthYSmallSnip <= sigy_limits(2);
% Flags_60mer = Flags_60mer & SpotInfo_60mer_50uW.GaussianWidthXSmallSnip+SpotInfo_60mer_50uW.GaussianWidthYSmallSnip >= sigsum_limits(1);
% Flags_60mer = Flags_60mer & SpotInfo_60mer_50uW.GaussianWidthXSmallSnip+SpotInfo_60mer_50uW.GaussianWidthYSmallSnip <= sigsum_limits(2);
% Flags_60mer = Flags_60mer & abs(SpotInfo_60mer_50uW.GaussianWidthXSmallSnip-SpotInfo_60mer_50uW.GaussianWidthYSmallSnip) >= sigdiff_limits(1);
% Flags_60mer = Flags_60mer & abs(SpotInfo_60mer_50uW.GaussianWidthXSmallSnip-SpotInfo_60mer_50uW.GaussianWidthYSmallSnip) <= sigdiff_limits(2);
Flags_60mer = Flags_60mer & SpotInfo_60mer_50uW.nearest_neighbors_Corrected > 10;
Flags_60mer = Flags_60mer & SpotInfo_60mer_50uW.nearest_neighbors_SS > 10;
Flags_60mer = Flags_60mer & SpotInfo_60mer_50uW.nearest_neighbors_kernel > 10;






%%
close all
Fig1 = figure(1);
set(Fig1,'units', 'normalized', 'position',[0.05, 0.05, 0.5, 0.5]);
set(gcf,'color','w');
ax1 = axes(Fig1);

h_120mer_50uW = histogram(SpotInfo_120mer_50uW.GaussianPeakSmallSnip(Flags_120mer), 'EdgeColor', 'none', 'Normalization', 'probability',...
    'BinEdges', 0:.1:max(SpotInfo_120mer_50uW.GaussianPeakSmallSnip(Flags_120mer)), 'FaceColor', 'blue', 'FaceAlpha', 0.6);
hold on
grid on
[f1_120mer_50uW, x1_120mer_50uW] = ksdensity(SpotInfo_120mer_50uW.GaussianPeakSmallSnip(Flags_120mer),0:0.1:max(SpotInfo_120mer_50uW.GaussianPeakSmallSnip(Flags_120mer)), 'Support', 'positive');
pl1 = plot(x1_120mer_50uW,f1_120mer_50uW*0.1, 'k-', 'LineWidth',2);
%xlim([0 300])
xlabel('Gaussian Peak Small Snip (AU)')
ylabel('p(I)')
[max1_120mer_50uW, idx1_120mer_50uW] = max(f1_120mer_50uW);
Fluo1_120mer_50uW = x1_120mer_50uW(idx1_120mer_50uW);
title('Gaussian Peak Small Snip 120mer 50uW 20ºC')
ax1.FontSize = 18;
hold off
xlim([ 0 100])
legend('Gaussian Peak Small Snip', ['Peak = ', num2str(round(Fluo1_120mer_50uW, 1)), ' AU']);
outpath = [plotdir, filesep, '120mer20C50uW_GaussianPeakSmallSnipDistWidth1.png'];
saveas(Fig1,outpath);

%%

Fig2 = figure(2);
set(Fig2,'units', 'normalized', 'position',[0.05, 0.05, 0.5, 0.5]);
set(gcf,'color','w');
ax2 = axes(Fig2);

h2_120mer_50uW = histogram(SpotInfo_120mer_50uW.GaussianIntensityCorrected(Flags_120mer), 'EdgeColor', 'none', 'Normalization', 'probability',...
    'BinEdges', 0:1:max(SpotInfo_120mer_50uW.GaussianIntensityCorrected(Flags_120mer)), 'FaceColor', 'blue', 'FaceAlpha', 0.6);
hold on
grid on
[f2_120mer_50uW, x2_120mer_50uW] = ksdensity(SpotInfo_120mer_50uW.GaussianIntensityCorrected(Flags_120mer),0:1:max(SpotInfo_120mer_50uW.GaussianIntensityCorrected(Flags_120mer)), 'Support', 'positive');
pl2 = plot(x2_120mer_50uW,f2_120mer_50uW, 'k-', 'LineWidth',2);
xlim([0 300])
xlabel('Gaussian Intensity Corrected (AU)')
ylabel('p(I)')
[max2_120mer_50uW, idx2_120mer_50uW] = max(f2_120mer_50uW);
Fluo2_120mer_50uW = x2_120mer_50uW(idx2_120mer_50uW);
title('Gaussian Intensity Corrected 120mer 50uW 20ºC')
ax2.FontSize = 18;
hold off
legend('Gaussian Intensity Corrected', ['Peak = ', num2str(round(Fluo2_120mer_50uW, 1)), ' AU']);
outpath = [plotdir, filesep, '120mer20C50uW_GaussianIntensityCorrectedDistWidth1.png'];
saveas(Fig2,outpath);

%%
Fig3 = figure(3);
set(Fig3,'units', 'normalized', 'position',[0.05, 0.05, 0.5, 0.5]);
set(gcf,'color','w');
ax3 = axes(Fig3);

h3_120mer_50uW = histogram(SpotInfo_120mer_50uW.GaussianIntensityCropped(Flags_120mer), 'EdgeColor', 'none', 'Normalization', 'probability',...
    'BinEdges', 0:1:max(SpotInfo_120mer_50uW.GaussianIntensityCropped(Flags_120mer)), 'FaceColor', 'blue', 'FaceAlpha', 0.6);
hold on
grid on
[f3_120mer_50uW, x3_120mer_50uW] = ksdensity(SpotInfo_120mer_50uW.GaussianIntensityCropped(Flags_120mer),0:1:max(SpotInfo_120mer_50uW.GaussianIntensityCropped(Flags_120mer)), 'Support', 'positive');
pl3 = plot(x3_120mer_50uW,f3_120mer_50uW, 'k-', 'LineWidth',2);
%xlim([0 300])
xlabel('Cropped Gaussian Intensity (AU)')
ylabel('p(I)')
[max3_120mer_50uW, idx3_120mer_50uW] = max(f3_120mer_50uW);
Fluo3_120mer_50uW = x3_120mer_50uW(idx3_120mer_50uW);
title('Cropped Gaussian Intensity 120mer 50uW 20ºC')
ax3.FontSize = 18;
hold off
xlim([0 200])
legend('Cropped Gaussian Intensity', ['Peak = ', num2str(round(Fluo3_120mer_50uW, 1)), ' AU']);
outpath = [plotdir, filesep, '120mer20C50uW_GaussianIntensityCroppedDistWidth1.png'];
saveas(Fig3,outpath);

%%
Fig4 = figure(4);
set(Fig4,'units', 'normalized', 'position',[0.05, 0.05, 0.5, 0.5]);
set(gcf,'color','w');
ax4 = axes(Fig4);

h4_120mer_50uW = histogram(SpotInfo_120mer_50uW.GaussianIntensitySmallSnip(Flags_120mer), 'EdgeColor', 'none', 'Normalization', 'probability',...
    'BinEdges', 0:1:max(SpotInfo_120mer_50uW.GaussianIntensitySmallSnip(Flags_120mer)), 'FaceColor', 'blue', 'FaceAlpha', 0.6);
hold on
grid on
[f4_120mer_50uW, x4_120mer_50uW] = ksdensity(SpotInfo_120mer_50uW.GaussianIntensitySmallSnip(~isnan(SpotInfo_120mer_50uW.GaussianIntensitySmallSnip) &SpotInfo_120mer_50uW.GaussianIntensitySmallSnip > 0 ),0:1:300, 'Support', 'positive');
pl4 = plot(x4_120mer_50uW,f4_120mer_50uW, 'k-', 'LineWidth',2);
%xlim([0 300])
xlabel('Small Snip Gaussian Intensity (AU)')
ylabel('p(I)')
[max4_120mer_50uW, idx4_120mer_50uW] = max(f4_120mer_50uW);
Fluo4_120mer_50uW = x4_120mer_50uW(idx4_120mer_50uW);
title('Small Snip Gaussian Intensity 120mer 50uW 20ºC')
ax4.FontSize = 18;
hold off
xlim([0 200])
legend('Small Snip Gaussian Intensity', ['Peak = ', num2str(round(Fluo4_120mer_50uW, 1)), ' AU']);
outpath = [plotdir, filesep, '120mer20C50uW_GaussianIntensitySmallSnipDistWidth1.png'];
saveas(Fig4,outpath);

%%
Fig5 = figure(5);
set(Fig5,'units', 'normalized', 'position',[0.05, 0.05, 0.5, 0.5]);
set(gcf,'color','w');
ax5 = axes(Fig5);

h5_120mer_50uW = histogram(SpotInfo_120mer_50uW.GaussianKernelIntensity(Flags_120mer), 'EdgeColor', 'none', 'Normalization', 'probability',...
    'BinEdges', 0:0.25:max(SpotInfo_120mer_50uW.GaussianKernelIntensity(Flags_120mer)), 'FaceColor', 'blue', 'FaceAlpha', 0.6);
hold on
grid on
[f5_120mer_50uW, x5_120mer_50uW] = ksdensity(SpotInfo_120mer_50uW.GaussianKernelIntensity(Flags_120mer), 0:0.25:max(SpotInfo_120mer_50uW.GaussianKernelIntensity(Flags_120mer)), 'Support', 'positive');
pl5 = plot(x5_120mer_50uW,f5_120mer_50uW*0.25, 'k-', 'LineWidth',2);
xlim([0 15])
xlabel('Gaussian Kernel Intensity (AU)')
ylabel('p(I)')
[max5_120mer_50uW, idx5_120mer_50uW] = max(f5_120mer_50uW);
Fluo5_120mer_50uW = x5_120mer_50uW(idx5_120mer_50uW);
title('Gaussian Kernel Intensity 120mer 50uW 20ºC')
ax5.FontSize = 18;
hold off
legend('Gaussian Kernel Intensity', ['Peak = ', num2str(round(Fluo5_120mer_50uW, 1)), ' AU']);
outpath = [plotdir, filesep, '120mer20C50uW_GaussianKernelIntensityDistWidth1.png'];
saveas(Fig5,outpath);

%%
Fig6 = figure(6);
set(Fig6,'units', 'normalized', 'position',[0.05, 0.05, 0.5, 0.5]);
set(gcf,'color','w');
ax6 = axes(Fig6);

h6_120mer_50uW = histogram(SpotInfo_120mer_50uW.IntegratedGaussianKernelIntensity(Flags_120mer), 'EdgeColor', 'none', 'Normalization', 'probability',...
    'BinEdges', 0:1:max(SpotInfo_120mer_50uW.IntegratedGaussianKernelIntensity(Flags_120mer)), 'FaceColor', 'blue', 'FaceAlpha', 0.6);
hold on
grid on
[f6_120mer_50uW, x6_120mer_50uW] = ksdensity(SpotInfo_120mer_50uW.IntegratedGaussianKernelIntensity(Flags_120mer),0:1:max(SpotInfo_120mer_50uW.IntegratedGaussianKernelIntensity(Flags_120mer)), 'Support', 'positive');
pl6 = plot(x6_120mer_50uW,f6_120mer_50uW, 'k-', 'LineWidth',2);
xlabel('Integrated Gaussian Kernel Intensity (AU)')
ylabel('p(I)')
[max6_120mer_50uW, idx6_120mer_50uW] = max(f6_120mer_50uW);
Fluo6_120mer_50uW = x6_120mer_50uW(idx6_120mer_50uW);
title('Integrated Gaussian Kernel Intensity 120mer 50uW 20ºC')
ax6.FontSize = 18;
hold off
xlim([0 200])
legend('Integrated Gaussian Kernel Intensity', ['Peak = ', num2str(round(Fluo6_120mer_50uW, 1)), ' AU']);
outpath = [plotdir, filesep, '120mer20C50uW_IntegratedGaussianKernelIntensityDistWidth1.png'];
saveas(Fig6,outpath);

%%
Fig7 = figure(7);
set(Fig7,'units', 'normalized', 'position',[0.05, 0.05, 0.5, 0.5]);
set(gcf,'color','w');
ax7 = axes(Fig7);

h_60mer_50uW = histogram(SpotInfo_60mer_50uW.GaussianPeakSmallSnip(Flags_60mer), 'EdgeColor', 'none', 'Normalization', 'probability',...
    'BinEdges', 0:.1:max(SpotInfo_60mer_50uW.GaussianPeakSmallSnip(Flags_60mer)), 'FaceColor', 'red', 'FaceAlpha', 0.6);
hold on
grid on
[f1_60mer_50uW, x1_60mer_50uW] = ksdensity(SpotInfo_60mer_50uW.GaussianPeakSmallSnip(Flags_60mer),0:.1:max(SpotInfo_60mer_50uW.GaussianPeakSmallSnip(Flags_60mer)), 'Support', 'positive');
pl1 = plot(x1_60mer_50uW,f1_60mer_50uW*.1, 'k-', 'LineWidth',2);
%xlim([0 300])
xlabel('Gaussian Peak Small Snip (AU)')
ylabel('p(I)')
[max1_60mer_50uW, idx1_60mer_50uW] = max(f1_60mer_50uW);
Fluo1_60mer_50uW = x1_60mer_50uW(idx1_60mer_50uW);
title('Gaussian Peak Small Snip  60mer 50uW 20ºC')
ax7.FontSize = 18;
hold off
legend('Gaussian Peak Small Snip ', ['Peak = ', num2str(round(Fluo1_60mer_50uW, 1)), ' AU']);
outpath = [plotdir, filesep, '60mer20C50uW_GaussianPeakSmallSnipDistWidth1.png'];
saveas(Fig7,outpath);

%%

Fig8 = figure(8);
set(Fig8,'units', 'normalized', 'position',[0.05, 0.05, 0.5, 0.5]);
set(gcf,'color','w');
ax8 = axes(Fig8);

h8_60mer_50uW = histogram(SpotInfo_60mer_50uW.GaussianIntensityCorrected(Flags_60mer), 'EdgeColor', 'none', 'Normalization', 'probability',...
    'BinEdges', 0:1:max(SpotInfo_60mer_50uW.GaussianIntensityCorrected(Flags_60mer)), 'FaceColor', 'red', 'FaceAlpha', 0.6);
hold on
grid on
[f8_60mer_50uW, x8_60mer_50uW] = ksdensity(SpotInfo_60mer_50uW.GaussianIntensityCorrected(Flags_60mer),0:1:max(SpotInfo_60mer_50uW.GaussianIntensityCorrected(Flags_60mer)), 'Support', 'positive');
pl8 = plot(x8_60mer_50uW,f8_60mer_50uW, 'k-', 'LineWidth',2);
%xlim([0 300])
xlabel('Gaussian Intensity Corrected (AU)')
ylabel('p(I)')
[max8_60mer_50uW, idx8_60mer_50uW] = max(f8_60mer_50uW);
Fluo8_60mer_50uW = x8_60mer_50uW(idx8_60mer_50uW);
title('Gaussian Intensity Corrected 60mer 50uW 20ºC')
ax8.FontSize = 18;
hold off
xlim([0 200])
legend('Gaussian Intensity Corrected', ['Peak = ', num2str(round(Fluo8_60mer_50uW, 1)), ' AU']);
outpath = [plotdir, filesep, '60mer20C50uW_GaussianIntensityCorrectedDistWidth1.png'];
saveas(Fig8,outpath);

%%
Fig9 = figure(9);
set(Fig9,'units', 'normalized', 'position',[0.05, 0.05, 0.5, 0.5]);
set(gcf,'color','w');
ax9 = axes(Fig9);

h9_60mer_50uW = histogram(SpotInfo_60mer_50uW.GaussianIntensityCropped(Flags_60mer), 'EdgeColor', 'none', 'Normalization', 'probability',...
    'BinEdges', 0:1:max(SpotInfo_60mer_50uW.GaussianIntensityCropped(Flags_60mer)), 'FaceColor', 'red', 'FaceAlpha', 0.6);
hold on
grid on
[f9_60mer_50uW, x9_60mer_50uW] = ksdensity(SpotInfo_60mer_50uW.GaussianIntensityCropped(Flags_60mer),0:1:max(SpotInfo_60mer_50uW.GaussianIntensityCropped(Flags_60mer)), 'Support', 'positive');
pl9 = plot(x9_60mer_50uW,f9_60mer_50uW, 'k-', 'LineWidth',2);
%xlim([0 300])
xlabel('Cropped Gaussian Intensity (AU)')
ylabel('p(I)')
[max9_60mer_50uW, idx9_60mer_50uW] = max(f9_60mer_50uW);
Fluo9_60mer_50uW = x9_60mer_50uW(idx9_60mer_50uW);
title('Cropped Gaussian Intensity 60mer 50uW 20ºC')
ax9.FontSize = 18;
hold off
xlim([0 200])
legend('Cropped Gaussian Intensity', ['Peak = ', num2str(round(Fluo9_60mer_50uW, 1)), ' AU']);
outpath = [plotdir, filesep, '60mer20C50uW_GaussianIntensityCroppedDistWidth1.png'];
saveas(Fig9,outpath);

%%
Fig10 = figure(10);
set(Fig10,'units', 'normalized', 'position',[0.05, 0.05, 0.5, 0.5]);
set(gcf,'color','w');
ax10 = axes(Fig10);

h10_60mer_50uW = histogram(SpotInfo_60mer_50uW.GaussianIntensitySmallSnip(Flags_60mer), 'EdgeColor', 'none', 'Normalization', 'probability',...
    'BinEdges', 0:1:max(SpotInfo_60mer_50uW.GaussianIntensitySmallSnip(Flags_60mer)), 'FaceColor', 'red', 'FaceAlpha', 0.6);
hold on
grid on
[f10_60mer_50uW, x10_60mer_50uW] = ksdensity(SpotInfo_60mer_50uW.GaussianIntensitySmallSnip(Flags_60mer),0:1:max(SpotInfo_60mer_50uW.GaussianIntensitySmallSnip(Flags_60mer)), 'Support', 'positive');
pl10 = plot(x10_60mer_50uW,f10_60mer_50uW, 'k-', 'LineWidth',2);
xlim([0 200])
xlabel('Small Snip Gaussian Intensity (AU)')
ylabel('p(I)')
[max10_60mer_50uW, idx10_60mer_50uW] = max(f10_60mer_50uW);
Fluo10_60mer_50uW = x10_60mer_50uW(idx10_60mer_50uW);
title('Small Snip Gaussian Intensity 120mer 50uW 20ºC')
ax10.FontSize = 18;
hold off
legend('Small Snip Gaussian Intensity', ['Peak = ', num2str(round(Fluo10_60mer_50uW, 1)), ' AU']);
outpath = [plotdir, filesep, '60mer20C50uW_GaussianIntensitySmallSnipDistWidth1.png'];
saveas(Fig10,outpath);

%%
Fig11 = figure(11);
set(Fig11,'units', 'normalized', 'position',[0.05, 0.05, 0.5, 0.5]);
set(gcf,'color','w');
ax11 = axes(Fig11);

h11_60mer_50uW = histogram(SpotInfo_60mer_50uW.GaussianKernelIntensity(Flags_60mer), 'EdgeColor', 'none', 'Normalization', 'probability',...
    'BinEdges', 0:0.1:max(SpotInfo_60mer_50uW.GaussianKernelIntensity(Flags_60mer)), 'FaceColor', 'red', 'FaceAlpha', 0.6);
hold on
grid on
[f11_60mer_50uW, x11_60mer_50uW] = ksdensity(SpotInfo_60mer_50uW.GaussianKernelIntensity(Flags_60mer),0:0.1:max(SpotInfo_60mer_50uW.GaussianKernelIntensity(Flags_60mer)), 'Support', 'positive');
pl11 = plot(x11_60mer_50uW,f11_60mer_50uW*.1, 'k-', 'LineWidth',2);
xlim([0 15])
xlabel('Gaussian Kernel Intensity (AU)')
ylabel('p(I)')
[max11_60mer_50uW, idx11_60mer_50uW] = max(f11_60mer_50uW);
Fluo11_60mer_50uW = x11_60mer_50uW(idx11_60mer_50uW);
title('Gaussian Kernel Intensity 60mer 50uW 20ºC')
ax11.FontSize = 18;
hold off
legend('Gaussian Kernel Intensity', ['Peak = ', num2str(round(Fluo11_60mer_50uW, 1)), ' AU']);
outpath = [plotdir, filesep, '60mer20C50uW_GaussianKernelIntensityDistWidth1.png'];
saveas(Fig11,outpath);

%%
Fig12 = figure(12);
set(Fig12,'units', 'normalized', 'position',[0.05, 0.05, 0.5, 0.5]);
set(gcf,'color','w');
ax12 = axes(Fig12);

h12_60mer_50uW = histogram(SpotInfo_60mer_50uW.IntegratedGaussianKernelIntensity(Flags_60mer), 'EdgeColor', 'none', 'Normalization', 'probability',...
    'BinEdges', 0:1:max(SpotInfo_60mer_50uW.IntegratedGaussianKernelIntensity(Flags_60mer)), 'FaceColor', 'red', 'FaceAlpha', 0.6);
hold on
grid on
[f12_60mer_50uW, x12_60mer_50uW] = ksdensity(SpotInfo_60mer_50uW.IntegratedGaussianKernelIntensity(Flags_60mer),0:1:max(SpotInfo_60mer_50uW.IntegratedGaussianKernelIntensity(Flags_60mer)), 'Support', 'positive');
pl12 = plot(x12_60mer_50uW,f12_60mer_50uW, 'k-', 'LineWidth',2);
xlim([0 200])
xlabel('Integrated Gaussian Kernel Intensity (AU)')
ylabel('p(I)')
[max12_60mer_50uW, idx12_60mer_50uW] = max(f12_60mer_50uW);
Fluo12_60mer_50uW = x12_60mer_50uW(idx12_60mer_50uW);
title('Integrated Gaussian Kernel Intensity 60mer 50uW 20ºC')
ax12.FontSize = 18;
hold off
legend('Integrated Gaussian Kernel Intensity', ['Peak = ', num2str(round(Fluo12_60mer_50uW, 1)), ' AU']);
outpath = [plotdir, filesep, '120mer20C50uW_IntegratedGaussianKernelIntensityDistWidth1.png'];
saveas(Fig12,outpath);



