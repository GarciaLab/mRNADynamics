% AddNewGaussianIntensityEstimates.m
% author: Gabriella Martini
% date last modified: 1/26/22
%

function AddNewGaussianIntensityEstimates(Prefix, varargin)
UseGoodSpots = true;


k=1;
while k<=length(varargin)
    if strcmpi(varargin{k},'useallspots')
        UseGoodSpots=false;
        k=k+1;
    end
end

warning('off','MATLAB:singularMatrix')

liveExperiment = LiveExperiment(Prefix);
if UseGoodSpots
outpath = [liveExperiment.resultsFolder, filesep, 'StoreSpotInfo.mat'];
load(outpath);%, 'GoodSpots');
Spots = GoodSpots;
% Spots = BestSpots;
else
    Spots = getSpots(liveExperiment);
end

% Add New Variables
Bgd_path = [liveExperiment.resultsFolder, filesep, 'BackgroundEstimateEllipses.mat'];
load(Bgd_path);

xDim = liveExperiment.xDim;
yDim = liveExperiment.yDim;
NumFrames = length(Spots);
TifList = dir([liveExperiment.preFolder, filesep, '*ch01.tif']);
ProbList = dir([liveExperiment.procFolder,'dogs',  filesep, 'prob*']);
MovieCells = cell(1, NumFrames);
ProbCells = cell(1, NumFrames);
NumSlices = 0;
for i = 1:NumFrames
    MovieCells{i} = imreadStack([TifList(i).folder, filesep,  TifList(i).name]);
    ProbCells{i} = imreadStack([ProbList(i).folder, filesep,  ProbList(i).name]);
    NumSlices = NumSlices + size(MovieCells{i}, 3);
end
spot_diameter = 1.5; % microns
crop_diameter = .5;
pixelsize = liveExperiment.pixelSize_um;
spot_pixel_diameter = spot_diameter/pixelsize;
crop_pixel_diameter = crop_diameter/pixelsize;
new_snippet_size = round((spot_pixel_diameter - 1)/2);
crop_snippet_size = max([round((crop_pixel_diameter-1)/2), 2]);
%%
%new_snippet_size = 3;
widthGuess = .05/pixelsize;
lsqOptions=optimset('Display','none');
snippet_size = Spots(1).Fits(1).snippet_size;
% persistent y x
% if isempty(y), 
[y_new,x_new] =meshgrid(1:2*new_snippet_size+1, 1:2*new_snippet_size+1);% end
y_new = single(y_new);
x_new = single(x_new);
dist_mat = sqrt(((y_new-single(new_snippet_size)-1).^2 +(x_new-single(new_snippet_size)-1).^2));
keep_index = dist_mat <= single(new_snippet_size);
crop_dist_mat = sqrt(((y_new-single(crop_snippet_size)-1).^2 +(x_new-single(crop_snippet_size)-1).^2));
crop_index = crop_dist_mat <= single(crop_snippet_size);
[y_old,x_old] = meshgrid(1:2*snippet_size+1, 1:2*snippet_size+1);% end
y_old = single(y_old);
x_old = single(x_old);
y_sub = y_old(snippet_size-new_snippet_size+1:snippet_size-new_snippet_size+2*new_snippet_size+1,...
    snippet_size-new_snippet_size+1:snippet_size-new_snippet_size+2*new_snippet_size+1);
x_sub = x_old(snippet_size-new_snippet_size+1:snippet_size-new_snippet_size+2*new_snippet_size+1,...
    snippet_size-new_snippet_size+1:snippet_size-new_snippet_size+2*new_snippet_size+1);
% [y,x] = meshgrid(1:xDim, 1:yDim);



%for FrameIndex = 1:NumFrames
FrameIndex = 1;

% New Vars
[Spots(FrameIndex).Fits.GaussianIntensityArea] = deal([]);
[Spots(FrameIndex).Fits.GaussianIntensityBackground] = deal([]);
[Spots(FrameIndex).Fits.GaussianIntensityZBackground] = deal([]);
[Spots(FrameIndex).Fits.GaussianIntensityCorrected] = deal([]);
[Spots(FrameIndex).Fits.GaussianIntensityCorrectedArea] = deal([]);
[Spots(FrameIndex).Fits.GaussianIntensityCorrectedZBackground] = deal([]);
[Spots(FrameIndex).Fits.GaussianIntensityCorrectedBackground] = deal([]);
[Spots(FrameIndex).Fits.GaussianIntensityCropped] = deal([]);
[Spots(FrameIndex).Fits.GaussianIntensityCroppedArea] = deal([]);
[Spots(FrameIndex).Fits.GaussianIntensityCroppedZBackground] = deal([]);
[Spots(FrameIndex).Fits.GaussianIntensityCroppedBackground] = deal([]);
[Spots(FrameIndex).Fits.CroppedGaussianIntensitySmallSnip] = deal([]);
[Spots(FrameIndex).Fits.GaussianIntensitySmallSnip] = deal([]);
[Spots(FrameIndex).Fits.GaussianZSmallSnip] = deal([]);
[Spots(FrameIndex).Fits.GaussianAsSmallSnip] = deal([]);
[Spots(FrameIndex).Fits.GaussianWidthSmallSnip] = deal([]);
[Spots(FrameIndex).Fits.GaussianSigmasSmallSnip] = deal([]);
[Spots(FrameIndex).Fits.GaussianSigmaXsSmallSnip] = deal([]);
[Spots(FrameIndex).Fits.GaussianSigmaYsSmallSnip] = deal([]);
[Spots(FrameIndex).Fits.GaussianSigmaXYsSmallSnip] = deal([]);
[Spots(FrameIndex).Fits.GaussianWidthXSmallSnip] = deal([]);
[Spots(FrameIndex).Fits.GaussianWidthYSmallSnip] = deal([]);
[Spots(FrameIndex).Fits.GaussianRhoSmallSnip] = deal([]);
[Spots(FrameIndex).Fits.GaussianPeakSmallSnip] = deal([]);
[Spots(FrameIndex).Fits.gaussParamsSmallSnip] = deal({});
[Spots(FrameIndex).Fits.GaussianErrorSmallSnip] = deal({});
[Spots(FrameIndex).Fits.GaussianResidualsSmallSnip] = deal({});
[Spots(FrameIndex).Fits.GaussianFitValuesSmallSnip] = deal({});
[Spots(FrameIndex).Fits.GaussianAreaSmallSnip] = deal([]);
[Spots(FrameIndex).Fits.GaussianZBackgroundSmallSnip] = deal([]);
[Spots(FrameIndex).Fits.GaussianBackgroundSmallSnip] = deal([]);
[Spots(FrameIndex).Fits.GaussianKernelValues] = deal({});
[Spots(FrameIndex).Fits.GaussianKernelIntensity] = deal([]);
[Spots(FrameIndex).Fits.GaussianKernelZ] = deal([]);
[Spots(FrameIndex).Fits.GaussianKernelArea] = deal([]);
[Spots(FrameIndex).Fits.GaussianKernelZBackground] = deal([]);
[Spots(FrameIndex).Fits.GaussianKernelBackground] = deal([]);
[Spots(FrameIndex).Fits.IntegratedGaussianKernelIntensity] = deal([]);
[Spots(FrameIndex).Fits.IntegratedGaussianKernelZ] = deal([]);
[Spots(FrameIndex).Fits.IntegratedGaussianKernelArea] = deal([]);
[Spots(FrameIndex).Fits.IntegratedGaussianKernelZBackground] = deal([]);
[Spots(FrameIndex).Fits.IntegratedGaussianKernelBackground] = deal([]);


NumSpots = length(Spots(FrameIndex).Fits);

%%


for SpotIndex = 1:NumSpots%1:NumSpots
    %disp(num2str(SpotIndex))
snippet_size = Spots(FrameIndex).Fits(SpotIndex).snippet_size;
NumZ = length(Spots(FrameIndex).Fits(SpotIndex).z);



for zIndex = 1:NumZ
PixelList = Spots(FrameIndex).Fits(SpotIndex).PixelList{zIndex};
z = Spots(FrameIndex).Fits(SpotIndex).z(zIndex);

xmin = min(PixelList(:,1));
ymin = min(PixelList(:,2));
xmax = max(PixelList(:,1));
ymax = max(PixelList(:,2));
centroid = uint16(round([mean(PixelList(:,1)), mean(PixelList(:,2))]));


snippet_old = MovieCells{FrameIndex}(centroid(2)-uint16(snippet_size):centroid(2)+uint16(snippet_size),...
    centroid(1)-uint16(snippet_size):centroid(1)+uint16(snippet_size), Spots(FrameIndex).Fits(SpotIndex).z(zIndex));
snippet_new = MovieCells{FrameIndex}(centroid(2)-uint16(new_snippet_size):centroid(2)+uint16(new_snippet_size),...
    centroid(1)-uint16(new_snippet_size):centroid(1)+uint16(new_snippet_size), Spots(FrameIndex).Fits(SpotIndex).z(zIndex));
gaussParams = Spots(FrameIndex).Fits(SpotIndex).gaussParams{zIndex};
%%
singleGaussian_large = gaussianForSpot(y_old, x_old, snippet_old);

Spots(FrameIndex).Fits(SpotIndex).GaussianIntensityCorrected(zIndex) = single(sum(sum(singleGaussian_large(gaussParams) +...
    snippet_old - gaussParams(7)-gaussParams(8)*double(x_old)-gaussParams(9)*double(y_old))));

Spots(FrameIndex).Fits(SpotIndex).GaussianIntensityArea = size(snippet_old,1)*size(snippet_old,2);
Spots(FrameIndex).Fits(SpotIndex).GaussianIntensityBackground = double(size(snippet_old,1)*size(snippet_old,2))*TotalBackgroundAverage;
Spots(FrameIndex).Fits(SpotIndex).GaussianIntensityZBackground(zIndex) = ...
    double(size(snippet_old,1)*size(snippet_old,2))*BackgroundZAverages(Spots(FrameIndex).Fits(SpotIndex).z(zIndex));

Spots(FrameIndex).Fits(SpotIndex).GaussianIntensityCorrectedArea = size(snippet_old,1)*size(snippet_old,2);
Spots(FrameIndex).Fits(SpotIndex).GaussianIntensityCorrectedBackground = double(size(snippet_old,1)*size(snippet_old,2))*TotalBackgroundAverage;
Spots(FrameIndex).Fits(SpotIndex).GaussianIntensityCorrectedZBackground(zIndex) = ...
    double(size(snippet_old,1)*size(snippet_old,2))*BackgroundZAverages(Spots(FrameIndex).Fits(SpotIndex).z(zIndex));


singleGaussian_sub = gaussianForSpot(y_sub(keep_index), x_sub(keep_index), snippet_new(keep_index));

Spots(FrameIndex).Fits(SpotIndex).GaussianIntensityCropped(zIndex) = single(sum(sum(singleGaussian_sub(gaussParams) +...
    snippet_new(keep_index) - gaussParams(7)-gaussParams(8)*double(x_sub(keep_index))-gaussParams(9)*double(y_sub(keep_index)))));

Spots(FrameIndex).Fits(SpotIndex).GaussianIntensityCroppedArea = length(y_sub(keep_index));
Spots(FrameIndex).Fits(SpotIndex).GaussianIntensityCroppedBackground =...
    double(length(y_sub(keep_index)))*TotalBackgroundAverage;
Spots(FrameIndex).Fits(SpotIndex).GaussianIntensityCroppedZBackground(zIndex) = ...
    double(length(y_sub(keep_index)))*BackgroundZAverages(Spots(FrameIndex).Fits(SpotIndex).z(zIndex));


snippet_new = double(snippet_new);
singleGaussian = gaussianForSpot(double(y_new(keep_index)), double(x_new(keep_index)), snippet_new(keep_index));
singleGaussianCropped = gaussianForSpot(double(y_new(crop_index)), double(x_new(crop_index)), snippet_new(crop_index));

%%

m = max(snippet_new(keep_index));

%let's cache this for efficiency
%realized caching does weird stuff in parpools
%gonna disable that til i figure out why. 

% [y,x] = meshgrid(1:xDim, 1:yDim);

med = min(snippet_new(keep_index));
%fits: [amplitude, x position, x width, y position, y width, offset, angle]


%Define some more initial parameters for fitting
initial_parameters = [m, round((2*new_snippet_size+1)/2), round((2*new_snippet_size+1)/2), ...
    0, widthGuess, widthGuess,med, 0, 0];

%fits: [amplitude, x position, x width, y position, y width, offset, angle]
lb_offset = 1/10; %this is empirical. corresponds to a weak background of 1 pixel per ten having a value of 1.
% lb = [0, 0, 0, 0, 0,lb_offset, 0, -med/2, -med/2];
% ub = [max(snippet(:))*1.5, size(snippet, 2), size(snippet, 1), size(snippet, 2), size(snippet, 1), max(snippet(:)), 2*pi, med/2, med/2];
%
% [single_fit, ~, residual, ~, ~, ~, ~] = lsqnonlin(singleGaussian, ...
%     initial_parameters,lb,ub, lsqOptions);

% @(A, x0, y0, rho, sigma_x, sigma_y, offset, offset_x, offset_y)
lb = [0, 0, 0, -1, 0, 0,BackgroundZAverages(z)*.9, -median(snippet_new(:))/2, -median(snippet_new(:))/2];
ub = [m*1.5, (2*new_snippet_size+1), (2*new_snippet_size+1), 1,...
    new_snippet_size*2+1, new_snippet_size*2+1, m, median(snippet_new(:))/2, median(snippet_new(:))/2];

%let's cache this for efficiency

try

[single_fit, ~, residual, ~, ~, ~, jacobian] = lsqnonlin(singleGaussian, ...
    initial_parameters,lb,ub, lsqOptions);

%quality control. probably just noisy background so amplitude
%should be close to offset in this case.
%
% if .1 > single_fit(3) | .1 > single_fit(5)...
%         | single_fit(3) > size(snippet)/2 | single_fit(5) > size(snippet)/2
%     initial_parameters = [lb_offset, round(length(snippet)/2), widthGuess, round(length(snippet)/2), ...
%         widthGuess,med, 0, 0, 0];
%     lb = [0, 0, 0, 0, 0,lb_offset, 0, -med/2, -med/2];
%     ub = [median(snippet(:))/2, size(snippet, 1), size(snippet, 1), size(snippet, 2), size(snippet, 2), max(snippet(:)), 2*pi, med/2, med/2];
%
%     [single_fit, res1, residual, exitflag, output, lambda, jacobian] = lsqnonlin(singleGaussian, ...
%         initial_parameters,lb,ub, lsqOptions);
% end

confidence_intervals = single(nlparci(single_fit,residual,'jacobian',jacobian));
errors = zeros(1, length(single_fit));
for i = 1:length(confidence_intervals)
    errors(i) = (1/2)* abs(...
        (abs( confidence_intervals(i, 1) )...
        - abs( confidence_intervals(i, 2)) )...
        );
end

relative_errors = single(abs(errors./single_fit));

residual = single(residual);

fits = single(single_fit);
fitvalues = singleGaussian(single_fit) + snippet_new(keep_index) - ...
    single_fit(7)-single_fit(8)*x_new((keep_index))-single_fit(9)*y_new((keep_index));
reconstructed_residual = zeros(2*new_snippet_size+1, 2*new_snippet_size+1, 'single');
reconstructed_fv =zeros(2*new_snippet_size+1, 2*new_snippet_size+1, 'single');
y_kept = y_new(keep_index);
x_kept = x_new(keep_index);
for k = 1:length(residual)
    reconstructed_residual(x_kept(k), y_kept(k)) = residual(k);
    reconstructed_fv(x_kept(k), y_kept(k)) = fitvalues(k);
end
GaussianIntensity = single(sum(sum(fitvalues)));
crop_fitvalues = singleGaussianCropped(single_fit) + snippet_new(crop_index) - ...
    single_fit(7)-single_fit(8)*x_new(crop_index)-single_fit(9)*y_new(crop_index);
CroppedGaussianIntensity = single(sum(sum(crop_fitvalues)));

%%
% Spots(FrameIndex).Fits(SpotIndex).gaussParams{zIndex}
% single_fit
% save_vals = [round(single_fit(5),2) round(single_fit(6),2) round(GaussianIntensity)]
% FigAx = cell(1, 5);
% close all
% TempFigure = figure(1);
% 
% set(TempFigure,'units', 'normalized', 'position',[0.05, 0.4, 0.9, 0.3]);
% set(gcf,'color','w');
% clim = [0, 15];
% FigAx{1} = subplot(1,5,1);
% imagesc(FigAx{1} , snippet_old,clim);   
% colormap(FigAx{1},'gray');
% FigAx{2} = subplot(1,5,2); 
% imagesc(FigAx{2} , singleGaussian_large(gaussParams) +snippet_old - gaussParams(7),clim); 
% colormap(FigAx{2},'gray');
% FigAx{3} = subplot(1,5,3);
% imagesc(FigAx{3} , singleGaussian_large(gaussParams) +snippet_old - gaussParams(7)-gaussParams(8)*double(x_old)-gaussParams(9)*double(y_old),clim);  
% colormap(FigAx{3},'gray');
% FigAx{4} = subplot(1,5,4);
% imagesc(FigAx{4} , snippet_new,clim);  
% colormap(FigAx{4},'gray');
% FigAx{5} = subplot(1,5,5);
% imagesc(FigAx{5} , reconstructed_fv,clim);  
% 
% colormap(FigAx{5},'gray');
% 
%     
% 
% plot_title = ['SpotIndex: ', num2str(SpotIndex),', z: ',...
%     num2str(Spots(FrameIndex).Fits(SpotIndex).z(zIndex))];
% plot_title = [plot_title, ', Old Gaussian Intensity: ',...
%     num2str(Spots(FrameIndex).Fits(SpotIndex).GaussianIntensity(zIndex))];
% OldGaussRecalc = single(sum(sum(singleGaussian_large(gaussParams) +snippet_old - gaussParams(7)-gaussParams(8)*double(x_old)-gaussParams(9)*double(y_old))));
% plot_title = [plot_title, ', Recalculated Old Gaussian Intensity: ', num2str(OldGaussRecalc)];
% plot_title = [plot_title, ', New Gaussian Intensity: ', num2str(GaussianIntensity)];
% 
% sgtitle(plot_title)

%%
Spots(FrameIndex).Fits(SpotIndex).CroppedGaussianIntensitySmallSnip(zIndex) = CroppedGaussianIntensity;
Spots(FrameIndex).Fits(SpotIndex).GaussianIntensitySmallSnip(zIndex) = GaussianIntensity;
Spots(FrameIndex).Fits(SpotIndex).GaussianAreaSmallSnip = length(y_new(keep_index));
Spots(FrameIndex).Fits(SpotIndex).GaussianBackgroundSmallSnip =...
    double(length(y_new(keep_index)))*TotalBackgroundAverage;
Spots(FrameIndex).Fits(SpotIndex).GaussianZBackgroundSmallSnip(zIndex) = ...
    double(length(y_new(keep_index)))*BackgroundZAverages(Spots(FrameIndex).Fits(SpotIndex).z(zIndex));

Spots(FrameIndex).Fits(SpotIndex).GaussianSigmasSmallSnip(zIndex) = mean([single_fit(5) single_fit(6)]);
Spots(FrameIndex).Fits(SpotIndex).GaussianSigmaXsSmallSnip(zIndex) =single_fit(5);
Spots(FrameIndex).Fits(SpotIndex).GaussianSigmaYsSmallSnip(zIndex) =single_fit(6);
Spots(FrameIndex).Fits(SpotIndex).GaussianSigmaXYsSmallSnip(zIndex) =single_fit(4);
Spots(FrameIndex).Fits(SpotIndex).GaussianAsSmallSnip(zIndex) =single_fit(1);
Spots(FrameIndex).Fits(SpotIndex).gaussParamsSmallSnip{zIndex} = single_fit;
Spots(FrameIndex).Fits(SpotIndex).GaussianErrorSmallSnip{zIndex} = errors;
Spots(FrameIndex).Fits(SpotIndex).GaussianResidualsSmallSnip{zIndex} = reconstructed_residual;
Spots(FrameIndex).Fits(SpotIndex).GaussianFitValuesSmallSnip{zIndex} = reconstructed_fv;
catch
Spots(FrameIndex).Fits(SpotIndex).CroppedGaussianIntensitySmallSnip(zIndex) = NaN;
Spots(FrameIndex).Fits(SpotIndex).GaussianIntensitySmallSnip(zIndex) = NaN;
Spots(FrameIndex).Fits(SpotIndex).GaussianAreaSmallSnip = 0;
Spots(FrameIndex).Fits(SpotIndex).GaussianBackgroundSmallSnip =...
    double(length(y_new(keep_index)))*TotalBackgroundAverage;
Spots(FrameIndex).Fits(SpotIndex).GaussianZBackgroundSmallSnip(zIndex) = ...
    double(length(y_new(keep_index)))*BackgroundZAverages(Spots(FrameIndex).Fits(SpotIndex).z(zIndex));

Spots(FrameIndex).Fits(SpotIndex).GaussianSigmasSmallSnip(zIndex) = NaN;
Spots(FrameIndex).Fits(SpotIndex).GaussianSigmaXsSmallSnip(zIndex) =NaN;
Spots(FrameIndex).Fits(SpotIndex).GaussianSigmaYsSmallSnip(zIndex) =NaN;
Spots(FrameIndex).Fits(SpotIndex).GaussianSigmaXYsSmallSnip(zIndex) =NaN;
Spots(FrameIndex).Fits(SpotIndex).GaussianAsSmallSnip(zIndex) =NaN;
Spots(FrameIndex).Fits(SpotIndex).gaussParamsSmallSnip{zIndex} = [NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN];
Spots(FrameIndex).Fits(SpotIndex).GaussianErrorSmallSnip{zIndex} = [NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN];
Spots(FrameIndex).Fits(SpotIndex).GaussianResidualsSmallSnip{zIndex} = NaN(size(snippet_new));
Spots(FrameIndex).Fits(SpotIndex).GaussianFitValuesSmallSnip{zIndex} =  NaN(size(snippet_new));  
end
end
[maxValue, brightZ] = nanmax(Spots(FrameIndex).Fits(SpotIndex).GaussianIntensitySmallSnip);
Spots(FrameIndex).Fits(SpotIndex).GaussianZSmallSnip = Spots(FrameIndex).Fits(SpotIndex).z( brightZ);
Spots(FrameIndex).Fits(SpotIndex).GaussianWidthSmallSnip = Spots(FrameIndex).Fits(SpotIndex).GaussianSigmasSmallSnip( brightZ);
Spots(FrameIndex).Fits(SpotIndex).GaussianWidthXSmallSnip = Spots(FrameIndex).Fits(SpotIndex).GaussianSigmaXsSmallSnip( brightZ);
Spots(FrameIndex).Fits(SpotIndex).GaussianWidthYSmallSnip = Spots(FrameIndex).Fits(SpotIndex).GaussianSigmaYsSmallSnip( brightZ);
Spots(FrameIndex).Fits(SpotIndex).GaussianRhoSmallSnip = Spots(FrameIndex).Fits(SpotIndex).GaussianSigmaXYsSmallSnip( brightZ);
Spots(FrameIndex).Fits(SpotIndex).GaussianPeakSmallSnip = Spots(FrameIndex).Fits(SpotIndex).GaussianAsSmallSnip( brightZ);
end

%%
mean_rho = mean([Spots(FrameIndex).Fits(:).GaussianRhoSmallSnip]);
mean_sigx = mean([Spots(FrameIndex).Fits(:).GaussianWidthXSmallSnip]);
mean_sigy = mean([Spots(FrameIndex).Fits(:).GaussianWidthYSmallSnip]);
filter_size = new_snippet_size;
middle = (filter_size+1)/2;
custom_filter = zeros(filter_size,filter_size,'double');
for i = 1:filter_size
    for j = 1:filter_size
        custom_filter(i,j) = exp(-1/(2*(1-mean_rho)^2)*(((j-middle)^2)/(2*mean_sigx^2)+...
            ((i-middle)^2)/(2*mean_sigy^2)-(2*mean_rho*(j-middle)*(i-middle))/(mean_sigx*mean_sigy)));
    end
end

custom_filter = custom_filter/sum(sum(custom_filter));
sum_filter = ones(filter_size,filter_size,'double');

%%
for SpotIndex =1:NumSpots
snippet_size = Spots(FrameIndex).Fits(SpotIndex).snippet_size;
NumZ = length(Spots(FrameIndex).Fits(SpotIndex).z);


for zIndex = 1:NumZ
PixelList = Spots(FrameIndex).Fits(SpotIndex).PixelList{zIndex};
z = Spots(FrameIndex).Fits(SpotIndex).z(zIndex);

xmin = min(PixelList(:,1));
ymin = min(PixelList(:,2));
xmax = max(PixelList(:,1));
ymax = max(PixelList(:,2));
centroid = [mean(PixelList(:,1)), mean(PixelList(:,2))];
new_gaussParams =  Spots(FrameIndex).Fits(SpotIndex).gaussParamsSmallSnip{zIndex};
fit_center = [new_gaussParams(2), new_gaussParams(3)];
diff_center=fit_center-[new_snippet_size+1, new_snippet_size+1];
try
centroid = uint16(round(double(centroid) + diff_center));
snippet_new = MovieCells{FrameIndex}(centroid(2)-uint16(new_snippet_size):centroid(2)+uint16(new_snippet_size),...
    centroid(1)-uint16(new_snippet_size):centroid(1)+uint16(new_snippet_size), Spots(FrameIndex).Fits(SpotIndex).z(zIndex));
catch
    centroid = [mean(PixelList(:,1)), mean(PixelList(:,2))];
  snippet_new = MovieCells{FrameIndex}(centroid(2)-uint16(new_snippet_size):centroid(2)+uint16(new_snippet_size),...
    centroid(1)-uint16(new_snippet_size):centroid(1)+uint16(new_snippet_size), Spots(FrameIndex).Fits(SpotIndex).z(zIndex));
end  
filtered_snippet = imfilter(snippet_new, custom_filter);
summed_snippet = imfilter(filtered_snippet, sum_filter);
Spots(FrameIndex).Fits(SpotIndex).GaussianKernelValues{zIndex} = filtered_snippet;
Spots(FrameIndex).Fits(SpotIndex).GaussianKernelIntensity(zIndex) = max(max(filtered_snippet));
Spots(FrameIndex).Fits(SpotIndex).IntegratedGaussianKernelIntensity(zIndex) = max(max(summed_snippet));

Spots(FrameIndex).Fits(SpotIndex).IntegratedGaussianKernelArea = double(new_snippet_size)^2;
Spots(FrameIndex).Fits(SpotIndex).IntegratedGaussianKernelBackground =...
    double(new_snippet_size)^2*TotalBackgroundAverage;
Spots(FrameIndex).Fits(SpotIndex).IntegratedGaussianKernelZBackground(zIndex) = ...
    double(new_snippet_size)^2*BackgroundZAverages(Spots(FrameIndex).Fits(SpotIndex).z(zIndex));

Spots(FrameIndex).Fits(SpotIndex).GaussianKernelArea = 1.0;
Spots(FrameIndex).Fits(SpotIndex).GaussianKernelBackground =...
    TotalBackgroundAverage;
Spots(FrameIndex).Fits(SpotIndex).GaussianKernelZBackground(zIndex) = ...
    BackgroundZAverages(Spots(FrameIndex).Fits(SpotIndex).z(zIndex));

end
[maxValue, brightZ] = max(Spots(FrameIndex).Fits(SpotIndex).GaussianKernelIntensity);
Spots(FrameIndex).Fits(SpotIndex).GaussianKernelZ = Spots(FrameIndex).Fits(SpotIndex).z( brightZ);
[maxValue, brightZ] = max(Spots(FrameIndex).Fits(SpotIndex).IntegratedGaussianKernelIntensity);
Spots(FrameIndex).Fits(SpotIndex).IntegratedGaussianKernelZ = Spots(FrameIndex).Fits(SpotIndex).z( brightZ);


end
%%
if UseGoodSpots
GoodSpots = Spots;
outpath = [liveExperiment.resultsFolder, filesep, 'StoreSpotInfo.mat'];
save(outpath, 'GoodSpots', 'BadSpots', 'MultiSpots', 'BestSpots', 'AggregateSpots',...
    'GoodIndices', 'BadIndices', 'MultiIndices', 'BestIndices', 'AggregateIndices');
else
   mkdir(liveExperiment.resultsFolder);
if whos(var2str(Spots)).bytes < 2E9
    save([liveExperiment.resultsFolder,...
        filesep, 'Spots.mat'], 'Spots', '-v6');
else
    save([liveExperiment.resultsFolder,...
        filesep, 'Spots.mat'], 'Spots', '-v7.3', '-nocompression');
end
 
end

