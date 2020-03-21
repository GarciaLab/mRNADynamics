clear
close all
% script to test alternative 3D fitting approaches

% make figure path 
FigPath = 'fit3D_validation_figs/';
mkdir(FigPath);

% define 3D snip dimensions
xDim = 25;
yDim = 25;
zDim = 12;
dim_vec = [yDim,xDim,zDim];
% num test samples
n_sim = 100;
min_sigma = 1;
sigma_scale = 3;
amp_scale = 50;
noise_scale = 10;
dxy_scale = noise_scale / xDim;
dz_scale = noise_scale / zDim;
PixelDims = [0.1,0.1,0.5];
% define positional reference arrays
[x_ref_vol, y_ref_vol, z_ref_vol] = meshgrid(1:dim_vec(1),1:dim_vec(2),1:dim_vec(3));
% define helper functions
cov_mat = @(params) [params(5)^2,    params(8)*params(5)*params(6),    params(9)*params(5)*params(7);
                   params(8)*params(5)*params(6),    params(6)^2,    params(10)*params(6)*params(7);
                   params(9)*params(5)*params(7),    params(10)*params(6)*params(7),    params(7)^2];
gauss_int = @(params) params(1)*(2*pi)^1.5 *sqrt(det(cov_mat(params)));
% initialize cell array to store spots
disp('generating simulated gaussians...')
gauss_cell = cell(1,n_sim);
true_param_mat = NaN(n_sim,14);
true_fluo_vec = NaN(1,n_sim);
for n = 1:n_sim    
    s_flag = 1;
    while s_flag 
        % define Gaussian parameters
        xPos = ceil(xDim/2) + rand(1,1)*sigma_scale;
        yPos = ceil(yDim/2) + rand(1,1)*sigma_scale;
        zPos = zDim/2 + rand(1,1)*sigma_scale;
        sigmaX = min_sigma+rand(1,1)*(sigma_scale-min_sigma);
        sigmaY = min_sigma+rand(1,1)*(sigma_scale-min_sigma);
        sigmaZ = min_sigma+rand(1,1)*(sigma_scale-min_sigma);
        rhoXY = 2*(rand(1,1)-.5);
        rhoXZ = 2*(rand(1,1)-.5);
        rhoYZ = 2*(rand(1,1)-.5);
        gAmp = amp_scale*rand(1,1);

        % check that cov matrix is of proper form
        params = [gAmp,yPos,xPos,zPos,sigmaY,sigmaX,sigmaZ,rhoXY,rhoYZ,rhoXZ];

        cov_mat = @(params) [params(5)^2,    params(8)*params(5)*params(6),    params(9)*params(5)*params(7);
                       params(8)*params(5)*params(6),    params(6)^2,    params(10)*params(6)*params(7);
                       params(9)*params(5)*params(7),    params(10)*params(6)*params(7),    params(7)^2];

        [L,s_flag] = chol(cov_mat(params));
    end
    % generate noisy background
    flat_bkg = rand(yDim,xDim,zDim)*rand()*noise_scale;
    dyn = 2*(rand()-.5)*nanmean(flat_bkg(:))/(yDim+xDim+zDim);
    dxn = 2*(rand()-.5)*nanmean(flat_bkg(:))/(yDim+xDim+zDim);
    dzn = 2*(rand()-.5)*nanmean(flat_bkg(:))/(yDim+xDim+zDim);
    snip3D = flat_bkg + dyn*x_ref_vol + dxn*y_ref_vol + dzn*z_ref_vol;
    true_param_mat(n,:) = [params nanmean(flat_bkg(:)) dyn dxn dzn];
    true_fluo_vec(n) = gauss_int(params);
    gauss_cell{n} = simulate3DGaussianRho(dim_vec, params) + snip3D;
end

% conduct fits
disp('conducting 3D fits...')
rho_inf_mat = NaN(n_sim,numel(params)+4);
np_flag_vec  = NaN(1,n_sim);
inf_fluo_vec  = NaN(1,n_sim);
raw_fluo_vec  = NaN(1,n_sim);
for n = 1:n_sim
    snip3D = gauss_cell{n};
    
    % rho first
    [rho_inf_mat(n,:), ~, inf_fluo_vec(n), ~, raw_fluo_vec(n) ,np_flag_vec(n)] = fit3DGaussianRho(snip3D,PixelDims);
%     % now cholesky method
%     [chol_inf_cell{n}, FitDeltas, GaussIntegral, GaussIntegralSE, GaussIntegralRaw] = ...
%                                                     fit3DGaussianCholesky(snip3D,PixelDims);
end                                              
    
    
% Check accuracy
disp('done.')
close all
% amplitude
amp_fig = figure;
cmap = flipud(brewermap([],'Spectral'));
colormap(cmap);
scatter(true_param_mat(:,1),rho_inf_mat(:,1),'filled')
xlabel('true amplitude')
ylabel('inferred amplitude')
set(gca,'Fontsize',14)
grid on
saveas(amp_fig,[FigPath 'amplitude.png'])

% background
bkg_fig = figure;
colormap(cmap);
scatter(true_param_mat(:,11),rho_inf_mat(:,11),'filled')
xlabel('true background')
ylabel('inferred background')
set(gca,'Fontsize',14)
grid on
saveas(bkg_fig,[FigPath 'bkg_flat.png'])


% background grad
bkg_grad = figure;
colormap(cmap);
hold on
scatter(true_param_mat(:,12),rho_inf_mat(:,12),'filled')
scatter(true_param_mat(:,13),rho_inf_mat(:,13),'filled')
scatter(true_param_mat(:,14),rho_inf_mat(:,14),'filled')
xlabel('true background gradient')
ylabel('inferred background gradient')
set(gca,'Fontsize',14)
legend('dy','dx','dz','Location','northwest')
grid on
saveas(bkg_grad,[FigPath 'bkg_grad.png'])

% position
pos_fig = figure;
colormap(cmap);
hold on
scatter(true_param_mat(:,2),rho_inf_mat(:,2),'filled')
scatter(true_param_mat(:,3),rho_inf_mat(:,3),'filled')
scatter(true_param_mat(:,4),rho_inf_mat(:,4),'filled')
xlabel('true position')
ylabel('inferred position')
legend('y','x','z','Location','northwest')
grid on
set(gca,'Fontsize',14)
saveas(pos_fig,[FigPath 'position.png'])

% sigmas
sigma_fig = figure;
colormap(cmap);
hold on
scatter(true_param_mat(:,5),rho_inf_mat(:,5),'filled')
scatter(true_param_mat(:,6),rho_inf_mat(:,6),'filled')
scatter(true_param_mat(:,7),rho_inf_mat(:,7),'filled')
xlabel('true size')
ylabel('inferred size')
legend('\sigma_y','\sigma_x','\sigma_z','Location','northwest')
set(gca,'Fontsize',14)
grid on
saveas(sigma_fig,[FigPath 'sigma.png'])

% rhos
rho_fig = figure;
colormap(cmap);
hold on
scatter(true_param_mat(:,8),rho_inf_mat(:,8),'filled')
scatter(true_param_mat(:,9),rho_inf_mat(:,9),'filled')
scatter(true_param_mat(:,10),rho_inf_mat(:,10),'filled')
xlabel('true cross correlation')
ylabel('inferred cross correlation')
legend('\rho_{xy}','\rho_{yz}','\rho_{xz}','Location','northwest')
set(gca,'Fontsize',14)
grid on
saveas(rho_fig,[FigPath 'rho.png'])
% Examine accuracy of raw and inferred spot fluo
% rhos
fluo_inf_fig = figure;
colormap(cmap);
hold on
scatter(true_fluo_vec,inf_fluo_vec,50,true_param_mat(:,end),'filled')
xlabel('true spot intensity')
ylabel('inferred spot intensity (integral)')
h = colorbar;
ylabel(h,'noise (au per pixel)')
set(gca,'Fontsize',14)
% legend('\rho_{xy}','\rho_{yz}','\rho_{xz}','Location','northwest')
grid on
saveas(fluo_inf_fig,[FigPath 'fluo_inf.png'])

fluo_raw_fig = figure;
colormap(cmap);
hold on
scatter(true_fluo_vec,raw_fluo_vec,50,true_param_mat(:,end),'filled')
xlabel('true spot intensity')
ylabel('inferred spot intensity (raw)')
h = colorbar;
ylabel(h,'noise (au per pixel)')
set(gca,'Fontsize',14)
% legend('\rho_{xy}','\rho_{yz}','\rho_{xz}','Location','northwest')
grid on
saveas(fluo_raw_fig,[FigPath 'fluo_raw.png'])