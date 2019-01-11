% Pipeline for Determining Dorsal Gaussian Gradient from embryo
% Paul Marchando
% 12-10-18

function [working_dat_final, fake_dv_comb, fake_dv_check, polyfit_data, polyfit_out,...
    polyfit_maximums, gauss_out, para_values] = DetermineDVGradient(Prefix)

tic
%Get the folders, including the default Dropbox one
[SourcePath, FISHPath, DefaultDropboxFolder, DropboxFolder, MS2CodePath, PreProcPath,...
configValues, movieDatabasePath] = DetermineAllLocalFolders(Prefix);

%Determine division times
%Load the information about the nc from moviedatabase file
[Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
nc9, nc10, nc11, nc12, nc13, nc14, CF] = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder);

% Calls functions that clean DV data, fit each nucleus's z-stack to a
% polynomial to determine fluorescence values, then fit fluorescence values
% for each frame to a Gaussian to determine an embryo's DV gradient. Plots
% of fitted Gaussians and each parameter in the fits are then generated.

load([DropboxFolder '/' Prefix '/' 'CompiledNuclei.mat'])
load([DropboxFolder '/' Prefix '/' Prefix '_lin.mat'])

% Cleans data from CompiledNuclei and schnitzcells.

[working_dat_final] = CleanDVData(schnitzcells,CompiledNuclei);

% Fits polynomials to each nucleus z-stack

[fake_dv_comb, fake_dv_check, polyfit_data, polyfit_out, ...
    polyfit_maximums] = DVpolyfit(working_dat_final,0.5);

% Fits a Gaussian to fluorescence values from each frame.

[gauss_out, para_values] = fitDVGaussian(fake_dv_comb);

% Generates plots.

test_ind = nc14:10:length(gauss_out);

if length(test_ind) > 7
    
    test_ind = test_ind(1:7);
    
end

colors = [0    0.4470    0.7410;
          0.8500    0.3250    0.0980;
          0.9290    0.6940    0.1250;
          0.4940    0.1840    0.5560;
          0.4660    0.6740    0.1880;
          0.3010    0.7450    0.9330;
          0.6350    0.0780    0.1840];
      
hold on
for i = 1:length(test_ind)
    
    if isempty(gauss_out{test_ind(i),1}) == 0
        
        temp_range = find(fake_dv_comb(:,3) == test_ind(i));
        xdata = fake_dv_comb(temp_range(1):temp_range(end),1);
        ydata = fake_dv_comb(temp_range(1):temp_range(end),2);
        scatter(xdata,ydata,[],colors(i,:),'.',...
            'DisplayName',['Frame: ' num2str(test_ind(i))])
        sort_xdata = sort(xdata);
        tempgauss = gauss_out{test_ind(i),1};
        plot(sort_xdata,tempgauss(sort_xdata),'color',colors(i,:),...
            'DisplayName',['Frame ' num2str(test_ind(i)) ' Fit'])
       
    end
end
hold off

legend
title('Dorsal gradient plots following NC14')
xlabel('"Fake" DV position')
ylabel('Nuclear fluorescence (AFU)')

figure
scatter(para_values(:,1),para_values(:,3),'.','k')
title('Center')
axis([ 1 length(gauss_out) -200 200])
xlabel('Frame number')
ylabel('Parameter value')

figure
figure
bar(para_values(:,1),para_values(:,2),'r')
title('Amplitude')
axis([ 1 length(gauss_out) 0 2500])
xlabel('Frame number')
ylabel('Parameter value')

figure
bar(para_values(:,1),para_values(:,4),'b')
title('Spread')
axis([ 1 length(gauss_out) 0 2000])
xlabel('Frame Number')
ylabel('Parameter value')
toc

end