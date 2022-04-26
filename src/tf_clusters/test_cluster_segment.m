close all; clear

preproc_dir = 'P:\DorsalClusters\Data\PreProcessedData\';
exp_name = '2022-03-29-Dl-mNeonGreen-MCP-mCh_snaBAC_settingsTest06_embryo04';
filename = [preproc_dir, exp_name, filesep, exp_name, '_020_ch02.tif'];

tiff_file = tiffreadVolume(filename);
im_stack = tiff_file;

im_3dgauss = imgaussfilt3(im_stack, 1);

snip_x = [575:675];
snip_y = [225:325];
center_slice = 15;
n_plots = [2, 3];
map = hot;
disp_range = [0, 7000];
subplot(n_plots(1),n_plots(2),1), imshow(im_stack(snip_y,snip_x,center_slice-1), DisplayRange=disp_range)%, ...
%                        Colormap=map)
subplot(n_plots(1),n_plots(2),2), imshow(im_stack(snip_y,snip_x,center_slice), DisplayRange=disp_range)%, ...
%                        Colormap=map)
subplot(n_plots(1),n_plots(2),3), imshow(im_stack(snip_y,snip_x,center_slice+1), DisplayRange=disp_range)%, ...
%                        Colormap=map)

subplot(n_plots(1),n_plots(2),4), imshow(im_3dgauss(snip_y,snip_x,center_slice-1), DisplayRange=disp_range)%, ...
%                        Colormap=map)
subplot(n_plots(1),n_plots(2),5), imshow(im_3dgauss(snip_y,snip_x,center_slice), DisplayRange=disp_range)%, ...
%                        Colormap=map)
subplot(n_plots(1),n_plots(2),6), imshow(im_3dgauss(snip_y,snip_x,center_slice+1), DisplayRange=disp_range)%, ...
%                        Colormap=map)