function testFilters(im)

supportedFilters = {'Gaussian_blur', 'Identity', 'Anisotropic_diffusion', 'bilateral',...
    'Edges', 'Sobel', 'Difference_of_Gaussian', 'Laplacian', 'Mean', 'Structure_smallest', 'Structure_largest',...
    'Median',  'Maximum', 'Minimum', 'Std', 'Variance', 'Hessian_smallest', 'Hessian_largest', 'Entropy', 'Range'};

s = 16;
% s = {16, 32};
    
fim = filterImage(im, supportedFilters{8}, s);

figure(1); imshowpair(im, fim, 'montage');

end