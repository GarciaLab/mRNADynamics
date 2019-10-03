function [GaussParams1, SpotInt1, GaussCI1, IntCI1, GaussParams2,SpotInt2,GaussCI2,IntCI2,offset, offsetCI]...
        = fit3DGaussian2spot(snip3D,PixelSize, varargin)
    % INPUT ARGUMENTS:
    % snip3D: 3D array containing spot to fit. Should contain only one spot
    % PixelSize: size of pixels (in nm)
    % OPTIONS:
    % initial_parameters: 1 x 11 vector containing values from which to
    %                     initialize inferemce
    % ub_vec: 1 x 11 vector specifying upper bounds for inference parameters
    % lb_vec: 1 x 11 vector specifying lower bounds
    
    % RETURNS:
    % GaussFit: 1 x 11 vector of inferred parameter values
    % Parameter identity is as follows: 
    %        (1) amplitude of gaussian (spot 1)    
    %        (2-4) y,x,and z center positions (spot 1)    
    %        (5-7) y,x,and z sigma values (spot 1)  
    %        (8) amplitude of gaussian (spot 1)
    %        (9-11) dy, dx, and dz distances from spot 1 to spot 2
    %        (12-14) y,x,and z sigma values (spot 2)  
    %        (15) inferred background offset    
    
    sigma_guess = 200/PixelSize;
    % define initial parameters
    xDim = size(snip3D,1);
    yDim = size(snip3D,2);
    zDim = size(snip3D,3);
    % initialize parameters
    initial_parameters =[max(snip3D(:)), floor(yDim/2)-1,floor(xDim/2)-1, floor(zDim/2)-1,sigma_guess,sigma_guess...
        ,.1,ceil(yDim/2)+1,ceil(xDim/2)+1, ceil(zDim/2)+1,1,1,prctile(snip3D(:),10)];
    
    % initialize upper and lower parameter bounds
    ub_vec = [3*max(snip3D(:)),yDim,xDim,zDim,3*sigma_guess,3*sigma_guess,...
              3*max(snip3D(:)),yDim,xDim,zDim,3*sigma_guess,3*sigma_guess,...
              Inf];
    lb_vec = [0,1,1,1,.5,.5,0,1,1,1,.1,.1,0];
    
    % check for additional arguments
    for i = 1:(numel(varargin)-1)  
        if ischar(varargin{i})
            eval([varargin{i} '= varargin{i+1};']);        
        end
    end
       
    % define objective function
    dim_vec = [yDim, xDim, zDim];        
    double3DObjective = @(params) simulate3DGaussSymmetric(dim_vec, params(1:6))...
        + simulate3DGaussSymmetric(dim_vec, params(7:12)) ...
        + params(end) - double(snip3D);     
    % attempt to fit
    options.Display = 'off';
    [GaussFit, ~, residual, ~, ~, ~, jacobian] = lsqnonlin(double3DObjective,...
        initial_parameters,lb_vec,ub_vec,options);
    % store parameters
    GaussParams1 = GaussFit(1:6);
    GaussParams2 = GaussFit(7:12);        
    offset = GaussFit(end);
    % estimate confidence intervals
    FitCI = nlparci(GaussFit,residual,'jacobian',jacobian);
    GaussCI1 = FitCI(1:6,:);
    Deltas1 = diff(GaussCI1');
    GaussCI2 = FitCI(7:12,:);
    Deltas2 = diff(GaussCI2');
    offsetCI = FitCI(end,:);
    % calculate spot integrals and propagate uncertainties
    % spot 1
    sXY1 = GaussParams1(5);
    sZ1 = GaussParams1(6);
    amp1 = GaussParams1(1);    
    SpotInt1 = amp1 * (2*pi)^1.5 * sXY1^2 * sZ1;   
    IntErr1 = (2*pi)^1.5 * sqrt((Deltas1(1)*sXY1^2 * sZ1)^2 + (Deltas1(6)*sXY1^2 * amp1)^2 ...
        + (2*Deltas1(5)*amp1 * sXY1^2 * sZ1)^2);    
    IntCI1 =  [SpotInt1 - IntErr1/2, SpotInt1 + IntErr1/2];
    % spot 2
    sXY2 = GaussParams2(5);
    sZ2 = GaussParams2(6);
    amp2 = GaussParams2(1);
    SpotInt2 = amp2 * (2*pi)^1.5 * sXY2^2 * sZ2;   
    IntErr2 = (2*pi)^1.5 * sqrt((Deltas2(1)*sXY2^2 * sZ2)^2 + (Deltas2(6)*sXY2^2 * amp2)^2 ...
        + (2*Deltas2(5)*amp2 * sXY2^2 * sZ2)^2);    
    IntCI2 =  [SpotInt2 - IntErr2/2, SpotInt2 + IntErr2/2];
    % perform error propagation to obtain error bounds on integral  
    
    
    
end