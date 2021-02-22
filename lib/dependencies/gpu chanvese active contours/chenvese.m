%==========================================================================
%
%   Active contour with Chen-Vese Method 
%   for image segementation
%
%   Implemented by Yue Wu (yue.wu@tufts.edu)
%   Tufts University
%   Feb 2009
%   http://sites.google.com/site/rexstribeofimageprocessing/
% 
%   all rights reserved 
%   Last update 02/26/2009
%--------------------------------------------------------------------------
%   Usage of varibles:
%   input: 
%       I           = any gray/double/RGB input image
%       mask        = initial mask, either customerlized or built-in
%       num_iter    = total number of iterations
%       mu          = weight of length term
%       method      = submethods pick from ('chen','vector','multiphase')
%
%   Types of built-in mask functions
%       'small'     = create a small circular mask
%       'medium'    = create a medium circular mask
%       'large'     = create a large circular mask
%       'whole'     = create a mask with holes around
%       'whole+small' = create a two layer mask with one layer small
%                       circular mask and the other layer with holes around
%                       (only work for method 'multiphase')
%   Types of methods
%       'chen'      = general CV method
%       'vector'    = CV method for vector image
%       'multiphase'= CV method for multiphase (2 phases applied here)
%
%   output: 
%       phi0        = updated level set function 
%
%--------------------------------------------------------------------------
%
% Description: This code implements the paper: "Active Contours Without
% Edges" by Chan and Vese for method 'chen', the paper:"Active Contours Without
% Edges for vector image" by Chan and Vese for method 'vector', and the paper
% "A Multiphase Level Set Framework for Image Segmentation Using the 
% Mumford and Shah Model" by Chan and Vese. 
%
%--------------------------------------------------------------------------
% Deomo: Please see HELP file for details
%==========================================================================

function seg = chenvese(I,mask,num_iter,mu,method)

%%
%-- Default settings
%   length term mu = 0.2 and default method = 'chan'
  if(~exist('mu','var')) 
    mu=0.2; 
  end
  
  if(~exist('method','var')) 
    method = 'chan'; 
  end

  %   auto mask settings
  if ischar(mask)
    switch lower (mask)
     case 'small'
      mask = maskcircle2(I,'small');
     case 'medium'
      mask = maskcircle2(I,'medium');
     case 'large'
      mask = maskcircle2(I,'large');              
     case 'whole'
      mask = maskcircle2(I,'whole'); 
     case 'whole+small'
      m1 = maskcircle2(I,'whole');
      m2 = maskcircle2(I,'small');
      mask = zeros(size(I,1),size(I,2),2);
      mask(:,:,1) = m1(:,:,1);
      mask(:,:,2) = m2(:,:,2);
     otherwise
      error('unrecognized mask shape name (MASK).');
    end
  else
    if size(mask,1)>size(I,1) || size(mask,2)>size(I,2)
      error('dimensions of mask unmathch those of the image.')
    end
    switch lower(method)
     case 'multiphase'
      if  (size(mask,3) == 1)  
        error('multiphase requires two masks but only gets one.')
      end
    end
    
  end       

  % accelereyes: changed hard-coded type-casts
  switch lower(method)
   case 'chan'
    if size(I,3)== 3
      P = rgb2gray(smartcast(I,'uint8'));
    elseif size(I,3) == 2
      P = smartcast(0.5.*(I(:,:,1)+I(:,:,2)),'double');
    else
      P = smartcast(I,'double');
    end
    layer = 1;
   otherwise
    error('!invalid method')
  end
  %-- End Initializations on input image I and mask

  %%
  %--   Core function
  switch lower(method)
   case {'chan','vector'}
    %-- SDF
    %   Get the distance map of the initial mask
    
    mask = mask(:,:,1);
    phi0 = bwdist(mask)-bwdist(1-mask)+im2double(mask)-.5; 
    %   initial force, set to eps to avoid division by zeros
    force = eps; 
    %-- End Initialization

    %-- Display settings
%     figure();
%     subplot(2,2,1); imshow(uint8(I)); title('Input Image');
%     subplot(2,2,2); contour(flipud(phi0), [0 0], 'LineWidth',1); title('initial contour');
%     subplot(2,2,3); title('Segmentation');
    %-- End Display original image and mask

    %-- Main loop
    for n=1:num_iter
      inidx = find(phi0>=0); % frontground index
      outidx = find(phi0<0); % background index
      force_image = 0; % initial image force for each layer 
      for i=1:layer
        L = im2double(P(:,:,i)); % get one image component
        c1 = sum(sum(L.*Heaviside(phi0)))/(length(inidx)+eps); % average inside of Phi0
        c2 = sum(sum(L.*(1-Heaviside(phi0))))/(length(outidx)+eps); % verage outside of Phi0
        force_image=-(L-c1).^2+(L-c2).^2+force_image; 
        % sum Image Force on all components (used for vector image)
        % if 'chan' is applied, this loop become one sigle code as a
        % result of layer = 1
      end

      % calculate the external force of the image 
      kappaphi0 = kappa(phi0);
      force = mu*kappaphi0./max(max(abs(kappaphi0)))+1/layer.*force_image;

      % normalized the force
      force = force./(max(max(abs(force))));

      % get stepsize dt
      dt=0.5;
      
      % get parameters for checking whether to stop
      old = phi0;
      phi0 = phi0+dt.*force;
      new = phi0;
      indicator = checkstop(old,new,dt);

%       % intermediate output
%       if(mod(n,20) == 0) 
%         showphi(I,phi0,n);  
%       end
      if indicator % decide to stop or continue 
%         showphi(I,phi0,n);

        %make mask from SDF
        seg = phi0<=0; %-- Get mask from levelset

%         subplot(2,2,4); imshow(seg); title('Global Region-Based Segmentation');

        return;
      end
    end
%     showphi(I,phi0,n);

    %make mask from SDF
    seg = phi0<=0; %-- Get mask from levelset

%     subplot(2,2,4); imshow(uint8(seg)); title('Global Region-Based Segmentation');
   case 'multiphase'
    error('multiphase case Unsupported. removed from the original code for conciseness.');
  end

