function KG = kappa(I)
% get curvature information of input image
% input: 2D image I
% output: curvature matrix KG

% Copyright (c) 2009, 
% Yue Wu @ ECE Department, Tufts University
% All Rights Reserved  

% accelereyes: I = double(I) was changed to be more flexible
% to Jacket data-types

I = smartcast(I,'double');
[m,n] = size(I);
P = padarray(I,[1,1],1,'pre');
P = padarray(P,[1,1],1,'post');

% accelereyes: consolidated multiple subsrefs for GPU-computing efficiency

%subindex the padded image
subP1 = P(3:end,2:n+1);
subP2 = P(1:m,2:n+1);
subP3 = P(2:m+1,3:end);
subP4 = P(2:m+1,1:n);
subP5 = P(3:end,3:end);
subP6 = P(1:m,3:end);
subP7 = P(3:end,1:n);
subP8 = P(1:m,1:n);
% central difference
fy = subP1-subP2;
fx = subP3-subP4;
fyy = subP1+subP2-2*I;
fxx = subP3+subP4-2*I;
fxy = 0.25.*(subP5-subP6+subP7-subP8);

fx2 = fx.^2;
fy2 = fy.^2;
G = (fx2+fy2).^(0.5);

K = (fxx.*fy2-2*fxy.*fx.*fy+fyy.*fx2)./((fx2+fy2+eps).^(1.5));
KG = K.*G;
KG(1,:) = eps;
KG(end,:) = eps;
KG(:,1) = eps;
KG(:,end) = eps;
KG = KG./max(max(abs(KG)));
