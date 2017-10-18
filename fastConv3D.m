function C = fastConv3D(A, B)

%Makes a final matrix the same size as B. Assumes matrices are 3
%dimensional

M = size(A,1); %A should be a cube.
L1= size(B, 1);
L2 = size(B,2);
L3 = size(B,3);

% C = ifftn(fftn(A, [M+L1-1, M+L2-1, M+L3-1]).*fftn(B,[M+L1-1, M+L2-1, M+L3-1]), [L1, L2, L3]);

C = ifftn(fftn(A, [L1,L2,L3]).*fftn(B,[L1,L2,L3]), [L1, L2, L3]);
