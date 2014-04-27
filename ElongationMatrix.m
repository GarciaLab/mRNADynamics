function A = ElongationMatrix(NF,NT,Memory)

if nargin<3
    Memory=8;
end

A=zeros(NF,NT);
for i=1:NT
    MaxIdx=min(i+Memory-1,NF);
    A(i:MaxIdx, i)=1;
end
