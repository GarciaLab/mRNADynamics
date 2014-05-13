function A = ElongationMatrix(NF,NT,Memory,ExcessF)

if nargin<3
    Memory=8;
end
if nargin<4
    ExcessF=0;
end

A=zeros(NF,NT);
for i=1:NT
    % Propagate through normal memory (i.e. multiple of TimeRes)
    MaxIdx=min(i+Memory-1,NF);
    A(i:MaxIdx, i)=1;
    
    % Get excess fluorescence
%     MaxIdx=min(i+Memory,NF);
%     A(MaxIdx,i)=ExcessF;

    % Get excess fluorescence
    MaxIdx=i+Memory;
    if MaxIdx<=NF
        A(MaxIdx,i)=ExcessF;
    end
end
