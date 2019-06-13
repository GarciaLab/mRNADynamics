function c = SumArrays(A,B)

%Maps each element of a larger array into a smaller array

a = [];
b = [];

if length(B) >= length(A)
    b = B;
    a = A;
else 
    b = A;
    a = B;
end

c = 0*a;
ratio = length(a) / length(b);

for i=1:length(b)
    c(ceil(i*ratio)) = c(ceil(i*ratio))+b(i);
end

c = c + a;

    