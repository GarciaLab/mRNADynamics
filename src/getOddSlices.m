function imStackOdd = getOddSlices(imStack)

imStackOdd = [];
n = 0;
for i  = 1:size(imStack, 3)
    if mod(i, 2)
        n = n+1;
        imStackOdd(:, :, n) = imStack(:, :, i);
    end
end