function result = concatArrayElement(A, B)
    A(end + 1) = B;
    result = A;
end
