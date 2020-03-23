function [precision, image] = getImagePrecision(image)

isInteger = sum(mod(image(:), 1));

if isInteger
    if max(image(:)) > 256
        precision = 'uint16';
    else
        precision = 'uint8';
        image = uint8(image);
    end
else
    precision = 'double';
end

end