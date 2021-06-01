function output = p_on(Bcd, Kd, KdRatio)
Numerator = (1 + Bcd/Kd).^6;
Denominator = (1 + Bcd/Kd).^6 + KdRatio;
output = double(Numerator)./double(Denominator);