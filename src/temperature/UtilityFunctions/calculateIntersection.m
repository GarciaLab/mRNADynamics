function [x0, y0] = calculateIntersection(a,b,c,d)
x0 = (d-b)/(a-c);
y0 = a*x0 + b;
end