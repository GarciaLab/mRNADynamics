function y =  leftHalfTrapezoidFitFunction(t, a, b, t1)
    y = zeros(size(t));
    a = abs(a);
    %b = -abs(b);
    %t1 = abs(t1);
    %t2 = max(t1, t2);
    y(t < t1) = a*t(t < t1) + b;
    y(( t >= t1) ) = a*t1+b;
    %y[:int(tau1)] =  a*x[:int(tau1)] + b
    %y[int(tau1):int(tau2)] =  a*tau1 + b 
    %y[int(tau2):] = c*(x[int(tau2):]-tau2) + (a*tau1 + b)
end