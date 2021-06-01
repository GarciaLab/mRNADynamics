function Bcd = BcdGradient(x, D, tau, BcdMax)
lambda = sqrt(D*tau);
Bcd = BcdMax*exp(-x./lambda);
