function gausslogL = gaussLogLFunction(dx,sigma)
  len = length(dx);
  sigma_inv = inv(sigma);
  seigma_det = det(sigma);
  gausslogL = 0.5*log((2*pi)^len * sigma_det) - 0.5*dx'*sigma_inv*dx;