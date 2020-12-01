function logL_vec = calculateIndividualProbs(KFTrack,z)

H = KFTrack.kalmanFilter.MeasurementModel;
P = KFTrack.kalmanFilter.StateCovariance;
R = KFTrack.kalmanFilter.MeasurementNoise;
x = KFTrack.kalmanFilter.State;

% perform calculations
% z_res = z - H*x;
% S = R + H*P*H';
% d2 = z_res'*inv(S)*z_res;
% d_n = d2 + log(det(S));

logL_vec = [];
for i = 1:length(z)   
    z_res = z(i) - H(i,:)*x;
    S = R(i,i) + P(2*(i-1)+1,2*(i-1)+1);
    d2 = z_res'*inv(S)*z_res;
    logL_vec(i) = d2 + log(det(S));
end
logL_vec = -logL_vec;