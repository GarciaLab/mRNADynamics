function ysample = GetMasterProfileForEmbryoTimes(xsample, MasterProfile, xfits)
%%
TFxfits = (sum(isnan(MasterProfile(:,:)), 2) < size(MasterProfile, 2)).';
xfits = xfits(TFxfits);
MasterProfile = MasterProfile(TFxfits,:);
ysample = interp1(xfits, MasterProfile, xsample);
f = @(b,x) b(1).*exp(b(2).*x) + b(3);
