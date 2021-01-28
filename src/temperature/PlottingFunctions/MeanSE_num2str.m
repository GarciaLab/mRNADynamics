function legstr = MeanSE_num2str(mean_param, se_param, Nsigfigs)
legstr = {};
mean_order = floor(log(abs(mean_param))./log(10));
se_order = floor(log(abs(se_param))./log(10));
legstr.m = num2str(mean_param, Nsigfigs);
if mean_order >= se_order
    
    if Nsigfigs >(mean_order-se_order)
        legstr.se = num2str(se_param, Nsigfigs-(mean_order-se_order));
    else
        legstr.se = num2str(0);
    end
else
    legstr.se = num2str(se_param, Nsigfigs);
end