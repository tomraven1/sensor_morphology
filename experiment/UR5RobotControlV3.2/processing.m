

for i=1:6
    
    y = voltage(1:100000,i);
    % Iterated median filter
    name{1} = 'Iterated medians';
    x_med(:,i) = pwc_medfiltit(y,15);
    
    % Robust total variation denoising
    name{2} = 'Robust TVD';
    x_tvd(:,i) = pwc_tvdrobust(y,4);
    
end

for i=1:6
    
    y_m(:,i)=x_med(2*interval/2:4*interval:end,i);
    y_v(:,i)=x_tvd(2*interval/2:4*interval:end,i);
    
    yo_m(:,i)=x_med(6*interval/2:4*interval:end,i);
    yo_v(:,i)=x_tvd(6*interval/2:4*interval:end,i);
    
    yn_m(:,i)=x_med(10*interval/2:4*interval:end,i);
    yn_v(:,i)=x_tvd(10*interval/2:4*interval:end,i);
end