for z=1
    rng(1)
    %K.object_size=1;
    %K.str_shp=4;
    sample_no=1024*1;
    K.grid_size=2^z;
    
    asd=linspace(1,K.grid_size,sqrt(sample_no));
    for i=1:sqrt(sample_no)
        
        if mod(i,2)==1
            pos_s(1,1+(i-1)*sqrt(sample_no):(i)*sqrt(sample_no))=linspace(1,K.grid_size,sqrt(sample_no));
            pos_s(2,1+(i-1)*sqrt(sample_no):(i)*sqrt(sample_no))=asd(i)*ones(sqrt(sample_no),1);
        else
            pos_s(1,1+(i-1)*sqrt(sample_no):(i)*sqrt(sample_no))=linspace(K.grid_size,1,sqrt(sample_no));
            pos_s(2,1+(i-1)*sqrt(sample_no):(i)*sqrt(sample_no))=asd(i)*ones(sqrt(sample_no),1);
        end
        
    end
    
   % para=K.grid_size*rand(1,2*K.grid_size*2);
    
    for i=1:sample_no
        pos.x=pos_s(1,i)+rand*5/loc;
        pos.y=pos_s(2,i)+rand*5/loc;
        sens_raw_all(i,:)=sensor_run(pos,para,K);
        loc_all(i,:)=[pos.x;pos.y];
    end
    
%     fname = sprintf('z_1para%d.mat', z);
%     save(fname)
%     clear sens_raw_all
%     clear loc_all
    
end