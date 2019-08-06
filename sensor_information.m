function entro=sensor_information(para,constants)

curve_points=constants.curve_points;
grid_size=constants.grid_size;
fineness=constants.fineness;
shp= constants.num_shapes;
loc= constants.num_location;
rng(constants.rng_seed);

para=reshape(para,curve_points,grid_size*4);

for i=1:grid_size
    sensor(i).loc_x=linspace(0,grid_size+1,curve_points+2)';
    sensor(i).loc_y=i*ones(curve_points+2,1);
    sensor(i+grid_size).loc_y=linspace(0,grid_size+1,curve_points+2)';
    sensor(i+grid_size).loc_x=i*ones(curve_points+2,1);
    
end

for i=1:grid_size
    
    sensor(i).loc_x(2:curve_points+1)=para(:,i);
    sensor(i).loc_y(2:curve_points+1)=para(:,grid_size+i);
    sensor(i+grid_size).loc_x(2:curve_points+1)=para(:,i+grid_size*2);
    sensor(i+grid_size).loc_y(2:curve_points+1)=para(:,i+grid_size*3);
    
    
    
    p = fit(sensor(i).loc_x,sensor(i).loc_y,'pchipinterp'); %nearestinterp%pchipinterp%smoothingspline
    x1 = linspace(0,grid_size+1,fineness);
    y1 = p(x1);
    sensor(i).curve_x=x1';
    sensor(i).curve_y=y1;
    p = fit(sensor(i+grid_size).loc_y,sensor(i+grid_size).loc_x,'pchipinterp'); %nearestinterp%pchipinterp
    x1 = linspace(0,grid_size+1,fineness);
    y1 = p(x1);
    sensor(i+grid_size).curve_x=y1;
    sensor(i+grid_size).curve_y=x1';

    
end
for i=2:grid_size
    var1 =([sensor(i:grid_size).curve_y]-[sensor(i-1).curve_y]);
    var2=([sensor(grid_size+i:end).curve_x]-[sensor(grid_size+i-1).curve_x]);
    dis(i-1)=abs(sum(var1(var1<0)));
    dis(i-2+grid_size)=abs(sum(var2(var2<0)));
end

asd=linspace(1,grid_size,sqrt(loc));
for i=1:sqrt(loc)
    
    if mod(i,2)==1
        pos_s(1,1+(i-1)*sqrt(loc):(i)*sqrt(loc))=linspace(1,grid_size,sqrt(loc));
        pos_s(2,1+(i-1)*sqrt(loc):(i)*sqrt(loc))=asd(i)*ones(sqrt(loc),1);
    else
        pos_s(1,1+(i-1)*sqrt(loc):(i)*sqrt(loc))=linspace(grid_size,1,sqrt(loc));
        pos_s(2,1+(i-1)*sqrt(loc):(i)*sqrt(loc))=asd(i)*ones(sqrt(loc),1);
    end
    
end


ct=1;
for k=1:shp
    for j=1:loc
        L = linspace(0,2*pi,constants.str_shp+k);
        %   xv = constants.object_size*cos(L)'+(2+grid_size)/2+randn*(2+grid_size)/4;
        %  yv = constants.object_size*sin(L)'+(2+grid_size)/2+randn*(2+grid_size)/4;
        %         posx(ct)=rand;
        %         posy(ct)=rand;
        %        posx(ct)=0.5;
        %       posy(ct)=0.5;
        %   xv = constants.object_size*cos(L)'+(1+grid_size)/2+(posx(ct)-0.5)*(1+grid_size);
        %  yv = constants.object_size*sin(L)'+(1+grid_size)/2+(posy(ct)-0.5)*(1+grid_size);
        
        posx(ct)= pos_s(1,ct)+rand*2/loc;
        posy(ct)= pos_s(2,ct)+rand*2/loc;
        xv = constants.object_size*cos(L)'+posx(ct);
        yv = constants.object_size*sin(L)'+posy(ct);
        
        for i=1:grid_size*2
            
            %             xv = constants.object_size*cos(L)'+(1+grid_size)/2;
            %             yv = constants.object_size*sin(L)'+(1+grid_size)/2;
            
            xq=sensor(i).curve_x;
            yq=sensor(i).curve_y;
            in = inpolygon(xq,yq,xv,yv);
            dx=abs(diff(xq));
            dy=abs(diff(yq));
           % dis_str=constants.object_size+(xq(in(1:end-1))-posx(ct)).^2+(yq(in(1:end-1))-posy(ct)).^2;
           % sensor(i).strain_x=sum(dx(in(1:end-1))./dis_str);
          %  sensor(i).strain_y=sum(dy(in(1:end-1))./dis_str);
                sensor(i).strain_x=sum(dx(in(1:end-1)));
                sensor(i).strain_y=sum(dy(in(1:end-1)));
            
            

        end

        mat_xy(ct,:)=[sensor.strain_x]+[sensor.strain_y];
        ct=ct+1;
    end
end


%mat_xy=mat_x+mat_y;
mat_xy=normalize(mat_xy,'range',[0 constants.resolution]);
mat_xy=ceil(mat_xy);

% posx=normalize(posx,'range',[0 constants.resolution*10]);
% posx=ceil(posx);
% posy=normalize(posy,'range',[0 constants.resolution*10]);
% posy=ceil(posy);

 smooth=sum(sum(abs(diff(diff(mat_xy)))));


entro=smooth/10000000+1/(1+h(mat_xy(:,1:grid_size))+h(mat_xy(:,grid_size+1:end)))+sum(dis)/200000;
%entro=smooth/10000000+1/(1+h(mat_xy(:,1:end)))+sum(dis)/200000;







%
% 
% entro1=0;
% for i=1:2*grid_size
%     entro1=entro1+h(mat_xy(:,i));
%     % entro1=condh(posx',mat_xy(:,i))+condh(posy',mat_xy(:,i));
% end
% % j_entro=h(mat_xy);
% % %entro=1/(0.01+entro1-j_entro);
% % entro=smooth/2000000+10/(0.01+entro1-j_entro)+sum(dis)/10000;
% 
% %entro=smooth/1000+10/(0.01+j_entro)+sum(dis)/2000;
%  entro=smooth/10000000+1/(1+entro1)+sum(dis)/2000000;
% % entro=1/j_entro+sum(dis)/200;
% %entro=1/(1+entro1);

% 
% entro1=0;
% entro2=0;
% for i=1:1*grid_size-1
% 
% %     for i=j+1:1*grid_size
% %         % dis(i-1)=min(min(abs([sensor(i:grid_size).curve_y]-[sensor(i-1).curve_y])));
% %         %  dis(i-2+grid_size)=min(min(abs([sensor(grid_size+i:end).curve_x]-[sensor(grid_size+i-1).curve_x])));
% % %         entro1=entro1+nvi(mat_xy(:,j),mat_xy(:,i));
% % %         entro2=entro2+nvi(mat_xy(:,j+grid_size),mat_xy(:,i+grid_size));
% % 
% %     end
%            entro1=entro1+entropy2(mat_xy(:,i));
%           entro2=entro2+entropy2(mat_xy(:,i+grid_size));
% 
% end


% entro=condh(posx',mat_xy)+condh(posy',mat_xy);

% entro=condh(posx',mat_xy)+condh(posy',mat_xy);
%entro=smooth/1000000+1/(1+entro1+entro2)+sum(dis)/200000;
%entro=smooth/100000000000+100/(1+h(mat_xy))+sum(dis)/2000;
end