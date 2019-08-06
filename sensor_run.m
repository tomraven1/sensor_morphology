function out=sensor_run(pos,para,constants)

curve_points=constants.curve_points;
grid_size=constants.grid_size;
fineness=constants.fineness;
shp= constants.num_shapes;
loc= constants.num_location;




for i=1:grid_size
    sensor(i).loc_x=linspace(0,grid_size+1,curve_points+2)';
    sensor(i).loc_y=i*ones(curve_points+2,1);
    sensor(i+grid_size).loc_y=linspace(0,grid_size+1,curve_points+2)';
    sensor(i+grid_size).loc_x=i*ones(curve_points+2,1);
    
end

for i=1:grid_size
    %
    sensor(i).loc_x(2:curve_points+1)=para(:,i);
    sensor(i).loc_y(2:curve_points+1)=para(:,grid_size+i);
    sensor(i+grid_size).loc_x(2:curve_points+1)=para(:,i+grid_size*2);
    sensor(i+grid_size).loc_y(2:curve_points+1)=para(:,i+grid_size*3);

    %
    
    
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
    % sensor(i).curve=fit( sensor(i).loc_x', sensor(i).loc_y,'poly2');
    %plot(sensor(i).curve_x,sensor(i).curve_y,'b')
    %hold on
    %plot(sensor(i+grid_size).curve_x,sensor(i+grid_size).curve_y,'r')
    
end



L = linspace(0,2*pi,constants.str_shp+1);
%   xv = constants.object_size*cos(L)'+(2+grid_size)/2+randn*(2+grid_size)/4;
%  yv = constants.object_size*sin(L)'+(2+grid_size)/2+randn*(2+grid_size)/4;
xv = constants.object_size*cos(L)'+pos.x;
yv = constants.object_size*sin(L)'+pos.y;
for i=1:grid_size*2
    
    %             xv = constants.object_size*cos(L)'+(1+grid_size)/2;
    %             yv = constants.object_size*sin(L)'+(1+grid_size)/2;
    
    xq=sensor(i).curve_x;
    yq=sensor(i).curve_y;
    in = inpolygon(xq,yq,xv,yv);
    dx=abs(diff(xq));
    dy=abs(diff(yq));
    dis_str=constants.object_size+(xq(in(1:end-1))-pos.x).^2+(yq(in(1:end-1))-pos.y).^2;
    sensor(i).strain_x=sum(dx(in(1:end-1))./dis_str);
    sensor(i).strain_y=sum(dy(in(1:end-1))./dis_str);
%     sensor(i).strain_x=sum(dx(in(1:end-1)));
%     sensor(i).strain_y=sum(dy(in(1:end-1)));
    
    %
    %     plot(sensor(i).curve_x,sensor(i).curve_y,'b')
    %     hold on
    %
    %     plot(xq(in),yq(in),'r+') % points inside
    %    hold on
    %    plot(xq(~in),yq(~in),'bo') % points outside
end
raw_values(1,:)=[sensor.strain_x]+[sensor.strain_y];



out=raw_values;
end