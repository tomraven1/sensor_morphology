
curve_points=K.curve_points;% parameters of the curve (~complexity)
grid_size=K.grid_size;% planar symmetric grid
fineness=K.fineness; %precision
shp= K.num_shapes;
loc= K.num_location;
rng(K.rng_seed);

para=reshape(val,curve_points,grid_size*4);

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
    
    %  sensor(i).loc_x(2:curve_points+1)=1+grid_size*rand(curve_points,1);
    % sensor(i).loc_y(2:curve_points+1)=1+grid_size*rand(curve_points,1);
    % sensor(i+grid_size).loc_x(2:curve_points+1)=1+grid_size*rand(curve_points,1);
    % sensor(i+grid_size).loc_y(2:curve_points+1)=1+grid_size*rand(curve_points,1);
    
    p = fit(sensor(i).loc_x,sensor(i).loc_y,'pchipinterp'); %nearestinterp%pchipinterp
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
    plot(sensor(i).curve_x,sensor(i).curve_y,'b','LineWidth',2)
    hold on
    plot(sensor(i+grid_size).curve_x,sensor(i+grid_size).curve_y,'r','LineWidth',2)
    
    scatter(sensor(i).loc_x(2:curve_points+1),sensor(i).loc_y(2:curve_points+1),50,'k','filled','h')
   scatter(  sensor(i+grid_size).loc_x(2:curve_points+1), sensor(i+grid_size).loc_y(2:curve_points+1),50,'k','filled','h')
end
ct=1;
for k=1:shp
    for j=1:loc
        L = linspace(0,2*pi,K.str_shp+k);
        %xv = K.object_size*cos(L)'+(2+grid_size)/2+randn*(2+grid_size)/4;
        %yv = K.object_size*sin(L)'+(2+grid_size)/2+randn*(2+grid_size)/4;
        %xv = K.object_size*cos(L)'+(1+grid_size)/2+(rand-0.5)*(1+grid_size);
        %yv = K.object_size*sin(L)'+(1+grid_size)/2+(rand-0.5)*(1+grid_size);
        
        xv = K.object_size*cos(L)'+(1+grid_size)/2+(0)*(1+grid_size);
        yv = K.object_size*sin(L)'+(1+grid_size)/2+(0)*(1+grid_size);
        for i=1:grid_size*2
            
            %  xv = K.object_size*cos(L)'+(1+grid_size)/2;
            % yv = K.object_size*sin(L)'+(1+grid_size)/2;
            
            
            xq=sensor(i).curve_x;
            yq=sensor(i).curve_y;
            in = inpolygon(xq,yq,xv,yv);
            dx=abs(diff(xq));
            dy=abs(diff(yq));
            sensor(i).strain_x=sum(dx(in(1:end-1)));
            sensor(i).strain_y=sum(dy(in(1:end-1)));
            
            
            
            % plot(sensor(i).curve_x,sensor(i).curve_y,'b')
            % hold on
           % plot(xq(in),yq(in),'g+') % points inside
            hold on
            %  plot(xq(~in),yq(~in),'bo') % points outside
        end
        
        mat_xy(ct,:)=[sensor.strain_x]+[sensor.strain_y];
        ct=ct+1;
    end
end
%entro=1/(entropy(normalize([sensor.strain_x]))+entropy(normalize([sensor.strain_y])));
%entro=(1+entropy(normalize(mat_x))+entropy(normalize(mat_y)))
%entro=1/(1+entropy(normalize(mat_x))+entropy(normalize(mat_y)))

% 
% 
mat_xy=normalize(mat_xy,'range',[0 K.resolution]);
mat_xy=ceil(mat_xy);

entro1=0;
entro2=0;
for j=1:1*grid_size-1
    
    for i=j+1:1*grid_size
        % dis(i-1)=min(min(abs([sensor(i:grid_size).curve_y]-[sensor(i-1).curve_y])));
        %  dis(i-2+grid_size)=min(min(abs([sensor(grid_size+i:end).curve_x]-[sensor(grid_size+i-1).curve_x])));
        entro1=entro1+nvi(mat_xy(:,j),mat_xy(:,i));
        entro2=entro2+nvi(mat_xy(:,j+grid_size),mat_xy(:,i+grid_size));
        
        
    end
    %         entro1=entro1+entropy2(mat_xy(:,i));
    %         entro2=entro2+entropy2(mat_xy(:,i+grid_size));
    
end
%entro=1/(1+entro1+entro2)+sum(dis)/1000;
entro=entro1+entro2
