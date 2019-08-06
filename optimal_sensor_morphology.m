addpath('InfoTheory')
addpath('MIToolbox-3.0.0\matlab')

iters=20;

K.resolution=4096;
K.curve_points=10;% parameters of the curve (~complexity)
K.grid_size=4;% planar symmetric grid
K.fineness=2000; %precision 
K.object_size=0.6;
K.num_shapes=1;
K.num_location=49; %number of points used for measuring entropy. square number
K.rng_seed=214;
K.str_shp=30;%3 is traingle 4 is square..etc

nvars=K.curve_points*2*2*K.grid_size;%number of variables

fun=@(x)sensor_information(x,K);
lb=1.001*ones(nvars,1);
ub=(K.grid_size-0.001)*ones(nvars,1);
lb(1:K.grid_size)=0.001;
ub(1:K.grid_size)=0.999+K.grid_size;
lb(1+K.grid_size*3:end)=0.001;
ub(1+K.grid_size*3:end)=0.999+K.grid_size;

options = optimoptions(@ga,'Display','iter','MaxGenerations',iters);
val = ga(fun,nvars,[],[],[],[],lb,ub,[],1,options);

% x0=rand(nvars,1);
% val = fmincon(fun,x0,[],[],[],[],lb,ub);
%entro=sensor_information(5,K)
