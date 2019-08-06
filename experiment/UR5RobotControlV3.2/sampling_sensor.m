clear all

rng(100)

ard=arduino();

% Connect to robot
Robot_IP = '169.254.35.190';
Socket_conn = tcpip(Robot_IP,30000,'NetworkRole','server');
fclose(Socket_conn);
disp('Press Play on Robot...')
fopen(Socket_conn);
disp('Connected!');
Orientation=[-3.14,0,0];
%Translation=[-472,-79,-80.2];%grid
Translation=[-318,-80,-80.2];
moverobot(Socket_conn,Translation,Orientation);
Translation_inner=Translation;
pause(5)

flag=1;
grid_length=15;
interval=10;
deformation=5;
k=1;
tic
for i=1:500000
    
%     
    for j=1:7
        str=sprintf('A%d',j-1);
        voltage(i,j) = readVoltage(ard,str);
    end
 
    
    if flag==interval
        pos(:,k)=(rand(2,1)-0.5)*grid_length*2*1.5;
       % asd=[15,15;10,5;-15,-10]';
        %pos(:,k)=asd(:,randi(3));%mm
        Translation_inner(1:2) = Translation(1:2)+pos(:,k)';
        Translation_inner(3) = Translation(3);
        moverobot(Socket_conn,Translation_inner,Orientation);
        flag=flag+1;
    elseif flag==interval*2
        Translation_inner(1:2) = Translation(1:2)+pos(:,k)';
        Translation_inner(3) = Translation(3)-deformation ;
        moverobot(Socket_conn,Translation_inner,Orientation);
        flag=flag+1;
      
    elseif flag==interval*4
        Translation_inner(1:2) = Translation(1:2)+pos(:,k)';
        Translation_inner(3) = Translation(3);
        moverobot(Socket_conn,Translation_inner,Orientation);
        k=k+1;
        flag=1;
    else
        
        flag=flag+1;
        
    end
    %pause(0.1)
    tim(i)=toc;
    
end

moverobot(Socket_conn,Translation,Orientation);
save para5_long_n.mat


% for i=1:8
%     plot(voltage(:,i)-voltage(1,i))
%     hold on
% end

clear ard
instrreset