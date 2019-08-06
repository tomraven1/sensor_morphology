clear all


% ard=arduino();

% Connect to robot
Robot_IP = '169.254.35.190';
Socket_conn = tcpip(Robot_IP,30000,'NetworkRole','server');
fclose(Socket_conn);
disp('Press Play on Robot...')
fopen(Socket_conn);
disp('Connected!');
Orientation=[-3.14,0,0];
Translation=[-468,-83,-83.2];
moverobot(Socket_conn,Translation,Orientation);

tic
for i=1:1000
    
%     for j=1:8
%         str=sprintf('A%d',j-1);
%         voltage(i,j) = readVoltage(ard,str);
%     end
    
    
    if mod(i,200)==0
        Translation(3) = Translation(3) + 1;
        moverobot(Socket_conn,Translation,Orientation);
    elseif mod(i,100)==0
        Translation(3) = Translation(3) - 1;
        moverobot(Socket_conn,Translation,Orientation);
    end
    pause(0.05)
    tim(i)=toc;
    
end



for i=1:8
     plot(voltage(:,i)-voltage(1,i))
     hold on
end
clear ard
instrreset