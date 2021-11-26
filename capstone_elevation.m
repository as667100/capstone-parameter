clear;close all;clc

% imshow(imread('elevation_profile.png'))
% roi = drawpolyline('Color','green');
% addlistener(roi,'MovingROI',@allevents);
% addlistener(roi,'ROIMoved',@allevents);
% x = roi.Position(:,1)
% y = roi.Position(:,2)
load('elevationData.mat', 'x');load('elevationData.mat', 'y');
figure(1);plot(x,y);title('elevation profile')
difY = [];timeInt = [];
for i = 1:length(x)-1
    
        difY =[difY (y(i+1)-y(i))/(x(i+1)-x(i))] %separate each slope into a different section
        timeInt = [timeInt (x(i+1)-x(i))/20000];%time in hr
    
end
smoothdifY = smooth(difY,'moving');
newslope = zeros(length(x));
figure(2);plot(difY);title('slope');xlabel('section');
figure(2);plot(smoothdifY);title('slope')
prefsum = 0.01
unitP = (0.002*cos(atan(difY))+sin(atan(difY))); %calculate rolling and gravity resistance

figure(2);plot(smoothdifY);title('slope');xlabel('section');ylabel('slope[m]');

% for i= 1:length(difY)-1
%     sum = timeInt(i)
%     for n = prefsum:sum
%         if sum<=20000*timeInt(i) 
%         newslope(n) = difY(i);
%         else
%         end
%         
%     end
%     prefsum=sum
%     newsum = sum+timeInt(i)
% end




%unitP = (0.002*cos(atan(difY))+sin(atan(difY)));

%Perameters from experiment
M=92 %kg of bike and average human
g=9.81 %m/s^2
Fres=M*g*unitP; % Calculate resistance

Cadance = 90; %rpm
wref= Cadance * 2*pi/60; %convert cadance to rad/s
Vel = 6.03504
CM=0.2 %Nm/A
T_motor = 24
VVat = 48 %V

%low gear force transfer
R4 = 0.675 %Rwheel
R3 = 0.05066612378891729 %Rgear
R2= 0.08093597919011615 %Rring
R1=0.1775 %Rcrankarm
pedal_res = Fres*(R2*R4)/(R1*R2); %resistance felt by user at pedal
torque_res=R4*Fres;%torque required to overcome


HR_ref = 128 %Pulled from mathmatical model

P_ref= 163; %reference power considered constant from mathmatical model
Torque_ref=P_ref/wref; %convert to torque assuming constant power output
F_ref=Torque_ref/R4; % covert to force applied by user
phi = (R1*R3)/(R2*R4) %constant for this particular bike
F_ped_ref = F_ref/phi; %convert to force at wheel
T_mot_ref= F_ped_ref*R4 
t_signal = 2*sin(1:113) + T_mot_ref;
P_F_req = pedal_res - F_ped_ref; %force to be applied by motor to keep user power and therefore HR constant




% Assume that motor will turn off when going down hill
for n = 1:length(unitP)
    if difY(n)<=0
        P_F_req(n)=0;
        timeInt(n)=0;
    end
    
end
F_wheel_req = phi*P_F_req;


% for n = 1:length(Fres)
%     if Fres(n)<0
%         Fres(n)=0;
%     
%     end
%     if Fres(n)==1
%         Fres(n)=0;
%          
%     end
% end

T_motor_req = R4*F_wheel_req;
P_motor = T_motor_req*wref;
C_used = 8.108108108*T_motor_req;
figure(2);plot(Fres);title('force applied');xlabel('section');ylabel('Force [N]');
figure(3);plot(T_motor_req);title('torque applied');xlabel('section');ylabel('Torque [Nm]');
figure(4);plot(C_used);title('current used');xlabel('section');ylabel('Current [A]');
total = x(length(x))-x(1);

% reconstruct motor 
totalP = trapz(P_motor);
totalA = trapz(C_used);
avgA = mean(C_used)
Capacity=20 %Amphours
% for n = 1:length
capacity_use =C_used.*timeInt;
motorTime = sum(timeInt)
capacity_used = sum(capacity_use)/60
%comp = 17-motorTime*20;
comp=Capacity-capacity_used;

if comp>0
    disp('travel is valid')
else
    disp('consider increase heart rate')
end

input = Torque_ref*ones(113);
disturbance = torque_res;