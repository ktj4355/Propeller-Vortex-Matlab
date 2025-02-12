clc
clear
close all
    R=@(angle)[cos(angle) -sin(angle) 0;sin(angle) cos(angle) 0;0 0 1];

Ini_coordinate=[1 0 0];

rpm=1000;
r=norm(Ini_coordinate);
n=rpm/60.0;
vth=[0,0,n*2*pi]; % Radial Coordinate rad/s
vFree=[0 0 -10];  % orthogonal Coordinate

dt=0.005;  %sec
r_history=[];
r_now=Ini_coordinate;
maxiter=300;
iter=0;
rold=[];

figure(1)
clf;
hold on

while iter<maxiter

    iter=iter+1;
    rold=r_now;
    r_history=[rold;r_history];
    dangle=(vth(3).*dt);
    rot=R(dangle);
    r_now=(rot*rold')';

    %radial to Orthogonal
  
    vtotal=vFree;
    for idx=1:size(r_history,1)
        r_history(idx,:)=r_history(idx,:)+vtotal.*dt;
    end
    figure(1)
    clf
    view([1,1,1])
    hold on
    plot3(r_history(:,1),r_history(:,2),r_history(:,3),'b-')
    plot3([0 r_now(1)],[0 r_now(2)],[0 r_now(3)],'ko-')
    axis equal
    
end
