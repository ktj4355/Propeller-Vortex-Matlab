%% Ini
clc
clear
close all
format longG

%% initial Value Setting
vFree=[0,0,-10];
Tilt_angle =0;
Tilt_Phi =0;
RPM=5000;
nAzimuth=36;
Blade=2;

R_p=@(angle)    [cosd(angle), 0 sind(angle);
                 0          ,0,0];
R_y


alt =0;
[T, P, den, D_vis, a]=STD_Atm(alt);

if vFree(1)==0; vFree(1)=eps;end
if vFree(2)==0; vFree(2)=eps;end
if vFree(3)==0; vFree(3)=eps;end

