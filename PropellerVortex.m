clc
clear
close all



%% Read Aerodata
load("InterpolatedModel.mat")
%inducedDragClac
Cd_Tot=@(Cd,Cl,Cl0,AR) ((Cl-Cl0).^2)./(pi.*AR);

%% Read Geometry

% Read Geomety FIle
% (1)r/R, (2)r(m), (3)Chord(m), (4)Beta(deg), (5)Thickness
inputGeom=readmatrix("Geometry.xlsx");
R=inputGeom(end,2);
bmean=median(inputGeom(:,3));
Rhub=0.0084;
ZeroGeom=[0,0,2*Rhub,0,0.008];
ZeroGeom(1)=0;
ZeroGeom(2)=0;
ZeroGeom(3)=2*Rhub;
%inputGeom=[ZeroGeom;inputGeom];
AR=2*(R)./bmean;
%r/R r Chord Beta(deg) TicknessRatio
BlendPosition=[2.85 4.65].*0.0254;
BlendPosition=BlendPosition./R;      % Blend Position should be Converted to Nondimension
BladeChord=inputGeom(:,3);

%% Set Operation Condition

% Setup Condition
alt     = 0 ;         %   m
D	    = R.*2;       %   m
%V      = eps;
RPM	    = 3000    ;  %   rev/min
Blade   = 2;
[Tmp, Pressure, rho, D_vis, a] = STD_Atm(alt);
n       =RPM./60;
%J	    =V./(n.*D);
V_tip	=2*pi*n.*R;
M_tip	=V_tip./340;
Afan    =pi*R*R;
vFree= [0,0,0];
nAzmuth=1;

%% Initial Gammma
gamma_Panel_ini=ones(1,size(now_Panel_Point_BD,1))
gamma_BoundVortex=ones(1,size(now_Panel_Point_BD,1))

%% set Defalt Geometry
Geom_rR=inputGeom(:,1)
Geom_r=inputGeom(:,2)
Geom_chord=inputGeom(:,3)
Geom_Beta=inputGeom(:,4)
Geom_Thick=inputGeom(:,5)

Geom_Point_Colocation  =   [Geom_r,    -0.5.*Geom_chord.*cosd(Geom_Beta),    -0.5.*Geom_chord.*sind(Geom_Beta)] 
Geom_Point_BD          =   [Geom_r,     Geom_r.*0,                             Geom_r.*0]
Geom_Point_LE          =   [Geom_r,     0.25.*Geom_chord.*cosd(Geom_Beta),    0.25.*Geom_chord.*sind(Geom_Beta)]
Geom_Point_TE          =   [Geom_r,    -0.75.*Geom_chord.*cosd(Geom_Beta),   -0.75.*Geom_chord.*sind(Geom_Beta)]

%% set panel Geometry
Panel_rR=0.5*(Geom_rR(1:end-1,1)+Geom_rR(2:end,1))
Panel_r=0.5*(Geom_r(1:end-1,1)+Geom_r(2:end,1))
Panel_chord=0.5*(Geom_chord(1:end-1,1)+Geom_chord(2:end,1))
Panel_Beta=0.5*(Geom_Beta(1:end-1,1)+Geom_Beta(2:end,1))
Panel_Thick=0.5*(Geom_Thick(1:end-1,1)+Geom_Thick(2:end,1))

% radius (x), Up_Chord(y) Up Axis(z)
Panel_Point_Colocation  =   [Panel_r,    -0.5.*Panel_chord.*cosd(Panel_Beta),    -0.5.*Panel_chord.*sind(Panel_Beta)] 
Panel_Point_BD          =   [Panel_r,     Panel_r.*0,                             Panel_r.*0]
Panel_Point_LE          =   [Panel_r,     0.25.*Panel_chord.*cosd(Panel_Beta),    0.25.*Panel_chord.*sind(Panel_Beta)]
Panel_Point_TE          =   [Panel_r,    -0.75.*Panel_chord.*cosd(Panel_Beta),   -0.75.*Panel_chord.*sind(Panel_Beta)]

rc_panel=Panel_chord.*rc_Ratio

%% Rotation Matrix
Rz=@(theta)[cosd(theta) -sind(theta), 0;sind(theta),cosd(theta), 0;0,0,1];

R_Pitch=@(Pitch)[cosd(Pitch),  0,  sind(Pitch);
                 0,            1,  0;
                -sind(Pitch),  0,  cosd(Pitch)];

%% Position Setting
Tilit_Angle=10
AzimuthAngle=45

now_Geom_Point_Colocation  = transpose(Rz(AzimuthAngle)*Geom_Point_Colocation') 
now_Geom_Point_BD         =transpose(Rz(AzimuthAngle)*Geom_Point_BD')
now_Geom_Point_LE       =transpose(Rz(AzimuthAngle)*Geom_Point_LE')
now_Geom_Point_TE       =transpose(Rz(AzimuthAngle)*Geom_Point_TE')

now_Panel_Point_Colocation  = transpose(Rz(AzimuthAngle)*Panel_Point_Colocation') 
now_Panel_Point_BD         =transpose(Rz(AzimuthAngle)*Panel_Point_BD')
now_Panel_Point_LE       =transpose(Rz(AzimuthAngle)*Panel_Point_LE')
now_Panel_Point_TE       =transpose(Rz(AzimuthAngle)*Panel_Point_TE')


%% Calculation Rotation 



%% Geom test
figure(1)
hold on
plot(Geom_rR,Geom_Thick,'ko')
plot(Panel_rR,Panel_Thick,'rx')
grid on 


figure(2)
hold on
plot3(Panel_Point_Colocation(:,1),Panel_Point_Colocation(:,2),Panel_Point_Colocation(:,3),'r:')
plot3(Panel_Point_BD(:,1),Panel_Point_BD(:,2),Panel_Point_BD(:,3),'b:')
plot3(Panel_Point_LE(:,1),Panel_Point_LE(:,2),Panel_Point_LE(:,3),'g:')
plot3(Panel_Point_TE(:,1),Panel_Point_TE(:,2),Panel_Point_TE(:,3),'k:')
axis equal

plot3(now_Panel_Point_Colocation(:,1),now_Panel_Point_Colocation(:,2),now_Panel_Point_Colocation(:,3),'ro-')
plot3(now_Panel_Point_BD(:,1),now_Panel_Point_BD(:,2),now_Panel_Point_BD(:,3),'bx-')
plot3(now_Panel_Point_LE(:,1),now_Panel_Point_LE(:,2),now_Panel_Point_LE(:,3),'g*-')
plot3(now_Panel_Point_TE(:,1),now_Panel_Point_TE(:,2),now_Panel_Point_TE(:,3),'k*-')
axis equal