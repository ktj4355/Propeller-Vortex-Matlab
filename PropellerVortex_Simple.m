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
RPM	    = 5000    ;  %   rev/min
vFree= [0,0,-0];
nAzmuth=1;
rc_Ratio=0.2;
vortex_n=1.06;
Blade   = 2;
[Tmp, Pressure, rho, D_vis, a] = STD_Atm(alt);
n       =RPM./60;
%J	    =V./(n.*D);
V_tip	=2*pi*n.*R;
M_tip	=V_tip./340;
Afan    =pi*R*R;





%% set Defalt Geometry
Geom_rR=inputGeom(:,1);
Geom_r=inputGeom(:,2);
Geom_chord=inputGeom(:,3);
Geom_Beta=inputGeom(:,4);
Geom_Thick=inputGeom(:,5);

Geom_Point_Colocation  =   [Geom_r,    -0.5.*Geom_chord.*cosd(Geom_Beta),    -0.5.*Geom_chord.*sind(Geom_Beta)];
Geom_Point_BD          =   [Geom_r,     Geom_r.*0,                             Geom_r.*0];
Geom_Point_LE          =   [Geom_r,     0.25.*Geom_chord.*cosd(Geom_Beta),    0.25.*Geom_chord.*sind(Geom_Beta)];
Geom_Point_TE          =   [Geom_r,    -0.75.*Geom_chord.*cosd(Geom_Beta),   -0.75.*Geom_chord.*sind(Geom_Beta)];

%% set panel Geometry
Panel_rR=0.5*(Geom_rR(1:end-1,1)+Geom_rR(2:end,1));
Panel_r=0.5*(Geom_r(1:end-1,1)+Geom_r(2:end,1));
Panel_chord=0.5*(Geom_chord(1:end-1,1)+Geom_chord(2:end,1));
Panel_Beta=0.5*(Geom_Beta(1:end-1,1)+Geom_Beta(2:end,1));
Panel_Thick=0.5*(Geom_Thick(1:end-1,1)+Geom_Thick(2:end,1));

Panel_dr=abs((Geom_r(1:end-1,1)-Geom_r(2:end,1)));
Panel_ds=Panel_dr.*Panel_chord;
Blade_Area=sum(Panel_ds);

% radius (x), Up_Chord(y) Up Axis(z)
Panel_Point_Colocation  =   [Panel_r,    -0.5.*Panel_chord.*cosd(Panel_Beta),    -0.5.*Panel_chord.*sind(Panel_Beta)];
Panel_Point_BD          =   [Panel_r,     Panel_r.*0,                             Panel_r.*0];
Panel_Point_LE          =   [Panel_r,     0.25.*Panel_chord.*cosd(Panel_Beta),    0.25.*Panel_chord.*sind(Panel_Beta)];
Panel_Point_TE          =   [Panel_r,    -0.75.*Panel_chord.*cosd(Panel_Beta),   -0.75.*Panel_chord.*sind(Panel_Beta)];

%% Initial Gammma /Rc

gamma_Panel_ini=ones(size(Panel_Point_BD,1),1);

gamma_BoundVortex=ones(size(Panel_Point_BD,1),1)*0;

rc_panel=Panel_chord.*rc_Ratio;
rc_Geom=Geom_chord.*rc_Ratio;

%% Rotation Matrix
Rz=@(theta)[cosd(theta) -sind(theta), 0;sind(theta),cosd(theta), 0;0,0,1];

R_Pitch=@(Pitch)[cosd(Pitch),  0,  sind(Pitch);
    0,            1,  0;
    -sind(Pitch),  0,  cosd(Pitch)];

%% Position Setting




dAngle=10; %deg
rotSpeed=2*pi*RPM/60;  %rad/s
rot_deg_speed=360*RPM/60; %deg/s
dt=dAngle./rot_deg_speed;
%initial Setting

Time=0;
Tilit_Angle=0;
AzimuthAngle=0;
Total_Rotate=0;
Wake_Geom_Position=[]; %x y z x y z
Wake_velocity=[];
Wake_Gamma=[];
globalData=[];
while Total_Rotate<2500


    now_Geom_Point_Colocation  = transpose(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle)*Geom_Point_Colocation');
    now_Geom_Point_BD         =transpose(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle)*Geom_Point_BD');
    now_Geom_Point_LE       =transpose(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle)*Geom_Point_LE');
    now_Geom_Point_TE       =transpose(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle)*Geom_Point_TE');

    now_Panel_Point_Colocation  = transpose(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle)*Panel_Point_Colocation');  %[x y z] coordinate
    now_Panel_Point_BD         =transpose(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle)*Panel_Point_BD');
    now_Panel_Point_LE       =transpose(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle)*Panel_Point_LE');
    now_Panel_Point_TE       =transpose(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle)*Panel_Point_TE');


    %% Calculation IJ Matrix


    for ColloIDX=1:size(now_Panel_Point_Colocation,1) %Collocation Point Loop
        ColoP=now_Panel_Point_Colocation(ColloIDX,:);
        for gamma_idx=1:length(gamma_Panel_ini) %Bound GaMMA Loop

            %Vortex Modeling


            gamma=gamma_Panel_ini(gamma_idx);

            % induced velocity by Bound Vortex
            vind_Bound(gamma_idx,:)     =Vortex_Vatistas(now_Geom_Point_BD(gamma_idx,:),now_Geom_Point_BD(gamma_idx+1,:),ColoP,gamma,rc_panel(gamma_idx),vortex_n);
            vind_trLeft(gamma_idx,:)    =Vortex_Vatistas(now_Geom_Point_TE(gamma_idx,:),now_Geom_Point_BD(gamma_idx,:),ColoP,gamma,rc_panel(gamma_idx),vortex_n);
            vind_trRight(gamma_idx,:)   =Vortex_Vatistas(now_Geom_Point_BD(gamma_idx+1,:),now_Geom_Point_TE(gamma_idx+1,:),ColoP,gamma,rc_panel(gamma_idx),vortex_n);
            vTotal=vind_Bound(gamma_idx,:)+ vind_trLeft(gamma_idx,:)+ vind_trRight(gamma_idx,:);

            I_ini_x(ColloIDX,gamma_idx)=vTotal(1);
            I_ini_y(ColloIDX,gamma_idx)=vTotal(2);
            I_ini_z(ColloIDX,gamma_idx)=vTotal(3);
        end
        Vind_collo(ColloIDX,:)=sum(vind_Bound)+sum(vind_trLeft)+sum(vind_trRight);


    end
    %% Wake Induced Velocity

    for ColloIDX=1:size(now_Panel_Point_Colocation,1) %Collocation Point Loop
        ColoP=now_Panel_Point_Colocation(ColloIDX,:);
        vind_Wake(ColloIDX,:)=[0,0,0];
        for idx=1:size(Wake_Gamma,2)
            xind=(idx-1)*3+1;
            yind=(idx-1)*3+2;
            zind=(idx-1)*3+3;
            if size(Wake_Geom_Position,1)<=2; break; end
            for Tidx=1:size(Wake_Geom_Position,1)-1
                
                A=Wake_Geom_Position(Tidx,xind:zind);
                B=Wake_Geom_Position(Tidx+1,xind:zind);
                gamma=Wake_Gamma(Tidx,idx);
                vind_Wake(ColloIDX,:)=vind_Wake(ColloIDX,:)+Vortex_Vatistas(A,B,ColoP,gamma,rc_panel(gamma_idx),vortex_n);
            end




        end
    end
     figure(33)
    clf
     hold on

    quiver3(now_Panel_Point_Colocation(:,1),now_Panel_Point_Colocation(:,2),now_Panel_Point_Colocation(:,3),vind_Wake(:,1),vind_Wake(:,2),vind_Wake(:,3))




    %% Blade Load Model
    Told=999;
    unit_ground_Span=[1,0,0]';
    unit_ground_Chord=[0,1,0]';
    unit_ground_Axis=[0,0,1]';

    % Axis Position Coordinate

    unit_conv_Span=(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle)*unit_ground_Span);
    unit_conv_Chord=(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle)*unit_ground_Chord);
    unit_conv_Axis=(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle)*unit_ground_Axis);
    vinduced=[];

    %gamma=0.5Veff*cl
    for iter=1:1000
        vinduced=[];
V_inflow=[];
        T=0;
        P=0;
        Q=0;
        F=0;
        v_infow_induced_x_set=I_ini_x*gamma_BoundVortex;
        v_infow_induced_y_set=I_ini_y*gamma_BoundVortex;
        v_infow_induced_z_set=I_ini_z*gamma_BoundVortex;

        Data_local_raidus=[];
        for ridx=1:size(now_Panel_Point_Colocation,1)



            V_inflow_free_Axis = dot(vFree,unit_conv_Axis);
            V_inflow_free_span = dot(vFree,unit_conv_Span);
            V_inflow_free_chord = dot(vFree,unit_conv_Chord);


            r=Panel_r(ridx);

            V_inflow_rot_Axis=0;
            V_inflow_rot_span=0;
            V_inflow_rot_chord=-(rotSpeed.*r);

            v_infow_induced_x=v_infow_induced_x_set(ridx);
            v_infow_induced_y=v_infow_induced_y_set(ridx);
            v_infow_induced_z=v_infow_induced_z_set(ridx);

            v_infow_induced=[v_infow_induced_x,v_infow_induced_y,v_infow_induced_z];
            V_inflow_induced_Axis = dot(v_infow_induced,unit_conv_Axis);
            V_inflow_induced_span = dot(v_infow_induced,unit_conv_Span);
            V_inflow_induced_chord = dot(v_infow_induced,unit_conv_Chord);

            v_infow_WakeInduced= vind_Wake(ridx,:);
            V_inflow_WakeInduced_Axis = dot(v_infow_WakeInduced,unit_conv_Axis);
            V_inflow_WakeInduced_Span = dot(v_infow_WakeInduced,unit_conv_Span);
            V_inflow_WakeInduced_Chord = dot(v_infow_WakeInduced,unit_conv_Chord);

            
            V_induced_sum_axis=V_inflow_rot_Axis+V_inflow_induced_Axis+V_inflow_WakeInduced_Axis;
            V_induced_sum_span=V_inflow_rot_span+V_inflow_induced_span+V_inflow_WakeInduced_Span;
            V_induced_sum_chord=V_inflow_rot_chord+V_inflow_induced_chord+V_inflow_WakeInduced_Chord;


            V_inflow_sum_axis=V_inflow_free_Axis+V_induced_sum_axis;
            V_inflow_sum_span=V_inflow_free_span+V_induced_sum_span;
            V_inflow_sum_chord=V_inflow_free_chord+V_induced_sum_chord;
            vinduced=[vinduced;V_inflow_sum_span,V_inflow_sum_chord,V_inflow_sum_axis];
            V_inflow=[ V_inflow;vFree];

            Flow_Angle=rad2deg(atan2(-V_inflow_sum_axis,-V_inflow_sum_chord));
            Veff=sqrt((V_inflow_sum_axis.^2+V_inflow_sum_chord.^2));
            beta=Panel_Beta(ridx);
            rR_local=r./R;
            c=Panel_chord(ridx);
            ds=Panel_ds(ridx);
            dr=Panel_dr(ridx);


            alpha=beta-Flow_Angle;


            Re=(rho*Veff*c)./D_vis;
            ThickRatio=Panel_Thick(ridx).*100;


            if BlendPosition(1)>=rR_local
                airfoil=1;
                CL0=Section1_CL(ThickRatio, Re ,0);
                inp_Cl=Section1_CL(ThickRatio, Re ,alpha);
                inp_Cd=Section1_CD(ThickRatio, Re ,alpha);
            elseif BlendPosition(2)<=rR_local
                airfoil=2;
                CL0=Section1_CL(ThickRatio, Re ,0);
                inp_Cl=Section2_CL(ThickRatio ,Re, alpha);
                inp_Cd=Section2_CD(ThickRatio, Re, alpha);
            else
                airfoil=0;
                BR=abs((rR_local-BlendPosition(1))./(BlendPosition(2)-BlendPosition(1)));
                CL0=BLD_CL(BR, Re ,0);
                inp_Cl=BLD_CL(BR, Re ,alpha);
                inp_Cd=BLD_CD(BR ,Re, alpha);
            end

            Cl=inp_Cl;
            CD=inp_Cl;

            LD=Cl/CD;
            %Calculate BET
            gamma = rad2deg(atan(1/LD));
            K=Cl.*c./((sin(deg2rad(Flow_Angle)).^2).*cos(deg2rad(gamma)));
            Tc=K.*cos(deg2rad(Flow_Angle+gamma));
            Qc=r.*K.*sin(deg2rad(Flow_Angle+gamma));
            Fc=K.*sin(deg2rad(Flow_Angle+gamma));

            Tdr=(0.5*rho*V_inflow_sum_axis.^2).*Tc.*dr;
            Qdr=(0.5*rho*V_inflow_sum_axis.^2).*Qc.*dr;
            Fdr=(0.5*rho*V_inflow_sum_axis.^2).*Fc.*dr;

            % figure(10)
            % hold on
            % quiver(-V_inflow_free_chord,-V_inflow_free_Axis,V_inflow_free_chord,V_inflow_free_Axis,1,'r')
            % quiver(-V_inflow_sum_chord,-V_inflow_sum_axis,V_inflow_sum_chord,V_inflow_sum_axis,1,'g')
            % quiver(-V_inflow_rot_chord,0,V_inflow_rot_chord,0,1,'b')
            % grid on
            % axis equal
            %
            boundVortex=0.5*Cl*Veff*c;
            kfac=0.05;
            gamma_BoundVortex(ridx)=gamma_BoundVortex(ridx)+(-gamma_BoundVortex(ridx)+boundVortex)*kfac;

            Data_local_raidus=[Data_local_raidus;r./R,r,Flow_Angle,alpha,Re,beta,boundVortex,Tdr];

            T=T+2*Tdr;
            Q=Q+2*dr;
            F=F+2*dr;
        end
        T;


        if abs(T-Told)<0.001
            disp("수렴!")
            break;
        else

        end
        Told=T;

    end


    %% Wake Structure (Semi-Prescribed Wake Sturcture)
    %TE -> Wake Structure 추가

    tmp=now_Geom_Point_TE';
    wakeForm=tmp(:)';
    vinduced_Ax=2*vinduced(:,3); %Blade 개수만큼 유도 전파속도 증가
    vinduced_Ax_set=transpose(vinduced_Ax.*unit_conv_Axis'+V_inflow);
    vinduced_Ax_GeomSet=[];
    vinduced_Ax_GeomSet=0.5.*(vinduced_Ax_set(:,2:end)+vinduced_Ax_set(:,1:end-1));
    vinduced_Ax_GeomSet=[vinduced_Ax_set(:,1), vinduced_Ax_GeomSet, vinduced_Ax_set(:,end)];
    vinduced_Ax_GeomSet=vinduced_Ax_GeomSet(:);
    WakeGammaset=[0;gamma_BoundVortex]-[gamma_BoundVortex;0];
    WakeGammaset(1)=0;
    Wake_Geom_Position=[wakeForm;Wake_Geom_Position]; %x y z x y z
    Wake_velocity=[vinduced_Ax_GeomSet';Wake_velocity];
    Wake_Gamma=[WakeGammaset';Wake_Gamma];
    Wake_Geom_Position=Wake_Geom_Position+Wake_velocity.*dt;

    %Trail Gamma 계산
    %shed Gamma 계산

    %% Rotate Next Structure
    Total_Rotate=Total_Rotate+dAngle;
    AzimuthAngle=AzimuthAngle+dAngle;
    if AzimuthAngle>360
        AzimuthAngle=AzimuthAngle-360;
    end

    Time=Time+dt;

    figure(100)
    clf
    hold on

    plot3(now_Panel_Point_Colocation(:,1),now_Panel_Point_Colocation(:,2),now_Panel_Point_Colocation(:,3),'ro-')
    plot3(now_Panel_Point_BD(:,1),now_Panel_Point_BD(:,2),now_Panel_Point_BD(:,3),'bx-')
    plot3(now_Panel_Point_LE(:,1),now_Panel_Point_LE(:,2),now_Panel_Point_LE(:,3),'g*-')
    plot3(now_Panel_Point_TE(:,1),now_Panel_Point_TE(:,2),now_Panel_Point_TE(:,3),'k*-')
    for idx=1:length(WakeGammaset)
        xind=(idx-1)*3+1;
        yind=(idx-1)*3+2;
        zind=(idx-1)*3+3;

        plot3(Wake_Geom_Position(1:end,xind),Wake_Geom_Position(1:end,yind),Wake_Geom_Position(1:end,zind),'k:')
    end
    view([1,1,1])
    axis equal
globalData=[globalData;Time,Total_Rotate,AzimuthAngle, T];
figure(6)
clf
hold on
plot(globalData(:,2),globalData(:,4))


end

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
figure(3)
hold on
quiver3(now_Panel_Point_Colocation(:,1),now_Panel_Point_Colocation(:,2),now_Panel_Point_Colocation(:,3),Vind_collo(:,1),Vind_collo(:,2),Vind_collo(:,3))


axis equal

figure(4)
clf
hold on
plot(Data_local_raidus(:,1),Data_local_raidus(:,3),'r') % Flow
plot(Data_local_raidus(:,1),Data_local_raidus(:,4),'g') % alpha
plot(Data_local_raidus(:,1),Data_local_raidus(:,6),'b') %beta
legend("Flow","Alpha","beta")


figure(5)
clf
hold on
plot(Data_local_raidus(:,1),Data_local_raidus(:,7),'r') % Flow
