clc
clear
close all


GlobalTimer=tic;
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
vFree= [0,5,-0];
nAzmuth=1;
rc_Ratio=1;
vortex_n=1.06;
Blade   = 2;


dAngle=10; %deg
rotSpeed=2*pi*RPM/60;  %rad/s
rot_deg_speed=360*RPM/60; %deg/s
dt=dAngle./rot_deg_speed;
%initial Setting


[Tmp, Pressure, rho, D_vis, a] = STD_Atm(alt);
n       =RPM./60;
%J	    =V./(n.*D);
V_tip	=2*pi*n.*R;
M_tip	=V_tip./340;
Afan    =pi*R*R;



figure(100)
set( gcf, 'Position', [0 0 1800 960] )

        mkdir 'fig1'

        mkdir 'fig2'

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
gamma_BoundVortex2=ones(size(Panel_Point_BD,1),1)*0;

rc_panel=Panel_chord.*rc_Ratio;
%rc_panel=ones(size(Panel_chord)).*mean(Panel_chord).*rc_Ratio

rc_Geom=Geom_chord.*rc_Ratio;
%rc_Geom=ones(size(Geom_chord)).*mean(Geom_chord).*rc_Ratio

%% Rotation Matrix
Rz=@(theta)[cosd(theta) -sind(theta), 0;sind(theta),cosd(theta), 0;0,0,1];

R_Pitch=@(Pitch)[cosd(Pitch),  0,  sind(Pitch);
    0,            1,  0;
    -sind(Pitch),  0,  cosd(Pitch)];

%% Position Setting

Time=0;
Tilit_Angle=0;
AzimuthAngle=0;
Total_Rotate=0;
Wake_Geom_Position=[]; %x y z x y z
Wake_velocity=[];
Wake_Gamma=[];

Wake2_Geom_Position=[]; %x y z x y z
Wake2_velocity=[];
Wake2_Gamma=[];
globalData=[];
rotateCount=0;
while Total_Rotate<5000


    now_Geom_Point_Colocation  = transpose(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle)*Geom_Point_Colocation');
    now_Geom_Point_BD         =transpose(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle)*Geom_Point_BD');
    now_Geom_Point_LE       =transpose(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle)*Geom_Point_LE');
    now_Geom_Point_TE       =transpose(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle)*Geom_Point_TE');

    now_Panel_Point_Colocation  = transpose(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle)*Panel_Point_Colocation');  %[x y z] coordinate
    now_Panel_Point_BD         =transpose(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle)*Panel_Point_BD');
    now_Panel_Point_LE       =transpose(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle)*Panel_Point_LE');
    now_Panel_Point_TE       =transpose(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle)*Panel_Point_TE');


    now_Geom2_Point_Colocation  = transpose(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle+180)*Geom_Point_Colocation');
    now_Geom2_Point_BD         =transpose(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle+180)*Geom_Point_BD');
    now_Geom2_Point_LE       =transpose(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle+180)*Geom_Point_LE');
    now_Geom2_Point_TE       =transpose(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle+180)*Geom_Point_TE');

    now_Panel2_Point_Colocation  = transpose(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle+180)*Panel_Point_Colocation');  %[x y z] coordinate
    now_Panel2_Point_BD         =transpose(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle+180)*Panel_Point_BD');
    now_Panel2_Point_LE       =transpose(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle+180)*Panel_Point_LE');
    now_Panel2_Point_TE       =transpose(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle+180)*Panel_Point_TE');


    %% Calculation IJ Matrix


    for ColloIDX=1:size(now_Panel_Point_Colocation,1) %Collocation Point Loop
        ColoP=now_Panel_Point_Colocation(ColloIDX,:);
        ColoP2=now_Panel2_Point_Colocation(ColloIDX,:);

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

             gamma=gamma_Panel_ini(gamma_idx);

            % induced velocity by Bound Vortex
            vind_Bound(gamma_idx,:)     =Vortex_Vatistas(now_Geom2_Point_BD(gamma_idx,:),now_Geom2_Point_BD(gamma_idx+1,:),ColoP2,gamma,rc_panel(gamma_idx),vortex_n);
            vind_trLeft(gamma_idx,:)    =Vortex_Vatistas(now_Geom2_Point_TE(gamma_idx,:),now_Geom2_Point_BD(gamma_idx,:),ColoP2,gamma,rc_panel(gamma_idx),vortex_n);
            vind_trRight(gamma_idx,:)   =Vortex_Vatistas(now_Geom2_Point_BD(gamma_idx+1,:),now_Geom2_Point_TE(gamma_idx+1,:),ColoP2,gamma,rc_panel(gamma_idx),vortex_n);
            vTotal=vind_Bound(gamma_idx,:)+ vind_trLeft(gamma_idx,:)+ vind_trRight(gamma_idx,:);

            I_ini2_x(ColloIDX,gamma_idx)=vTotal(1);
            I_ini2_y(ColloIDX,gamma_idx)=vTotal(2);
            I_ini2_z(ColloIDX,gamma_idx)=vTotal(3);
        end


    end
    %% Wake Induced Velocity

    for ColloIDX=1:size(now_Panel_Point_Colocation,1) %Collocation Point Loop
        ColoP=now_Panel_Point_Colocation(ColloIDX,:);
        ColoP2=now_Panel2_Point_Colocation(ColloIDX,:);
        vind_Wake(ColloIDX,:)=[0,0,0];
        vind_Wake2(ColloIDX,:)=[0,0,0];
        for idx=1:size(Wake_Gamma,2)
            xind=(idx-1)*3+1;
            yind=(idx-1)*3+2;
            zind=(idx-1)*3+3;
            if size(Wake_Geom_Position,1)<=2; break; end
            for Tidx=1:size(Wake_Geom_Position,1)-1
                %Blade [Wake1->Blade1]
                A=Wake_Geom_Position(Tidx,xind:zind);
                B=Wake_Geom_Position(Tidx+1,xind:zind);
                gamma=Wake_Gamma(Tidx,idx);
                vind_Wake(ColloIDX,:)=vind_Wake(ColloIDX,:)+Vortex_Vatistas(A,B,ColoP,gamma,rc_panel(gamma_idx),vortex_n);
                %Blade [Wake2->Blade1]
                A=Wake2_Geom_Position(Tidx,xind:zind);
                B=Wake2_Geom_Position(Tidx+1,xind:zind);
                gamma=Wake2_Gamma(Tidx,idx);
                vind_Wake(ColloIDX,:)=vind_Wake(ColloIDX,:)+Vortex_Vatistas(A,B,ColoP,gamma,rc_panel(gamma_idx),vortex_n);

                %Blade [Wake1->Blade2]
                A=Wake_Geom_Position(Tidx,xind:zind);
                B=Wake_Geom_Position(Tidx+1,xind:zind);
                gamma=Wake_Gamma(Tidx,idx);
                vind_Wake2(ColloIDX,:)=vind_Wake2(ColloIDX,:)+Vortex_Vatistas(A,B,ColoP2,gamma,rc_panel(gamma_idx),vortex_n);
                %Blade [Wake2->Blade2]
                A=Wake2_Geom_Position(Tidx,xind:zind);
                B=Wake2_Geom_Position(Tidx+1,xind:zind);
                gamma=Wake2_Gamma(Tidx,idx);
                vind_Wake2(ColloIDX,:)=vind_Wake2(ColloIDX,:)+Vortex_Vatistas(A,B,ColoP2,gamma,rc_panel(gamma_idx),vortex_n);

            end
        end
    end


    %% Blade Load Model
    Told=999;
    T1old=999;
    T2old=999;
    unit_ground_Span=[1,0,0]';
    unit_ground_Chord=[0,1,0]';
    unit_ground_Axis=[0,0,1]';

    % Axis Position Coordinate

    unit_conv_Span=(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle)*unit_ground_Span);
    unit_conv_Chord=(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle)*unit_ground_Chord);
    unit_conv_Axis=(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle)*unit_ground_Axis);

    unit_conv2_Span=(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle+180)*unit_ground_Span);
    unit_conv2_Chord=(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle+180)*unit_ground_Chord);
    unit_conv2_Axis=(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle+180)*unit_ground_Axis);

    vinduced=[];

    %gamma=0.5Veff*cl
    for iter=1:1000
        %% Blade 1
        vinduced1=[];
        V_inflow1=[];
        T1=0;
        P1=0;
        Q1=0;
        F1=0;
        v_infow_induced_x_set=I_ini_x*gamma_BoundVortex;
        v_infow_induced_y_set=I_ini_y*gamma_BoundVortex;
        v_infow_induced_z_set=I_ini_z*gamma_BoundVortex;

        Data_local_raidus1=[];
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
            vinduced1=[vinduced1;V_inflow_sum_span,V_inflow_sum_chord,V_inflow_sum_axis];
            V_inflow1=[ V_inflow1;vFree];

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
            CD=inp_Cd;
            

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
            kfac=0.01;
            gamma_BoundVortex(ridx)=gamma_BoundVortex(ridx)+(-gamma_BoundVortex(ridx)+boundVortex)*kfac;

            Data_local_raidus1=[Data_local_raidus1;r./R,r,Flow_Angle,alpha,Re,beta,boundVortex,Tdr];

            T1=T1+Tdr;
            Q1=Q1+Qdr;
            F1=F1+Fdr;
        end
        %% Blade 2
        vinduced2=[];
        V_inflow2=[];
        T2=0;
        P2=0;
        Q2=0;
        F2=0;
        v_infow_induced_x_set=I_ini2_x*gamma_BoundVortex2;
        v_infow_induced_y_set=I_ini2_y*gamma_BoundVortex2;
        v_infow_induced_z_set=I_ini2_z*gamma_BoundVortex2;

        Data_local_raidus2=[];
        for ridx=1:size(now_Panel_Point_Colocation,1)


            % Blade Setting
            V_inflow_free_Axis = dot(vFree,unit_conv2_Axis);
            V_inflow_free_span = dot(vFree,unit_conv2_Span);
            V_inflow_free_chord = dot(vFree,unit_conv2_Chord);
            
            % Geomtry Setting
            
            r=Panel_r(ridx);

            V_inflow_rot_Axis=0;
            V_inflow_rot_span=0;
            V_inflow_rot_chord=-(rotSpeed.*r);

            v_infow_induced_x=v_infow_induced_x_set(ridx);
            v_infow_induced_y=v_infow_induced_y_set(ridx);
            v_infow_induced_z=v_infow_induced_z_set(ridx);

            v_infow_induced=[v_infow_induced_x,v_infow_induced_y,v_infow_induced_z];
            V_inflow_induced_Axis = dot(v_infow_induced,unit_conv2_Axis);
            V_inflow_induced_span = dot(v_infow_induced,unit_conv2_Span);
            V_inflow_induced_chord = dot(v_infow_induced,unit_conv2_Chord);

            v_infow_WakeInduced= vind_Wake2(ridx,:);
            V_inflow_WakeInduced_Axis = dot(v_infow_WakeInduced,unit_conv2_Axis);
            V_inflow_WakeInduced_Span = dot(v_infow_WakeInduced,unit_conv2_Span);
            V_inflow_WakeInduced_Chord = dot(v_infow_WakeInduced,unit_conv2_Chord);


            V_induced_sum_axis=V_inflow_rot_Axis+V_inflow_induced_Axis+V_inflow_WakeInduced_Axis;
            V_induced_sum_span=V_inflow_rot_span+V_inflow_induced_span+V_inflow_WakeInduced_Span;
            V_induced_sum_chord=V_inflow_rot_chord+V_inflow_induced_chord+V_inflow_WakeInduced_Chord;


            V_inflow_sum_axis=V_inflow_free_Axis+V_induced_sum_axis;
            V_inflow_sum_span=V_inflow_free_span+V_induced_sum_span;
            V_inflow_sum_chord=V_inflow_free_chord+V_induced_sum_chord;
            vinduced2=[vinduced2;V_inflow_sum_span,V_inflow_sum_chord,V_inflow_sum_axis];
            V_inflow2= [V_inflow2;vFree];

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
            CD=inp_Cd;

            
            LD = Cl / CD;
            Gam=rad2deg(atan(1 / LD));
            Ks=0.5*rho*(V_inflow_sum_axis.^2)*Cl/((sind(Flow_Angle).^2)*cosd(Gam));

            gamma = rad2deg(atan(1/LD));
            K=Cl.*c./((sin(deg2rad(Flow_Angle)).^2).*cos(deg2rad(gamma)));
            Tc=K.*cos(deg2rad(Flow_Angle+gamma));
            Qc=r.*K.*sin(deg2rad(Flow_Angle+gamma));
            Fc=K.*sin(deg2rad(Flow_Angle+gamma));


            Tds=Ks*cosd(Gam+Flow_Angle)*ds;
            Qds=r*Ks*sind(Gam+Flow_Angle)*ds;
            Fds=-Ks*sind(Gam+Flow_Angle)*ds;

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
            kfac=0.01;
            gamma_BoundVortex2(ridx)=gamma_BoundVortex2(ridx)+(-gamma_BoundVortex2(ridx)+boundVortex)*kfac;

            Data_local_raidus2=[Data_local_raidus2;r./R,r,Flow_Angle,alpha,Re,beta,boundVortex,Tdr];

            T2=T2+Tds;
            Q2=Q2+Qds;
            F2=F2+Fds;
        end

        T=T1+T2;
        Q=Q1+Q2;


        if abs(T1old-T1)<0.0001 & abs(T2old-T2)<0.0001 
            fprintf("수렴! T= %.3fN  Q= %.3fNm\n",T,Q)
            break;
        else

        end
        Told=T;
        T1old=T1;
        T2old=T2;

    end


    %% Wake Structure (Semi-Prescribed Wake Sturcture)






    %TE -> Wake Structure 추가
        %Blade 1
    tmp=now_Geom_Point_TE';
    wakeForm=tmp(:)';
    vinduced_Ax=vinduced1(:,3);
    vinduced_Ax_set=transpose(vinduced_Ax.*unit_conv_Axis'+V_inflow1);
    vinduced_Ax_GeomSet=[];
    vinduced_Ax_GeomSet=0.5.*(vinduced_Ax_set(:,2:end)+vinduced_Ax_set(:,1:end-1));
    vinduced_Ax_GeomSet=[vinduced_Ax_set(:,1), vinduced_Ax_GeomSet, vinduced_Ax_set(:,end)];
    vinduced_Ax_GeomSet=vinduced_Ax_GeomSet(:);
    WakeGammaset=[0;gamma_BoundVortex]-[gamma_BoundVortex;0];
    %WakeGammaset(1:5)=0;
    Wake_Geom_Position=[wakeForm;Wake_Geom_Position]; %x y z x y z
    Wake_velocity=[vinduced_Ax_GeomSet';Wake_velocity];
    

    %Blade 2
    tmp=now_Geom2_Point_TE';
    wakeForm2=tmp(:)';
    vinduced2_Ax=vinduced2(:,3);
    vinduced2_Ax_set=transpose(vinduced2_Ax.*unit_conv2_Axis'+V_inflow2);
    vinduced2_Ax_GeomSet=[];
    vinduced2_Ax_GeomSet=0.5.*(vinduced2_Ax_set(:,2:end)+vinduced2_Ax_set(:,1:end-1));
    vinduced2_Ax_GeomSet=[vinduced2_Ax_set(:,1), vinduced2_Ax_GeomSet, vinduced2_Ax_set(:,end)];
    vinduced2_Ax_GeomSet=vinduced2_Ax_GeomSet(:);
    WakeGammaset2=[0;gamma_BoundVortex2]-[gamma_BoundVortex2;0];
    %WakeGammaset2(1:5)=0;
    Wake2_Geom_Position=[wakeForm2;Wake2_Geom_Position]; %x y z x y z
    Wake2_velocity=[vinduced2_Ax_GeomSet';Wake2_velocity];
    
    Wake_Gamma=[WakeGammaset';Wake_Gamma];
    Wake2_Gamma=[WakeGammaset2';Wake2_Gamma];


    %% Wake - Wake Interaction
    V_wake_1=zeros(size(Wake_Geom_Position));
    V_wake_2=zeros(size(Wake_Geom_Position));

    if (Total_Rotate>180)
    vOutVel = WakeVortexCalculator_mex_Fast([Wake_Geom_Position Wake2_Geom_Position], [Wake_Gamma Wake2_Gamma],[rc_Geom rc_Geom]);
    V_wake_1=vOutVel(:,1:size(vOutVel,2)/2);
    V_wake_2=vOutVel(:,1+size(vOutVel,2)/2:end);
    else
        Wake_Gamma=Wake_Gamma.*0;
        Wake2_Gamma=Wake2_Gamma.*0;
    end
        





%% Blade-Wake Interaction
% sizeofWakemetrx=size(Wake_Geom_Position);
% Idx=1;
%     while Idx<sizeofWakemetrx(1)*sizeofWakemetrx(2)
%         %Wake, geom에 따라 하나씩 계산
%         wakeIdx=fix(Idx/sizeofWakemetrx(2))+1;
%         geomIdx=(mod(Idx,sizeofWakemetrx(2))-1)/3+1;
%         xind=3*(geomIdx-1)+1;
%         yind=3*(geomIdx-1)+2;
%         zind=3*(geomIdx-1)+3;
%         WakeMarker=Wake_Geom_Position(wakeIdx,xind:zind);
%         WakeMarker2=Wake2_Geom_Position(wakeIdx,xind:zind);
% 
%         for gamma_idx=1:size(gamma_BoundVortex,1)
%             %Blade Vortex로인한 후류 후퇴속도 계산
%             gamma=gamma_BoundVortex(gamma_idx);
%             gamma2=gamma_BoundVortex2(gamma_idx);
%              % induced velocity by Bound Vortex
%             vind1_W1_Bound(gamma_idx,:)     =Vortex_Vatistas(now_Geom_Point_BD(gamma_idx,:),now_Geom_Point_BD(gamma_idx+1,:),WakeMarker,gamma,rc_panel(gamma_idx),vortex_n);
%             vind1_W1_trLeft(gamma_idx,:)    =Vortex_Vatistas(now_Geom_Point_TE(gamma_idx,:),now_Geom_Point_BD(gamma_idx,:),WakeMarker,gamma,rc_panel(gamma_idx),vortex_n);
%             vind1_W1_trRight(gamma_idx,:)   =Vortex_Vatistas(now_Geom_Point_BD(gamma_idx+1,:),now_Geom_Point_TE(gamma_idx+1,:),WakeMarker,gamma,rc_panel(gamma_idx),vortex_n);
%             vind1_W2_Bound(gamma_idx,:)     =Vortex_Vatistas(now_Geom_Point_BD(gamma_idx,:),now_Geom_Point_BD(gamma_idx+1,:),WakeMarker2,gamma2,rc_panel(gamma_idx),vortex_n);
%             vind1_W2_trLeft(gamma_idx,:)    =Vortex_Vatistas(now_Geom_Point_TE(gamma_idx,:),now_Geom_Point_BD(gamma_idx,:),WakeMarker2,gamma2,rc_panel(gamma_idx),vortex_n);
%             vind1_W2_trRight(gamma_idx,:)   =Vortex_Vatistas(now_Geom_Point_BD(gamma_idx+1,:),now_Geom_Point_TE(gamma_idx+1,:),WakeMarker2,gamma2,rc_panel(gamma_idx),vortex_n);
% 
% 
% 
%             % induced velocity by Bound Vortex
%             vind2_W1_Bound(gamma_idx,:)     =Vortex_Vatistas(now_Geom2_Point_BD(gamma_idx,:),now_Geom2_Point_BD(gamma_idx+1,:),WakeMarker,gamma,rc_panel(gamma_idx),vortex_n);
%             vind2_W1_trLeft(gamma_idx,:)    =Vortex_Vatistas(now_Geom2_Point_TE(gamma_idx,:),now_Geom2_Point_BD(gamma_idx,:),WakeMarker,gamma,rc_panel(gamma_idx),vortex_n);
%             vind2_W1_trRight(gamma_idx,:)   =Vortex_Vatistas(now_Geom2_Point_BD(gamma_idx+1,:),now_Geom2_Point_TE(gamma_idx+1,:),WakeMarker,gamma,rc_panel(gamma_idx),vortex_n);
%             vind2_W2_Bound(gamma_idx,:)     =Vortex_Vatistas(now_Geom2_Point_BD(gamma_idx,:),now_Geom2_Point_BD(gamma_idx+1,:),WakeMarker2,gamma2,rc_panel(gamma_idx),vortex_n);
%             vind2_W2_trLeft(gamma_idx,:)    =Vortex_Vatistas(now_Geom2_Point_TE(gamma_idx,:),now_Geom2_Point_BD(gamma_idx,:),WakeMarker2,gamma2,rc_panel(gamma_idx),vortex_n);
%             vind2_W2_trRight(gamma_idx,:)   =Vortex_Vatistas(now_Geom2_Point_BD(gamma_idx+1,:),now_Geom2_Point_TE(gamma_idx+1,:),WakeMarker2,gamma2,rc_panel(gamma_idx),vortex_n);
% 
% 
%         end
% 
%         Vind_w1=sum([vind1_W1_Bound;vind2_W1_Bound])...
%             +sum([vind1_W1_trLeft;vind2_W1_trLeft])...
%             +sum([vind1_W1_trRight;vind2_W1_trRight]);
% 
%         Vind_w2=sum([vind1_W2_Bound;vind2_W2_Bound])...
%             +sum([vind1_W2_trLeft;vind2_W2_trLeft])...
%             +sum([vind1_W2_trRight;vind2_W2_trRight]);
% 
%         V_wake_1(wakeIdx,xind:zind)= V_wake_1(wakeIdx,xind:zind)+Vind_w1;
%         V_wake_2(wakeIdx,xind:zind)= V_wake_2(wakeIdx,xind:zind)+Vind_w2;
% 
%         Idx=Idx+3;
% 
% 
%     end




    Wake_Geom_Position=Wake_Geom_Position+(Wake_velocity+V_wake_1).*dt;
    Wake2_Geom_Position=Wake2_Geom_Position+(Wake2_velocity+V_wake_2).*dt;

    


    %Trail Gamma 계산
    %shed Gamma 계산

    %% Rotate Next Structure
    Total_Rotate=Total_Rotate+dAngle;
    AzimuthAngle=AzimuthAngle+dAngle;
    if AzimuthAngle>360
        AzimuthAngle=AzimuthAngle-360;
        rotateCount=rotateCount+1;
    end

    Time=Time+dt;
    if rotateCount==1
     Wake_Geom_Position(end,:)=[]; %x y z x y z
     Wake_velocity(end,:)=[];
     Wake_Gamma(end,:)=[];

     Wake2_Geom_Position(end,:)=[]; %x y z x y z
     Wake2_velocity(end,:)=[];
     Wake2_Gamma(end,:)=[];
     end


     % if rotateCount==360
     %     size(Wake_Geom_Position,1)-int16(180/dAngle);
     % Wake_Geom_Position(end-int16(180/dAngle):end,:)=[]; %x y z x y z
     % Wake_velocity(end-int16(180/dAngle):end,:)=[];
     % Wake_Gamma(end-int16(180/dAngle):end,:)=[];
     % 
     % Wake2_Geom_Position(end-int16(180/dAngle):end,:)=[]; %x y z x y z
     % Wake2_velocity(end-int16(180/dAngle):end,:)=[];
     % Wake2_Gamma(end-int16(180/dAngle):end,:)=[];
     % end
     % 
     % if Total_Rotate==720
     %     size(Wake_Geom_Position,1)-int16(180/dAngle);
     %    Wake_Geom_Position(end-int16(180/dAngle):end,:)=[]; %x y z x y z
     % Wake_velocity(end-int16(180/dAngle):end,:)=[];
     % Wake_Gamma(end-int16(180/dAngle):end,:)=[];
     % 
     % Wake2_Geom_Position(end-int16(180/dAngle):end,:)=[]; %x y z x y z
     % Wake2_velocity(end-int16(180/dAngle):end,:)=[];
     % Wake2_Gamma(end-int16(180/dAngle):end,:)=[];
     % end

% 
%     figure(100)
%     clf
%     hold on
% 
%     plot3(now_Panel_Point_Colocation(:,1),now_Panel_Point_Colocation(:,2),now_Panel_Point_Colocation(:,3),'ro-')
%     plot3(now_Panel_Point_BD(:,1),now_Panel_Point_BD(:,2),now_Panel_Point_BD(:,3),'bx-')
%     plot3(now_Panel_Point_LE(:,1),now_Panel_Point_LE(:,2),now_Panel_Point_LE(:,3),'g*-')
%     plot3(now_Panel_Point_TE(:,1),now_Panel_Point_TE(:,2),now_Panel_Point_TE(:,3),'k*-')
%     for idx=int16(length(WakeGammaset)/4):length(WakeGammaset)
%         xind=(idx-1)*3+1;
%         yind=(idx-1)*3+2;
%         zind=(idx-1)*3+3;
% 
%         plot3(Wake_Geom_Position(1:end,xind),Wake_Geom_Position(1:end,yind),Wake_Geom_Position(1:end,zind),'k:')
%     end
% 
%     plot3(now_Panel2_Point_Colocation(:,1),now_Panel2_Point_Colocation(:,2),now_Panel2_Point_Colocation(:,3),'ro-')
%     plot3(now_Panel2_Point_BD(:,1),now_Panel2_Point_BD(:,2),now_Panel2_Point_BD(:,3),'bx-')
%     plot3(now_Panel2_Point_LE(:,1),now_Panel2_Point_LE(:,2),now_Panel2_Point_LE(:,3),'g*-')
%     plot3(now_Panel2_Point_TE(:,1),now_Panel2_Point_TE(:,2),now_Panel2_Point_TE(:,3),'k*-')
%     for idx=int16(length(WakeGammaset)/4):length(WakeGammaset)
%         xind=(idx-1)*3+1;
%         yind=(idx-1)*3+2;
%         zind=(idx-1)*3+3;
% 
%         plot3(Wake2_Geom_Position(1:end,xind),Wake2_Geom_Position(1:end,yind),Wake2_Geom_Position(1:end,zind),'k:')
%     end
%             %plot3(Wake2_Geom_Position(1:end,xind),Wake2_Geom_Position(1:end,yind),Wake2_Geom_Position(1:end,zind),'r:')
%         %plot3(Wake_Geom_Position(1:end,xind),Wake_Geom_Position(1:end,yind),Wake_Geom_Position(1:end,zind),'r:')
% 
%         xlabel("x (m)")
%         ylabel("y (m)")
%         zlabel("z (m)")
%         title(["Propeller & Wake Structure ";"(semi-Free Wake method)"])
%         view([1,1,1])
%         axis equal
% 
%         filename=sprintf("FIG_%05dRPM_%04dANGLE.png",RPM,Total_Rotate)
%         exportgraphics(figure(100),'fig1/'+filename,'Resolution',300)
%         view(unit_conv2_Chord')
%         axis equal 
% 
% 
%         filename=sprintf("FIG_%05dRPM_%04dANGLE.png",RPM,Total_Rotate)
%         exportgraphics(figure(100),'fig2/'+filename,'Resolution',300)
% view([1,1,1])
         globalData=[globalData;Time,Total_Rotate,AzimuthAngle, T,T1,T2];

        figure(6)
    clf
    hold on
    plot(globalData(:,2),globalData(:,4),'k')
     plot(globalData(:,2),globalData(:,5),'r')

      plot(globalData(:,2),globalData(:,6),'b')
      xlabel("Total Rotation Angle (deg)")
      ylabel("Thrust(N)")
      title("Thrust Calculation")
      legend("Total Thrust", "Blade 1 Thrust", "Blade 2 Thrust")


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

figure(4)
clf
hold on
plot(Data_local_raidus1(:,1),Data_local_raidus1(:,3),'r') % Flow
plot(Data_local_raidus1(:,1),Data_local_raidus1(:,4),'g') % alpha
plot(Data_local_raidus1(:,1),Data_local_raidus1(:,6),'b') %beta
legend("Flow","Alpha","beta")


figure(5)
clf
hold on
plot(Data_local_raidus1(:,1),Data_local_raidus1(:,8),'r') % Flow
clc
Clac_time=toc(GlobalTimer);
fprintf("Thrust : %.3f N\n",T)
fprintf("Calculate Time : %.3f s\n",Clac_time)
