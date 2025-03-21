

clc
clear
close all
load("InterpolatedModel_ncrit7.mat")
hub=0;


RPM	        =   8000    ;  %   rev/min
method      =   0; % method 0 is local, 1 is mean Local
rc_Ratio    =  0.6;
vortex_n    =   1.06;


alt     = 0 ;         %   m
[Tmp, Pressure, rho, D_vis, a] = STD_Atm(alt);
rho=1.225;

caseName    =   "A1_7_1.06_Local";


dAngle      =   10; %deg
Full_Logfilename=sprintf("[LOG]_%s.txt",caseName);
fid=fopen(Full_Logfilename,'w');
fprintf(fid,"\n");
fclose(fid);

Calc_case=[
    RPM,0,0,0;
    RPM,0,0,-3;
    RPM,0,0,-5;
    RPM,0,0,-7;
    RPM,0,0,-10;
    RPM,0,0,-12;
    RPM,0,0,-15;
    RPM,0,0,-17;
    RPM,0,0,-20;];

Calc_data=[];
mkdir 'OutputData'
for CaseInd=1:size(Calc_case,1)
    clc
    close all
    % set propeller
    RPM	        =   Calc_case(CaseInd,1);  %   rev/min
    vFree       =   Calc_case(CaseInd,2:4);
    CaseTitle =   sprintf("%s_%d_%d_%d_%d",caseName,RPM,vFree);
    GlobalTimer =   tic;
    Logfilename=CaseTitle+"_log.txt";
    fid=fopen(Logfilename,'w');
    fprintf(fid,"%s Log\n",CaseTitle);
    fprintf(fid,"Total_Rotate, CaseInd, method, vortex_n, rc_Ratio, Calc_case(CaseInd), Tav,Qav,Pav,FxAver,FyAver,FzAver \n");

    %% Read Aerodata
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
    Avarge_globalData=[];
    %% Set Operation Condition

    % Setup Condition

    D	    = R.*2;       %   m
    %V      = eps;



    % Calc variable

    nAzmuth=2;

    Blade   = 2;
    DiskArea=pi*R*R;

    rotSpeed=2*pi*RPM/60;  %rad/s
    rot_deg_speed=360*RPM/60; %deg/s
    dt=dAngle./rot_deg_speed;
    %initial Setting
    RotateTotalAngle=4000;
    initialAngle=0;


    n       =RPM./60;
    %J	    =V./(n.*D);
    V_tip	=2*pi*n.*R;
    M_tip	=V_tip./340;
    Afan    =pi*R*R;

    OperationData={vFree,alt,Tmp Pressure, rho, D_vis, a, RotateTotalAngle,dAngle,R,D,Afan};

    figure(100)
    set( gcf, 'Position', [0 0 1800 960] )

    mkdir 'fig1'
    mkdir 'fig2'

    Blade1Alpha =[];
    Blade1Re    =[];
    Blade1Gamma =[];
    Blade1dT    =[];
    Blade1dF    =[];
    Blade1dQ    =[];
    Blade1Cl    =[];
    Blade1Cd    =[];

    Blade2Alpha =[];
    Blade2Re    =[];
    Blade2Gamma =[];
    Blade2dT    =[];
    Blade2dF    =[];
    Blade2dQ    =[];
    Blade2Cl    =[];
    Blade2Cd    =[];
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
    Panel_dr2=abs((Geom_r(1:end-1,1).^2-Geom_r(2:end,1).^2));

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

    switch method
        case 0
            rc_panel=Panel_chord.*rc_Ratio;
            rc_Geom=Geom_chord.*rc_Ratio;
        case 1
            rc_panel=ones(size(Panel_chord)).*mean(Panel_chord).*rc_Ratio
            rc_Geom=ones(size(Geom_chord)).*mean(Geom_chord).*rc_Ratio
    end

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
    ini_Wake_velocity=[];
    Wake_Gamma=[];

    Wake2_Geom_Position=[]; %x y z x y z
    ini_Wake2_velocity=[];
    Wake2_Gamma=[];
    globalData=[];
    rotateCount=0;

    TotalWake1=zeros(size(Wake_Geom_Position));
    TotalWake2=zeros(size(Wake2_Geom_Position));
    dT_R1_old=999;
    dT_R1=999;
    dT_R2_old=999;
    dT_R2=999;
    while Total_Rotate<RotateTotalAngle


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
        timer_1=tic;

        for ColloIDX=1:size(now_Panel_Point_Colocation,1) %Collocation Point Loop
            ColoP=now_Panel_Point_Colocation(ColloIDX,:);
            ColoP2=now_Panel2_Point_Colocation(ColloIDX,:);

            for gamma_idx=1:length(gamma_Panel_ini) %Bound GaMMA Loop

                %Vortex Modeling
                rc_panel_IJ=rc_panel;

                gamma=gamma_Panel_ini(gamma_idx);
                %Vortex_bio;
                %Vortex_Vatistas;
                % induced velocity by Bound Vortex
                vind_Bound(gamma_idx,:)     =Vortex_Vatistas(now_Geom_Point_BD(gamma_idx,:),now_Geom_Point_BD(gamma_idx+1,:),ColoP,gamma,rc_panel_IJ(gamma_idx),vortex_n);
                vind_trLeft(gamma_idx,:)    =Vortex_Vatistas(now_Geom_Point_TE(gamma_idx,:),now_Geom_Point_BD(gamma_idx,:),ColoP,gamma,rc_panel_IJ(gamma_idx),vortex_n);
                vind_trRight(gamma_idx,:)   =Vortex_Vatistas(now_Geom_Point_BD(gamma_idx+1,:),now_Geom_Point_TE(gamma_idx+1,:),ColoP,gamma,rc_panel_IJ(gamma_idx),vortex_n);
                vTotal=vind_Bound(gamma_idx,:)+ vind_trLeft(gamma_idx,:)+ vind_trRight(gamma_idx,:);

                gamma=gamma_Panel_ini(gamma_idx);

                % induced velocity by Bound Vortex
                vind_Bound(gamma_idx,:)     =Vortex_Vatistas(now_Geom2_Point_BD(gamma_idx,:),now_Geom2_Point_BD(gamma_idx+1,:),ColoP,gamma,rc_panel_IJ(gamma_idx),vortex_n);
                vind_trLeft(gamma_idx,:)    =Vortex_Vatistas(now_Geom2_Point_TE(gamma_idx,:),now_Geom2_Point_BD(gamma_idx,:),ColoP,gamma,rc_panel_IJ(gamma_idx),vortex_n);
                vind_trRight(gamma_idx,:)   =Vortex_Vatistas(now_Geom2_Point_BD(gamma_idx+1,:),now_Geom2_Point_TE(gamma_idx+1,:),ColoP,gamma,rc_panel_IJ(gamma_idx),vortex_n);
                vTotal=vTotal+vind_Bound(gamma_idx,:)+ vind_trLeft(gamma_idx,:)+ vind_trRight(gamma_idx,:);

                I_ini_x(ColloIDX,gamma_idx)=vTotal(1);
                I_ini_y(ColloIDX,gamma_idx)=vTotal(2);
                I_ini_z(ColloIDX,gamma_idx)=vTotal(3);


                gamma=gamma_Panel_ini(gamma_idx);

                % induced velocity by Bound Vortex
                vind_Bound(gamma_idx,:)     =Vortex_Vatistas(now_Geom_Point_BD(gamma_idx,:),now_Geom_Point_BD(gamma_idx+1,:),ColoP2,gamma,rc_panel_IJ(gamma_idx),vortex_n);
                vind_trLeft(gamma_idx,:)    =Vortex_Vatistas(now_Geom_Point_TE(gamma_idx,:),now_Geom_Point_BD(gamma_idx,:),ColoP2,gamma,rc_panel_IJ(gamma_idx),vortex_n);
                vind_trRight(gamma_idx,:)   =Vortex_Vatistas(now_Geom_Point_BD(gamma_idx+1,:),now_Geom_Point_TE(gamma_idx+1,:),ColoP2,gamma,rc_panel_IJ(gamma_idx),vortex_n);
                vTotal=vind_Bound(gamma_idx,:)+ vind_trLeft(gamma_idx,:)+ vind_trRight(gamma_idx,:);


                gamma=gamma_Panel_ini(gamma_idx);

                % induced velocity by Bound Vortex
                vind_Bound(gamma_idx,:)     =Vortex_Vatistas(now_Geom2_Point_BD(gamma_idx,:),now_Geom2_Point_BD(gamma_idx+1,:),ColoP2,gamma,rc_panel_IJ(gamma_idx),vortex_n);
                vind_trLeft(gamma_idx,:)    =Vortex_Vatistas(now_Geom2_Point_TE(gamma_idx,:),now_Geom2_Point_BD(gamma_idx,:),ColoP2,gamma,rc_panel_IJ(gamma_idx),vortex_n);
                vind_trRight(gamma_idx,:)   =Vortex_Vatistas(now_Geom2_Point_BD(gamma_idx+1,:),now_Geom2_Point_TE(gamma_idx+1,:),ColoP2,gamma,rc_panel_IJ(gamma_idx),vortex_n);
                vTotal=vTotal+vind_Bound(gamma_idx,:)+ vind_trLeft(gamma_idx,:)+ vind_trRight(gamma_idx,:);

                I_ini2_x(ColloIDX,gamma_idx)=vTotal(1);
                I_ini2_y(ColloIDX,gamma_idx)=vTotal(2);
                I_ini2_z(ColloIDX,gamma_idx)=vTotal(3);
            end


        end
        time_IJcalc=toc(timer_1);
        timer_1=tic;

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
                gamma_idx=idx;
                if size(Wake_Geom_Position,1)<=2; break; end
                for Tidx=1:size(Wake_Geom_Position,1)-1
                    %Blade [Wake1->Blade1]
                    A=Wake_Geom_Position(Tidx,xind:zind);
                    B=Wake_Geom_Position(Tidx+1,xind:zind);
                    gamma=Wake_Gamma(Tidx,idx);
                    vind_Wake(ColloIDX,:)=vind_Wake(ColloIDX,:)+Vortex_Vatistas(A,B,ColoP,gamma,rc_Geom(gamma_idx),vortex_n);
                    %Blade [Wake2->Blade1]
                    A=Wake2_Geom_Position(Tidx,xind:zind);
                    B=Wake2_Geom_Position(Tidx+1,xind:zind);
                    gamma=Wake2_Gamma(Tidx,idx);
                    vind_Wake(ColloIDX,:)=vind_Wake(ColloIDX,:)+Vortex_Vatistas(A,B,ColoP,gamma,rc_Geom(gamma_idx),vortex_n);

                    %Blade [Wake1->Blade2]
                    A=Wake_Geom_Position(Tidx,xind:zind);
                    B=Wake_Geom_Position(Tidx+1,xind:zind);
                    gamma=Wake_Gamma(Tidx,idx);
                    vind_Wake2(ColloIDX,:)=vind_Wake2(ColloIDX,:)+Vortex_Vatistas(A,B,ColoP2,gamma,rc_Geom(gamma_idx),vortex_n);
                    %Blade [Wake2->Blade2]
                    A=Wake2_Geom_Position(Tidx,xind:zind);
                    B=Wake2_Geom_Position(Tidx+1,xind:zind);
                    gamma=Wake2_Gamma(Tidx,idx);
                    vind_Wake2(ColloIDX,:)=vind_Wake2(ColloIDX,:)+Vortex_Vatistas(A,B,ColoP2,gamma,rc_Geom(gamma_idx),vortex_n);

                end
            end
        end

        time_WakeInduced=toc(timer_1);
        timer_1=tic;
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


        dT_R1_old=999;
        dT_R2_old=999;
        %gamma=0.5Veff*cl
        for iter=1:1000
            %% Blade 1
            dT_R1=[];
            dT_R2=[];


            vinduced1=[];
            V_inflow1=[];
            T1=0;
            P1=0;
            Q1=0;
            F1=0;
            momentumTdr=0;
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


                % Calculate Alpha_Wake
                V_FlowInduced_Axis=V_inflow_free_Axis+V_inflow_rot_Axis+V_inflow_WakeInduced_Axis;
                V_FlowInduced_Span=V_inflow_free_span+V_inflow_rot_span+V_inflow_WakeInduced_Span;
                V_FlowInduced_Chord=V_inflow_free_chord+V_inflow_rot_chord+V_inflow_WakeInduced_Chord;


                % Calculate Alpha_ind
                V_BladeInduced_Axis=V_inflow_induced_Axis;
                V_BladeInduced_Span=V_inflow_induced_span;
                V_BladeInduced_Chord=V_inflow_induced_chord;

                Flow_Angle=rad2deg(atan2(-V_FlowInduced_Axis,-V_FlowInduced_Chord));
                Flow_Angle_ind=rad2deg(atan2(-V_inflow_sum_axis,-V_inflow_sum_chord));




                V_inflow1=[ V_inflow1;vFree];

                Veff=sqrt((V_FlowInduced_Axis.^2+V_FlowInduced_Chord.^2));
                beta=Panel_Beta(ridx);
                rR_local=r./R;
                c=Panel_chord(ridx);
                ds=Panel_ds(ridx);
                dr=Panel_dr(ridx);
                dr2=Panel_dr2(ridx);


                alpha=beta-Flow_Angle;
                alpha_ind=beta-Flow_Angle_ind;

                Re=(rho*Veff*c)./D_vis;
                ThickRatio=Panel_Thick(ridx).*100;


                if BlendPosition(1)>=rR_local
                    airfoil=1;
                    inp_Cl=BLD_CL(0,ThickRatio, Re ,alpha);
                    inp_Cleff=BLD_CL(0,ThickRatio ,Re, alpha_ind);

                    inp_Cd=BLD_CD(0,ThickRatio, Re ,alpha_ind);
                elseif BlendPosition(2)<=rR_local
                    airfoil=2;
                    inp_Cl=BLD_CL(1,ThickRatio ,Re, alpha);
                    inp_Cleff=BLD_CL(1,ThickRatio ,Re, alpha_ind);

                    inp_Cd=BLD_CD(1,ThickRatio, Re, alpha_ind);
                else
                    airfoil=0;
                    BR=abs((rR_local-BlendPosition(1))./(BlendPosition(2)-BlendPosition(1)));
                    %CL0=BLD_CL(BR, Re ,0);
                    inp_Cl=BLD_CL(BR,ThickRatio, Re ,alpha);
                    inp_Cleff=BLD_CL(BR,ThickRatio, Re ,alpha_ind);

                    inp_Cd=BLD_CD(BR ,ThickRatio,Re, alpha_ind);
                end

                Cl=inp_Cleff;
                ai=Flow_Angle_ind-Flow_Angle;
                CD=inp_Cd+inp_Cl.*sind(ai);


                LD=Cl/CD;
                %Calculate BET
                gamma = rad2deg(atan(1/LD));
                Angle=Flow_Angle;
                if Angle==0; Angle=eps;                end

                K=Cl.*c./((sin(deg2rad(Angle)).^2).*cos(deg2rad(gamma)));
                Tc=K.*cos(deg2rad(Angle+gamma));
                Qc=r.*K.*sin(deg2rad(Angle+gamma));
                Fc=K.*sin(deg2rad(Angle+gamma));

                Tdr=(0.5*rho*V_FlowInduced_Axis.^2).*Tc.*dr;
                Qdr=(0.5*rho*V_FlowInduced_Axis.^2).*Qc.*dr;
                Fdr=(0.5*rho*V_FlowInduced_Axis.^2).*Fc.*dr;

                momentumTdr=momentumTdr+(2*rho*DiskArea*V_induced_sum_axis*V_inflow_sum_axis*(dr2./(R.^2)));


                % figure(10)
                % hold on
                % quiver(-V_inflow_free_chord,-V_inflow_free_Axis,V_inflow_free_chord,V_inflow_free_Axis,1,'r')
                % quiver(-V_inflow_sum_chord,-V_inflow_sum_axis,V_inflow_sum_chord,V_inflow_sum_axis,1,'g')
                % quiver(-V_inflow_rot_chord,0,V_inflow_rot_chord,0,1,'b')
                % grid on
                % axis equal
                %
                boundVortex=0.5*Cl*Veff*c./Blade;
                kfac=0.05;
                gamma_BoundVortex(ridx)=gamma_BoundVortex(ridx)+(-gamma_BoundVortex(ridx)+boundVortex)*kfac;

                Data_local_raidus1=[Data_local_raidus1;r./R,r,Flow_Angle,alpha,Re,beta,boundVortex,Tdr,Fdr,Qdr,Cl,CD];

                T1=T1+Tdr;
                Q1=Q1+Qdr;
                F1=F1+Fdr;

                dT_R1=[dT_R1;Tdr];
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

                % Calculate Alpha_Wake
                V_FlowInduced_Axis=V_inflow_free_Axis+V_inflow_rot_Axis+V_inflow_WakeInduced_Axis;
                V_FlowInduced_Span=V_inflow_free_span+V_inflow_rot_span+V_inflow_WakeInduced_Span;
                V_FlowInduced_Chord=V_inflow_free_chord+V_inflow_rot_chord+V_inflow_WakeInduced_Chord;


                % Calculate Alpha_ind
                V_BladeInduced_Axis=V_inflow_induced_Axis;
                V_BladeInduced_Span=V_inflow_induced_span;
                V_BladeInduced_Chord=V_inflow_induced_chord;

                Flow_Angle=rad2deg(atan2(-V_FlowInduced_Axis,-V_FlowInduced_Chord));
                Flow_Angle_ind=rad2deg(atan2(-V_inflow_sum_axis,-V_inflow_sum_chord));


                Veff=sqrt((V_FlowInduced_Axis.^2+V_FlowInduced_Chord.^2));
                %Veff=sqrt((V_inflow_sum_axis.^2+V_inflow_sum_chord.^2));

                beta=Panel_Beta(ridx);
                rR_local=r./R;
                c=Panel_chord(ridx);
                ds=Panel_ds(ridx);
                dr=Panel_dr(ridx);


                alpha=beta-Flow_Angle;
                alpha_ind=beta-Flow_Angle_ind;


                Re=(rho*Veff*c)./D_vis;
                ThickRatio=Panel_Thick(ridx).*100;


                if BlendPosition(1)>=rR_local
                    airfoil=1;
                    inp_Cl=BLD_CL(0,ThickRatio, Re ,alpha);
                    inp_Cleff=BLD_CL(0,ThickRatio ,Re, alpha_ind);

                    inp_Cd=BLD_CD(0,ThickRatio, Re ,alpha_ind);
                elseif BlendPosition(2)<=rR_local
                    airfoil=2;
                    inp_Cl=BLD_CL(1,ThickRatio ,Re, alpha);
                    inp_Cleff=BLD_CL(1,ThickRatio ,Re, alpha_ind);

                    inp_Cd=BLD_CD(1,ThickRatio, Re, alpha_ind);
                else
                    airfoil=0;
                    BR=abs((rR_local-BlendPosition(1))./(BlendPosition(2)-BlendPosition(1)));
                    %CL0=BLD_CL(BR, Re ,0);
                    inp_Cl=BLD_CL(BR,ThickRatio, Re ,alpha);
                    inp_Cleff=BLD_CL(BR,ThickRatio, Re ,alpha_ind);

                    inp_Cd=BLD_CD(BR ,ThickRatio,Re, alpha_ind);
                end

                Cl=inp_Cleff;
                ai=Flow_Angle_ind-Flow_Angle;

                CD=inp_Cd+inp_Cl.*sind(ai);
                LD = Cl / CD;
                gamma = rad2deg(atan(1/LD));


                Angle=Flow_Angle;
                if Angle==0; Angle=eps;                end
                K=Cl.*c./((sin(deg2rad(Angle)).^2).*cos(deg2rad(gamma)));
                Tc=K.*cos(deg2rad(Angle+gamma));
                Qc=r.*K.*sin(deg2rad(Angle+gamma));
                Fc=K.*sin(deg2rad(Angle+gamma));

                Tdr=(0.5*rho*V_FlowInduced_Axis.^2).*Tc.*dr;
                Qdr=(0.5*rho*V_FlowInduced_Axis.^2).*Qc.*dr;
                Fdr=(0.5*rho*V_FlowInduced_Axis.^2).*Fc.*dr;
                % figure(10)
                % hold on
                % quiver(-V_inflow_free_chord,-V_inflow_free_Axis,V_inflow_free_chord,V_inflow_free_Axis,1,'r')
                % quiver(-V_inflow_sum_chord,-V_inflow_sum_axis,V_inflow_sum_chord,V_inflow_sum_axis,1,'g')
                % quiver(-V_inflow_rot_chord,0,V_inflow_rot_chord,0,1,'b')
                % grid on
                % axis equal
                %
                boundVortex=0.5*Cl*Veff*c./Blade;
                kfac=0.05;
                gamma_BoundVortex2(ridx)=gamma_BoundVortex2(ridx)+(-gamma_BoundVortex2(ridx)+boundVortex)*kfac;

                Data_local_raidus2=[Data_local_raidus2;r./R,r,Flow_Angle,alpha,Re,beta,gamma_BoundVortex2(ridx),Tdr,Fdr,Qdr,Cl,CD];
                T2=T2+Tdr;
                Q2=Q2+Qdr;
                F2=F2+Fdr;
                dT_R2=[dT_R2;Tdr];

            end
            % unit_conv_Span=(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle)*unit_ground_Span);
            % unit_conv_Chord=(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle)*unit_ground_Chord);
            % unit_conv_Axis=(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle)*unit_ground_Axis);
            %
            % unit_conv2_Span=(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle+180)*unit_ground_Span);
            % unit_conv2_Chord=(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle+180)*unit_ground_Chord);
            % unit_conv2_Axis=(R_Pitch(Tilit_Angle)*Rz(AzimuthAngle+180)*unit_ground_Axis);
            T=T1+T2;
            Q=Q1+Q2;

            Total_forceVec=T1.*unit_conv_Axis+T2.*unit_conv2_Axis-F1.*unit_conv_Chord-F2*unit_conv2_Chord;
            %
            %         if abs(T1old-T1)<0.00001 & abs(T2old-T2)<0.00001
            %             %
            %             break;
            %
            %
            %
            %
            %
            %         else
            %             %T
            %
            %         end
            cri1=max(abs(dT_R1_old-dT_R1));
            cri2=max(abs(dT_R2_old-dT_R2)) ;
            if cri1<0.00001 & cri2<0.00001

                %
                break;





            else
                dT_R1_old=dT_R1;
                dT_R2_old=dT_R2;
                %disp("수렴하지 않았습니다.");

            end


            Told=T;
            T1old=T1;
            T2old=T2;

        end
        time_BladeModel=toc(timer_1);
        timer_1=tic;

        %% Wake Structure (Semi-Prescribed Wake Sturcture)
        % Data_local_raidus1=[Data_local_raidus1; r./R ,r ,Flow_Angle, alpha, Re, beta, boundVortex, Tdr,Fdr,Qdr];

        Blade1Alpha =[Blade1Alpha ;Data_local_raidus1(:,4)'];
        Blade1Re    =[Blade1Re    ;Data_local_raidus1(:,5)'];
        Blade1Gamma =[Blade1Gamma ;Data_local_raidus1(:,7)'];
        Blade1dT    =[Blade1dT    ;Data_local_raidus1(:,8)'];
        Blade1dF    =[Blade1dF    ;Data_local_raidus1(:,9)'];
        Blade1dQ    =[Blade1dQ    ;Data_local_raidus1(:,10)'];
        Blade1Cl    =[Blade1Cl    ;Data_local_raidus1(:,11)'];
        Blade1Cd    =[Blade1Cd    ;Data_local_raidus1(:,12)'];

        Blade2Alpha =[Blade2Alpha ;Data_local_raidus2(:,4)'];
        Blade2Re    =[Blade2Re    ;Data_local_raidus2(:,5)'];
        Blade2Gamma =[Blade2Gamma ;Data_local_raidus2(:,7)'];
        Blade2dT    =[Blade2dT    ;Data_local_raidus2(:,8)'];
        Blade2dF    =[Blade2dF    ;Data_local_raidus2(:,9)'];
        Blade2dQ    =[Blade2dQ    ;Data_local_raidus2(:,10)'];
        Blade2Cl    =[Blade2Cl    ;Data_local_raidus2(:,11)'];
        Blade2Cd    =[Blade2Cd    ;Data_local_raidus2(:,12)'];



        %TE -> Wake Structure 추가
        %Blade 1
        tmp=now_Geom_Point_TE';
        wakeForm=tmp(:)';
        vinduced_Ax=vinduced1(:,3);
        vinduced_Ax_set=transpose(0.*vinduced_Ax.*unit_conv_Axis'+V_inflow1);
        vinduced_Ax_GeomSet=[];
        vinduced_Ax_GeomSet=0.5.*(vinduced_Ax_set(:,2:end)+vinduced_Ax_set(:,1:end-1));
        vinduced_Ax_GeomSet=[vinduced_Ax_set(:,1), vinduced_Ax_GeomSet, vinduced_Ax_set(:,end)];
        vinduced_Ax_GeomSet=vinduced_Ax_GeomSet(:);
        WakeGammaset=[0;gamma_BoundVortex]-[gamma_BoundVortex;0];
        WakeGammaset(find(Geom_rR<hub))=0;
        Wake_Geom_Position=[wakeForm;Wake_Geom_Position]; %x y z x y z


        %Blade 2
        tmp=now_Geom2_Point_TE';
        wakeForm2=tmp(:)';
        vinduced2_Ax=vinduced2(:,3);
        vinduced2_Ax_set=transpose(0.*vinduced2_Ax.*unit_conv2_Axis'+V_inflow2);
        vinduced2_Ax_GeomSet=[];
        vinduced2_Ax_GeomSet=0.5.*(vinduced2_Ax_set(:,2:end)+vinduced2_Ax_set(:,1:end-1));
        vinduced2_Ax_GeomSet=[vinduced2_Ax_set(:,1), vinduced2_Ax_GeomSet, vinduced2_Ax_set(:,end)];
        vinduced2_Ax_GeomSet=vinduced2_Ax_GeomSet(:);
        WakeGammaset2=[0;gamma_BoundVortex2]-[gamma_BoundVortex2;0];
        WakeGammaset2(find(Geom_rR<hub))=0;




        % WakeGammaset2(1)=0;
        Wake2_Geom_Position=[wakeForm2;Wake2_Geom_Position]; %x y z x y z


        Wake_Gamma=[WakeGammaset';Wake_Gamma];
        Wake2_Gamma=[WakeGammaset2';Wake2_Gamma];

        time_Wake_inducedV=toc(timer_1);
        timer_1=tic;


        %% Wake - Wake Interaction

        TotalWake1_old=TotalWake1;
        TotalWake2_old=TotalWake2;
        V_wake_1=zeros(size(Wake_Geom_Position));
        V_wake_2=zeros(size(Wake_Geom_Position));



        %
        % Wake_velocity=[vinduced_Ax_GeomSet';Wake_velocity];
        % Wake2_velocity=[vinduced2_Ax_GeomSet';Wake2_velocity];
        %


        time_Wake_Wake=toc(timer_1);
        timer_1=tic;



        %% Blade-Wake Interaction
        sizeofWakemetrx=size(Wake_Geom_Position);
        Idx=1;
        if (Total_Rotate>initialAngle)
            vOutVel = WakeVortexCalculator_mex_Fast([Wake_Geom_Position Wake2_Geom_Position], [Wake_Gamma Wake2_Gamma],[rc_Geom rc_Geom]);
            V_wake_1=vOutVel(:,1:size(vOutVel,2)/2);
            V_wake_2=vOutVel(:,1+size(vOutVel,2)/2:end);
            ini_Wake_velocity=[];
            ini_Wake_velocity2=[];
            rc_panel_local=rc_panel;
            while Idx<sizeofWakemetrx(1)*sizeofWakemetrx(2)
                %Wake, geom에 따라 하나씩 계산
                wakeIdx=fix(Idx/sizeofWakemetrx(2))+1;
                geomIdx=(mod(Idx,sizeofWakemetrx(2))-1)/3+1;
                xind=3*(geomIdx-1)+1;
                yind=3*(geomIdx-1)+2;
                zind=3*(geomIdx-1)+3;
                WakeMarker=Wake_Geom_Position(wakeIdx,xind:zind);
                WakeMarker2=Wake2_Geom_Position(wakeIdx,xind:zind);

                for gamma_idx=1:size(gamma_BoundVortex,1)
                    %Blade Vortex로인한 후류 후퇴속도 계산
                    gamma=gamma_BoundVortex(gamma_idx);
                    gamma2=gamma_BoundVortex2(gamma_idx);
                    % induced velocity by Bound Vortex
                    vind1_W1_Bound(gamma_idx,:)     =Vortex_Vatistas(now_Geom_Point_BD(gamma_idx,:),now_Geom_Point_BD(gamma_idx+1,:),WakeMarker,gamma,rc_panel_local(gamma_idx),vortex_n);
                    vind1_W1_trLeft(gamma_idx,:)    =Vortex_Vatistas(now_Geom_Point_TE(gamma_idx,:),now_Geom_Point_BD(gamma_idx,:),WakeMarker,gamma,rc_panel_local(gamma_idx),vortex_n);
                    vind1_W1_trRight(gamma_idx,:)   =Vortex_Vatistas(now_Geom_Point_BD(gamma_idx+1,:),now_Geom_Point_TE(gamma_idx+1,:),WakeMarker,gamma,rc_panel_local(gamma_idx),vortex_n);
                    vind1_W2_Bound(gamma_idx,:)     =Vortex_Vatistas(now_Geom_Point_BD(gamma_idx,:),now_Geom_Point_BD(gamma_idx+1,:),WakeMarker2,gamma2,rc_panel_local(gamma_idx),vortex_n);
                    vind1_W2_trLeft(gamma_idx,:)    =Vortex_Vatistas(now_Geom_Point_TE(gamma_idx,:),now_Geom_Point_BD(gamma_idx,:),WakeMarker2,gamma2,rc_panel_local(gamma_idx),vortex_n);
                    vind1_W2_trRight(gamma_idx,:)   =Vortex_Vatistas(now_Geom_Point_BD(gamma_idx+1,:),now_Geom_Point_TE(gamma_idx+1,:),WakeMarker2,gamma2,rc_panel_local(gamma_idx),vortex_n);



                    % induced velocity by Bound Vortex
                    vind2_W1_Bound(gamma_idx,:)     =Vortex_Vatistas(now_Geom2_Point_BD(gamma_idx,:),now_Geom2_Point_BD(gamma_idx+1,:),WakeMarker,gamma,rc_panel_local(gamma_idx),vortex_n);
                    vind2_W1_trLeft(gamma_idx,:)    =Vortex_Vatistas(now_Geom2_Point_TE(gamma_idx,:),now_Geom2_Point_BD(gamma_idx,:),WakeMarker,gamma,rc_panel_local(gamma_idx),vortex_n);
                    vind2_W1_trRight(gamma_idx,:)   =Vortex_Vatistas(now_Geom2_Point_BD(gamma_idx+1,:),now_Geom2_Point_TE(gamma_idx+1,:),WakeMarker,gamma,rc_panel_local(gamma_idx),vortex_n);
                    vind2_W2_Bound(gamma_idx,:)     =Vortex_Vatistas(now_Geom2_Point_BD(gamma_idx,:),now_Geom2_Point_BD(gamma_idx+1,:),WakeMarker2,gamma2,rc_panel_local(gamma_idx),vortex_n);
                    vind2_W2_trLeft(gamma_idx,:)    =Vortex_Vatistas(now_Geom2_Point_TE(gamma_idx,:),now_Geom2_Point_BD(gamma_idx,:),WakeMarker2,gamma2,rc_panel_local(gamma_idx),vortex_n);
                    vind2_W2_trRight(gamma_idx,:)   =Vortex_Vatistas(now_Geom2_Point_BD(gamma_idx+1,:),now_Geom2_Point_TE(gamma_idx+1,:),WakeMarker2,gamma2,rc_panel_local(gamma_idx),vortex_n);


                end

                Vind_w1=sum([vind1_W1_Bound;vind2_W1_Bound])...
                    +sum([vind1_W1_trLeft;vind2_W1_trLeft])...
                    +sum([vind1_W1_trRight;vind2_W1_trRight]);

                Vind_w2=sum([vind1_W2_Bound;vind2_W2_Bound])...
                    +sum([vind1_W2_trLeft;vind2_W2_trLeft])...
                    +sum([vind1_W2_trRight;vind2_W2_trRight]);

                V_wake_1(wakeIdx,xind:zind)= V_wake_1(wakeIdx,xind:zind)+Vind_w1;
                V_wake_2(wakeIdx,xind:zind)= V_wake_2(wakeIdx,xind:zind)+Vind_w2;

                Idx=Idx+3;
                TotalWake1=V_wake_1+vinduced_Ax_GeomSet';
                TotalWake2=V_wake_2+vinduced2_Ax_GeomSet';
            end

        else
            vinf=sqrt(sum(V_inflow2.^2,2));
            vinf=[vinf(1);vinf];
            C=2.*0.5.*(T+momentumTdr)./(rho.*DiskArea);
            V_wake=-transpose((0.5.*(-vinf+sqrt(vinf.^2+C))).*unit_conv_Axis');
            IniAx_Vel=V_wake(:)+vinduced_Ax_GeomSet;
            % IniAx_Vel=vinduced_Ax_GeomSet;
            ini_Wake_velocity=[IniAx_Vel';ini_Wake_velocity];
            ini_Wake2_velocity=[IniAx_Vel';ini_Wake2_velocity];
            TotalWake1=ini_Wake_velocity;
            TotalWake2=ini_Wake_velocity;
        end




        TotalWake1_old=[TotalWake1(1,:);TotalWake1_old];
        TotalWake2_old=[TotalWake2(1,:);TotalWake2_old];

        Wake_Geom_Position=Wake_Geom_Position+0.5.*(TotalWake1+TotalWake1_old).*dt;
        Wake2_Geom_Position=Wake2_Geom_Position+0.5*(TotalWake2+TotalWake2_old).*dt;



        time_Wake_Blade=toc(timer_1);
        timer_1=tic;

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
        if (Total_Rotate>160 & Total_Rotate<225)|(Total_Rotate>5000)
            Wake_Geom_Position(end,:)=[]; %x y z x y z

            Wake_Gamma(end,:)=[];
            TotalWake1(end,:)=[];
            Wake2_Geom_Position(end,:)=[]; %x y z x y z

            Wake2_Gamma(end,:)=[];
            TotalWake2(end,:)=[];
            try
                ini_Wake2_velocity(end,:)=[];
                ini_Wake_velocity(end,:)=[];
            end
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

        if (Total_Rotate>0)

            % Wake_Geom_Position(end,:)=[]; %x y z x y z
            % Wake_velocity(end,:)=[];
            % Wake_Gamma(end,:)=[];
            % Wake2_Geom_Position(end,:)=[]; %x y z x y z
            % Wake2_velocity(end,:)=[];
            % Wake2_Gamma(end,:)=[];


        end

        %% Post Calculate Data output
        Fx=Total_forceVec(1);
        Fy=Total_forceVec(2);
        Fz=Total_forceVec(3);
        P=Q.*2.*pi.*n;
        globalData=[globalData;Time,Total_Rotate,AzimuthAngle,T,Q,P,T1,T2,Q1,Q2,Fx,Fy,Fz,momentumTdr];
        if size(globalData,1)>int16(360/dAngle)
            Avarge_Global=mean(globalData(end-int16(360/dAngle)+1:end,4:end),1);
        else
            Avarge_Global=globalData(end,4:end).*0;
        end

        Avarge_globalData=[Avarge_globalData;Time,Total_Rotate,AzimuthAngle,Avarge_Global];
        Tav=Avarge_Global(1);
        Qav=Avarge_Global(2);
        Pav=Avarge_Global(3);

        FxAver=Avarge_Global(8);
        FyAver=Avarge_Global(9);
        FzAver=Avarge_Global(10);

        MTav=Avarge_Global(end);
        err=(abs(MTav-Tav)./MTav)*1000;


        time_Postcalc=toc(timer_1);

        fprintf("Iter Calculated!\n")
        fprintf("Time= %.3fsec\n",Time)
        fprintf("Criteria1 : %f\n",cri1)
        fprintf("Criteria2 : %f\n",cri2)

        fprintf("T= %.3fN  Q= %.3fNm\n",T,Q)
        fprintf("T(avarage)= %.3fN | Q(avarge)= %.3fNm | P(avearge)= %.3f‰\n",Tav,Qav,err)
        fprintf("Excell Form : %f, %f, %f, %f, %f, %f\n",Tav,Qav,Pav,FxAver,FyAver,FzAver)

        fprintf("momentumT= %.3fN\n",momentumTdr)

        fprintf("Now Angle= %.3fdeg  TotalAngle= %.3fdeg\n",AzimuthAngle,Total_Rotate)

        fprintf("IJ Matrix Calculate Time    = %.3fsec\n",time_IJcalc)
        fprintf("Wake Induced Calculate Time = %.3fsec\n",time_WakeInduced)
        fprintf("VLM Calc                    = %.3fsec\n",time_BladeModel)

        fprintf("Induced Velocity Prepare    = %.3fsec\n",time_Wake_inducedV)
        fprintf("Wake Wake Interaction Calc. = %.3fsec\n",time_Wake_Wake)
        fprintf("Blade Wake Interaction Calc.= %.3fsec\n",time_Wake_Blade)
        fprintf("Post Calcualte Time         = %.3fsec\n",time_Postcalc)
        fprintf("\n\n")


        OutputBlade1=...
            {
            globalData
            Avarge_globalData
            Blade1Alpha
            Blade1Re
            Blade1Gamma
            Blade1dT
            Blade1dF
            Blade1dQ
            Blade1Cl
            Blade1Cd


            };

        OutputBlade2=...
            {
            globalData
            Avarge_globalData
            Blade2Alpha
            Blade2Re
            Blade2Gamma
            Blade2dT
            Blade2dF
            Blade2dQ
            Blade2Cl
            Blade2Cd


            };

        OutputVortexSturcture=...
            {
            globalData
            Avarge_globalData
            Wake_Geom_Position
            Wake2_Geom_Position
            Wake_Gamma
            Wake2_Gamma
            rc_panel
            };


        figure(100)
        clf
        tiledlayout(3,3)

        nexttile([2,2])

        %% Geometry Calculate

        hold on; grid off

        plot3(now_Panel_Point_Colocation(:,1),now_Panel_Point_Colocation(:,2),now_Panel_Point_Colocation(:,3),'ro-')
        plot3(now_Panel_Point_BD(:,1),now_Panel_Point_BD(:,2),now_Panel_Point_BD(:,3),'rx-')
        plot3(now_Panel_Point_LE(:,1),now_Panel_Point_LE(:,2),now_Panel_Point_LE(:,3),'r*-')
        plot3(now_Panel_Point_TE(:,1),now_Panel_Point_TE(:,2),now_Panel_Point_TE(:,3),'r*-')
        for idx=int16(length(WakeGammaset)/4):length(WakeGammaset)
            %for idx=length(WakeGammaset)
            xind=(idx-1)*3+1;
            yind=(idx-1)*3+2;
            zind=(idx-1)*3+3;

            plot3(Wake_Geom_Position(1:end,xind),Wake_Geom_Position(1:end,yind),Wake_Geom_Position(1:end,zind),'k:')
        end

        plot3(now_Panel2_Point_Colocation(:,1),now_Panel2_Point_Colocation(:,2),now_Panel2_Point_Colocation(:,3),'bo-')
        plot3(now_Panel2_Point_BD(:,1),now_Panel2_Point_BD(:,2),now_Panel2_Point_BD(:,3),'bx-')
        plot3(now_Panel2_Point_LE(:,1),now_Panel2_Point_LE(:,2),now_Panel2_Point_LE(:,3),'b*-')
        plot3(now_Panel2_Point_TE(:,1),now_Panel2_Point_TE(:,2),now_Panel2_Point_TE(:,3),'b*-')
        for idx=int16(length(WakeGammaset)/4):length(WakeGammaset)
            %for idx=length(WakeGammaset)
            xind=(idx-1)*3+1;
            yind=(idx-1)*3+2;
            zind=(idx-1)*3+3;

            plot3(Wake2_Geom_Position(1:end,xind),Wake2_Geom_Position(1:end,yind),Wake2_Geom_Position(1:end,zind),'k:')
        end
        %plot3(Wake2_Geom_Position(1:end,xind),Wake2_Geom_Position(1:end,yind),Wake2_Geom_Position(1:end,zind),'r:')
        %plot3(Wake_Geom_Position(1:end,xind),Wake_Geom_Position(1:end,yind),Wake_Geom_Position(1:end,zind),'r:')

        xlabel("x (m)")
        ylabel("y (m)")
        zlabel("z (m)")
        title(["Propeller & Wake Structure ";"(semi-Free Wake method)"])


        view([1,1,1])
        axis equal
        %
        % filename=sprintf("FIG_%05dRPM_%04dANGLE.png",RPM,Total_Rotate)
        % exportgraphics(figure(100),'fig1/'+filename,'Resolution',300)
        % view(unit_conv2_Chord')
        % axis equal
        % filename=sprintf("FIG_%05dRPM_%04dANGLE.png",RPM,Total_Rotate)
        % exportgraphics(figure(100),'fig2/'+filename,'Resolution',300)
        % view([1,1,1])



        nexttile

        hold on; grid on
        plot(globalData(:,2),Avarge_globalData(:,4),'k')
        plot(globalData(:,2),Avarge_globalData(:,7),'r')

        plot(globalData(:,2),Avarge_globalData(:,8),'b')
        plot(globalData(:,2),Avarge_globalData(:,end),'g')

        plot(globalData(:,2),globalData(:,4),'k:')
        plot(globalData(:,2),globalData(:,7),'r:')

        plot(globalData(:,2),globalData(:,8),'b:')
        plot(globalData(:,2),globalData(:,end),'g:')

        xlabel("Total Rotation Angle (deg)")
        ylabel("Thrust(N)")
        title("Thrust Calculation")
        legend("Total Thrust", "Blade 1 Thrust", "Blade 2 Thrust","Momentum Thrust",'Location',"northwest")

        nexttile

        hold on; grid on
        plot(Data_local_raidus1(:,1),Data_local_raidus1(:,8),'rx-') % Flow
        plot(Data_local_raidus2(:,1),Data_local_raidus2(:,8),'bo-') % Flow

        legend("Blade 1 Thrust", "Blade 2 Thrust",'Location','northwest')
        title("Local Blade Thrust ")
        xlabel("r/R")
        ylabel("N")

        nexttile


        hold on; grid on
        plot(Data_local_raidus1(:,1),Data_local_raidus1(:,3),'r') % Flow
        plot(Data_local_raidus1(:,1),Data_local_raidus1(:,4),'g') % alpha
        plot(Data_local_raidus1(:,1),Data_local_raidus1(:,6),'b') %beta
        legend("Flow","Alpha","beta",'Location','northwest')
        title("Inflow Angle")
        xlabel("r/R")
        ylabel("Deg")

        nexttile
        hold on; grid on
        plot(Data_local_raidus1(:,1),Blade1Re(end,:),'rx-') % Flow
        plot(Data_local_raidus2(:,1),Blade2Re(end,:),'bo-') % Flow
        legend("Blade 1", "Blade 2",'Location','northwest')
        grid on
        title("Reynolds Number")
        xlabel("r/R")
        ylabel("Reynolds Number")

        nexttile
        hold on; grid on
        plot(Data_local_raidus1(:,1),V_inflow1(:,3),'r')
        plot(Data_local_raidus1(:,1),vind_Wake(:,3),'g') % Flow
        plot(Data_local_raidus1(:,1),v_infow_induced_z_set,'b') % alpha
        plot(Data_local_raidus1(:,1),V_inflow1(:,3)+vind_Wake(:,3)+v_infow_induced_z_set,'k') % alpha

        legend("Flow","Wake Induced","Blade Induced","Total",'Location','best')
        title("Induced Velocity (blade1)")


        xlabel("r/R")
        ylabel("Reynolds Number")
        Logfilename=CaseTitle+"_log.txt";
        fid=fopen(Logfilename,'a');
        CalcData=[CaseInd, method vortex_n rc_Ratio Calc_case(CaseInd), Tav,Qav,Pav,FxAver,FyAver,FzAver];

        fprintf(fid,"%d,%d,%d,%f,%f,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f \n",Total_Rotate ,CaseInd, method, vortex_n, rc_Ratio, Calc_case(CaseInd), Tav,Qav,Pav,FxAver,FyAver,FzAver);
                fprintf(fid,"\n");

        fclose(fid);
    end


    Logfilename=CaseTitle+"_log.txt";
    fid=fopen(Logfilename,'a');
    LocalCalcData=[CaseInd, method vortex_n rc_Ratio Calc_case(CaseInd), Tav,Qav,Pav,FxAver,FyAver,FzAver];
    Calc_data=[Calc_data;CaseInd, method vortex_n rc_Ratio Calc_case(CaseInd), Tav,Qav,Pav,FxAver,FyAver,FzAver];

    %% Local log
    Logfilename=CaseTitle+"_log.txt";
    fid=fopen(Logfilename,'a');

    Clac_time=toc(GlobalTimer);
    fprintf(fid,"Calculate Done!!\n");

    fprintf(fid,"Calculate Time : %.3f s\n",Clac_time);

    fprintf(fid,"Time= %.3fsec\n",Time);

    fprintf(fid,"T= %.3fN  Q= %.3fNm\n",T,Q);
    fprintf(fid,"T(avarage)= %.3fN | Q(avarge)= %.3fNm | P(avearge)= %.3f‰\n",Tav,Qav,err);
    fprintf(fid,"Excell Form : %f, %f, %f, %f, %f, %f\n",Tav,Qav,Pav,FxAver,FyAver,FzAver);

    fprintf(fid,"momentumT= %.3fN\n",momentumTdr);

    fprintf(fid,"Now Angle= %.3fdeg  TotalAngle= %.3fdeg\n",AzimuthAngle,Total_Rotate);

    fprintf(fid,"IJ Matrix Calculate Time    = %.3fsec\n",time_IJcalc);
    fprintf(fid,"Wake Induced Calculate Time = %.3fsec\n",time_WakeInduced);
    fprintf(fid,"VLM Calc                    = %.3fsec\n",time_BladeModel);

    fprintf(fid,"Induced Velocity Prepare    = %.3fsec\n",time_Wake_inducedV);
    fprintf(fid,"Wake Wake Interaction Calc. = %.3fsec\n",time_Wake_Wake);
    fprintf(fid,"Blade Wake Interaction Calc.= %.3fsec\n",time_Wake_Blade);
    fprintf(fid,"Post Calcualte Time         = %.3fsec\n",time_Postcalc);
    fprintf(fid,"\n\n");

    fclose(fid);



    %% Full log


    Full_Logfilename=sprintf("[LOG]_%s.txt",caseName);
    fid=fopen(Full_Logfilename,'a');
    fprintf(fid,"%f, ",LocalCalcData);
    fprintf(fid,"\n");
    fclose(fid);

    save("OutputData\"+CaseTitle+".mat",'OutputVortexSturcture','OutputBlade1','OutputBlade2','OutputVortexSturcture');

end
