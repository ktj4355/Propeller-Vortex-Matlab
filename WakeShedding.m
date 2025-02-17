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

%inputGeom(:,4)=0

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

%% Calcualte Colllocation Point of Blade Coordinate
PannelGeom=[];
PropGeom_r=inputGeom(:,1:2);
%Clollocation Point Data Matrix
PannelGeom_r=[]; %rR ,r
PannelGeom_Chorld=[];
PannelGeom_Beta=[];
PannelGeom_Thickness=[];
collocation_point_local=[];
BDVortex_local=[];

% prop Geom IDX  1---2---3---4---5---6---7---8
% Panel Geom IDX |-1-|-2-|-3-|-4-|-5-|-6-|-7-|
% Coloation Geom IDX is Same to Panel Geom
%r/R r Chord Beta(deg) TicknessRatio

for idx=1:size(inputGeom,1)-1
    Section1=inputGeom(idx,:);
    Section2=inputGeom(idx+1,:);
    PannelGeom(idx,:)=(Section1+Section2)./2;
    PannelGeom_r(idx,:)=[PannelGeom(idx,1), PannelGeom(idx,2)];
    PannelGeom_Chorld(idx)=PannelGeom(idx,3);
    PannelGeom_Beta(idx)=PannelGeom(idx,4);
    PannelGeom_Thickness(idx)=PannelGeom(idx,5);


    %LE:y=0

    % Twist center is 1/4 point = y=0
    collocation_point_local(idx,:)=[PannelGeom_r(idx,2),-PannelGeom_Chorld(idx).*0.5,0];
    %Collocation Point at 3/4 chord

    BDVortex_local(idx,:)=[PannelGeom_r(idx,2),-PannelGeom_Chorld(idx).*0.25,0];
    %Bound Vortex lIne at 1/4 chord -> x=0, Twist Center line

    %le -------- BD --------- | ------- Collocation ------ te


end

%% Calculate Real Vector
% Calculate Geom Real Vector (Bound Vortex Line is Rotate Axis
% Geom Point Vector : Real Blade Positon -> VLM Position

GeomPointVector_BD_local=[inputGeom(:,2),-inputGeom(:,3).*0.25,0.*inputGeom(:,2)];
% 0.25C를 기준으로하여 블레이드 Leading Edge point 계산
% r, 반경좌표
% n, Disk 수직 벡터
% t, 디스크 접선벡터

r=GeomPointVector_BD_local(:,1);                        % local section radial position
n=GeomPointVector_BD_local(:,2).*sind(inputGeom(:,4));  %
t=GeomPointVector_BD_local(:,2).*cosd(inputGeom(:,4));
GeomPointVector_BD_Global_ini=[r,r.*0,r.*0];            %[r,0,0]
GeomPointVector_LE_Global_ini=[r,-t,-n];
GeomPointVector_TE_Global_ini=[r,3.*t,4.*n];


%36point

% Calculate Collocation Real Vector (Bound Vortex Line is Rotate Axis

r=BDVortex_local(:,1);
n=BDVortex_local(:,2).*sind(PannelGeom_Beta)';
t=BDVortex_local(:,2).*cosd(PannelGeom_Beta)';

BD_point_Global_ini=[r,r.*0,r.*0];
LE_point_Global_ini=[r,-t,-n];
TE_point_Global_ini=[r,3.*t,4.*n];
collocation_point_Global_ini=[r,2.*t,n];



%% Calculate Azimuth Angle

nGeom=size(inputGeom,1);
nColo=nGeom-1;

BladeAngle=linspace(0,360-360/nAzmuth,nAzmuth);
Rz=@(theta)[cosd(theta) -sind(theta), 0;sind(theta),cosd(theta), 0;0,0,1];
for angleIdx=1:size(BladeAngle,2)
    angle=BladeAngle(angleIdx);

    point=[];
    LEpoint=[];
    TEpoint=[];
    BDpoint=[];

    angle=BladeAngle(angleIdx);

    for ridx=1:size(collocation_point_Global_ini,1)
        point(ridx,:)=(Rz(angle)*collocation_point_Global_ini(ridx,:)')';
        LEpoint(ridx,:)=(Rz(angle)*LE_point_Global_ini(ridx,:)')';
        TEpoint(ridx,:)=(Rz(angle)*TE_point_Global_ini(ridx,:)')';
        BDpoint(ridx,:)=(Rz(angle)*BD_point_Global_ini(ridx,:)')';


        Geom_LEpoint(ridx,:)=(Rz(angle)*GeomPointVector_LE_Global_ini(ridx,:)')';
        Geom_TEpoint(ridx,:)=(Rz(angle)*GeomPointVector_TE_Global_ini(ridx,:)')';
        Geom_BDpoint(ridx,:)=(Rz(angle)*GeomPointVector_BD_Global_ini(ridx,:)')';
        RLine=[LEpoint(ridx,:);TEpoint(ridx,:)];
    end
    Geom_LEpoint(ridx+1,:)=(Rz(angle)*GeomPointVector_LE_Global_ini(end,:)')';
    Geom_TEpoint(ridx+1,:)=(Rz(angle)*GeomPointVector_TE_Global_ini(end,:)')';
    Geom_BDpoint(ridx+1,:)=(Rz(angle)*GeomPointVector_BD_Global_ini(end,:)')';

    % Debug plot
    % figure(2)
    % hold on
    % axis equal
    % scatter3(point(:,1), point(:,2), point(:,3),10+angle./10,'red','filled');
    % scatter3(BDpoint(:,1), BDpoint(:,2), BDpoint(:,3),10+angle./10,'go');
    % scatter3(Geom_LEpoint(:,1), Geom_LEpoint(:,2), Geom_LEpoint(:,3),10+angle./10,'cx');
    % scatter3(Geom_TEpoint(:,1), Geom_TEpoint(:,2), Geom_TEpoint(:,3),10+angle./10,'mx');
    % scatter3(Geom_BDpoint(:,1), Geom_BDpoint(:,2), Geom_BDpoint(:,3),10+angle./10,'kx');
    % pause(0.1);

    % update Global Cell
    collocation_point_Global{angleIdx}=point;
    BD_point_Global{angleIdx}=BDpoint;
    LE_point_Global{angleIdx}=LEpoint;
    TE_point_Global{angleIdx}=TEpoint;

    Geom_BD_point_Global{angleIdx}=Geom_BDpoint;
    Geom_LE_point_Global{angleIdx}=Geom_LEpoint;
    Geom_TE_point_Global{angleIdx}=Geom_TEpoint;
end

%% Calculate Iij Matrix
I_global={};
for angleIdx=1:size(BladeAngle,2)
    % set Point
    CollocPoint=collocation_point_Global{angleIdx};
    % Set Boundary Vortex, Trailing Vortex (in Blade)
    % Trailing Vortex IDX   1---2---3---4---5---6---7---8
    % Boundary Vortex IDX   |-1-|-2-|-3-|-4-|-5-|-6-|-7-|

    gammaSet_BoundVortex=ones(size(CollocPoint,1),1);
    %gammaSet_BoundVortex=linspace(0.1,2,40)'

    BDpoint=BD_point_Global{angleIdx};
    LEpoint=LE_point_Global{angleIdx};
    TEpoint=TE_point_Global{angleIdx};

    Geom_BDpoint=Geom_BD_point_Global{angleIdx};
    Geom_LEpoint=Geom_LE_point_Global{angleIdx};
    Geom_TEpoint=Geom_TE_point_Global{angleIdx};
    rc=PannelGeom_Chorld.*0.1;
    for ColloIDX=1:size(CollocPoint,1) %Collocation Point Loop
        ColoP=CollocPoint(ColloIDX,:);
        for gamma_idx=1:size(gammaSet_BoundVortex,1) %Bound GaMMA Loop

            %Vortex Modeling

            n=2;
            gamma=gammaSet_BoundVortex(gamma_idx);

            % induced velocity by Bound Vortex
            vind_Bound(gamma_idx,:)=Vortex_Scully(Geom_BDpoint(gamma_idx,:),Geom_BDpoint(gamma_idx+1,:),ColoP,gamma,rc(gamma_idx));
            vind_trLeft(gamma_idx,:)=Vortex_Scully(Geom_TEpoint(gamma_idx,:),Geom_BDpoint(gamma_idx,:),ColoP,gamma,rc(gamma_idx));
            vind_trRight(gamma_idx,:)=Vortex_Scully(Geom_BDpoint(gamma_idx+1,:),Geom_TEpoint(gamma_idx+1,:),ColoP,gamma,rc(gamma_idx));
            vTotal=vind_Bound(gamma_idx,:)+ vind_trLeft(gamma_idx,:)+ vind_trRight(gamma_idx,:);

            I_ini_x(ColloIDX,gamma_idx)=vTotal(1);
            I_ini_y(ColloIDX,gamma_idx)=vTotal(2);
            I_ini_z(ColloIDX,gamma_idx)=vTotal(3);
        end

        %    figure(3);  axis equal;    hold on;      quiver3(Geom_BDpoint(1:end-1,1),Geom_BDpoint(1:end-1,2),Geom_BDpoint(1:end-1,3),r1data(:,1),r1data(:,2),r1data(:,3))
        %    figure(3);   axis equal;   hold on;      quiver3(Geom_BDpoint(2:end,1),Geom_BDpoint(2:end,2),Geom_BDpoint(2:end,3),r2data(:,1),r2data(:,2),r2data(:,3))
        Vind_collo(ColloIDX,:)=sum(vind_Bound)+sum(vind_trLeft)+sum(vind_trRight);
        I_global{angleIdx,1}=I_ini_x;
        I_global{angleIdx,2}=I_ini_y;
        I_global{angleIdx,3}=I_ini_z;
        I_global{angleIdx,4}=Vind_collo;

    end



    % figure(4)
    % hold on
    % angleIdx;
    % %scatter3(r1(1),r1(2),r1(3),angleIdx.*10);
    % %scatter3(r2(1),r2(2),r2(3),angleIdx.*10);
    % figure(3);     hold on;  axis equal;    quiver3(CollocPoint(:,1),CollocPoint(:,2),CollocPoint(:,3),Vind_collo(:,1),Vind_collo(:,2),Vind_collo(:,3))

end


%% Prop Geometric Plot
%
% for idx=1:size(inputGeom,1)
%     % xb yb zb : BladeLocalAxis
%     % (yc) (Chord Axis)
%     % |
%     % |
%     % ---------------------------- (xc) (Span Axis)
%     xb_LETE=inputGeom(idx,2);
%     b_local=inputGeom(idx,4);
%     yb_LE=(inputGeom(idx,3)./2).*cosd(b_local);
%     yb_TE=(-inputGeom(idx,3)./2).*cosd(b_local);
%     zb_LE=(inputGeom(idx,3)./2).*sind(b_local);
%     zb_TE=(-inputGeom(idx,3)./2).*sind(b_local);
%
%     LEPosition(idx,:)=[xb_LETE yb_LE zb_LE];
%     TEPosition(idx,:)=[xb_LETE yb_TE zb_TE];
%     MeanPosition(idx,:)=[xb_LETE,0,0];
% end
% theta=linspace(0,2.*pi,36);
% Hubposition=[Rhub.*cos(theta)',Rhub.*sin(theta)'];
% figure(1)
% hold on
%
% plot3(LEPosition(:,1),LEPosition(:,2),LEPosition(:,3),'r-')
% plot3(TEPosition(:,1),TEPosition(:,2),TEPosition(:,3),'b-')
% plot3(MeanPosition(:,1),MeanPosition(:,2),MeanPosition(:,3),'k-')
% plot3(-LEPosition(:,1),-LEPosition(:,2),LEPosition(:,3),'r-')
% plot3(-TEPosition(:,1),-TEPosition(:,2),TEPosition(:,3),'b-')
% plot3(-MeanPosition(:,1),-MeanPosition(:,2),-MeanPosition(:,3),'k-')
% plot(Hubposition(:,1),Hubposition(:,2),'k-')
% for idx=1:size(inputGeom,1)
%     LEP=LEPosition(idx,:);
%     TEP=TEPosition(idx,:);
%     lineP=[LEP;TEP];
%     plot3(lineP(:,1),lineP(:,2),lineP(:,3),'k-')
%     plot3(-lineP(:,1),-lineP(:,2),lineP(:,3),'k-')
%
% end
% scatter(collocation_point_Global(:,1),collocation_point_Global(:,2),'k')
%
% axis equal
%

%% time Marching Structure
% trailing Vortex를 Vind, Vfree, WakeVind에 맞춰서 프로펠러 전진시킴(Vortex를 뒤로 후퇴시킴)

%t=0, -1st, -2nd, -3rd, -4th........
positon_Traling_Wake=[];
gammaSet_Traling_Wake=[];
now_Traling_Wake=[];
Pre_Traling_Wake=[];
rc_trail=BladeChord.*0.1;
dAngle_deg=360/nAzmuth;
dAngle_rad=2*pi/nAzmuth;
v_rot_rad=(RPM/60).*(2*pi);  %rad/s
dt=dAngle_rad/v_rot_rad;
t_global=0;
angle_count=0;
%gammaSet_BoundVortex=gammaSet_BoundVortex.*(linspace(0.05,0.7,length(gammaSet_BoundVortex)).*linspace(0.7,0.05,length(gammaSet_BoundVortex)))';
loadlibrary('VortexDLL.dll','Vortex_DLL.h')

    WakeGeom_ptr = libpointer('doublePtr', positon_Traling_Wake);
    GammaMatrix_ptr = libpointer('doublePtr', gammaSet_Traling_Wake);
    rc_ptr = libpointer('doublePtr', rc_trail);

while t_global<0.1

    t_global;
    size(positon_Traling_Wake,1)
    angle_ind=mod(angle_count,nAzmuth)+1;
    CollocPoint=collocation_point_Global{angle_ind};
    BDpoint=BD_point_Global{angle_ind};
    LEpoint=LE_point_Global{angle_ind};
    TEpoint=TE_point_Global{angle_ind};
    Geom_BDpoint=Geom_BD_point_Global{angle_ind};
    Geom_LEpoint=Geom_LE_point_Global{angle_ind};
    Geom_TEpoint=Geom_TE_point_Global{angle_ind};
    I_ini_x=I_global{angle_ind,1};
    I_ini_y=I_global{angle_ind,2};
    I_ini_z=I_global{angle_ind,3};


    gammaSet_BoundVortex=gammaSet_BoundVortex;
    vx_induced_bladeVortex=I_ini_x*gammaSet_BoundVortex;

    vy_induced_bladeVortex=I_ini_y*gammaSet_BoundVortex;
    vz_induced_bladeVortex=I_ini_z*gammaSet_BoundVortex;
    % Blade Vortex로 인한 유도속도
    v_induced_bladeVortex=[vx_induced_bladeVortex vy_induced_bladeVortex vz_induced_bladeVortex];
    v_induced_bladeVortex_TE=interp1(PannelGeom_r(:,1),v_induced_bladeVortex,PropGeom_r(:,1),"linear","extrap");




    % Local BoundGamma Setting
    % TE gamma를 적분하는 필라멘트 벡터 방향은 LE -> Te 방향 (t-1 -> t-2 방향)
    % 그래서 왼쪽께 -로시작
    gammaSet_TralingVortex=-gammaSet_BoundVortex(1);
    for ind=1:length(gammaSet_BoundVortex)-1
        %ind = bound vortex gamma set index, [1 ~ end-1]  (1, 1-2, 2-3 ...)
        gammaSet_TralingVortex(ind+1,1)=gammaSet_BoundVortex(ind)-gammaSet_BoundVortex(ind+1);
    end
    gammaSet_TralingVortex(length(gammaSet_BoundVortex)+1)=gammaSet_BoundVortex(end);


    sizeofWakemetrx=size(positon_Traling_Wake);
    rev_positon_Traling_Wake=positon_Traling_Wake';
    % Wake 전파
    Idx=1;
    %{
    while Idx<sizeofWakemetrx(1)*sizeofWakemetrx(2)
        %Wake, geom에 따라 하나씩 계산
        wakeIdx=fix(Idx/sizeofWakemetrx(2))+1;
        geomIdx=(mod(Idx,sizeofWakemetrx(2))-1)/3+1;
        xind=3*(geomIdx-1)+1;
        yind=3*(geomIdx-1)+2;
        zind=3*(geomIdx-1)+3;
        WakeMarker=positon_Traling_Wake(wakeIdx,xind:zind);
        for gamma_idx=1:nGeom-1
            %Blade Vortex로인한 후류 후퇴속도 계산
            gamma=gammaSet_BoundVortex(gamma_idx);
            vind_Bound(gamma_idx,:)=Vortex_Scully(Geom_BDpoint(gamma_idx,:),Geom_BDpoint(gamma_idx+1,:),WakeMarker,gamma,rc(gamma_idx));
            vind_trLeft(gamma_idx,:)=Vortex_Scully(Geom_TEpoint(gamma_idx,:),Geom_BDpoint(gamma_idx,:),WakeMarker,gamma,rc(gamma_idx));
            vind_trRight(gamma_idx,:)=Vortex_Scully(Geom_BDpoint(gamma_idx+1,:),Geom_TEpoint(gamma_idx+1,:),WakeMarker,gamma,rc(gamma_idx));
            vTotal=vind_Bound(gamma_idx,:)+ vind_trLeft(gamma_idx,:)+ vind_trRight(gamma_idx,:);
        end

        Vind_blade_position=sum(vind_Bound)+sum(vind_trLeft)+sum(vind_trRight);
        v_induced_Axis=[0 0 v_induced_bladeVortex_TE(geomIdx,3)];
        vtotal=Vind_blade_position+vFree+v_induced_Axis;

        positon_Traling_Wake(wakeIdx,xind)= positon_Traling_Wake(wakeIdx,xind)+vtotal(1)*dt;
        positon_Traling_Wake(wakeIdx,yind)= positon_Traling_Wake(wakeIdx,yind)+vtotal(2)*dt;
        positon_Traling_Wake(wakeIdx,zind)= positon_Traling_Wake(wakeIdx,zind)+vtotal(3)*dt;


        Idx=Idx+3;


    end
    %}
    %save("WAKE_DLLTEST","positon_Traling_Wake","gammaSet_Traling_Wake","rc_trail")
%    positon_Traling_Wake_1D=WakeGeom(:);
   % gammaSet_Traling_Wake_1D=GammaMatrix(:);
    WakeGeom_ptr = libpointer('doublePtr', positon_Traling_Wake);
    GammaMatrix_ptr = libpointer('doublePtr', gammaSet_Traling_Wake);
    




    if angle_count>0

        DLLname='VortexDLL';
        vOutVel=WakeVortexCalculator_C(DLLname,WakeGeom_ptr,GammaMatrix_ptr,rc_ptr);
        % vOutVel2=WakeVortexCalculator(positon_Traling_Wake,gammaSet_Traling_Wake,rc_trail);
        %disp(max(abs(vOutVel2-vOutVel),[],'all'))
        positon_Traling_Wake=positon_Traling_Wake+vOutVel.*dt;

        %pause(0.01)

    end
    % Wake Matrix Update
    %n번째 Wake Position 의 위치 =[3(n-1)+1,3(n-1)+2,3(n-1)+3]
    % [x1 y1 z1 x2 y2 z2 ........]
    % [gamma1 gamma2 gamma3 gamma4...]
    Pre_Traling_Wake=now_Traling_Wake;
    now_Traling_Wake=reshape(Geom_TEpoint',1,[]);

    
    gammaSet_Traling_Wake=[gammaSet_TralingVortex';gammaSet_Traling_Wake];   %Wake position Matrix
    positon_Traling_Wake=[now_Traling_Wake;positon_Traling_Wake];            %Wake gamma Matrix

    if (angle_count/nAzmuth)>10
        positon_Traling_Wake(end,:)=[];

        break;
    end
    % Debug show




    figure(2)
    clf
    hold on
    axis equal
    scatter3(CollocPoint(:,1), CollocPoint(:,2), CollocPoint(:,3),10+angle./10,'red','filled');
    scatter3(BDpoint(:,1), BDpoint(:,2), BDpoint(:,3),10+angle./10,'go');
    scatter3(Geom_LEpoint(:,1), Geom_LEpoint(:,2), Geom_LEpoint(:,3),10+angle./10,'cx');
    scatter3(Geom_TEpoint(:,1), Geom_TEpoint(:,2), Geom_TEpoint(:,3),10+angle./10,'mx');
    scatter3(Geom_BDpoint(:,1), Geom_BDpoint(:,2), Geom_BDpoint(:,3),10+angle./10,'kx');
    view ([1,1,1])
    xlim([-0.2,0.2])
    ylim([-0.2,0.2])
    %pause(0.1);
    for geomIdx=1:nGeom
        xind=3*(geomIdx-1)+1;
        yind=3*(geomIdx-1)+2;
        zind=3*(geomIdx-1)+3;
        grid on
        plot3(positon_Traling_Wake(:,xind),positon_Traling_Wake(:,yind),positon_Traling_Wake(:,zind),'k-')
    end



    %count on
    angle_count=angle_count+1;
    t_global=t_global+dt;


end
unloadlibrary('VortexDLL')
%% data plot
figure(3)
clf
hold on
axis equal
scatter3(CollocPoint(:,1), CollocPoint(:,2), CollocPoint(:,3),10+angle./10,'red','filled');
scatter3(BDpoint(:,1), BDpoint(:,2), BDpoint(:,3),10+angle./10,'go');
scatter3(Geom_LEpoint(:,1), Geom_LEpoint(:,2), Geom_LEpoint(:,3),10+angle./10,'cx');
scatter3(Geom_TEpoint(:,1), Geom_TEpoint(:,2), Geom_TEpoint(:,3),10+angle./10,'mx');
scatter3(Geom_BDpoint(:,1), Geom_BDpoint(:,2), Geom_BDpoint(:,3),10+angle./10,'kx');
view ([1,0,0])
xlim([-0.2,0.2])
ylim([-0.2,0.2])
%pause(0.1);
for geomIdx=1:nGeom-1
    xind=3*(geomIdx-1)+1;
    yind=3*(geomIdx-1)+2;
    zind=3*(geomIdx-1)+3;
    grid on
    plot3(positon_Traling_Wake(:,xind),positon_Traling_Wake(:,yind),positon_Traling_Wake(:,zind),"k:");

end

xind=3*(nGeom-1)+1;
yind=3*(nGeom-1)+2;
zind=3*(nGeom-1)+3;
grid on
plot3(positon_Traling_Wake(:,xind),positon_Traling_Wake(:,yind),positon_Traling_Wake(:,zind),"r-");
%    count on
save('FullData',"positon_Traling_Wake","CollocPoint","BDpoint","Geom_LEpoint","Geom_TEpoint","Geom_BDpoint","nGeom");



%% Rotation Vector Calculate


%% Calculate Each Blade Section Lift


%% Calculate Influence Vector


%% Calculate Influence Vector


