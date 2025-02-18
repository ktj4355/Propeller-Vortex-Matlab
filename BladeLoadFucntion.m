function [outpuTtArg1,outputArg2] = BladeLoadFucntion(gamma_BoundVortex,inputArg2)
%BLADELOADFUCNTION 이 함수의 요약 설명 위치
%   자세한 설명 위치
%% Rotation Matrix
Rz=@(theta)[cosd(theta) -sind(theta), 0;sind(theta),cosd(theta), 0;0,0,1];

R_Pitch=@(Pitch)[cosd(Pitch),  0,  sind(Pitch);
    0,            1,  0;
    -sind(Pitch),  0,  cosd(Pitch)];

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

    T=T+Tdr;
    Q=Q+dr;
    F=F+dr;
end
T;
outputArg1 = inputArg1;
outputArg2 = inputArg2;
end

