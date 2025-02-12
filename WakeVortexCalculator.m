function [VelocityMatrix] = WakeVortexCalculator(WakeGeom,GammaMatrix,rc)
% [VelocityMatrix] = WakeVortexCalculator(WakeGeom,GammaMatrix,rc)
% GammaMatrix: Trailing Vortex Strength Array
% Sturcture Description 
%
% <Wake Geom Structure nx3m>
% Wake t   | [x1 y1 z1] [x2 y2 z2] - - - [xm ym zm]  
% Wake t-1 | [x1 y1 z1] [x2 y2 z2] - - - [xm ym zm]  
% Wake t-2 | [x1 y1 z1] [x2 y2 z2] - - - [xm ym zm]  
% Wake t-n | [x1 y1 z1] [x2 y2 z2] - - - [xm ym zm]  
%
% <Gamma Matrix Structure nx3m>
% Wake t   | [g1] [g2] - - - [gm]  
% Wake t-1 | [g1] [g2] - - - [gm]  
% Wake t-2 | [g1] [g2] - - - [gm]  
% Wake t-n | [g1] [g2] - - - [gm]
%
% rc : [rc1 rc2 rc3 ... rcm]
%
% <Output Matrix nx3m >
% Wake t   | [xVel1 yVel1 zVel1] [xVel2 yVel2 zVel2] - - - [xVelm yVelm zVelm]  
% Wake t-1 | [xVel1 yVel1 zVel1] [xVel2 yVel2 zVel2] - - - [xVelm yVelm zVelm]  
% Wake t-2 | [xVel1 yVel1 zVel1] [xVel2 yVel2 zVel2] - - - [xVelm yVelm zVelm]  
% Wake t-n | [xVel1 yVel1 zVel1] [xVel2 yVel2 zVel2] - - - [xVelm yVelm zVelm]  

WakeGeom_size=size(WakeGeom);
WakeSize=size(WakeGeom,1);
GeomPointSize=size(WakeGeom,2)/3;
VelocityMatrix=zeros(WakeGeom_size(1),WakeGeom_size(2));

for idx=1:3:WakeGeom_size(1)*WakeGeom_size(2) % MATLAB에서 다중 For문을 최소화 하기위해 Index로 탐색하기위함
%C코드에서는 다중For문으로 바꾸는게 좋을듯
    clc
    wakeIdx=fix(idx/WakeGeom_size(2))+1;  %n번째 Wake
    geomIdx=(mod(idx,WakeGeom_size(2))-1)/3+1; % n번째 Wake에서 m번쨰 x인 xm의 위치
        % Wake에서 각 x y z좌표
      geomind_X=(3*(geomIdx-1)+1); 
      geomind_Y=(3*(geomIdx-1)+2);
      geomind_Z=(3*(geomIdx-1)+3);
        % 현재 계산중인 Wake Marker의 좌표
      Vind_Point=WakeGeom(wakeIdx,geomind_X:geomind_Z);

     %trailing Wake Vortex Calculation by wake and gamma
     vind=[0 0 0]; % 다른 Wake Marker로부터 만들어지는 Wake에 유도되는 속도
     
     %Trailing Wake
     for idx_LocalWake=1:WakeSize-1
        
        for idx_LocalGamma= 1:GeomPointSize
        LocalGamma=GammaMatrix(idx_LocalWake,idx_LocalGamma);
        wakeind_X=3*(idx_LocalGamma-1)+1;
        wakeind_Y=3*(idx_LocalGamma-1)+2;
        wakeind_Z=3*(idx_LocalGamma-1)+3;
        wakePoint1=WakeGeom(idx_LocalWake,wakeind_X:wakeind_Z);
        wakePoint2=WakeGeom(idx_LocalWake+1,wakeind_X:wakeind_Z);
        %integrate Velocity
        vind=vind+Vortex_Scully(wakePoint1,wakePoint2,Vind_Point,LocalGamma,rc(idx_LocalGamma));
        
        end
    end

    VelocityMatrix(wakeIdx,geomind_X:geomind_Z)=vind;


end

end

