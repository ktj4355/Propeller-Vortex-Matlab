function [VelocityMatrix] = WakeVortexCalculator_C(DLLname,WakeGeom_ptr,GammaMatrix_ptr,rc_ptr)
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



    vOutVel = zeros(size(WakeGeom_ptr.Value));
    vOutVel_1D=vOutVel(:);
    VelocityMatrix_ptr = libpointer('doublePtr', vOutVel_1D);
    
    calllib(DLLname, 'WakeVortexCalculator',WakeGeom_ptr,GammaMatrix_ptr,rc_ptr,size(GammaMatrix_ptr.Value,1),size(GammaMatrix_ptr.Value,2),VelocityMatrix_ptr);
    vOutVel=VelocityMatrix_ptr.Value;
    VelocityMatrix=reshape(vOutVel,size(WakeGeom_ptr.Value));


end

