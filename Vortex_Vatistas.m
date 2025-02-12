function [Vout] = Vortex_Vatistas(A,B,ColocationPoint,vortexStrength,rc,n)
% Vortex_Vatistas(A,B,ColocationPoint,vortexStrength,rc,n)
% A, B, Colocation Point 좌표를 이용하여 Colocation Point의 유도속도 계산
%           
%     A-----B
%  r1->\   /
%       \ / -> r2
%        P(colocation Point
r1=ColocationPoint-A;
r2=ColocationPoint-B;

r1s=norm(r1);
r2s=norm(r2);
dL=abs(norm(A-B));

cosA=dot(r1,(B-A))/(r1s*dL);
cosB=dot(r2,(B-A))/(r2s*dL);
r=norm(cross(r1,r2))./dL;

if r<0.000001
    Vout=[0 0 0];
    return
end


r_norm=r/rc;
   
    vortexFactor=(vortexStrength./(4*pi*r))*(          (r^2)/( (rc^(2*n)+r^(2*n))^(1/n) )        );
    Vout=vortexFactor*(cosA-cosB)*cross(r1,r2)./(norm(cross(r1,r2)));

end

