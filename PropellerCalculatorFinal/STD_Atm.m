function [T, P, den, D_vis, a] = STD_Atm(Height)
%STD_ATM 이 함수의 요약 설명 위치
%   자세한 설명 위치


if Height <11,000;
    T=288.15-0.0065.*Height;
    P=101.325.*(T./288.15).^5.2559;
elseif Height<25000
    T=273.15-56.46;
    P=22.65*exp(1.73-0.000157.*Height);
else Height>=25000;
    T=273.15-131.21+0.00299*Height;
    P=2.488.*(T./216.15).^(-11.388);
end

den=P./(0.2869.*T);

D_vis=((0.000001458).*(T.^(1.5)))/(T+110.4);
a=sqrt(1.4.*287.*T);

end