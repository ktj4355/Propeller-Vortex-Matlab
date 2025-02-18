clc
clear
close all



% 예제 데이터 (크기 및 값은 실제 문제에 맞게 설정)
n = 200; m = 40;
WakeGeom = rand(n,3*m);
GammaMatrix = rand(n,m);
rc = rand(1,m);

% MEX 파일 호출 (컴파일된 파일 이름이 WakeVortexCalculator_mex)
timer1=tic;
VelocityMatrix1 = WakeVortexCalculator_mex(WakeGeom, GammaMatrix, rc);

T1=toc(timer1);
disp(T1)
timer2=tic;
VelocityMatrix2= WakeVortexCalculator_mex_Fast(WakeGeom, GammaMatrix, rc);
T2=toc(timer2);
disp(T2)
diffff=VelocityMatrix1-VelocityMatrix2;
max(max(diffff))