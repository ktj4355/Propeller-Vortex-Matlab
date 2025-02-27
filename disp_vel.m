

clc
clear
close all
%Validation(100)

%%
% 
alt=0;
[Tmp, Pressure, rho, D_vis, a] = STD_Atm(alt);

vFree= [0,0,0];
PtFree=Pressure.*1000+0.5*rho.*(norm(vFree)).^2
inputGeom=readmatrix("Geometry.xlsx");
Geom_chord=inputGeom(:,3);
rc_Ratio=0.5
rc_Geom=Geom_chord.*rc_Ratio;
vortex_n=1.06

[filename, path] =uigetfile(".mat")

fullname=[path,filename]

load(fullname)


globalData=OutputVortexSturcture{1}
Wake_Geom_Position=OutputVortexSturcture{2}
Wake2_Geom_Position=OutputVortexSturcture{3}
Wake_Gamma=OutputVortexSturcture{4}
Wake2_Gamma=OutputVortexSturcture{5}
rc_panel=OutputVortexSturcture{6}


%% 최적화된 Vortex 유도속도 계산 및 시각화 코드

% 1. 계산할 좌표 그리드 생성 (y는 0으로 고정)
% xy plane
xpos = linspace(-0.2, 0.2, 100);
ypos = linspace(-0.2, 0.2, 100);
[X, Y] = meshgrid(xpos, ypos);
zpos=-0.1
Z = zpos+zeros(size(X));  % ypos = 0
% xz plane
xpos = linspace(-0.2, 0.2, 100);
zpos = linspace(-0.2, 0.2, 100);
[X, Z] = meshgrid(xpos, zpos);
ypos=0
Y = ypos+zeros(size(X));  % ypos = 0

% 모든 점을 N×3 행렬로 정리 (각 행: 하나의 colocation point)
points = [X(:), Y(:), Z(:)];

% 2. 유도속도 누적 변수 초기화 (각 colocation point에 대해 Nx3)
vWake_total = zeros(size(points));

% 3. Wake 계산 (여기서는 wake 데이터가 여러 blade segment에 대해 주어진다고 가정)
% Wake_Gamma, Wake_Geom_Position, Wake2_Gamma, Wake2_Geom_Position,
% rc_panel, vortex_n, gamma_idx 등은 미리 정의되어 있어야 함
cnt=0;
clc
for idx = 1:size(Wake_Gamma, 2)
    cnt=cnt+1;
    progress=double(idx/size(Wake_Gamma, 2));
    if cnt>2
        clc
        fprintf("Processing....   %.2f%% \n",progress*100)
        cnt=0;
    end
    % 각 blade의 좌표 인덱스 (3개씩 묶음)
    xind = (idx-1)*3 + 1;
    yind = (idx-1)*3 + 2;
    zind = (idx-1)*3 + 3;

    % 유효한 wake 데이터가 있는지 확인
    if size(Wake_Geom_Position, 1) <= 2
        break;
    end

    % 각 wake segment에 대해 계산 (Tidx: 행 인덱스)
    for Tidx = 1:(size(Wake_Geom_Position, 1) - 1)
        % Blade 1: Wake_Geom_Position과 Wake_Gamma 이용
        A = Wake_Geom_Position(Tidx, xind:zind);
        B = Wake_Geom_Position(Tidx+1, xind:zind);
        gamma_val =Wake_Gamma(Tidx, idx);
        % 벡터화된 Vortex_Vatistas 호출 (모든 colocation point 한 번에 계산)
        vWake_total = vWake_total + Vortex_Vatistas_multiColoc(A, B, points, gamma_val, rc_Geom(idx), vortex_n);

        % Blade 2: Wake2_Geom_Position과 Wake2_Gamma 이용
        A = Wake2_Geom_Position(Tidx, xind:zind);
        B = Wake2_Geom_Position(Tidx+1, xind:zind);
        gamma_val = Wake2_Gamma(Tidx, idx);
        vWake_total = vWake_total + Vortex_Vatistas_multiColoc(A, B, points, gamma_val, rc_Geom(idx), vortex_n);
    end
end
vWake_total=vWake_total+vFree;
% 4. 결과 정리: 각 colocation point에 대해 [x, y, z, Vx, Vy, Vz, norm(V)]
vNorm = sqrt(sum(vWake_total(:,[1,3]).^2, 2));
vout = [points, vWake_total, vNorm];
%% dispaly

% 5. 시각화 (quiver plot: x-z 평면 상에 벡터의 x, z 성분 표시)
%  : X, Z는 meshgrid 형태, Vx는 4번째 열, Vz는 6번째 열
X_grid = reshape(vout(:,1), size(X));
Y_grid = reshape(vout(:,2), size(Y));
Z_grid = reshape(vout(:,3), size(Z));
Vx_grid = reshape(vout(:,4), size(X));
Vy_grid = reshape(vout(:,5), size(Y));

Vz_grid = reshape(vout(:,6), size(Z));
Vmag_grid = reshape(vout(:,7), size(X));  % 속도 크기
Vswarl_grid=sqrt(Vx_grid.^2+Vy_grid.^2)

VS_max=max(Vswarl_grid,[],'all')
VS_min=min(Vswarl_grid,[],'all')

P_grid=(PtFree-0.5*rho.*(Vmag_grid).^2)-Pressure.*1000;
Pmax=max(P_grid,[],'all')
Pmin=min(P_grid,[],'all')


figure(11)
clf
hold on
contourf(X_grid, Z_grid, Vmag_grid,1000,  "lineColor",'none');  % 50개의 contour level
colormap(turbo)
    clim("auto")

%clim([VS_min,VS_max])
axis equal
colorbar
title("Velocity Magnitude Contour(m/s)")


figure(12)
clf
hold on
quiver(X_grid, Z_grid, Vx_grid, Vz_grid,0.75,'k')
xlabel('x'); ylabel('z');
title('Vortex에 의한 유도속도 (x-z 평면)');
axis equal


figure(13)
clf
hold on

contourf(X_grid, Z_grid, P_grid,100,  "lineColor",'none');  % 50개의 contour level
colormap(turbo)
clim([Pmin,Pmax])
axis equal
colorbar
title("Gauge Pressure Contour(Pa)")

figure(14)
clf
hold on
contourf(X_grid, Z_grid, abs(Vz_grid),1000,  "lineColor",'none');  % 50개의 contour level
colormap(turbo)
    clim("auto")

%clim([VS_min,VS_max])
axis equal
colorbar
title("Axial Velocity Contour(m/s)")

%% test_vortex.m
% 원래 함수(Vortex_Vatistas_orig)와 벡터화 함수(Vortex_Vatistas)의 결과 비교
function [max_error]=Validation(seed)
rng(seed); % 재현성을 위한 난수 시드 설정
N = 100;  % 테스트할 colocation point 개수

% 임의의 A, B 좌표 (1x3)
A = rand(1,3);
B = rand(1,3);

% 임의의 Colocation Point (Nx3)
ColocationPoints = rand(N,3);

% 임의의 기타 스칼라 파라미터
vortexStrength = rand();
rc = 0.1 + rand();  % rc > 0
n = 2;

% 원래 함수를 loop로 각 점에 대해 계산
Vout_loop = zeros(N,3);
for i = 1:N
    Vout_loop(i,:) = Vortex_Vatistas(A, B, ColocationPoints(i,:), vortexStrength, rc, n);
end

% 벡터화된 함수로 한 번에 계산
Vout_vec = Vortex_Vatistas_multiColoc(A, B, ColocationPoints, vortexStrength, rc, n);

% 두 결과의 차이 계산
diff = Vout_loop - Vout_vec;
max_error = max(abs(diff(:)));
fprintf('최대 절대 오차: %e\n', max_error);

% 각 성분 비교 플롯 (원래 vs 벡터화)
figure;
subplot(1,3,1);
plot(Vout_loop(:,1), Vout_vec(:,1), 'bo');
xlabel('Original Vx'); ylabel('Vectorized Vx');
title('Vx 비교');

subplot(1,3,2);
plot(Vout_loop(:,2), Vout_vec(:,2), 'ro');
xlabel('Original Vy'); ylabel('Vectorized Vy');
title('Vy 비교');

subplot(1,3,3);
plot(Vout_loop(:,3), Vout_vec(:,3), 'go');
xlabel('Original Vz'); ylabel('Vectorized Vz');
title('Vz 비교');







end
function [Vout] = Vortex_Vatistas_multiColoc(A, B, ColocationPoint, vortexStrength, rc, n)
% Vortex_Vatistas 벡터화 버전
%   A, B: 1x3 벡터 (와이막 선분의 양끝)
%   ColocationPoint: Nx3 행렬 (여러 점)
%   vortexStrength, rc, n: 스칼라 값
%
%   각 ColocationPoint에서의 유도속도 벡터를 Nx3 행렬 Vout로 반환

% 각 점에 대한 r1, r2 (Nx3)
r1 = ColocationPoint - A;
r2 = ColocationPoint - B;

% 각 점에 대한 크기 (Nx1)
r1s = sqrt(sum(r1.^2, 2));
r2s = sqrt(sum(r2.^2, 2));
dL = norm(A - B);  % 스칼라

% (B-A) 벡터, 모든 점에 대해 동일
BA = B - A;

% 각 점에서의 cosA, cosB (Nx1)
cosA = sum(r1 .* repmat(BA, size(ColocationPoint,1), 1), 2) ./ (r1s * dL);
cosB = sum(r2 .* repmat(BA, size(ColocationPoint,1), 1), 2) ./ (r2s * dL);

% 행별로 cross product (Nx3) 및 그 크기 (Nx1)
cross_r1r2 = cross(r1, r2, 2);
norm_cross = sqrt(sum(cross_r1r2.^2, 2));

% r는 cross 크기를 dL로 나눈 값 (Nx1)
r = norm_cross / dL;

% 결과 초기화
Vout = zeros(size(ColocationPoint));  % Nx3

% 임계값보다 큰 r에 대해서만 계산 (0에 가까운 값은 0으로 처리)
valid = r >= 1e-6;

% 유효한 점에 대해 벡터 인자 계산
r_valid = r(valid);
vortexFactor = (vortexStrength ./ (4*pi*r_valid)) .* ...
    ((r_valid.^2) ./ ((rc^(2*n) + r_valid.^(2*n)).^(1/n)));

% 단위 벡터 (유효한 점만 계산)
unit_vec = zeros(size(cross_r1r2));
unit_vec(valid,:) = cross_r1r2(valid,:) ./ repmat(norm_cross(valid),1,3);

% 최종 유도속도 계산 (각 유효점마다)
Vout(valid,:) = (vortexFactor .* (cosA(valid) - cosB(valid))) .* unit_vec(valid,:);
end
