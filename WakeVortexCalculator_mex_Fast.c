#include "mex.h"
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* 
   Vortex_Vatistas 함수 (MATLAB 코드와 동일한 계산식)
   입력:
     A[3]              : 세그먼트 시작점 좌표
     B[3]              : 세그먼트 종료점 좌표
     ColocationPoint[3]: 유도속도를 계산할 대상점 좌표
     vortexStrength    : 해당 세그먼트의 소용돌이 강도
     rc                : 코어 반경
     n                 : Vatistas 모델의 지수
   출력:
     Vout[3]           : 대상점에서 유도된 속도 (3성분 벡터)
*/
static inline void VortexVatistas(const double A[3], const double B[3], const double ColocationPoint[3],
                                  double vortexStrength, double rc, double n, double Vout[3])
{
    double r1[3], r2[3], d[3];
    int k;
    for (k = 0; k < 3; k++) {
        r1[k] = ColocationPoint[k] - A[k];
        r2[k] = ColocationPoint[k] - B[k];
        d[k]  = B[k] - A[k];
    }
    
    double r1s = sqrt(r1[0]*r1[0] + r1[1]*r1[1] + r1[2]*r1[2]);
    double r2s = sqrt(r2[0]*r2[0] + r2[1]*r2[1] + r2[2]*r2[2]);
    double dL  = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
    if (dL < 1e-10 || r1s < 1e-10 || r2s < 1e-10) {
        Vout[0] = Vout[1] = Vout[2] = 0.0;
        return;
    }
    
    double dotA = r1[0]*d[0] + r1[1]*d[1] + r1[2]*d[2];
    double dotB = r2[0]*d[0] + r2[1]*d[1] + r2[2]*d[2];
    double cosA = dotA/(r1s*dL);
    double cosB = dotB/(r2s*dL);
    
    double cross[3];
    cross[0] = r1[1]*r2[2] - r1[2]*r2[1];
    cross[1] = r1[2]*r2[0] - r1[0]*r2[2];
    cross[2] = r1[0]*r2[1] - r1[1]*r2[0];
    double cross_norm = sqrt(cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]);
    double r_val = cross_norm / dL;
    if (r_val < 1e-6) {
        Vout[0] = Vout[1] = Vout[2] = 0.0;
        return;
    }
    
    double numerator = r_val * r_val;  /* r^2 */
    double denom_inner = pow(rc, 2*n) + pow(r_val, 2*n);
    double denom = pow(denom_inner, 1.0/n);
    double vortexFactor = (vortexStrength/(4.0 * M_PI * r_val)) * (numerator/denom);
    double factor = vortexFactor * (cosA - cosB);
    
    double cross_unit[3];
    cross_unit[0] = cross[0] / cross_norm;
    cross_unit[1] = cross[1] / cross_norm;
    cross_unit[2] = cross[2] / cross_norm;
    for (k = 0; k < 3; k++) {
        Vout[k] = factor * cross_unit[k];
    }
}

/*
   mexFunction: MEX 파일의 진입점
   입력:
     prhs[0]: WakeGeom (n x (3*m))
     prhs[1]: GammaMatrix (n x m)
     prhs[2]: rc (벡터 길이 m)
   출력:
     plhs[0]: VelocityMatrix (n x (3*m))
*/
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* 입력 개수 검사 */
    if (nrhs != 3) {
        mexErrMsgIdAndTxt("WakeVortexCalculator_mex:invalidNumInputs",
                          "세 개의 입력 (WakeGeom, GammaMatrix, rc)이 필요합니다.");
    }
    if (nlhs > 1) {
        mexErrMsgIdAndTxt("WakeVortexCalculator_mex:invalidNumOutputs",
                          "하나의 출력만 허용됩니다.");
    }
    
    /* 입력 0: WakeGeom */
    double *wakeGeom = mxGetPr(prhs[0]);
    mwSize wakeSize = mxGetM(prhs[0]);   /* n */
    mwSize geomCols = mxGetN(prhs[0]);     /* 3*m */
    if (geomCols % 3 != 0) {
        mexErrMsgIdAndTxt("WakeVortexCalculator_mex:invalidInput",
                          "WakeGeom의 열 수는 3의 배수여야 합니다.");
    }
    int nWake = (int)wakeSize;
    int nGeom = (int)(geomCols / 3);  /* m */
    
    /* 입력 1: GammaMatrix */
    double *gammaMatrix = mxGetPr(prhs[1]);
    if (mxGetM(prhs[1]) != wakeSize || mxGetN(prhs[1]) != (mwSize)nGeom) {
        mexErrMsgIdAndTxt("WakeVortexCalculator_mex:invalidInput",
                          "GammaMatrix는 WakeGeom과 동일한 행 수와 기하점 수를 가져야 합니다.");
    }
    
    /* 입력 2: rc */
    double *rc_ptr = mxGetPr(prhs[2]);
    if (mxGetNumberOfElements(prhs[2]) != (mwSize)nGeom) {
        mexErrMsgIdAndTxt("WakeVortexCalculator_mex:invalidInput",
                          "rc의 길이는 기하점 수와 같아야 합니다.");
    }
    
    /* 출력: VelocityMatrix (크기: nWake x (3*nGeom)) */
    plhs[0] = mxCreateDoubleMatrix(wakeSize, geomCols, mxREAL);
    double *velocityMatrix = mxGetPr(plhs[0]);
    
    int i, j, lw, g;
    int idx_x, idx_y, idx_z;
    
    /* OpenMP를 이용해 외부 루프(i)를 병렬화합니다.
       (내부 루프(j)는 일반 반복문으로 처리)
    */
    #pragma omp parallel for default(none) private(i, j, lw, g, idx_x, idx_y, idx_z) shared(nWake, nGeom, wakeGeom, gammaMatrix, rc_ptr, velocityMatrix)
    for (i = 0; i < nWake; i++) {
        for (j = 0; j < nGeom; j++) {
            double Vind_Point[3];
            idx_x = i + (3 * j + 0) * nWake;
            idx_y = i + (3 * j + 1) * nWake;
            idx_z = i + (3 * j + 2) * nWake;
            Vind_Point[0] = wakeGeom[idx_x];
            Vind_Point[1] = wakeGeom[idx_y];
            Vind_Point[2] = wakeGeom[idx_z];
            
            double vind[3] = {0.0, 0.0, 0.0};
            double v_temp[3];
            double wakePoint1[3], wakePoint2[3];
            for (lw = 0; lw < nWake - 1; lw++) {
                for (g = 0; g < nGeom; g++) {
                    double LocalGamma = gammaMatrix[lw + g * nWake];
                    
                    int idx1 = lw + (3 * g + 0) * nWake;
                    int idx2 = lw + (3 * g + 1) * nWake;
                    int idx3 = lw + (3 * g + 2) * nWake;
                    wakePoint1[0] = wakeGeom[idx1];
                    wakePoint1[1] = wakeGeom[idx2];
                    wakePoint1[2] = wakeGeom[idx3];
                    
                    int idx4 = (lw + 1) + (3 * g + 0) * nWake;
                    int idx5 = (lw + 1) + (3 * g + 1) * nWake;
                    int idx6 = (lw + 1) + (3 * g + 2) * nWake;
                    wakePoint2[0] = wakeGeom[idx4];
                    wakePoint2[1] = wakeGeom[idx5];
                    wakePoint2[2] = wakeGeom[idx6];
                    
                    VortexVatistas(wakePoint1, wakePoint2, Vind_Point,
                                   LocalGamma, rc_ptr[g], 1.06, v_temp);
                    vind[0] += v_temp[0];
                    vind[1] += v_temp[1];
                    vind[2] += v_temp[2];
                }
            }
            velocityMatrix[idx_x] = vind[0];
            velocityMatrix[idx_y] = vind[1];
            velocityMatrix[idx_z] = vind[2];
        }
    }
}