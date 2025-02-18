#include "mex.h"
#include <math.h>

/* M_PI가 정의되어 있지 않은 경우 정의 */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* 
   Vortex_Vatistas 함수 (MATLAB 코드와 동일한 계산식)
   입력:
      A[3]             : 세그먼트 시작점 좌표
      B[3]             : 세그먼트 종료점 좌표
      ColocationPoint[3]: 유도속도를 계산할 대상점 좌표
      vortexStrength   : 해당 세그먼트의 소용돌이 강도
      rc               : 코어 반경
      n                : Vatistas 모델의 지수
   출력:
      Vout[3]          : 대상점에서 유도된 속도 (3성분 벡터)
*/
static void VortexVatistas(const double A[3], const double B[3], const double ColocationPoint[3],
                           double vortexStrength, double rc, double n, double Vout[3])
{
    double r1[3], r2[3];
    int i;
    for (i = 0; i < 3; i++) {
        r1[i] = ColocationPoint[i] - A[i];
        r2[i] = ColocationPoint[i] - B[i];
    }
    
    double r1s = sqrt(r1[0]*r1[0] + r1[1]*r1[1] + r1[2]*r1[2]);
    double r2s = sqrt(r2[0]*r2[0] + r2[1]*r2[1] + r2[2]*r2[2]);
    
    double d[3];
    for (i = 0; i < 3; i++) {
        d[i] = B[i] - A[i];
    }
    double dL = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
    if (dL < 1e-10) {
        Vout[0] = Vout[1] = Vout[2] = 0.0;
        return;
    }
    
    /* 계산: cosA = dot(r1, B-A)/(norm(r1)*dL) */
    if (r1s < 1e-10 || r2s < 1e-10) {
        Vout[0] = Vout[1] = Vout[2] = 0.0;
        return;
    }
    double dotA = r1[0]*d[0] + r1[1]*d[1] + r1[2]*d[2];
    double dotB = r2[0]*d[0] + r2[1]*d[1] + r2[2]*d[2];
    double cosA = dotA/(r1s*dL);
    double cosB = dotB/(r2s*dL);
    
    /* r = norm(cross(r1, r2))/dL */
    double cross[3];
    cross[0] = r1[1]*r2[2] - r1[2]*r2[1];
    cross[1] = r1[2]*r2[0] - r1[0]*r2[2];
    cross[2] = r1[0]*r2[1] - r1[1]*r2[0];
    double cross_norm = sqrt(cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]);
    double r_val = cross_norm/dL;
    
    if (r_val < 1e-6) {
        Vout[0] = Vout[1] = Vout[2] = 0.0;
        return;
    }
    
    /* r_norm는 계산되지만 이후 사용되지는 않음: double r_norm = r_val/rc; */
    double numerator = r_val * r_val;  /* r^2 */
    double denom_inner = pow(rc, 2*n) + pow(r_val, 2*n);
    double denom = pow(denom_inner, 1.0/n);
    double vortexFactor = (vortexStrength/(4.0 * M_PI * r_val)) * (numerator/denom);
    
    double factor = vortexFactor * (cosA - cosB);
    /* 단위 벡터 (cross) */
    double cross_unit[3];
    cross_unit[0] = cross[0] / cross_norm;
    cross_unit[1] = cross[1] / cross_norm;
    cross_unit[2] = cross[2] / cross_norm;
    
    for (i = 0; i < 3; i++) {
        Vout[i] = factor * cross_unit[i];
    }
}

/* mexFunction: MEX 파일의 진입점 
   입력:
      prhs[0] : WakeGeom (n×(3·m))
      prhs[1] : GammaMatrix (n×m)
      prhs[2] : rc (1×m 또는 m×1)
   출력:
      plhs[0]: VelocityMatrix (n×(3·m))
*/
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* 입력 개수 체크 */
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
    mwSize wakeSize = mxGetM(prhs[0]);   /* 행 수: wake의 개수 */
    mwSize geomCols = mxGetN(prhs[0]);     /* 열 수: 3*m */
    if (geomCols % 3 != 0) {
        mexErrMsgIdAndTxt("WakeVortexCalculator_mex:invalidInput",
                          "WakeGeom의 열 수는 3의 배수여야 합니다.");
    }
    int geomPoints = (int)(geomCols / 3);
    
    /* 입력 1: GammaMatrix */
    double *gammaMatrix = mxGetPr(prhs[1]);
    mwSize gammaRows = mxGetM(prhs[1]);
    mwSize gammaCols = mxGetN(prhs[1]);
    if (gammaRows != wakeSize || gammaCols != (mwSize)geomPoints) {
        mexErrMsgIdAndTxt("WakeVortexCalculator_mex:invalidInput",
                          "GammaMatrix는 WakeGeom과 동일한 행 수 및 (#열 = 기하점 수)여야 합니다.");
    }
    
    /* 입력 2: rc */
    double *rc_ptr = mxGetPr(prhs[2]);
    mwSize rc_num = mxGetNumberOfElements(prhs[2]);
    if (rc_num != (mwSize)geomPoints) {
        mexErrMsgIdAndTxt("WakeVortexCalculator_mex:invalidInput",
                          "rc의 길이는 기하점 수와 같아야 합니다.");
    }
    
    /* 출력: VelocityMatrix (크기: wakeSize x (3*geomPoints)) */
    plhs[0] = mxCreateDoubleMatrix(wakeSize, geomCols, mxREAL);
    double *velocityMatrix = mxGetPr(plhs[0]);
    
    int i, j, lw, g;
    /* 각 Wake (행)과 기하점 (열 그룹)에 대해 */
    for (i = 0; i < (int)wakeSize; i++) {
        for (j = 0; j < geomPoints; j++) {
            double Vind_Point[3];
            /* 현재 wake 행 i의 j번째 기하점 좌표 추출 */
            Vind_Point[0] = wakeGeom[i + (3 * j + 0) * wakeSize];
            Vind_Point[1] = wakeGeom[i + (3 * j + 1) * wakeSize];
            Vind_Point[2] = wakeGeom[i + (3 * j + 2) * wakeSize];
            
            double vind[3] = {0.0, 0.0, 0.0};
            
            /* trailing wake: 각 세그먼트(lw)와 각 기하점(g)에 대해 */
            for (lw = 0; lw < (int)wakeSize - 1; lw++) {
                for (g = 0; g < geomPoints; g++) {
                    double LocalGamma = gammaMatrix[lw + g * wakeSize];
                    
                    double wakePoint1[3], wakePoint2[3];
                    wakePoint1[0] = wakeGeom[lw + (3 * g + 0) * wakeSize];
                    wakePoint1[1] = wakeGeom[lw + (3 * g + 1) * wakeSize];
                    wakePoint1[2] = wakeGeom[lw + (3 * g + 2) * wakeSize];
                    
                    wakePoint2[0] = wakeGeom[(lw + 1) + (3 * g + 0) * wakeSize];
                    wakePoint2[1] = wakeGeom[(lw + 1) + (3 * g + 1) * wakeSize];
                    wakePoint2[2] = wakeGeom[(lw + 1) + (3 * g + 2) * wakeSize];
                    
                    double v_temp[3];
                    VortexVatistas(wakePoint1, wakePoint2, Vind_Point,
                                   LocalGamma, rc_ptr[g], 1.06, v_temp);
                    
                    /* 누적 */
                    vind[0] += v_temp[0];
                    vind[1] += v_temp[1];
                    vind[2] += v_temp[2];
                }
            }
            
            /* 결과를 출력 행렬에 저장 */
            velocityMatrix[i + (3 * j + 0) * wakeSize] = vind[0];
            velocityMatrix[i + (3 * j + 1) * wakeSize] = vind[1];
            velocityMatrix[i + (3 * j + 2) * wakeSize] = vind[2];
        }
    }
}
