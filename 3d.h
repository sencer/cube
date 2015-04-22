#ifndef THREED_SIMPLE_MATH
#define THREED_SIMPLE_MATH 123456

#define MIN(a, b) ((a<b)?a:b)
#define MAX(a, b) ((a>b)?a:b)
#define DEG 57.2957795786
#define B2A 0.529177249

#include <math.h>
double Determinant(double matrix[3][3]);
void Inv3D(double matrix[3][3], double inverse[3][3]);
double VecDot(double v1[3], double v2[3]);
double VecLen(double v[3]);
void VecNormalize(double v[3]);
void VecShift(double pnt[3], double vec[3], double c);
void VecScale(double v[3], double c);
double WeightedAverage(double a, double b, double weight);
#endif
