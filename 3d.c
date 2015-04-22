#include "3d.h"

double Determinant(double matrix[3][3])
{
  double det = 0;

  for(int i=0;i<3;++i)
  {
    det += ( matrix[1][(i + 1) % 3]  *  matrix[2][(i + 2) % 3]
           - matrix[1][(i + 2) % 3]  *  matrix[2][(i + 1) % 3] )  *  matrix[0][i];
  }

  return det;
}

void Inv3D(double matrix[3][3], double inverse[3][3])
{
  double det = Determinant(matrix);

  for(int i=1;i<4;++i)
  {
    for(int j=1;j<4;++j)
    {
      inverse[i-1][j-1] = (
          matrix[i % 3][j % 3]  *  matrix[(i + 1) % 3][(j + 1) % 3]
        - matrix[i % 3][(j + 1) % 3]  *  matrix[(i + 1) % 3][j % 3]) / det;
    }
  }
}

double VecDot(double v1[3], double v2[3])
{
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

double VecLen(double v[3])
{
  return sqrt(VecDot(v, v));
}

void VecNormalize(double v[3])
{
  double l = VecLen(v);

  for (int i = 0; i < 3; ++i)
  {
    v[i] /= l;
  }
}

void VecShift(double pnt[3], double vec[3], double c)
{
  for(int i = 0; i < 3; ++i)
  {
    pnt[i] += vec[i] * c;
  }
}

void VecScale(double v[3], double c)
{
  for (int i = 0; i < 3; ++i)
  {
    v[i] *= c;
  }
}

// TODO Not really 3D.
double WeightedAverage(double a, double b, double weight)
{
  return (1 - weight) * a + weight * b;
}

// vi: foldmethod=syntax
