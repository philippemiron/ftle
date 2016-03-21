// Construction and destruction of 2D, 3D and 4D matrices.
#include "matrix.h"

// Vector related functions
void vecResize(vector2d<double> &vec, int N1, int N2) {
  vec.resize(N1);
  for (int i = 0; i < N1; i++)
    vec[i].resize(N2);
};

void vecResize(vector3d<double> &vec, int N1, int N2, int N3) {
  vec.resize(N1);
  for (int i(0); i < N1; i++) {
    vec[i].resize(N2);
    for (int j(0); j < N2; j++)
      vec[i][j].resize(N3);
  }
};

void vecResize(vector4d<double> &vec, int N1, int N2, int N3, int N4) {
  vec.resize(N1);
  for (int i(0); i < N1; i++) {
    vec[i].resize(N2);
    for (int j(0); j < N2; j++) {
      vec[i][j].resize(N3);
      for (int k(0); k < N3; k++)
        vec[i][j][k].resize(N4);
    }
  }
};

void vecResize(vector5d<double> &vec, int N1, int N2, int N3, int N4, int N5) {
  vec.resize(N1);
  for (int i(0); i < N1; i++) {
    vec[i].resize(N2);
    for (int j(0); j < N2; j++) {
      vec[i][j].resize(N3);
      for (int k(0); k < N3; k++) {
        vec[i][j][k].resize(N4);
        for (int l(0); l < N4; l++)
          vec[i][j][k][l].resize(N5);
      }
    }
  }
};

void vecConvert(const vector<double> &vec, vector2d<double> &vec2) {
  int s1 = vec2.size();
  int s2 = vec2[0].size();
  for (int i(0); i < s1; i++)
    for (int j(0); j < s2; j++)
      vec2[i][j] = vec[i * s2 + j];
};

void vecConvert(const vector<double> &vec, vector3d<double> &vec3) {
  int s1 = vec3.size();
  int s2 = vec3[0].size();
  int s3 = vec3[0][0].size();
  for (int i(0); i < s1; i++)
    for (int j(0); j < s2; j++)
      for (int k(0); k < s3; k++)
        vec3[i][j][k] = vec[i * s2 * s3 + j * s3 + k];
};

// vector operations
// dot product
double dot(const vector<double> &a, const vector<double> &b) {
  double out(0.0);
  for (size_t i(0); i < a.size(); i++)
    out += a[i] * b[i];
  return out;
};

// distance between two points
double distance(const vector<double> &a, const vector<double> &b) {
  double length(0.0);
  for (size_t i(0); i < a.size(); i++)
    length += pow(a[i] - b[i], 2);
  return sqrt(length);
};

// norm
double norm(const vector<double> &a) {
  double value(0.0);
  for (size_t i(0); i < a.size(); i++)
    value += pow(a[i], 2);
  return sqrt(value);
};

// cross product
vector<double> cross(const vector<double> &a, const vector<double> &b) {
  vector<double> out(a.size(), 0.0);

  // vector component
  out[0] = a[1] * b[2] - a[2] * b[1];
  out[1] = a[2] * b[0] - a[0] * b[2];
  out[2] = a[0] * b[1] - a[1] * b[0];

  return out;
};

// Reverse vector orientation
void vecReverse(vector<double> &a) {
  for (size_t i(0); i < a.size(); i++)
    a[i] = -a[i];
};