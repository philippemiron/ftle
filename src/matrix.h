#ifndef _matrix_
#define _matrix_

using namespace std;
#include <math.h>
#include <vector>

template<typename T> using vector2d = std::vector<std::vector<T>>;
template<typename T> using vector3d = std::vector<std::vector<std::vector<T>>>;
template<typename T> using vector4d = std::vector<std::vector<std::vector<std::vector<T>>>>;
template<typename T> using vector5d = std::vector<std::vector<std::vector<std::vector<std::vector<T>>>>>;

void vecResize(vector2d<double> &vec, int N1, int N2);
void vecResize(vector3d<double> &vec, int N1, int N2, int N3);
void vecResize(vector4d<double> &vec, int N1, int N2, int N3, int N4);
void vecResize(vector5d<double> &vec, int N1, int N2, int N3, int N4, int N5);

void vecConvert(const vector<double> &vec, vector3d<double> &vec3);
void vecConvert(const vector<double> &vec, vector2d<double> &vec2);

double dot(const vector<double> &a, const vector<double> &b);
double distance(const vector<double> &a, const vector<double> &b);
double norm(const vector<double> &a);
vector<double> cross(const vector<double> &a, const vector<double> &b);
void vecReverse(vector<double> &a);

#endif