#include "parameter.h"

parameter::parameter(const char *file) {

  FILE *fp = fopen(file, "r");
  if (fp == NULL) {
    cerr << "Can't open config file " << file << ". Look at the example for help." << endl;
    exit(-1);
  }

  char token[1000];
  nx_ = 0;
  ny_ = 0;
  nz_ = 0;

  do {
    fscanf(fp, "%s", token);
    if (strcmp(token, "T0") == 0) fscanf(fp, "%lf \n", &t0_);
    if (strcmp(token, "T") == 0) fscanf(fp, "%lf \n", &t_);
    if (strcmp(token, "tolerance") == 0) fscanf(fp, "%lf \n", &tolerance_);
    if (strcmp(token, "hmin") == 0) fscanf(fp, "%lf \n", &hmin_);
    if (strcmp(token, "hmax") == 0) fscanf(fp, "%lf \n", &hmax_);
    if (strcmp(token, "nx") == 0) fscanf(fp, "%d \n", &nx_);
    if (strcmp(token, "ny") == 0) fscanf(fp, "%d \n", &ny_);
    if (strcmp(token, "nz") == 0) fscanf(fp, "%d \n", &nz_);
    if (strcmp(token, "xmin") == 0) fscanf(fp, "%lf \n", &xmin_);
    if (strcmp(token, "xmax") == 0) fscanf(fp, "%lf \n", &xmax_);
    if (strcmp(token, "ymin") == 0) fscanf(fp, "%lf \n", &ymin_);
    if (strcmp(token, "ymax") == 0) fscanf(fp, "%lf \n", &ymax_);
    if (strcmp(token, "zmin") == 0) fscanf(fp, "%lf \n", &zmin_);
    if (strcmp(token, "zmax") == 0) fscanf(fp, "%lf \n", &zmax_);

    if (strcmp(token, "function_u") == 0) {
      fscanf(fp, "%s \n", token);
      fu_ = token;
    }
    if (strcmp(token, "function_v") == 0) {
      fscanf(fp, "%s \n", token);
      fv_ = token;
    }
    if (strcmp(token, "function_w") == 0) {
      fscanf(fp, "%s \n", token);
      fw_ = token;
    }
    if (strcmp(token, "function_dudx") == 0) {
      fscanf(fp, "%s \n", token);
      fdudx_ = token;
    }
    if (strcmp(token, "function_dudy") == 0) {
      fscanf(fp, "%s \n", token);
      fdudy_ = token;
    }
    if (strcmp(token, "function_dudz") == 0) {
      fscanf(fp, "%s \n", token);
      fdudz_ = token;
    }
    if (strcmp(token, "function_dvdx") == 0) {
      fscanf(fp, "%s \n", token);
      fdvdx_ = token;
    }
    if (strcmp(token, "function_dvdy") == 0) {
      fscanf(fp, "%s \n", token);
      fdvdy_ = token;
    }
    if (strcmp(token, "function_dvdz") == 0) {
      fscanf(fp, "%s \n", token);
      fdvdz_ = token;
    }
    if (strcmp(token, "function_dwdx") == 0) {
      fscanf(fp, "%s \n", token);
      fdwdx_ = token;
    }
    if (strcmp(token, "function_dwdy") == 0) {
      fscanf(fp, "%s \n", token);
      fdwdy_ = token;
    }
    if (strcmp(token, "function_dwdz") == 0) {
      fscanf(fp, "%s \n", token);
      fdwdz_ = token;
    }
  } while (strcmp(token, "End") != 0);

};
