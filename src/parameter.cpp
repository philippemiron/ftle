#include "parameter.h"

parameter::parameter(const char* fichier)
{
FILE* fp = fopen(fichier,"r");
if(fp == NULL) {
	printf("Le fichier de configuration, %s, est absent.\n", fichier);
	exit(0);
}
 
char token [200];

// Valeur par défaut
Nx = 0;
Ny = 0;
Nz = 0;
do {
    fscanf(fp,"%s",token);
    if(strcmp(token,"T"          )==0) fscanf(fp,"%lf \n",&T);
    if(strcmp(token,"Npt"        )==0) fscanf(fp,"%d \n",&Npt);
    if(strcmp(token,"T0"         )==0) fscanf(fp,"%d \n",&T0);
    if(strcmp(token,"Nx"	       )==0) fscanf(fp,"%d \n",&Nx);
    if(strcmp(token,"Ny"	       )==0) fscanf(fp,"%d \n",&Ny);
    if(strcmp(token,"Nz"         )==0) fscanf(fp,"%d \n",&Nz);
    if(strcmp(token,"xmin")==0) fscanf(fp,"%lf \n", &xmin);
	if(strcmp(token,"xmax")==0) fscanf(fp,"%lf \n", &xmax);
	if(strcmp(token,"ymin")==0) fscanf(fp,"%lf \n", &ymin);
	if(strcmp(token,"ymax")==0) fscanf(fp,"%lf \n", &ymax);
	if(strcmp(token,"zmin")==0) fscanf(fp,"%lf \n", &zmin);
	if(strcmp(token,"zmax")==0) fscanf(fp,"%lf \n", &zmax);
	
    if(strcmp(token,"MthOde"   )==0) 
    {
     fscanf(fp,"%s",token);
     if(strcmp(token,"Euler")==0)  MthOde=euler;
     if(strcmp(token,"RK5")==0)  MthOde=rk5;
    }
    
    if(strcmp(token,"File_CG")==0) {
        fscanf(fp,"%s \n", token);
        output_file_cg = token; 
    }
    
    if(strcmp(token,"function_u")==0) {
        fscanf(fp,"%s \n", token);
        fu = token; 
    }
    if(strcmp(token,"function_v")==0) {
        fscanf(fp,"%s \n", token);
        fv = token; 
    }
    if(strcmp(token,"function_w")==0) {
        fscanf(fp,"%s \n", token);
        fw = token; 
    }
    if(strcmp(token,"function_dudx")==0) {
        fscanf(fp,"%s \n", token);
        fdudx = token; 
    }
    if(strcmp(token,"function_dudy")==0) {
        fscanf(fp,"%s \n", token);
        fdudy = token; 
    }
    if(strcmp(token,"function_dudz")==0) {
        fscanf(fp,"%s \n", token);
        fdudz = token; 
    }
    if(strcmp(token,"function_dvdx")==0) {
        fscanf(fp,"%s \n", token);
        fdvdx = token; 
    }
    if(strcmp(token,"function_dvdy")==0) {
        fscanf(fp,"%s \n", token);
        fdvdy = token; 
    }
    if(strcmp(token,"function_dvdz")==0) {
        fscanf(fp,"%s \n", token);
        fdvdz = token; 
    }
    if(strcmp(token,"function_dwdx")==0) {
        fscanf(fp,"%s \n", token);
        fdwdx = token; 
    }
    if(strcmp(token,"function_dwdy")==0) {
        fscanf(fp,"%s \n", token);
        fdwdy = token; 
    }
    if(strcmp(token,"function_dwdz")==0) {
        fscanf(fp,"%s \n", token);
        fdwdz = token; 
    }
} while(strcmp(token,"End")!=0  );

};

parameter::~parameter(void)
{
};
