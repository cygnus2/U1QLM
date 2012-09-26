/* Cluster Algorithm for the U(1) quantum link model */

#include<stdio.h>
#include<math.h>
#include "mt19937ar.h"

/* Parameters                      *
 * (nx,ny,nt) : lattice dimension  *
 * j          : coupling constant  *
 * lam        : potential energy   *
 * beta       : inverse temperature*
 * eps        : Trotter dscrt      *
 * ieq        : equilibrium iter   *
 * imeas      : # of measurements  *
 * iskp       : # of iter/meas     */

#define nx 30
#define ny 30
#define nt 30
const int nlink=2*nx*ny*nt;
const int nlink2=nx*ny*nt;
const int ncube=nx*ny*nt/2;

double p1,p2;
double j,lam,beta,eps;
int mark[ncube],markl[nlink];
double hist[31];

/* pointer variables for neighbours */
int neigh[nlink][7][2], icube[nlink][2];
int ixc[nlink],iyc[nlink],itc[nlink],idir[nlink];
int il[nlink2][2];

/* Key for Mersenne Twister */
unsigned long init[4]={0x123, 0x234, 0x345, 0x456};
int length=4;

/* Field variables */
int ising[nlink];

int main(){
    
   extern void neighbours();

   int ieq,imeas,iskp;
   /* Set parameters */
   j=1.0;
   lam = 0.0;
   beta=3.6;
   eps=2.0/((double)nt);
   ieq=100;
   imeas=100;
   iskp=1;
   printf("Cluster Algorithm for the U(1) quantum link model\n");
   printf("Nx=%d, Ny=%d, Nt=%d\n",nx,ny,nt);
   printf("beta=%2.3f; J=%2.3f; lam=%2.3f\n",beta,j,lam);

   /* Define the probabilities */
   double x,p1,p2,coshx,sinhx;
   x     = eps*beta*j
   coshx = (exp(x)+exp(-x))/2.0;
   sinhx = (exp(x)-exp(-x))/2.0;
   p1    = exp(-eps*lam)*coshx;
   p2    = exp(-eps*lam)*sinhx;

   /* initialize random no */
   init_by_array(init,length);

   /* initialize neighbours */
   neighbours();

   return 0;
 }

 int convert(int x,int y,int t){
   int n;
   n = t*ny*nx + y*nx + x;
   return n;
 }

void neighbours(){
 extern int convert(int,int,int);
 int i;
 int id,ix,iy,it,n;
 int ixp1,ixm1,iyp1,iym1,itp1,itm1;
 int ilxyt1,ilxyt2;
 /* relate cartesian and lattice co-ordinates */
 i=0;
 for(id=0;id<2;id++){
 for(it=0;it<nt;it++){
 for(iy=0;iy<ny;iy++){
 for(ix=0;ix<nx;ix++){
  il[convert(ix,iy,it)][id]=i;
  ixc[i]=ix;
  iyc[i]=iy;
  itc[i]=it;
  idir[i]=id;
  i++;
 } } } } 
 /* define the neighbours */
 for(ix=0;ix<nx;ix++){
  ixp1=(ix+1)%nx;
  ixm1=(ix-1+nx)%nx;
  for(iy=0;iy<ny;iy++){
   iyp1=(iy+1)%ny;
   iym1=(iy-1+ny)%ny;
   for(it=0;it<nt;it++){
    itp1=(it+1)%nt;
    itm1=(it-1+nt)%nt;

    ilxyt1=il[convert(ix,iy,it)][0];
    ilxyt2=il[convert(ix,iy,it)][1];

/* first if starts */
     if(it%2==0){ /* select the even time slices */
/* second if starts */
     if((ix+iy)%2==0){
      /* fwd neighbours */
      neigh[ilxyt1][0][0] = il[convert(ixp1,iy,it)][1];
      neigh[ilxyt1][1][0] = il[convert(ix,iyp1,it)][0];
      neigh[ilxyt1][2][0] = il[convert(ix,iy,it)][1];
      neigh[ilxyt1][3][0] = il[convert(ix,iy,itp1)][0];
      neigh[ilxyt1][4][0] = il[convert(ixp1,iy,itp1)][1];
      neigh[ilxyt1][5][0] = il[convert(ix,iyp1,itp1)][0];
      neigh[ilxyt1][6][0] = il[convert(ix,iy,itp1)][1];
      neigh[ilxyt2][0][0] = il[convert(ix,iyp1,it)][0];
      neigh[ilxyt2][1][0] = il[convert(ixp1,iy,it)][1];
      neigh[ilxyt2][2][0] = il[convert(ix,iy,it)][0];
      neigh[ilxyt2][3][0] = il[convert(ix,iy,itp1)][1];
      neigh[ilxyt2][4][0] = il[convert(ix,iyp1,itp1)][0];
      neigh[ilxyt2][5][0] = il[convert(ixp1,iy,itp1)][1];
      neigh[ilxyt2][6][0] = il[convert(ix,iy,itp1)][0];
      /* bwd neighbours */
      neigh[ilxyt1][0][1] = il[convert(ix,iym1,it)][1];
      neigh[ilxyt1][1][1] = il[convert(ix,iym1,it)][0];
      neigh[ilxyt1][2][1] = il[convert(ixp1,iym1,it)][1];
      neigh[ilxyt1][3][1] = il[convert(ix,iy,itm1)][0];
      neigh[ilxyt1][4][1] = il[convert(ix,iym1,itm1)][1];
      neigh[ilxyt1][5][1] = il[convert(ix,iym1,itm1)][0];
      neigh[ilxyt1][6][1] = il[convert(ixp1,iym1,itm1)][1];
      neigh[ilxyt2][0][1] = il[convert(ixm1,iy,it)][0];
      neigh[ilxyt2][1][1] = il[convert(ixm1,iy,it)][1];
      neigh[ilxyt2][2][1] = il[convert(ixm1,iyp1,it)][0];
      neigh[ilxyt2][3][1] = il[convert(ix,iy,itm1)][1];
      neigh[ilxyt2][4][1] = il[convert(ixm1,iy,itm1)][0];
      neigh[ilxyt2][5][1] = il[convert(ixm1,iy,itm1)][1];
      neigh[ilxyt2][6][1] = il[convert(ixm1,iyp1,itm1)][0];
    }
    else{
      /* fwd neighbours */
      neigh[ilxyt1][0][0] = il[convert(ix,iym1,it)][1];
      neigh[ilxyt1][1][0] = il[convert(ix,iym1,it)][0];
      neigh[ilxyt1][2][0] = il[convert(ixp1,iym1,it)][1];
      neigh[ilxyt1][3][0] = il[convert(ix,iy,itp1)][0];
      neigh[ilxyt1][4][0] = il[convert(ix,iym1,itp1)][1];
      neigh[ilxyt1][5][0] = il[convert(ix,iym1,itp1)][0];
      neigh[ilxyt1][6][0] = il[convert(ixp1,iym1,itp1)][1];
      neigh[ilxyt2][0][0] = il[convert(ixm1,iy,it)][0];
      neigh[ilxyt2][1][0] = il[convert(ixm1,iy,it)][1];
      neigh[ilxyt2][2][0] = il[convert(ixm1,iyp1,it)][0];
      neigh[ilxyt2][3][0] = il[convert(ix,iy,itp1)][1];
      neigh[ilxyt2][4][0] = il[convert(ixm1,iy,itp1)][0];
      neigh[ilxyt2][5][0] = il[convert(ixm1,iy,itp1)][1];
      neigh[ilxyt2][6][0] = il[convert(ixm1,iyp1,itp1)][0];
      /* bwd neighbours */
      neigh[ilxyt1][0][1] = il[convert(ixp1,iy,it)][1];
      neigh[ilxyt1][1][1] = il[convert(ix,iyp1,it)][0];
      neigh[ilxyt1][2][1] = il[convert(ix,iy,it)][1];
      neigh[ilxyt1][3][1] = il[convert(ix,iy,itm1)][0];
      neigh[ilxyt1][4][1] = il[convert(ixp1,iy,itm1)][1];
      neigh[ilxyt1][5][1] = il[convert(ix,iyp1,itm1)][0];
      neigh[ilxyt1][6][1] = il[convert(ix,iy,itm1)][1];
      neigh[ilxyt2][0][1] = il[convert(ix,iyp1,it)][0];
      neigh[ilxyt2][1][1] = il[convert(ixp1,iy,it)][1];
      neigh[ilxyt2][2][1] = il[convert(ix,iy,it)][0];
      neigh[ilxyt2][3][1] = il[convert(ix,iy,itm1)][1];
      neigh[ilxyt2][4][1] = il[convert(ix,iyp1,itm1)][0];
      neigh[ilxyt2][5][1] = il[convert(ixp1,iy,itm1)][1];
      neigh[ilxyt2][6][1] = il[convert(ix,iy,itm1)][0];
    }
/* second if over */
   }
/* else for first if */
   else{
/* third if starts */
    if((ix+iy)%2==1){
      /* fwd neighbours */
      neigh[ilxyt1][0][0] = il[convert(ixp1,iy,it)][1];
      neigh[ilxyt1][1][0] = il[convert(ix,iyp1,it)][0];
      neigh[ilxyt1][2][0] = il[convert(ix,iy,it)][1];
      neigh[ilxyt1][3][0] = il[convert(ix,iy,itp1)][0];
      neigh[ilxyt1][4][0] = il[convert(ixp1,iy,itp1)][1];
      neigh[ilxyt1][5][0] = il[convert(ix,iyp1,itp1)][0];
      neigh[ilxyt1][6][0] = il[convert(ix,iy,itp1)][1];
      neigh[ilxyt2][0][0] = il[convert(ix,iyp1,it)][0];
      neigh[ilxyt2][1][0] = il[convert(ixp1,iy,it)][1];
      neigh[ilxyt2][2][0] = il[convert(ix,iy,it)][0];
      neigh[ilxyt2][3][0] = il[convert(ix,iy,itp1)][1];
      neigh[ilxyt2][4][0] = il[convert(ix,iyp1,itp1)][0];
      neigh[ilxyt2][5][0] = il[convert(ixp1,iy,itp1)][1];
      neigh[ilxyt2][6][0] = il[convert(ix,iy,itp1)][0];
      /* bwd neighbours */
      neigh[ilxyt1][0][1] = il[convert(ix,iym1,it)][1];
      neigh[ilxyt1][1][1] = il[convert(ix,iym1,it)][0];
      neigh[ilxyt1][2][1] = il[convert(ixp1,iym1,it)][1];
      neigh[ilxyt1][3][1] = il[convert(ix,iy,itm1)][0];
      neigh[ilxyt1][4][1] = il[convert(ix,iym1,itm1)][1];
      neigh[ilxyt1][5][1] = il[convert(ix,iym1,itm1)][0];
      neigh[ilxyt1][6][1] = il[convert(ixp1,iym1,itm1)][1];
      neigh[ilxyt2][0][1] = il[convert(ixm1,iy,it)][0];
      neigh[ilxyt2][1][1] = il[convert(ixm1,iy,it)][1];
      neigh[ilxyt2][2][1] = il[convert(ixm1,iyp1,it)][0];
      neigh[ilxyt2][3][1] = il[convert(ix,iy,itm1)][1];
      neigh[ilxyt2][4][1] = il[convert(ixm1,iy,itm1)][0];
      neigh[ilxyt2][5][1] = il[convert(ixm1,iy,itm1)][1];
      neigh[ilxyt2][6][1] = il[convert(ixm1,iyp1,itm1)][0];
    }
    else{
      /* fwd neighbours */
      neigh[ilxyt1][0][0] = il[convert(ix,iym1,it)][1];
      neigh[ilxyt1][1][0] = il[convert(ix,iym1,it)][0];
      neigh[ilxyt1][2][0] = il[convert(ixp1,iym1,it)][1];
      neigh[ilxyt1][3][0] = il[convert(ix,iy,itp1)][0];
      neigh[ilxyt1][4][0] = il[convert(ix,iym1,itp1)][1];
      neigh[ilxyt1][5][0] = il[convert(ix,iym1,itp1)][0];
      neigh[ilxyt1][6][0] = il[convert(ixp1,iym1,itp1)][1];
      neigh[ilxyt2][0][0] = il[convert(ixm1,iy,it)][0];
      neigh[ilxyt2][1][0] = il[convert(ixm1,iy,it)][1];
      neigh[ilxyt2][2][0] = il[convert(ixm1,iyp1,it)][0];
      neigh[ilxyt2][3][0] = il[convert(ix,iy,itp1)][1];
      neigh[ilxyt2][4][0] = il[convert(ixm1,iy,itp1)][0];
      neigh[ilxyt2][5][0] = il[convert(ixm1,iy,itp1)][1];
      neigh[ilxyt2][6][0] = il[convert(ixm1,iyp1,itp1)][0];
      /* bwd neighbours */
      neigh[ilxyt1][0][1] = il[convert(ixp1,iy,it)][1];
      neigh[ilxyt1][1][1] = il[convert(ix,iyp1,it)][0];
      neigh[ilxyt1][2][1] = il[convert(ix,iy,it)][1];
      neigh[ilxyt1][3][1] = il[convert(ix,iy,itm1)][0];
      neigh[ilxyt1][4][1] = il[convert(ixp1,iy,itm1)][1];
      neigh[ilxyt1][5][1] = il[convert(ix,iyp1,itm1)][0];
      neigh[ilxyt1][6][1] = il[convert(ix,iy,itm1)][1];
      neigh[ilxyt2][0][1] = il[convert(ix,iyp1,it)][0];
      neigh[ilxyt2][1][1] = il[convert(ixp1,iy,it)][1];
      neigh[ilxyt2][2][1] = il[convert(ix,iy,it)][0];
      neigh[ilxyt2][3][1] = il[convert(ix,iy,itm1)][1];
      neigh[ilxyt2][4][1] = il[convert(ix,iyp1,itm1)][0];
      neigh[ilxyt2][5][1] = il[convert(ixp1,iy,itm1)][1];
      neigh[ilxyt2][6][1] = il[convert(ix,iy,itm1)][0];
    }
/* third if over */
   }
/* closing for else corresponding to first if */
/* next 3 braces closing for the three for loops */
   }
  }
 }

  /* label cubes */
  i=0;
  for(it=0;it<nt;it++){
   itp1 = (it+1)%nt;
   for(iy=0;iy<ny;iy++){
    iyp1 = (iy+1)%ny;
    for(ix=0;ix<nx;ix++){
     ixp1 = (ix+1)%nx;
     if((ix+iy+it)%2==0){
      icube[il[convert(ix,iy,it)][0]][0] = i;
      icube[il[convert(ix,iy,it)][1]][0] = i;
      icube[il[convert(ixp1,iy,it)][1]][0]=i;
      icube[il[convert(ix,iyp1,it)][0]][0]=i;
      icube[il[convert(ix,iy,itp1)][0]][1]=i;
      icube[il[convert(ix,iy,itp1)][1]][1]=i;
      icube[il[convert(ixp1,iy,itp1)][1]][1]=i;
      icube[il[convert(ix,iyp1,itp1)][0]][1]=i;
      i++;
     }
    }
   }
  }
 }

/* This makes a random start for all links */
void ranstrt(){
  int ix,iy,id,it;
  int ilink;
  for(ix=0;ix<nx;ix++){
  for(iy=0;iy<ny;iy++){
  for(id=0;id<2;id++){
    ilink=0;
    if(genrand_real3()>0.5) ilink=1;
    for(it=0;it<nt;it++){
     ising[il[convert(ix,iy,it)][id]]=ilink;
    }
  } } }
}

/* This makes a twist start for starting configuration */
void twiststrt(){
  int ix,iy,it;
  for(ix=0;ix<nx;ix++){
  for(iy=0;iy<ny;iy++){
  for(it=0;it<nt;it++){
   if((ix+iy+it)%2==0){
    ising[il[convert(ix,iy,it)][0]]=1;
    ising[il[convert(ix,iy,it)][1]]=0;
   }
   else{
    ising[il[convert(ix,iy,it)][0]]=0;
    ising[il[convert(ix,iy,it)][1]]=1;
   }
  } } }
}

/* Check for forbidden configurations */
void check(int *nontriv, int *nwind){
  int il,ic;
  int il1,il2,il3,il4,il5,il6,il7;
  nontriv=0;
  nwind=0;
  for(il=0;il<nlink;il++){
  for(ic=0;ic<2;ic++){
   il1=neigh[il][1][ic];
   il2=neigh[il][2][ic];
   il3=neigh[il][3][ic];
   il4=neigh[il][4][ic];
   il5=neigh[il][5][ic];
   il6=neigh[il][6][ic];
   il7=neigh[il][7][ic];
   if((ising[il]==ising[il4])&&(ising[il1]==ising[il5])&&(ising[il2]==ising[il6])&&(ising[il3]==ising[il7]))
    continue;
   else{ /* count the non-trivial cubes */
    if((ising[il]==ising[il1])&&(ising[il]!=ising[il2])&&(ising[il]!=ising[il3])&&(ising[il]!=ising[il4])&&
       (ising[il]!=ising[il5])&&(ising[il]==ising[il6])&&(ising[il]==ising[il7])){
        nontriv++;
	/* calculate the winding number */
        if((ixc[il]==0)&&(iyc[il]==0)&&(itc[il]%2==0)&&(idir[il]==0)&&(ic==1))
          nwind += 2*ising[il]-1;
       }
     else{
       printf("Forbidden configuration\n");
       exit(0);
     }
   }
  }
  }
}
