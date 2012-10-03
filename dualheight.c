/* Cluster Algorithm for the U(1) quantum link model */
/* using the representation of the dual height model */


#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>
#include "mt19937ar.h"

/* Parameters                      *
 * (LX,LY,LT) : lattice dimension  *
 * j          : coupling constant  *
 * lam        : potential energy   *
 * beta       : inverse temperature*
 * eps        : Trotter dscrt      *
 * ieq        : equilibrium iter   *
 * imeas      : # of measurements  *
 * iskp       : # of iter/meas     */

#define DIM 3
int LX,LY,LT,VOL,VOL2,VOL4;
double p1,p2;
double j,lam,beta,eps;
int SEED;
int nclus;

/* pointer variables for neighbours */
#define NNBR 10
int *neigh[NNBR];
/* checkerboard labeling for sites */
int *ixc,*iyc,*itc;
/* Field variables */
int *ising;

int main(){

   int ieq,imeas,iskp;
   int i,readval;    
   extern void neighbours();
   extern void neighchk();
   extern double ran();
   extern void cluster();

   FILE *fptr;
   char st[20];
   /* read file */
   fptr=fopen("QUEUE","r");
   if(fptr == NULL) {printf("QUEUE error.\n"); exit(1);}
   readval = fscanf(fptr,"%s %d\n",st,&LX);
   if(readval == -1) printf("Error\n");
   readval = fscanf(fptr,"%s %d\n",st,&LY);
   if(readval == -1) printf("Error\n");
   readval = fscanf(fptr,"%s %d\n",st,&LT);
   if(readval == -1) printf("Error\n");
   fclose(fptr);

   VOL = LX*LY*LT;
   VOL2= VOL/2;
   VOL4= VOL/4;

   /* allocate memory */
   ixc   = (int *)malloc(VOL*sizeof(int));
   iyc   = (int *)malloc(VOL*sizeof(int));
   itc   = (int *)malloc(VOL*sizeof(int));
   ising = (int *)malloc(VOL*sizeof(int));
   for(i=0;i<NNBR;i++)
     neigh[i] = (int *)malloc(VOL*sizeof(int));

   /* Set parameters */
   j=1.0;
   lam = 0.0;
   beta=3.6;
   eps=2.0/((double)LT);
   ieq=100;
   imeas=100;
   iskp=1;
   SEED=95540;

  /* Initialize MT */
  init_genrand(SEED);

  /* initialize neighbours */
  neighbours();

  printf("Cluster Algorithm for the U(1) quantum link model\n");
  printf("Nx=%d, Ny=%d, Nt=%d\n",LX,LY,LT);
  printf("beta=%2.3f; J=%2.3f; lam=%2.3f\n",beta,j,lam);
  printf("Starting seed=%d\n",SEED);

  neighchk();
  /* Define the probabilities */
  double x,coshx,sinhx;
  x     = eps*j;
  coshx = (exp(x)+exp(-x))/2.0;
  sinhx = (exp(x)-exp(-x))/2.0;
  p1    = exp(-2*x);

  /* update */
  for(i=0;i<ieq;i++) clust();
   

  /* free memory */
  free(ixc); free(iyc); free(itc);
  for(i=0;i<NNBR;i++) free(neigh[i]);
  return 0;
 }

/* Random no generator using Mersanne Twister */
double ran(){
  return genrand_real2();
}

/* returns the checkerboard pointer for site (x,y,z) */
int convert(int x,int y,int t){
  int n,parity;
  parity = (x+y+t)%2;
  n = t*LY*LX + y*LX + x;
  n = n/2;
  if(parity == 1) n = n + VOL2;
  return n;
}

/* sets the neighbors */
void neighbours(){
 extern int convert(int,int,int);
 int i,ieven,iodd,parity;
 int ixp1,ixm1,iyp1,iym1,itp1,itm1;
 int ix,iy,it;
 /* serial to checkerboard list */ 
 ieven=iodd=0;
 for(it=0;it<LT;it++){
 itp1=(it+1)%LT;
 itm1=(it-1+LT)%LT;
 for(iy=0;iy<LY;iy++){
 iyp1=(iy+1)%LY;
 iym1=(iy-1+LY)%LY;
 for(ix=0;ix<LX;ix++){
 ixp1=(ix+1)%LX;
 ixm1=(ix-1+LX)%LX;

 parity = (ix+iy+it)%2;
 /* create checkerboard neighbors */
 /*  7   2   6
  *
  *  3   x   1
  *
  *  8   4   9
  * neigh is a double indexed array. The second index tracks the
  * dual site x, while the first index tracks the 6-spin interaction.
  * 1,2,3,4 show the spins on time-slice t as seen above.
  * 0 denotes the spin at x in time t-1
  * 5 denotes the spin at x in time t+1
  */
 if(parity==0){
  ixc[ieven]=ix;
  iyc[ieven]=iy;
  itc[ieven]=it;
  neigh[0][ieven]=convert(ix,iy,itm1);
  neigh[1][ieven]=convert(ixp1,iy,it);
  neigh[2][ieven]=convert(ix,iyp1,it);
  neigh[3][ieven]=convert(ixm1,iy,it);
  neigh[4][ieven]=convert(ix,iym1,it);
  neigh[5][ieven]=convert(ix,iy,itp1);
  neigh[6][ieven]=convert(ixp1,iyp1,it);
  neigh[7][ieven]=convert(ixm1,iyp1,it);
  neigh[8][ieven]=convert(ixm1,iym1,it);
  neigh[9][ieven]=convert(ixp1,iym1,it);
  ieven++;
 }
 else{
  ixc[VOL2+iodd]=ix;
  iyc[VOL2+iodd]=iy;
  itc[VOL2+iodd]=it;
  neigh[0][VOL2+iodd]=convert(ix,iy,itm1);
  neigh[1][VOL2+iodd]=convert(ixp1,iy,it);
  neigh[2][VOL2+iodd]=convert(ix,iyp1,it);
  neigh[3][VOL2+iodd]=convert(ixm1,iy,it);
  neigh[4][VOL2+iodd]=convert(ix,iym1,it);
  neigh[5][VOL2+iodd]=convert(ix,iy,itp1);
  neigh[6][VOL2+iodd]=convert(ixp1,iyp1,it);
  neigh[7][VOL2+iodd]=convert(ixm1,iyp1,it);
  neigh[8][VOL2+iodd]=convert(ixm1,iym1,it);
  neigh[9][VOL2+iodd]=convert(ixp1,iym1,it);
  iodd++;
 }
 }}}
 if((ieven!=VOL2)||(iodd!=VOL2)) printf("ERROR\n");
 }

void neighchk(){
 int p;
 for(p=0;p<VOL;p++){
  printf("% d %d %d %d %d %d %d %d %d %d\n",p,ixc[p],iyc[p],itc[p],neigh[0][p],
    neigh[1][p],neigh[2][p],neigh[3][p],neigh[4][p],neigh[5][p]);
 }
 }

/* checks the occurance of forbidden configurations */
void chkconfig(){
  int p;
  int flag1,flag2;
  int s0,s1,s2,s3,s4,s5;
  for(p=0;p<VOL;p++){
    /* at site p, check if the spins 0 and 5 are the same */
    s0 = ising[neigh[0][p]];
    s1 = ising[neigh[5][p]];
    if(s0 == s1) flag1 = 0; else flag1 = 1;
    flag2=0;
    if(flag1==1){
     if((ising[neigh[2][p]]==ising[neigh[3][p]])&&
        (ising[neigh[4][p]]==ising[neigh[1][p]])&&
        (ising[neigh[1][p]]!=ising[neigh[2][p]])) flag2=1;
     if((flag1==1)&&(flag2==0)) printf("Forbidden config encountered.\n");
    }
  }
}

/* cluster update */
void clust(){
 int i,p;
 int cflag[VOL4];
 /* note that the variables ising[VOL2+i] with i in [0,VOL2-1] are the */
 /* ones that carry the flag for the reference configurations          */
 for(i=VOL2;i<VOL;i++){
  if((ising[neigh[2][i]]==ising[neigh[3][i]])&&
     (ising[neigh[4][i]]==ising[neigh[1][i]])&&
     (ising[neigh[1][i]]!=ising[neigh[2][i]])) ising[i]=1;
  else ising[i]=0;
 }
 /* mark spins on even time slices for growing clusters */
 for(p=0;p<VOL4;p++) cflag[p]=1;
 /* grow clusters on the even time slices */
 nclus = 0; 
 for(i=VOL2;i<VOL;i++){
   /* first check if the tracking is done on odd slice */
   if(itc[i]%2==0) continue;
   /* then check if ref config exists */
   if(ising[i]==0) {
     /* if it doesnt, then bind the spins at t+1 and t-1 */
     if(cflag[neigh[0][i]]==1) list
   }
 }
 
}
