/* Cluster Algorithm for the U(1) quantum link model */
/* using the representation of the dual height model */


#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>
#include "ranlxd.h"

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
int nclusevn,nclusodd,nclus;
double nclusevsq,nclusodsq;
int refA,refB,mA,mB;

/* pointer variables for neighbours */
#define NNBR 10
int *neigh[NNBR];
/* checkerboard labeling for sites */
int *ixc,*iyc,*itc;
/* Field variables */
int *ising;
int *list;

int main(){

   int ieq,imeas,iskp;
   int i,readval;    
   extern void neighbours();
   extern void neighchk();
   extern double ran();
   extern void clusteven();
   extern void clustodd();
   extern void chkconf();
   extern void initconf();

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
   readval = fscanf(fptr,"%s %d\n",st,&SEED);
   if(readval == -1) printf("Error\n");
   readval = fscanf(fptr,"%s %d\n",st,&ieq);
   if(readval == -1) printf("Error\n");
   readval = fscanf(fptr,"%s %d\n",st,&imeas);
   if(readval == -1) printf("Error\n");
   readval = fscanf(fptr,"%s %lf\n",st,&beta);
   if(readval == -1) printf("Error\n");
   readval = fscanf(fptr,"%s %lf\n",st,&j);
   if(readval == -1) printf("Error\n");
   readval = fscanf(fptr,"%s %lf\n",st,&lam);
   if(readval == -1) printf("Error\n");
   fclose(fptr);
  
  printf("Multi-Cluster Algorithm for the U(1) quantum link model\n");
  printf("Nx=%d, Ny=%d, Nt=%d\n",LX,LY,LT);
  printf("beta=%2.4f; J=%2.3f; lam=%2.3f\n",beta,j,lam);
  printf("Starting seed=%d\n",SEED);

   LT  = 2*LT; /* the dof are spread over twice actual length */
   VOL = LX*LY*LT;
   VOL2= VOL/2;
   VOL4= VOL/4;

   /* allocate memory */
   ixc   = (int *)malloc(VOL*sizeof(int));
   iyc   = (int *)malloc(VOL*sizeof(int));
   itc   = (int *)malloc(VOL*sizeof(int));
   ising = (int *)malloc(VOL*sizeof(int));
   list  = (int *)malloc(VOL*sizeof(int));
   for(i=0;i<NNBR;i++)
     neigh[i] = (int *)malloc(VOL*sizeof(int));

   /* Set parameters */
   eps=1.0*beta/((double)LT);
   iskp=1;

  /* Initialize ranlux */
  rlxd_init(1,SEED);

  /* initialize neighbours */
  neighbours();
  //neighchk();


  /* Define the probabilities */
  double x,coshx,sinhx;
  x     = eps*j;
  coshx = (exp(x)+exp(-x))/2.0;
  sinhx = (exp(x)-exp(-x))/2.0;
//  p1    = exp(-2*x);
  p1    = exp(-x)/coshx;
  p2    = 1.0 - exp(eps*lam)/coshx;
  printf("Prob p1: %f;  Prob p2: %f\n",p1,p2);

  initconf();
  chkconf();

  /* update */
  for(i=0;i<ieq;i++){
     nclusevn = 0; nclusevsq=0; mA=0;
     clusteven();
     nclusodd = 0; nclusodsq=0; mB=0;
     clustodd();
     chkconf();
  }
  
  /* measure */ 
  fptr=fopen("multi.dat","w");
  for(i=0;i<imeas;i++){
   nclusevn = 0; nclusevsq=0; mA=0;
   clusteven();
   nclusevsq = nclusevsq/VOL4;
   nclusodd = 0; nclusodsq=0; mB=0;
   clustodd();
   nclusodsq = nclusodsq/VOL4;
   nclus = nclusevn + nclusodd;
   chkconf();
   //fprintf(fptr,"%d %d %d %lf %lf %d %d\n",nclusevn,nclusodd,nclus,nclusevsq,nclusodsq,mA,mB);
   fprintf(fptr,"%e %e\n",(double)mA,(double)mB);
  }
  fclose(fptr);

  /* free memory */
  free(ixc); free(iyc); free(itc);
  for(i=0;i<NNBR;i++) free(neigh[i]);
  return 0;
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

void initconf(){
  int p;
  /* start with a uniform configuration */
  for(p=0;p<VOL2;p++) ising[p]=1;
  /* since none of the spins are now in a ref config *
   * assign 0 to all flag variables */
  for(p=VOL2;p<VOL;p++) ising[p]=0;
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
  printf("% d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",p,ixc[p],iyc[p],itc[p],neigh[0][p],
    neigh[1][p],neigh[2][p],neigh[3][p],neigh[4][p],neigh[5][p],neigh[6][p],neigh[7][p],neigh[8][p],neigh[9][p],
    neigh[5][neigh[5][p]]);
 }
 }

/* checks the occurance of forbidden configurations */
void chkconf(){
  int p,i,m;
  int flag1;
  int s0,s1,s2,s3,s4,s5;
  for(p=VOL2;p<VOL;p++){
    if(ising[p]<0) printf("Wrong cluster has been grown.\n");
    /* at site p, check if the spins 0 and 5 are the same */
    s0 = ising[neigh[0][p]];
    s5 = ising[neigh[5][p]];
    s1 = ising[neigh[1][p]];
    s2 = ising[neigh[2][p]];
    s3 = ising[neigh[3][p]];
    s4 = ising[neigh[4][p]];
    flag1 = 0;
    if((s0 != -1) && (s0 != 1)) printf("Error at %d\n",s0);
    if((s1 != -1) && (s1 != 1)) printf("Error at %d\n",s1);
    if((s2 != -1) && (s2 != 1)) printf("Error at %d\n",s2);
    if((s3 != -1) && (s3 != 1)) printf("Error at %d\n",s3);
    if((s4 != -1) && (s4 != 1)) printf("Error at %d\n",s4);
    if((s5 != -1) && (s5 != 1)) printf("Error at %d\n",s5);
    if(s0 != s5) {
     if((s2 == s3) && (s4 == s1) && (s1 != s2)) flag1=1;
     if(flag1==0) printf("Forbidden config encountered. Flag at = %d. %d %d %d %d\n",p,s2,s3,s4,s1);
    }
  }
}

/* cluster update for even time-slices */
void clusteven(){
 int i,p,d,m,imf6,imf7;
 int im,imf0,imf1,imf2,imf3,imf4,imf5,fwd,bwd;
 int r1,r2,r3,l1,l2,l3,u1,u2,u3,d1,d2,d3;
 int A,B,C,D,aux;
 int cflag[VOL];
 int bondflag;
 int sx,sy;
 double ran[1];
 /* note that the variables ising[VOL2+i] with i in [0,VOL2-1] are the */
 /* ones that carry the flag for the reference configurations          */
 /* reference configuration flags are 2 */
 /* if there is no ref config, then ising[i] is set to zero */
 for(i=VOL2;i<VOL;i++){
  if((ising[neigh[2][i]]==ising[neigh[3][i]])&&
     (ising[neigh[4][i]]==ising[neigh[1][i]])&&
     (ising[neigh[1][i]]!=ising[neigh[2][i]])) ising[i]=2;
  else ising[i]=0;
 }
 /* mark spins on even time slices for growing clusters */
 /* spins on odd-time slices are marked as 0; so that it will never join to a cluster */
 for(p=0;p<VOL2;p++) {
  if((itc[p]%2)==0) cflag[p]=1;
  if((itc[p]%2)==1) cflag[p]=0;
 }
 /* serve as flags for checking interactions while cluster building */
 for(p=VOL2;p<VOL;p++) cflag[p]=1;

 /* grow clusters on the even time slices */
 for(p=VOL2;p<VOL;p++){
   /* first check if the tracking is done on odd slice */
   if(itc[p]%2==0) continue;

   /* check for any inconsistency: if any of the flag variables is -1 */
   /* then the cluster building is flawed */
   if((ising[p]==-1) || (ising[p]==1)) printf("Wrong cluster grown\n");

   /* skip if the site already belongs to a cluster */
   if(cflag[neigh[0][p]]==0) continue;
   /* otherwise, start building a new cluster */
   m=0; i=0; list[i]=neigh[0][p]; cflag[neigh[0][p]]=0; nclusevn++;
   do{
    im=list[m]; /* m is the new or the starting site */
    /* first check the spin on time-slice t+1 wants to bind */
    /* remember that you are on even time-slice t-1*/
    imf0=neigh[5][im];
    imf1=neigh[5][imf0];
    if(cflag[imf0]==1){
      bondflag=0;
//    if(ising[imf0]==2) { ranlxd(ran,1); if(ran[0] < p1) bondflag=1; }
    if((ising[imf0]==2)&&(ising[im]==ising[imf1])) { ranlxd(ran,1); if(ran[0] < p1) bondflag=1; }
    else if(ising[imf0]==0) bondflag=1; 
    if((bondflag)&&(cflag[imf1]==1)){
     i++; list[i]=imf1; /* increase list*/
     cflag[imf1]=0;  /* unmark spins belonging to cluster */ }
    cflag[imf0]=0;   /* unmark the interaction */
   }
   /* ============================================== */
   /* Also check if the spin in time-slice t-3 wants to bind */
   imf6=neigh[0][im];
   imf7=neigh[0][imf6];
   if(cflag[imf6]==1){
    bondflag=0;
//  if(ising[imf6]==2){ ranlxd(ran,1); if(ran[0] < p1) bondflag=1; }
    if((ising[imf6]==2)&&(ising[im]==ising[imf7])) { ranlxd(ran,1); if(ran[0] < p1) bondflag=1; }
    else if(ising[imf6]==0) bondflag=1; 
    if((bondflag)&&(cflag[imf7]==1)){
     i++; list[i]=imf7; /* increase list*/
     cflag[imf7]=0;  /* unmark spins belonging to cluster */ }
    cflag[imf6]=0;   /* unmark the interaction */
   }
   /* ============================================== */
   /* Next check if other spins in the time-slice t-1 want to bind */
   /* To see if the spins to the right-side of neigh[0] want to bind */
   imf2=neigh[1][im];
   if(cflag[imf2]==1){
     fwd=neigh[0][imf2]; bwd=neigh[5][imf2];
     bondflag=0;
     if(ising[fwd]!=ising[bwd]) bondflag=1;
     else { if(ising[imf2]==2) { ranlxd(ran,1); if(ran[0] < p2) bondflag=1;}}
     r1=neigh[1][imf2];  /* x   r2     x */
     r2=neigh[2][imf2];  /* im imf2   r1 */
     r3=neigh[4][imf2];  /* x   r3     x */
     if((bondflag)&&(cflag[r1]==1)){
       i++; list[i]=r1; /* increase list*/
       cflag[r1]=0;  /* unmark spins belonging to cluster */ }
     if((bondflag)&&(cflag[r2]==1)){
       i++; list[i]=r2; /* increase list*/
       cflag[r2]=0;  /* unmark spins belonging to cluster */ }
     if((bondflag)&&(cflag[r3]==1)){
       i++; list[i]=r3; /* increase list*/
       cflag[r3]=0;  /* unmark spins belonging to cluster */ }
     cflag[imf2]=0;  /* unmark the interaction */
   }
   /* ============================================== */
   /* To see if the spins to the left-side of neigh[0] want to bind */
   imf3=neigh[3][im];
   if(cflag[imf3]==1){
     fwd=neigh[0][imf3]; bwd=neigh[5][imf3];
     bondflag=0;
     if(ising[fwd]!=ising[bwd]) bondflag=1;
     else { if(ising[imf3]==2) { ranlxd(ran,1); if(ran[0] < p2) bondflag=1;}}
     l1=neigh[2][imf3];  /* x   l1    x */
     l2=neigh[3][imf3];  /* l2 imf3  im */
     l3=neigh[4][imf3];  /* x   l3    x */
     if((bondflag)&&(cflag[l1]==1)){
       i++; list[i]=l1; /* increase list*/
       cflag[l1]=0;  /* unmark spins belonging to cluster */ }
     if((bondflag)&&(cflag[l2]==1)){
       i++; list[i]=l2; /* increase list*/
       cflag[l2]=0;  /* unmark spins belonging to cluster */ }
     if((bondflag)&&(cflag[l3]==1)){
       i++; list[i]=l3; /* increase list*/
       cflag[l3]=0;  /* unmark spins belonging to cluster */ }
     cflag[imf3]=0;  /* unmark the interaction */
   }
   /* ============================================== */
   /* To see if the spins to the top of neigh[0] want to bind */
   imf4=neigh[2][im];
   if(cflag[imf4]==1){
     fwd=neigh[0][imf4]; bwd=neigh[5][imf4];
     bondflag=0;
     if(ising[fwd]!=ising[bwd]) bondflag=1;
     else { if(ising[imf4]==2) { ranlxd(ran,1); if(ran[0] < p2) bondflag=1;}}
     u1=neigh[1][imf4];  /* x   u2    x */
     u2=neigh[2][imf4];  /* u3 imf4  u1 */
     u3=neigh[3][imf4];  /* x   im    x */
     if((bondflag)&&(cflag[u1]==1)){
       i++; list[i]=u1; /* increase list*/
       cflag[u1]=0;  /* unmark spins belonging to cluster */ }
     if((bondflag)&&(cflag[u2]==1)){
       i++; list[i]=u2; /* increase list*/
       cflag[u2]=0;  /* unmark spins belonging to cluster */ }
     if((bondflag)&&(cflag[u3]==1)){
       i++; list[i]=u3; /* increase list*/
       cflag[u3]=0;  /* unmark spins belonging to cluster */ }
     cflag[imf4]=0;  /* unmark the interaction */
   }
   /* ============================================== */
   /* To see if the spins to the down of neigh[0] want to bind */
   imf5=neigh[4][im];
   if(cflag[imf5]==1){
     fwd=neigh[0][imf5]; bwd=neigh[5][imf5];
     bondflag=0;
     if(ising[fwd]!=ising[bwd]) bondflag=1;
     else { if(ising[imf5]==2) { ranlxd(ran,1); if(ran[0] < p2) bondflag=1;}}
     d1=neigh[1][imf5];  /* x   im    x */
     d2=neigh[3][imf5];  /* d2 imf5  d1 */
     d3=neigh[4][imf5];  /* x   d3    x */
     if((bondflag)&&(cflag[d1]==1)){
       i++; list[i]=d1; /* increase list*/
       cflag[d1]=0;  /* unmark spins belonging to cluster */ }
     if((bondflag)&&(cflag[d2]==1)){
       i++; list[i]=d2; /* increase list*/
       cflag[d2]=0;  /* unmark spins belonging to cluster */ }
     if((bondflag)&&(cflag[d3]==1)){
       i++; list[i]=d3; /* increase list*/
       cflag[d3]=0;  /* unmark spins belonging to cluster */ }
     cflag[imf5]=0;  /* unmark the interaction */
   }
   /* Implement the Gauss law if this is the last time slice */
   if(itc[im]==(LT-2)){
       /* This involves checking of 4-spins in the same time-slice */
       /*     if they want to get added in the cluster             */
       /*       o-----x-----o      1-----A-----2                   */
       /*       |     |     |      |     |     |                   */
       /*       x-----o-----x      B-----0-----C                   */
       /*       |     |     |      |     |     |                   */
       /*       o-----x-----o      3-----D-----4                   */
       /* The central o decides to connect the diagonal o depending*/
       /* the configuration at the crosses. For example, 0 will    */
       /* decide to connect with 1 depending on whether AB are in  */
       /* the reference configuration                              */
       A=neigh[5][neigh[2][im]]; B=neigh[5][neigh[3][im]];
       C=neigh[5][neigh[1][im]]; D=neigh[5][neigh[4][im]];
       if(A > VOL2) printf("Error. Gauss Law counter\n");
       if(B > VOL2) printf("Error. Gauss Law counter\n");
       if(C > VOL2) printf("Error. Gauss Law counter\n");
       if(D > VOL2) printf("Error. Gauss Law counter\n");
       /* decide whether to join spin 1 to the cluster */
       if(ising[A] != ising[B]) {
         aux=neigh[7][im];
         if(cflag[aux]==1) { 
         i++; list[i]=aux; /* increase list */
         cflag[aux]=0; /* unmark spin belonging to cluster */ }
       }
       /* decide whether to join spin 2 to the cluster */
       if(ising[A] == ising[C]) {
         aux=neigh[6][im];
         if(cflag[aux]==1) { 
         i++; list[i]=aux; /* increase list */
         cflag[aux]=0; /* unmark spin belonging to cluster */ }
       }
       /* decide whether to join spin 3 to the cluster */
       if(ising[B] == ising[D]) {
         aux=neigh[8][im];
         if(cflag[aux]==1) { 
         i++; list[i]=aux; /* increase list */
         cflag[aux]=0; /* unmark spin belonging to cluster */ }
       }
       /* decide whether to join spin 4 to the cluster */
       if(ising[C] != ising[D]) {
         aux=neigh[9][im];
         if(cflag[aux]==1) { 
         i++; list[i]=aux; /* increase list */
         cflag[aux]=0; /* unmark spin belonging to cluster */ }
       }
   }
   m++;
   } while(m<=i);
   /* check if the list only contains genuine spins */
   for(d=0;d<=i;d++) if(list[d]>=VOL2) printf("Cluster grown in Flag.\n");
   /* decide orientation wrt the reference config */
   sx=ixc[list[0]]; sy=iyc[list[0]];
   if((sx-sy)%4==0) refA=1;
   else refA=-1;
   if(ising[list[0]]==refA) refA=1;
   else refA=-1;
   mA = mA + refA*(i+1);

   /* size */
   nclusevsq += (i+1)*(i+1);
   /* flip the cluster with a 50% probability */
   ranlxd(ran,1);
   if(ran[0]<0.5){
     for(d=0;d<=i;d++) ising[list[d]] = -ising[list[d]];
   }
 }
}

/* cluster update for odd  time-slices */
void clustodd(){
 int i,p,d,m,imf6,imf7;
 int im,imf0,imf1,imf2,imf3,imf4,imf5,fwd,bwd;
 int r1,r2,r3,l1,l2,l3,u1,u2,u3,d1,d2,d3;
 int A,B,C,D,aux;
 int cflag[VOL];
 int bondflag;
 int sx,sy,val;
 double ran[1];
 /* note that the variables ising[VOL2+i] with i in [0,VOL2-1] are the */
 /* ones that carry the flag for the reference configurations          */
 /* reference configuration flags are 2 */
 /* if there is no ref config, then ising[i] is set to zero */
 for(i=VOL2;i<VOL;i++){
  if((ising[neigh[2][i]]==ising[neigh[3][i]])&&
     (ising[neigh[4][i]]==ising[neigh[1][i]])&&
     (ising[neigh[1][i]]!=ising[neigh[2][i]])) ising[i]=2; 
  else ising[i]=0;
 } 

 /* mark spins on odd time slices for growing clusters */
 /* spins on even-time slices are marked as 0; so that it will never join to a cluster */
 for(p=0;p<VOL2;p++) {
  if((itc[p]%2)==0) cflag[p]=0;
  if((itc[p]%2)==1) cflag[p]=1;
 }
 /* serve as flags for joining spins on the same time-slice */
 for(p=VOL2;p<VOL;p++) cflag[p]=1;

 /* grow clusters on the odd time slices */
 for(p=VOL2;p<VOL;p++){
   /* first check if the tracking is done on even slice */
   if(itc[p]%2==1) continue;

   /* check for any inconsistency: if any of the flag variables is -1 */
   /* then the cluster building is flawed */
   if((ising[p]==-1) ||(ising[p]==1)) printf("Wrong cluster grown\n");

   /* skip if the site already belongs to a cluster */
   if(cflag[neigh[0][p]]==0) continue;
   /* otherwise, start building a new cluster */
   m=0; i=0; list[i]=neigh[0][p]; cflag[neigh[0][p]]=0; nclusodd++;
   do{
    im=list[m]; /* m is the new or the starting site */
    /* first check the spin on time-slice t+1 wants to bind */
    /* remember that you are on odd time-slice t-1*/
    imf0=neigh[5][im];
    imf1=neigh[5][imf0];
    if(cflag[imf0]==1) {
      bondflag=0;
//   if(ising[imf0]==2) {ranlxd(ran,1); if(ran[0] < p1) bondflag=1; }
    if((ising[imf0]==2)&&(ising[im]==ising[imf1])) { ranlxd(ran,1); if(ran[0] < p1) bondflag=1; }
      else if(ising[imf0]==0) bondflag=1; 
      if((bondflag)&&(cflag[imf1]==1)){
      i++; list[i]=imf1; /* increase list*/
      cflag[imf1]=0;  /* unmark spins belonging to cluster */ }
    cflag[imf0]=0;   /* unmark the interaction */
   }
   /* ============================================== */
   /* Also check if the spin in time-slice t-3 wants to bind */
   imf6=neigh[0][im];
   imf7=neigh[0][imf6];
   if(cflag[imf6]==1){
     bondflag=0;
//  if(ising[imf6]==2){ ranlxd(ran,1); if(ran[0] < p1) bondflag=1; }
    if((ising[imf6]==2)&&(ising[im]==ising[imf7])) { ranlxd(ran,1); if(ran[0] < p1) bondflag=1; }
     else if(ising[imf6]==0) bondflag=1; 
     if((bondflag)&&(cflag[imf7]==1)){
      i++; list[i]=imf7; /* increase list*/
      cflag[imf7]=0;  /* unmark spins belonging to cluster */ }
     cflag[imf6]=0;   /* unmark the interaction */
   }
   /* ============================================== */
   /* Next check if other spins in the time-slice t-1 want to bind */
   /* To see if the spins to the right-side of neigh[0] want to bind */
   imf2=neigh[1][im];
   if(cflag[imf2]==1){
     fwd=neigh[0][imf2]; bwd=neigh[5][imf2];
     bondflag=0;
     if(ising[fwd]!=ising[bwd]) bondflag=1;
     else { if(ising[imf2]==2) { ranlxd(ran,1); if(ran[0] < p2) bondflag=1;}}
     r1=neigh[1][imf2];  /* x   r2     x */
     r2=neigh[2][imf2];  /* im imf2   r1 */
     r3=neigh[4][imf2];  /* x   r3     x */
     if((bondflag)&&(cflag[r1]==1)){
       i++; list[i]=r1; /* increase list*/
       cflag[r1]=0;  /* unmark spins belonging to cluster */ }
     if((bondflag)&&(cflag[r2]==1)){
       i++; list[i]=r2; /* increase list*/
       cflag[r2]=0;  /* unmark spins belonging to cluster */ }
     if((bondflag)&&(cflag[r3]==1)){
       i++; list[i]=r3; /* increase list*/
       cflag[r3]=0;  /* unmark spins belonging to cluster */ }
     cflag[imf2]=0;  /* unmark the interaction */
   }
   /* ============================================== */
   /* To see if the spins to the left-side of neigh[0] want to bind */
   imf3=neigh[3][im];
   if(cflag[imf3]==1){
     fwd=neigh[0][imf3]; bwd=neigh[5][imf3];
     bondflag=0;
     if(ising[fwd]!=ising[bwd]) bondflag=1;
     else { if(ising[imf3]==2) { ranlxd(ran,1); if(ran[0] < p2) bondflag=1;}}
     l1=neigh[2][imf3];  /* x   l1    x */
     l2=neigh[3][imf3];  /* l2 imf3  im */
     l3=neigh[4][imf3];  /* x   l3    x */
     if((bondflag)&&(cflag[l1]==1)){
       i++; list[i]=l1; /* increase list*/
       cflag[l1]=0;  /* unmark spins belonging to cluster */ }
     if((bondflag)&&(cflag[l2]==1)){
       i++; list[i]=l2; /* increase list*/
       cflag[l2]=0;  /* unmark spins belonging to cluster */ }
     if((bondflag)&&(cflag[l3]==1)){
       i++; list[i]=l3; /* increase list*/
       cflag[l3]=0;  /* unmark spins belonging to cluster */ }
     cflag[imf3]=0;  /* unmark the interaction */
   }
   /* ============================================== */
   /* To see if the spins to the top of neigh[0] want to bind */
   imf4=neigh[2][im];
   if(cflag[imf4]==1){
     fwd=neigh[0][imf4]; bwd=neigh[5][imf4];
     bondflag=0;
     if(ising[fwd]!=ising[bwd]) bondflag=1;
     else { if(ising[imf4]==2) { ranlxd(ran,1); if(ran[0] < p2) bondflag=1;}}
     u1=neigh[1][imf4];  /* x   u2    x */
     u2=neigh[2][imf4];  /* u3 imf4  u1 */
     u3=neigh[3][imf4];  /* x   im    x */
     if((bondflag)&&(cflag[u1]==1)){
       i++; list[i]=u1; /* increase list*/
       cflag[u1]=0;  /* unmark spins belonging to cluster */ }
     if((bondflag)&&(cflag[u2]==1)){
       i++; list[i]=u2; /* increase list*/
       cflag[u2]=0;  /* unmark spins belonging to cluster */ }
     if((bondflag)&&(cflag[u3]==1)){
       i++; list[i]=u3; /* increase list*/
       cflag[u3]=0;  /* unmark spins belonging to cluster */ }
     cflag[imf4]=0;  /* unmark the interaction */
   }
   /* ============================================== */
   /* To see if the spins to the down of neigh[0] want to bind */
   imf5=neigh[4][im];
   if(cflag[imf5]==1){
     fwd=neigh[0][imf5]; bwd=neigh[5][imf5];
     bondflag=0;
     if(ising[fwd]!=ising[bwd]) bondflag=1;
     else { if(ising[imf5]==2) { ranlxd(ran,1); if(ran[0] < p2) bondflag=1;}}
     d1=neigh[1][imf5];  /* x   im    x */
     d2=neigh[3][imf5];  /* d2 imf5  d1 */
     d3=neigh[4][imf5];  /* x   d3    x */
     if((bondflag)&&(cflag[d1]==1)){
       i++; list[i]=d1; /* increase list*/
       cflag[d1]=0;  /* unmark spins belonging to cluster */ }
     if((bondflag)&&(cflag[d2]==1)){
       i++; list[i]=d2; /* increase list*/
       cflag[d2]=0;  /* unmark spins belonging to cluster */ }
     if((bondflag)&&(cflag[d3]==1)){
       i++; list[i]=d3; /* increase list*/
       cflag[d3]=0;  /* unmark spins belonging to cluster */ }
     cflag[imf5]=0;  /* unmark the interaction */
   }
   /* ============================================== */
   /* Implement the Gauss law if this is the last time slice */
   if(itc[im]==(LT-1)){
       /* This involves checking of 4-spins in the same time-slice */
       /*     if they want to get added in the cluster             */
       /*       o-----x-----o      1-----A-----2                   */
       /*       |     |     |      |     |     |                   */
       /*       x-----o-----x      B-----0-----C                   */
       /*       |     |     |      |     |     |                   */
       /*       o-----x-----o      3-----D-----4                   */
       /* The central o decides to connect the diagonal o depending*/
       /* the configuration at the crosses. For example, 0 will    */
       /* decide to connect with 1 depending on whether AB are in  */
       /* the reference configuration                              */
       A=neigh[0][neigh[2][im]]; B=neigh[0][neigh[3][im]];
       C=neigh[0][neigh[1][im]]; D=neigh[0][neigh[4][im]];
       if(A > VOL2) printf("Error. Gauss Law counter\n");
       if(B > VOL2) printf("Error. Gauss Law counter\n");
       if(C > VOL2) printf("Error. Gauss Law counter\n");
       if(D > VOL2) printf("Error. Gauss Law counter\n");
       /* decide whether to join spin 1 to the cluster */
       if(ising[A] != ising[B]) {
         aux=neigh[7][im];
         if(cflag[aux]==1) { 
         i++; list[i]=aux; /* increase list */
         cflag[aux]=0; /* unmark spin belonging to cluster */ }
       }
       /* decide whether to join spin 2 to the cluster */
       if(ising[A] == ising[C]) {
         aux=neigh[6][im];
         if(cflag[aux]==1) { 
         i++; list[i]=aux; /* increase list */
         cflag[aux]=0; /* unmark spin belonging to cluster */ }
       }
       /* decide whether to join spin 3 to the cluster */
       if(ising[B] == ising[D]) {
         aux=neigh[8][im];
         if(cflag[aux]==1) { 
         i++; list[i]=aux; /* increase list */
         cflag[aux]=0; /* unmark spin belonging to cluster */ }
       }
       /* decide whether to join spin 4 to the cluster */
       if(ising[C] != ising[D]) {
         aux=neigh[9][im];
         if(cflag[aux]==1) { 
         i++; list[i]=aux; /* increase list */
         cflag[aux]=0; /* unmark spin belonging to cluster */ }
       }
    }
    m++;
   } while(m<=i);
   /* check if the list only contains genuine spins */
   for(d=0;d<=i;d++) if(list[d]>=VOL2) printf("Cluster grown in Flag.\n");

   /* decide orientation wrt to the reference config */
   sx=ixc[list[0]]; sy=iyc[list[0]];
   val=(sx-sy)%4;
   if((val==-1)||(val==3)) refB=-1;
   else refB=1;
   if(ising[list[0]]==refB) refB=1;
   else refB=-1;
   mB = mB + refB*(i+1);

   /* size */
   nclusodsq += (i+1)*(i+1);

   /* flip the cluster with a 50% probability */
   ranlxd(ran,1);
   if(ran[0]<0.5){
     for(d=0;d<=i;d++) ising[list[d]] = -ising[list[d]];
   }
 }
}
