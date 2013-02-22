/* Cluster Algorithm for the U(1) quantum link model */
/* using the representation of the dual height model */
/* codes in the presence of static charges and anti-charges (+2Q,-2Q) */
/* This version creates the initial configurations and stores them in */
/* initconf.dat */


#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>
#include "ranlxd.h"

/* Parameters                      *
 * (LX,LY,LT) : lattice dimension  *
 * J          : coupling constant  *
 * lam        : potential energy   *
 * beta       : inverse temperature*
 * eps        : Trotter dscrt      *
 * ieq        : equilibrium iter   *
 * imeas      : # of measurements  *
 * iskp       : # of iter/meas     */

#define DIM 3
int LX,LY,LT2,LT,VOL,VOL2,VOL4;
int SVOL,SIZE;
double p1,p2;
double J,lam,beta,eps;
double tanhx,cothx;
int SEED;
int nclusevn,nclusodd,nclus;
int thermflag,diagflag;
/* flags to monitor the flux */
int flagxx,flagyy;
/* counters to count the superselection sectors whether the 
 * pass through the boundary or the bulk  */
int flxcnt1,flxcnt2,flxcnt3,flxcnt4;
double nclusevsq,nclusodsq;
double inten,intpe;
int refA,refB,mA,mB;
int cp,cm,cpos; /* records the position of the charges */
int minMA,minMB,maxMA,maxMB;
int *MA,*MB,**refC;
double **pMAB;

/* pointer variables for neighbours */
#define NNBR 10
int *neigh[NNBR];
/* checkerboard labeling for sites */
int *ixc,*iyc,*itc;
/* Field variables */
int *ising;
int *list;
int *chptr;
/* for the flux configuration */
int **fx,**fy,*next[2*DIM+1];
int SPV;
int flt1,flt2;
double flx,fly;
/* avfl(x,y)(1,2) tracks flux in bulk */
/* avfl(x,y)(3,4) tracks flux in boundary */
double *avflx1,*avfly1,*avflx2,*avfly2;
double *avflx3,*avfly3,*avflx4,*avfly4;
/* check variables for Gauss Law on the expectation values */
float glchk1,glchk2;
/* measure of "a energy density" */
double eden; 


int main(){

   int ieq,imeas,iskp;
   int i,j,q,readval,parity,val;   
   int ix,iy;
   double ex,ey;
   double norm; 
   extern void neighbours();
   extern void neighchk();
   extern double ran();
   extern void clusteven();
   extern void clustodd();
   extern void chkconf();
   extern void measureMAB();
   extern void constflux();
   extern void energy();
   extern void initconf();
   extern double **allocatedouble2d(int,int);
   extern void deallocatedouble2d(double**,int,int);
   extern int **allocate2d(int,int);
   extern void deallocate2d(int**,int,int);
   extern void creatediagconf(int,int); 
   extern void createaxisconf(int,int);

   FILE *fptr;
   char st[20];
   /* read file */
   fptr=fopen("QUEUECH","r");
   if(fptr == NULL) {printf("QUEUE error.\n"); exit(1);}
   readval = fscanf(fptr,"%s %d\n",st,&LX);
   if(readval == -1) printf("Error\n");
   readval = fscanf(fptr,"%s %d\n",st,&LY);
   if(readval == -1) printf("Error\n");
   readval = fscanf(fptr,"%s %d\n",st,&LT2);
   if(readval == -1) printf("Error\n");
   readval = fscanf(fptr,"%s %d\n",st,&SEED);
   if(readval == -1) printf("Error\n");
   readval = fscanf(fptr,"%s %d\n",st,&ieq);
   if(readval == -1) printf("Error\n");
   readval = fscanf(fptr,"%s %d\n",st,&imeas);
   if(readval == -1) printf("Error\n");
   readval = fscanf(fptr,"%s %lf\n",st,&beta);
   if(readval == -1) printf("Error\n");
   readval = fscanf(fptr,"%s %lf\n",st,&J);
   if(readval == -1) printf("Error\n");
   readval = fscanf(fptr,"%s %lf\n",st,&lam);
   if(readval == -1) printf("Error\n");
   readval = fscanf(fptr,"%s %d\n",st,&diagflag);
   if(readval == -1) printf("Error\n");
   readval = fscanf(fptr,"%s %d\n",st,&cp);
   if(readval == -1) printf("Error\n");
   readval = fscanf(fptr,"%s %d\n",st,&cm);
   if(readval == -1) printf("Error\n");
   fclose(fptr);
  
  printf("Multi-Cluster Algorithm for the U(1) quantum link model\n");
  printf("Nx=%d, Ny=%d, Nt=%d\n",LX,LY,LT2);
  printf("beta=%2.4f; J=%2.3f; lam=%2.5f\n",beta,J,lam);
  printf("Starting seed=%d\n",SEED);
  if(diagflag==1) {
  printf("Positive charge placed at=(%d,%d)\n",cp,cp);
  printf("Negative charge placed at=(%d,%d)\n",cm,cm);
  }
  else if(diagflag==0){
  printf("Positive charge placed at=(%d,%d)\n",cp,LY/2-1);
  printf("Negative charge placed at=(%d,%d)\n",cm,LY/2-1);
  }

  /* calculate the position of the +ve charge in lattice coordinates */
  if(diagflag==1) cpos = cp*LX + cp;
  else if(diagflag==0) cpos = (LY/2-1)*LX + cp;

   flagxx=flagyy=0;
   LT  = 2*LT2; /* the dof are spread over twice actual length */
   VOL = LX*LY*LT;
   SVOL= LX*LY/2;
   SPV = 2*SVOL;
   SIZE= 2*SVOL+1;
   VOL2= VOL/2;
   VOL4= VOL/4;
   minMA = minMB = -LX*LY/2;
   maxMA = maxMB =  LX*LY/2;
   
   /* allocate memory */
   ixc   = (int *)malloc(VOL*sizeof(int));
   iyc   = (int *)malloc(VOL*sizeof(int));
   itc   = (int *)malloc(VOL*sizeof(int));
   ising = (int *)malloc(VOL*sizeof(int));
   list  = (int *)malloc(VOL*sizeof(int));
   avflx1= (double *)malloc(SPV*sizeof(double));
   avfly1= (double *)malloc(SPV*sizeof(double));
   avflx2= (double *)malloc(SPV*sizeof(double));
   avfly2= (double *)malloc(SPV*sizeof(double));
   avflx3= (double *)malloc(SPV*sizeof(double));
   avfly3= (double *)malloc(SPV*sizeof(double));
   avflx4= (double *)malloc(SPV*sizeof(double));
   avfly4= (double *)malloc(SPV*sizeof(double));
   chptr = (int *)malloc(VOL2*sizeof(double)); 
   /* MA and MB will be measured for each of the LT/2 timeslices and avg-d */
   MA    = (int *)malloc((LT2)*sizeof(int));
   MB    = (int *)malloc((LT2)*sizeof(int));
   pMAB  = allocatedouble2d(SIZE,SIZE);
   refC  = allocate2d(LX,LY);
   fx    = allocate2d(LT,SPV);
   fy    = allocate2d(LT,SPV);
   for(i=0;i<NNBR;i++) neigh[i] = (int *)malloc(VOL*sizeof(int));
   for(i=0;i<2*DIM+1;i++) next[i] = (int *)malloc(VOL*sizeof(int));

   /* Set parameters */
   eps=1.0*beta/((double)LT2);
   iskp=1;

  /* Initialize ranlux */
  rlxd_init(1,SEED);

  /* initialize neighbours */
  neighbours();

  /* Define the probabilities */
  double x,coshx,sinhx;
  x     = eps*J;
  coshx = (exp(x)+exp(-x))/2.0;
  sinhx = (exp(x)-exp(-x))/2.0;
  tanhx = sinhx/coshx;
  cothx = coshx/sinhx;
  p1    = exp(-x)/coshx;
  p2    = 1.0 - exp(eps*lam)/coshx;
  printf("eps = %f\n",eps);
  printf("Prob p1: %f;  Prob p2: %f\n",p1,p2);

  /* create initial configuration */
  if(diagflag==1)
    creatediagconf(cm,cp);
  else
    createaxisconf(cm,cp);
  /* initialize configration with charges and check how the flux is channeled */
  initconf();
  chkconf();
  constflux();
  
  /* initialize the reference configuration */
  for(i=0;i<LX;i++) for(j=0;j<LY;j++){
       parity=(i+j)%2;
       val=(i-j)%4;
       if(parity==0){ /* for the even-time slices */
          if(val==0) refC[i][j]=1;
          else refC[i][j]=-1;
       }
       else{  /* for the odd-time slices */
          if((val==-1)||(val==3)) refC[i][j]=-1;
          else refC[i][j]=1;
       }
  }

  /* initialize average flux variable */
  for(i=0;i<SPV;i++) 
  avflx1[i]=avfly1[i]=avflx2[i]=avflx2[i]=avflx3[i]=avfly3[i]=avflx4[i]=avfly4[i]=0.0;

  /* update */
  thermflag=1;
  for(i=0;i<ieq;i++){
     nclusevn = 0; nclusevsq=0; mA=0;
     nclusodd = 0; nclusodsq=0; mB=0;
     clusteven();
     clustodd();
     chkconf();
  }
  thermflag=0;
  /* measure */
  flxcnt1=flxcnt2=flxcnt3=flxcnt4=0;
  for(i=0;i<SIZE;i++) for(j=0;j<SIZE;j++) pMAB[i][j]=0.0;
  fptr=fopen("multi.dat","w");
  for(i=0;i<imeas;i++){
     nclusevn = 0; nclusevsq=0; mA=0; 
     nclusodd = 0; nclusodsq=0; mB=0;
     //inten=intpe=0.0; 
     clusteven();
     nclusevsq = nclusevsq/VOL4;
     clustodd();
     nclusodsq = nclusodsq/VOL4;
     nclus = nclusevn + nclusodd;
     measureMAB();
     constflux();
     energy();
     chkconf();
     //fprintf(fptr,"%d %d %d %lf %lf %d %d %lf %lf\n",nclusevn,nclusodd,nclus,nclusevsq,nclusodsq,mA,mB,inten,intpe);
     fprintf(fptr,"%d %d %lf %lf\n",mA,mB,inten,intpe);
  }
  fclose(fptr);

  /* normalize and print histogram */
  norm=0.0;
  for(i=0;i<SIZE;i++) for(j=0;j<SIZE;j++) norm += pMAB[i][j];
  fptr=fopen("magdist.dat","w");
  for(i=0;i<SIZE;i++){
  for(j=0;j<SIZE;j++){
   pMAB[i][j] /= norm;
   fprintf(fptr,"%d %d %le\n",i,j,pMAB[i][j]);}
  fprintf(fptr,"\n");}
  fclose(fptr);

  /* average and normalize the flux profile and print it */
  printf("flux in bulk: %d + %d = %d\n",flxcnt1,flxcnt2,flxcnt1+flxcnt2);
  printf("# meas: %d + %d + %d = %d\n",flxcnt1+flxcnt2,flxcnt3,flxcnt4,
           flxcnt1+flxcnt2+flxcnt3+flxcnt4);
  printf("boundary flux: %d %d\n",flagxx,flagyy);
  for(i=0;i<SPV;i++){
    if((flxcnt1 != 0)&&(flxcnt2 != 0)){
    avflx1[i] /= (2.0*LT2*flxcnt1);
    avfly1[i] /= (2.0*LT2*flxcnt1);
    avflx2[i] /= (2.0*LT2*flxcnt2);
    avfly2[i] /= (2.0*LT2*flxcnt2);}
    if((flxcnt3 != 0)||(flxcnt4 != 0)){
    avflx3[i] /= (2.0*LT2*flxcnt3);
    avfly3[i] /= (2.0*LT2*flxcnt3);
    avflx4[i] /= (2.0*LT2*flxcnt4);
    avfly4[i] /= (2.0*LT2*flxcnt4);}
  }
  fptr=fopen("fprof_blk.dat","w");
  for(iy=0;iy<LY;iy++){
  for(ix=0;ix<LX;ix++){
    i=iy*LX+ix;
    /* also calculate and print the Gauss Law at the vertices */
    glchk1 = avflx1[i] - avflx1[next[DIM-1][i]] + avfly1[i] - avfly1[next[DIM-2][i]];
    glchk2 = avflx2[i] - avflx2[next[DIM-1][i]] + avfly2[i] - avfly2[next[DIM-2][i]];
    eden = avflx1[i]*avflx1[i] + avflx1[next[DIM-1][i]]*avflx1[next[DIM-1][i]] +
          avfly1[i]*avfly1[i] + avfly1[next[DIM-2][i]]*avfly1[next[DIM-2][i]];
    fprintf(fptr,"%d %d %d % .5f % .5f % .2f % .5f % .5f % .2f %.5f\n",i,ix,iy,avflx1[i],avfly1[i],
       glchk1,avflx2[i],avfly2[i],glchk2,eden);
  }
  fprintf(fptr,"\n");
  }
  fclose(fptr);
  fptr=fopen("fprof_bdy.dat","w");
  for(iy=0;iy<LY;iy++){
  for(ix=0;ix<LX;ix++){
    i=iy*LX+ix;
    /* also calculate and print the Gauss Law at the vertices */
    glchk1 = avflx3[i] - avflx3[next[DIM-1][i]] + avfly3[i] - avfly3[next[DIM-2][i]];
    glchk2 = avflx4[i] - avflx4[next[DIM-1][i]] + avfly4[i] - avfly4[next[DIM-2][i]];
    eden = avflx3[i]*avflx3[i] + avflx3[next[DIM-1][i]]*avflx3[next[DIM-1][i]] +
          avfly3[i]*avfly3[i] + avfly3[next[DIM-2][i]]*avfly3[next[DIM-2][i]];
    fprintf(fptr,"%d %d %d % .5f % .5f % .2f % .5f % .5f % .2f %.5f\n",i,ix,iy,avflx3[i],avfly3[i],
       glchk1,avflx4[i],avfly4[i],glchk2,eden);
  }
  fprintf(fptr,"\n");
  }
  fclose(fptr);



 /* check flux carried by average flux profile */
 fptr=fopen("CONFCHK.dat","a");
 fprintf(fptr,"==== CHK FOR FLUX CARRIED BY AVG BLK FLUX PROFILE ===\n"); 
 fprintf(fptr,"FLUX IN X-DIR \n");
 for(i=0;i<LX;i++){ 
    ex=0.0;
    for(j=0;j<LY;j++){
     q=j*LX+i; ex += avflx2[q]; }
    ex=ex/2.0; fprintf(fptr,"x=%d, EX=%f\n",i,ex);
   }
  fprintf(fptr,"FLUX IN Y-DIR \n");
  for(i=0;i<LY;i++){ 
    ey=0.0;
    for(j=0;j<LX;j++){
     q=i*LX+j; ey += avfly2[q]; }
    ey=ey/2.0; fprintf(fptr,"y=%d, EY=%f\n",i,ey);
   }
 fprintf(fptr,"==== CHK FOR FLUX CARRIED BY AVG BDY FLUX PROFILE ===\n"); 
 fprintf(fptr,"FLUX IN X-DIR \n");
 for(i=0;i<LX;i++){ 
    ex=0.0;
    for(j=0;j<LY;j++){
     q=j*LX+i; ex += avflx4[q]; }
    ex=ex/2.0; fprintf(fptr,"x=%d, EX=%f\n",i,ex);
   }
  fprintf(fptr,"FLUX IN Y-DIR \n");
  for(i=0;i<LY;i++){ 
    ey=0.0;
    for(j=0;j<LX;j++){
     q=i*LX+j; ey += avfly4[q]; }
    ey=ey/2.0; fprintf(fptr,"y=%d, EY=%f\n",i,ey);
   }
 fclose(fptr);

  /* free memory */
  free(ixc); free(iyc); free(itc);
  free(MA);  free(MB);
  free(avflx1); free(avfly1);
  free(avflx2); free(avfly2);
  free(avflx3); free(avfly3);
  free(avflx4); free(avfly4);
  free(chptr); 
  for(i=0;i<NNBR;i++) free(neigh[i]);
  for(i=0;i<2*DIM+1;i++) free(next[i]);
  deallocatedouble2d(pMAB,SIZE,SIZE);
  deallocate2d(refC,LX,LY);
  deallocate2d(fx,LT,SPV);
  deallocate2d(fy,LT,SPV);
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
  extern int **allocate2d(int,int);
  extern void deallocate2d(int**,int,int);
  int i,j,k,l,m,p,q,t;
  int p1,p2,sp1,sp2;
  double ex,ey;
    int id,jd;
  int ix,iy,it;
  int A,B,C,D;
  int val;
  int **height;
  FILE *iconf,*fout;

 /* read in initial configuration from initconf.dat */
 height=allocate2d(LX,LY);
 iconf = fopen("initconf.dat","r");
 for(i=0;i<LX;i++){
 for(j=0;j<LY;j++){
   fscanf(iconf,"%d %d %d",&id,&jd,&height[i][j]);
 }}

 fclose(iconf);

 for(p=0;p<VOL2;p++){
   ix=ixc[p]; 
   iy=iyc[p]; 
   it=itc[p];
   ising[p]=height[ix][iy];
 }

 deallocate2d(height,LX,LY);

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
 /* set pointers to the charge */
 for(i=0;i<VOL2;i++){
   A = neigh[5][neigh[2][i]]; B = neigh[5][neigh[3][i]];
   C = neigh[5][neigh[1][i]]; D = neigh[5][neigh[4][i]];
   j = neigh[6][i]; k = neigh[7][i];
   l = neigh[8][i]; m = neigh[9][i];
   chptr[i]=0;
   /* top right vertex */
   if((ising[i]!=ising[j])&&(ising[A]==ising[C])) chptr[i]=1;
   /* top left vertex */
   if((ising[i]==ising[k])&&(ising[A]!=ising[B])) chptr[i]=1;
   /* bottom left vertex */
   if((ising[i]!=ising[l])&&(ising[B]==ising[D])) chptr[i]=1;
   /* bottom right vertex */
   if((ising[i]==ising[m])&&(ising[C]!=ising[D])) chptr[i]=1;
 }


 fout = fopen("CONFCHK.dat","w");
  /* check where the flux is channeled through */
  /* For an explanation of the notation here, please refer to the routine constflux() */
  for(p=0;p<VOL2;p++){
     t=itc[p];
     q=iyc[p]*LX + ixc[p];
     sp1=ising[p];
     /* link-1 */
     p1=neigh[3][p]; p2=neigh[5][p1]; sp2=ising[p2];
     if(sp1==sp2) fy[t][q]=-1; else fy[t][q]=1;
     /* link-2 */
     p1=neigh[4][p]; p2=neigh[5][p1]; sp2=ising[p2];
     if(sp1==sp2) fx[t][q]=-1; else fx[t][q]=1;
     /* link-3 */
     p1=neigh[1][p]; p2=neigh[5][p1]; sp2=ising[p2];
     if(sp1==sp2) fy[t][next[DIM+1][q]]=-1; else fy[t][next[DIM+1][q]]=1;
     /* link-4 */
     p1=neigh[2][p]; p2=neigh[5][p1]; sp2=ising[p2];
     if(sp1==sp2) fx[t][next[DIM+2][q]]=-1; else fx[t][next[DIM+2][q]]=1;
  }
  /* print out the flux config */
/*  fprintf(fout,"====FLUX====\n");
  for(j=0;j<LY;j++){
  for(i=0;i<LX;i++){
  q=j*LY+i;
  fprintf(fout,"%d %d % d % d\n",i,j,fx[0][q],fy[0][q]);
  }}
  fprintf(fout,"====END FLUX====\n");*/


  /* Do this for the first timeslice */
  t=0;
  fprintf(fout,"FLUX IN X-DIR at t=0 \n");
  for(i=0;i<LX;i++){ 
    ex=0.0;
    for(j=0;j<LY;j++){
     q=j*LX+i;
     ex += fx[0][q]; }
    ex=ex/2.0; fprintf(fout,"x=%d, EX=%f\n",i,ex);
   }
  fprintf(fout,"FLUX IN Y-DIR at t=0 \n");
  for(i=0;i<LY;i++){ 
    ey=0.0;
    for(j=0;j<LX;j++){
     q=i*LX+j;
     ey += fy[0][q]; }
    ey=ey/2.0; fprintf(fout,"y=%d, EY=%f\n",i,ey);
   }
  fclose(fout);
}

/* sets the neighbors */
void neighbours(){
 extern int convert(int,int,int);
 int i,p,ieven,iodd,parity;
 int ixp1,ixm1,iyp1,iym1,itp1,itm1;
 int ix,iy,it;
 /* serial to checkerboard list */ 
 ieven=iodd=0;
 p=0;
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
 /* create neighbor index for flux config */
 next[DIM][p] = p;
 next[DIM+1][p] = it*SPV + iy*LX + ixp1;
 next[DIM+2][p] = it*SPV + iyp1*LX + ix;
 next[DIM+3][p] = itp1*SPV + iy*LX + ix;
 next[DIM-1][p] = it*SPV + iy*LX + ixm1;
 next[DIM-2][p] = it*SPV + iym1*LX + ix;
 next[DIM-3][p] = itm1*SPV + iy*LX + ix;
 p++;
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
  int p,q,i,j,t,m;
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
 int bondflag,chargeflag;
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
  else cflag[p]=0;
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

   /* initialize the charge-flag */
   chargeflag=0;

   /* otherwise, start building a new cluster */
   m=0; i=0; list[i]=neigh[0][p]; cflag[neigh[0][p]]=0; nclusevn++;

   do{
    im=list[m]; /* m is the new or the starting site */
    if(chptr[im]==1) chargeflag=1; /* mark the charge-carrying cluster */
    /* first check the spin on time-slice t+1 wants to bind */
    /* remember that you are on even time-slice t-1*/
    imf0=neigh[5][im];
    imf1=neigh[5][imf0];
    if(cflag[imf0]==1){
      bondflag=0;
    if((ising[imf0]==2)&&(ising[im]==ising[imf1])) { ranlxd(ran,1); if(ran[0] < p1) bondflag=1; }
    else if(ising[imf0]==0) bondflag=1; 
    if((bondflag)&&(cflag[imf1]==1)){
     i++; list[i]=imf1; /* increase list*/
    if(chptr[imf1]==1) chargeflag=1; /* mark the charge-carrying cluster */
     cflag[imf1]=0;  /* unmark spins belonging to cluster */ }
    cflag[imf0]=0;   /* unmark the interaction */
   }
   /* ============================================== */
   /* Also check if the spin in time-slice t-3 wants to bind */
   imf6=neigh[0][im];
   imf7=neigh[0][imf6];
   if(cflag[imf6]==1){
    bondflag=0;
    if((ising[imf6]==2)&&(ising[im]==ising[imf7])) { ranlxd(ran,1); if(ran[0] < p1) bondflag=1; }
    else if(ising[imf6]==0) bondflag=1; 
    if((bondflag)&&(cflag[imf7]==1)){
     i++; list[i]=imf7; /* increase list*/
    if(chptr[imf7]==1) chargeflag=1; /* mark the charge-carrying cluster */
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
     if(chptr[r1]==1) chargeflag=1; /* mark the charge-carrying cluster */
       cflag[r1]=0;  /* unmark spins belonging to cluster */ }
     if((bondflag)&&(cflag[r2]==1)){
       i++; list[i]=r2; /* increase list*/
     if(chptr[r2]==1) chargeflag=1; /* mark the charge-carrying cluster */
       cflag[r2]=0;  /* unmark spins belonging to cluster */ }
     if((bondflag)&&(cflag[r3]==1)){
       i++; list[i]=r3; /* increase list*/
     if(chptr[r3]==1) chargeflag=1; /* mark the charge-carrying cluster */
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
     if(chptr[l1]==1) chargeflag=1; /* mark the charge-carrying cluster */
       cflag[l1]=0;  /* unmark spins belonging to cluster */ }
     if((bondflag)&&(cflag[l2]==1)){
       i++; list[i]=l2; /* increase list*/
     if(chptr[l2]==1) chargeflag=1; /* mark the charge-carrying cluster */
       cflag[l2]=0;  /* unmark spins belonging to cluster */ }
     if((bondflag)&&(cflag[l3]==1)){
       i++; list[i]=l3; /* increase list*/
     if(chptr[l3]==1) chargeflag=1; /* mark the charge-carrying cluster */
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
     if(chptr[u1]==1) chargeflag=1; /* mark the charge-carrying cluster */
       cflag[u1]=0;  /* unmark spins belonging to cluster */ }
     if((bondflag)&&(cflag[u2]==1)){
       i++; list[i]=u2; /* increase list*/
     if(chptr[u2]==1) chargeflag=1; /* mark the charge-carrying cluster */
       cflag[u2]=0;  /* unmark spins belonging to cluster */ }
     if((bondflag)&&(cflag[u3]==1)){
       i++; list[i]=u3; /* increase list*/
     if(chptr[u3]==1) chargeflag=1; /* mark the charge-carrying cluster */
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
     if(chptr[d1]==1) chargeflag=1; /* mark the charge-carrying cluster */
       cflag[d1]=0;  /* unmark spins belonging to cluster */ }
     if((bondflag)&&(cflag[d2]==1)){
       i++; list[i]=d2; /* increase list*/
     if(chptr[d2]==1) chargeflag=1; /* mark the charge-carrying cluster */
       cflag[d2]=0;  /* unmark spins belonging to cluster */ }
     if((bondflag)&&(cflag[d3]==1)){
       i++; list[i]=d3; /* increase list*/
     if(chptr[d3]==1) chargeflag=1; /* mark the charge-carrying cluster */
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
      // if(A > VOL2) printf("Error. Gauss Law counter\n");
      // if(B > VOL2) printf("Error. Gauss Law counter\n");
      // if(C > VOL2) printf("Error. Gauss Law counter\n");
      // if(D > VOL2) printf("Error. Gauss Law counter\n");
       /* decide whether to join spin 1 to the cluster */
       if(ising[A] != ising[B]) {
         aux=neigh[7][im];
         if(cflag[aux]==1) { 
         i++; list[i]=aux; /* increase list */
         cflag[aux]=0; /* unmark spin belonging to cluster */ }
         /* check for the presence of the charge */
         if((chptr[A]==1)&&(chptr[B]==1)&&(chptr[im]==1)&&(chptr[aux]==1)){
           if(ising[im] != ising[aux]) printf("Error for charge\n"); }
       }
       /* decide whether to join spin 2 to the cluster */
       if(ising[A] == ising[C]) {
         aux=neigh[6][im];
         if(cflag[aux]==1) { 
         i++; list[i]=aux; /* increase list */
         cflag[aux]=0; /* unmark spin belonging to cluster */ }
         /* check for the presence of the charge */
         if((chptr[A]==1)&&(chptr[C]==1)&&(chptr[im]==1)&&(chptr[aux]==1)){
           if(ising[im] == ising[aux]) printf("Error for charge\n"); }
       }
       /* decide whether to join spin 3 to the cluster */
       if(ising[B] == ising[D]) {
         aux=neigh[8][im];
         if(cflag[aux]==1) { 
         i++; list[i]=aux; /* increase list */
         cflag[aux]=0; /* unmark spin belonging to cluster */ }
         /* check for the presence of the charge */
         if((chptr[B]==1)&&(chptr[D]==1)&&(chptr[im]==1)&&(chptr[aux]==1)){
           if(ising[im] == ising[aux]) printf("Error for charge\n"); }
       }
       /* decide whether to join spin 4 to the cluster */
       if(ising[C] != ising[D]) {
         aux=neigh[9][im];
         if(cflag[aux]==1) { 
         i++; list[i]=aux; /* increase list */
         cflag[aux]=0; /* unmark spin belonging to cluster */ }
         /* check for the presence of the charge */
         if((chptr[C]==1)&&(chptr[D]==1)&&(chptr[im]==1)&&(chptr[aux]==1)){
           if(ising[im] != ising[aux]) printf("Error for charge\n"); }
       }
   }
   m++;
   } while(m<=i);
   /* check if the list only contains genuine spins */
   for(d=0;d<=i;d++) if(list[d]>=VOL2) printf("Cluster grown in Flag.\n");
   /* check if the cluster touches the charge, and calculate the profile */
   //if(thermflag==0){
   //if(chargeflag==1) measureflux(i+1);}

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
 int chargeflag;
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

   /* initialize the charge-flag */
   chargeflag=0;

   /* otherwise, start building a new cluster */
   m=0; i=0; list[i]=neigh[0][p]; cflag[neigh[0][p]]=0; nclusodd++;
   do{
    im=list[m]; /* m is the new or the starting site */
    if(chptr[im]==1) chargeflag=1; /* mark the charge-carrying cluster */
    /* first check the spin on time-slice t+1 wants to bind */
    /* remember that you are on odd time-slice t-1*/
    imf0=neigh[5][im];
    imf1=neigh[5][imf0];
    if(cflag[imf0]==1) {
      bondflag=0;
    if((ising[imf0]==2)&&(ising[im]==ising[imf1])) { ranlxd(ran,1); if(ran[0] < p1) bondflag=1; }
      else if(ising[imf0]==0) bondflag=1; 
      if((bondflag)&&(cflag[imf1]==1)){
      i++; list[i]=imf1; /* increase list*/
      if(chptr[imf1]==1) chargeflag=1; /* mark the charge-carrying cluster */
      cflag[imf1]=0;  /* unmark spins belonging to cluster */ }
    cflag[imf0]=0;   /* unmark the interaction */
   }
   /* ============================================== */
   /* Also check if the spin in time-slice t-3 wants to bind */
   imf6=neigh[0][im];
   imf7=neigh[0][imf6];
   if(cflag[imf6]==1){
     bondflag=0;
    if((ising[imf6]==2)&&(ising[im]==ising[imf7])) { ranlxd(ran,1); if(ran[0] < p1) bondflag=1; }
     else if(ising[imf6]==0) bondflag=1; 
     if((bondflag)&&(cflag[imf7]==1)){
      i++; list[i]=imf7; /* increase list*/
     if(chptr[imf7]==1) chargeflag=1; /* mark the charge-carrying cluster */
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
     if(chptr[r1]==1) chargeflag=1; /* mark the charge-carrying cluster */
       cflag[r1]=0;  /* unmark spins belonging to cluster */ }
     if((bondflag)&&(cflag[r2]==1)){
       i++; list[i]=r2; /* increase list*/
     if(chptr[r2]==1) chargeflag=1; /* mark the charge-carrying cluster */
       cflag[r2]=0;  /* unmark spins belonging to cluster */ }
     if((bondflag)&&(cflag[r3]==1)){
       i++; list[i]=r3; /* increase list*/
     if(chptr[r3]==1) chargeflag=1; /* mark the charge-carrying cluster */
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
     if(chptr[l1]==1) chargeflag=1; /* mark the charge-carrying cluster */
       cflag[l1]=0;  /* unmark spins belonging to cluster */ }
     if((bondflag)&&(cflag[l2]==1)){
       i++; list[i]=l2; /* increase list*/
     if(chptr[l2]==1) chargeflag=1; /* mark the charge-carrying cluster */
       cflag[l2]=0;  /* unmark spins belonging to cluster */ }
     if((bondflag)&&(cflag[l3]==1)){
       i++; list[i]=l3; /* increase list*/
     if(chptr[l3]==1) chargeflag=1; /* mark the charge-carrying cluster */
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
     if(chptr[u1]==1) chargeflag=1; /* mark the charge-carrying cluster */
       cflag[u1]=0;  /* unmark spins belonging to cluster */ }
     if((bondflag)&&(cflag[u2]==1)){
       i++; list[i]=u2; /* increase list*/
     if(chptr[u2]==1) chargeflag=1; /* mark the charge-carrying cluster */
       cflag[u2]=0;  /* unmark spins belonging to cluster */ }
     if((bondflag)&&(cflag[u3]==1)){
       i++; list[i]=u3; /* increase list*/
     if(chptr[u3]==1) chargeflag=1; /* mark the charge-carrying cluster */
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
     if(chptr[d1]==1) chargeflag=1; /* mark the charge-carrying cluster */
       cflag[d1]=0;  /* unmark spins belonging to cluster */ }
     if((bondflag)&&(cflag[d2]==1)){
       i++; list[i]=d2; /* increase list*/
     if(chptr[d2]==1) chargeflag=1; /* mark the charge-carrying cluster */
       cflag[d2]=0;  /* unmark spins belonging to cluster */ }
     if((bondflag)&&(cflag[d3]==1)){
       i++; list[i]=d3; /* increase list*/
     if(chptr[d3]==1) chargeflag=1; /* mark the charge-carrying cluster */
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
       //if(A > VOL2) printf("Error. Gauss Law counter\n");
       //if(B > VOL2) printf("Error. Gauss Law counter\n");
       //if(C > VOL2) printf("Error. Gauss Law counter\n");
       //if(D > VOL2) printf("Error. Gauss Law counter\n");
       /* decide whether to join spin 1 to the cluster */
       if(ising[A] != ising[B]) {
         aux=neigh[7][im];
         if(cflag[aux]==1) { 
         i++; list[i]=aux; /* increase list */
         cflag[aux]=0; /* unmark spin belonging to cluster */ }
         /* check for the presence of the charge */
         if((chptr[A]==1)&&(chptr[B]==1)&&(chptr[im]==1)&&(chptr[aux]==1)){
           if(ising[im] != ising[aux]) printf("Error for charge\n"); }
       }
       /* decide whether to join spin 2 to the cluster */
       if(ising[A] == ising[C]) {
         aux=neigh[6][im];
         if(cflag[aux]==1) { 
         i++; list[i]=aux; /* increase list */
         cflag[aux]=0; /* unmark spin belonging to cluster */ }
         /* check for the presence of the charge */
         if((chptr[A]==1)&&(chptr[C]==1)&&(chptr[im]==1)&&(chptr[aux]==1)){
           if(ising[im] == ising[aux]) printf("Error for charge\n"); }
       }
       /* decide whether to join spin 3 to the cluster */
       if(ising[B] == ising[D]) {
         aux=neigh[8][im];
         if(cflag[aux]==1) { 
         i++; list[i]=aux; /* increase list */
         cflag[aux]=0; /* unmark spin belonging to cluster */ }
         /* check for the presence of the charge */
         if((chptr[B]==1)&&(chptr[D]==1)&&(chptr[im]==1)&&(chptr[aux]==1)){
           if(ising[im] == ising[aux]) printf("Error for charge\n"); }
       }
       /* decide whether to join spin 4 to the cluster */
       if(ising[C] != ising[D]) {
         aux=neigh[9][im];
         if(cflag[aux]==1) { 
         i++; list[i]=aux; /* increase list */
         cflag[aux]=0; /* unmark spin belonging to cluster */ }
         /* check for the presence of the charge */
         if((chptr[C]==1)&&(chptr[D]==1)&&(chptr[im]==1)&&(chptr[aux]==1)){
           if(ising[im] != ising[aux]) printf("Error for charge\n"); }
       }
    }
    m++;
   } while(m<=i);
   /* check if the list only contains genuine spins */
   for(d=0;d<=i;d++) if(list[d]>=VOL2) printf("Cluster grown in Flag.\n");

   /* check if the cluster touches the charge, and calculate the profile */
   //if(thermflag==0){
   //if(chargeflag==1) measureflux(i+1);}

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

double **allocatedouble2d(int row, int col)
{
  int i,j;
  double **mat;
  mat = (double **)malloc(row*sizeof(double*));
  if(mat==NULL) {printf("Out of memory\n"); exit(0);}

  for(i=0;i<row;i++){
   mat[i]=(double *)malloc(col*sizeof(double));
   if(mat[i]==NULL)  {printf("Out of memory\n"); exit(0);}}

 for(i=0;i<row;i++) for(j=0;j<col;j++) mat[i][j]=0.0;
 return mat;
}

void deallocatedouble2d(double **mat, int row, int col){
  int i;
  for(i=0;i<row;i++) free(mat[i]);
  free(mat);
}

int **allocate2d(int row, int col)
{
  int i,j,**mat;
  mat = (int **)malloc(row*sizeof(int*));
  if(mat==NULL) {printf("Out of memory\n"); exit(0);}

  for(i=0;i<row;i++){
   mat[i]=(int *)malloc(col*sizeof(int));
   if(mat[i]==NULL)  {printf("Out of memory\n"); exit(0);}
  }

 for(i=0;i<row;i++) for(j=0;j<col;j++)
  mat[i][j]=0;

 return mat;
}

void deallocate2d(int **mat, int row, int col)
{
  int i;
  for(i=0;i<row;i++) free(mat[i]);
  free(mat);
}

void measureMAB(){
 int ix,iy,it,it2,p;
 /* initialize MA,MB */
 for(it=0;it<LT2;it++)
   MA[it]=MB[it]=0;
 for(p=0;p<VOL2;p++){
  ix=ixc[p]; iy=iyc[p]; it=itc[p]; it2=it/2;
  if((it%2)==0){ /* even-time slice; measure MA*/
  if(ising[p]==refC[ix][iy]) MA[it2]++;
  else MA[it2]--;
  }
  else{ /* odd-time slice; measure MB */
  if(ising[p]==refC[ix][iy]) MB[it2]++;
  else MB[it2]--;
  } }
  /* get averaged probability dist over all time-slices */
  for(it=0;it<LT2;it++){
    pMAB[MA[it]-minMA][MB[it]-minMB]++;
  }
}

void constflux(){
  int p,p1,p2,q,x,y,t,tf;
  int x1,y1,x2,y2;
  int sp1,sp2,f1,f2;
  int ptA,ptB;
  int nc1,nac1,nc2,nac2;
  int Q;
  int flagx,flagy;

  /* initialize the flux-config */
  for(t=0;t<LT;t++) for(p=0;p<SPV;p++) fx[t][p]=fy[t][p]=0;

  /* Construct the flux-configuration from the height configurations 
     The fluxes in the x-and y-direction are obtained by considering 
     the fwd and backward height variables for each individual height variable
      O  |  X  |  O       Note that at for each height variable 
         |     |          (denoted by O), the four height variables on the
    -----o--4--o-----     surrounding sites (denoted by X) are used to construct
         |     |          the forward and backward links surrounding O (denoted by
      X  1  p  3  X       1,2,3,4). The current index is denoted by p.
         |     |
    -----o--2--o-----
         |     |
      O  |  X  |  O
  */

  for(p=0;p<VOL2;p++){
     t=itc[p];  
     q=iyc[p]*LX + ixc[p];
     sp1=ising[p]; 
     /* link-1 */ 
     p1=neigh[3][p]; p2=neigh[5][p1];  
     sp2=ising[p2];
     if(sp1==sp2) fy[t][q]=-1; else fy[t][q]=1;
     /* link-2 */
     p1=neigh[4][p]; p2=neigh[5][p1];
     sp2=ising[p2];
     if(sp1==sp2) fx[t][q]=-1; else fx[t][q]=1;
     /* link-3 */
     p1=neigh[1][p]; p2=neigh[5][p1];
     sp2=ising[p2];
     if(sp1==sp2) fy[t][next[DIM+1][q]]=-1; else fy[t][next[DIM+1][q]]=1;
     /* link-4 */
     p1=neigh[2][p]; p2=neigh[5][p1];
     sp2=ising[p2];
     if(sp1==sp2) fx[t][next[DIM+2][q]]=-1; else fx[t][next[DIM+2][q]]=1;
  }

  /* Measure charge profile */
  /* First, detect flux across the boundary at timeslice t=0 */
  tf=0;
  flagx=flagy=0;
  for(x=0;x<LX;x++){
     p = (LY-1)*LX + x;
     flagy = flagy + fy[tf][p]; }
  for(y=0;y<LY;y++){
     p = y*LX + LX-1;
     flagx = flagx + fx[tf][p]; }
  if(abs(flagx)==4) flagxx++;
  if(abs(flagy)==4) flagyy++;


  /* check the charge at site cpos at timeslice=0 */
  p = cpos;
  Q = fx[tf][p] - fx[tf][next[DIM-1][p]] + fy[tf][p] - fy[tf][next[DIM-2][p]];
  if(Q==4) { 
   if((flagx==0)&&(flagy==0)){   /* positive charge and flux in bulk */
     flxcnt1++;
     for(t=0;t<LT;t++){ for(p=0;p<SPV;p++){
          avflx1[p] += fx[t][p]; avfly1[p] += fy[t][p]; } }
   }
   else{                     /* positive charge and flux in boundary */
     flxcnt3++;
     for(t=0;t<LT;t++){ for(p=0;p<SPV;p++){
          avflx3[p] += fx[t][p]; avfly3[p] += fy[t][p]; } }
   } 
  }
  else if(Q==-4){ 
   if((flagx==0)&&(flagy==0)){   /* negative charge and flux in bulk */
     flxcnt2++;
     for(t=0;t<LT;t++){ for(p=0;p<SPV;p++){
          avflx2[p] += fx[t][p]; avfly2[p] += fy[t][p]; } }
   }
   else{                     /* negative charge and flux in boundary */
     flxcnt4++;
     for(t=0;t<LT;t++){ for(p=0;p<SPV;p++){
          avflx4[p] += fx[t][p]; avfly4[p] += fy[t][p]; } }
   } 
  }

   /* print height conf */
   //for(p=0;p<VOL;p++) printf("%d %d %d\n",p,itc[p],ising[p]);
   //printf("Flux conf:\n");
   //for(t=0;t<LT;t++)
   //for(p=0;p<SPV;p++) printf("%d %d %d %d\n",t,p,fx[t][p],fy[t][p]);
}

 void energy(){
  int p,imf,imb;
  inten=intpe=0.0;
  /* go over all the interactions */
  for(p=VOL2;p<VOL;p++){
     if(ising[p]==2){
       intpe--;
       imf=neigh[5][p]; imb=neigh[0][p];
       if(ising[imf]==ising[imb]) inten += -lam + tanhx;
       else if(ising[imf]!=ising[imb]) inten += -lam + cothx;
     }
  }
  inten = -inten/((double)(LT2)); 
  intpe = -intpe/((double)(LT2)); 
 }

 /* This routine provides an initial height configuration corresponding to */
 /* a negative charge Q=-2 placed at (ixm,ixm) and a positive charge Q=+2  */
 /* placed at (ixp,ixp), with ixm < ixp                                    */
 void creatediagconf(int ixm,int ixp){
    int **iex,**iey,**ih;
    int L;
    int ix,iy,ixmm1,ixpp1,ix1,iy1,ix1m1,iy1m1;
    int ix2,iy2,ix2m1,iy2m1;
    double avx,avy;
    FILE *fptr;

    iex=allocate2d(LX,LY); iey=allocate2d(LX,LY);
    ih =allocate2d(LX,LY);

    /* assuming square box */
    L=LX;

    /* initialize a reference configuration of height variables */
    for(ix=0;ix<L;ix++) for(iy=0;iy<L;iy++) {
      iex[ix][iy] = 2*((ix+iy+1)%2) - 1;
      iey[ix][iy] = 2*((ix+iy)%2) - 1;
    }

    /* insert a charge -2 at (ixm,ixm) and a charge +2 at (ixp,ixp) */
    /* by modifying fluxes accordingly                              */
    ixmm1 = ((ixm-1+L)%L); ixpp1 = ((ixp+1)%L);
    iex[ixmm1][ixm] = 1;  iey[ixmm1][ixm] =-1; 
    iex[ixm][ixmm1] =-1;  iey[ixm][ixmm1] = 1;
    iex[ixp-1][ixpp1]=-1; iey[ixpp1][ixp-1]=-1; 
    iex[ixp][ixp]   = 1;  iey[ixp][ixp]  = 1;
    for(ix=ixm;ix<ixp;ix++){
      ix1 = ix; ix1m1 = (ix1-1+L)%L;
      iy1 = ix + 1; iy1m1 = iy1 - 1;
      iex[ix1][iy1]  = -1; iey[ix1][iy1] = -1; 
      iex[ix1m1][iy1]= -1; iey[ix1] [iy1m1]=-1;
      ix2 = ix+1;  ix2m1 = ix2-1;
      iy2 = ix;  iy2m1=(iy2-1+L)%L;
      iex[ix2][iy2] = -1; iey[ix2][iy2] = -1;
      iex[ix2m1][iy2]=-1; iey[ix2][iy2m1]=-1;
    }

    /* calculate flux carried */
    printf("Flux in the x-dir: \n");
    for(ix=0;ix<L;ix++){
    avx=0.0;
    for(iy=0;iy<L;iy++){
    avx += iex[ix][iy];
    }
    printf("%d % f\n",ix,avx/2.0);
    }
    printf("Flux in the y-dir: \n");
    for(iy=0;iy<L;iy++){
    avy=0.0;
    for(ix=0;ix<L;ix++){
    avy += iey[ix][iy];
    }
    printf("%d % f\n",iy,avy/2.0);
    }

    /* calculate the height configuration */
    ih[0][0]=1;
    for(ix=1;ix<L;ix++){
     if(iey[ix][0]==1) ih[ix][0] = -ih[ix-1][0];
     if(iey[ix][0]==-1)ih[ix][0] =  ih[ix-1][0];
    }
    for(iy=1;iy<L;iy++){
    for(ix=0;ix<L;ix++){
       if(iex[ix][iy]==1) ih[ix][iy] = -ih[ix][iy-1];
       if(iex[ix][iy]==-1)ih[ix][iy] =  ih[ix][iy-1];
    }
    }

    /* write out the height variables */
    fptr = fopen("initconf.dat","w");
    for(ix=0;ix<L;ix++){
    for(iy=0;iy<L;iy++){
     fprintf(fptr,"%d %d %d\n",ix,iy,ih[ix][iy]);
    }}
    fclose(fptr);

    deallocate2d(iex,LX,LY);
    deallocate2d(iey,LX,LY);
    deallocate2d(ih,LX,LY);
}

  
 /* This routine provides an initial height configuration corresponding to   */
 /* a negative charge Q=-2 placed at (ixm,L/2-1) and a positive charge Q=+2  */
 /* placed at (ixp,L/2+1), with ixm < ixp, and ixm%2==0                      */
 void createaxisconf(int ixm,int ixp){

    int **iex,**iey,**ih;
    int L,l2,l2p2,l2p1,l2m1,l2m2,ixpp1,ixp1;
    int ix,iy;
    double avx,avy;
    FILE *fptr;

    iex=allocate2d(LX,LY); iey=allocate2d(LX,LY);
    ih =allocate2d(LX,LY);

    /* assuming square box */
    L=LX;

    /* initialize a reference configuration of height variables */
    for(ix=0;ix<L;ix++) for(iy=0;iy<L;iy++) {
      iex[ix][iy] = 2*((ix+iy+1)%2) - 1;
      iey[ix][iy] = 2*((ix+iy)%2) - 1;
    }
    /* insert a charge -2 at (ixm,L/2-1) and a charge +2 at (ixp,L/2+1) */
    /* by modifying fluxes accordingly                                  */
    l2 = L/2-1;
    l2p2 = (l2+2)%L;
    l2p1 = (l2+1)%L;
    l2m1 = (l2-1+L)%L;
    l2m2 = (l2-2+L)%L;
    ixpp1= ((ixp+1)%L);
    iey[ixm][l2m1] =  1;
    iey[ixm][l2]   = -1;
    for(ix=ixm;ix<ixp;ix++){
      ixp1 = ((ix+1)%L);
      iex[ix][l2m1] = -1;
      iex[ix][l2] = -1;
      iex[ix][l2p1] = -1;
      iex[ix][l2p2] = -1;
      iey[ixp1][l2m2] = 2*(ix%2)-1;
      iey[ixp1][l2m1] = 2*(ix%2)-1;
      iey[ixp1][l2] = 2*(ix%2)-1;
      iey[ixp1][l2p1] = 2*(ix%2)-1;
    }
    if((ixp%2)==0){
      iey[ixp][l2m1] = -1;
      iex[ixp][l2] = 1;
      iex[ixp][l2p1] = -1;
      iey[ixpp1][l2] = 1;
    }
    else{
      iey[ixp][l2p1] = 1;
      iey[ixp][l2] = 1;
      iex[ixp][l2] = 1;
      iex[ixp][l2m1] = -1;
      iey[ixpp1][l2m1] = -1;
    }
    /* print the electric fluxes */
    for(ix=0;ix<L;ix++) for(iy=0;iy<L;iy++)
    printf("%d %d %d %d\n",ix+1,iy+1,iex[ix][iy],iey[ix][iy]);

    /* calculate flux carried */
    printf("Flux in the x-dir: \n");
    for(ix=0;ix<L;ix++){
    avx=0.0;
    for(iy=0;iy<L;iy++){
    avx += iex[ix][iy];
    }
    printf("%d % f\n",ix,avx/2.0);
    }
    printf("Flux in the y-dir: \n");
    for(iy=0;iy<L;iy++){
    avy=0.0;
    for(ix=0;ix<L;ix++){
    avy += iey[ix][iy];
    }
    printf("%d % f\n",iy,avy/2.0);
    }

    /* calculate the height configuration */
    ih[0][0]=1;
    for(ix=1;ix<L;ix++){
     if(iey[ix][0]==1) ih[ix][0] = -ih[ix-1][0];
     if(iey[ix][0]==-1)ih[ix][0] =  ih[ix-1][0];
    }
    for(iy=1;iy<L;iy++){
    for(ix=0;ix<L;ix++){
       if(iex[ix][iy]==1) ih[ix][iy] = -ih[ix][iy-1];
       if(iex[ix][iy]==-1)ih[ix][iy] =  ih[ix][iy-1];
    }
    }

    /* write out the height variables */
    fptr = fopen("initconf.dat","w");
    for(ix=0;ix<L;ix++){
    for(iy=0;iy<L;iy++){
     fprintf(fptr,"%d %d %d\n",ix,iy,ih[ix][iy]);
    }}
    fclose(fptr);

    deallocate2d(iex,LX,LY);
    deallocate2d(iey,LX,LY);
    deallocate2d(ih,LX,LY);
 }
