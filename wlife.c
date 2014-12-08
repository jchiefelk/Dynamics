#include	<stdio.h>
#include <stdlib.h>
#include	<math.h>
#define NW 4559 /* number of water must be less than this number*/
#define MaxT 73061 /*length of trajectory in ps*/
#define MaxNT 146122 /*MaxT/ts: number of time points*/
/*****Input file structure***********
 0  5 SOL 157 OW  0.224664 (nm)00
 0  24 SOL 214 OW  0.222265 (nm)
 0  266 SOL 940 OW  0.418863 (nm)
 0  276 SOL 970 OW  0.197458 (nm)
 0  316 SOL 1090 OW  0.487795 (nm)
 0  483 SOL 1591 OW  0.244727 (nm)
 0  526 SOL 1720 OW  0.471246 (nm)
 0  540 SOL 1762 OW  0.322147 (nm)
 0  707 SOL 2263 OW  0.364377 (nm)
 0  826 SOL 2620 OW  0.480705 (nm)
 0  898 SOL 2836 OW  0.212354 (nm)
 0  950 SOL 2992 OW  0.225353 (nm)
 0  954 SOL 3004 OW  0.460464 (nm)
 0.05  5 SOL 157 OW  0.221062 (nm)
....
....
*****/

main(argc, argv)
int argc;
char *argv[];
{
char inFile[80],buf[180],s1[3],s2[2];
int nw, nO, i,l,j,k, Td;
double t, r, avT, ts, atof();
int vec[NW] = {0};
int pvec[NW] = {0};
int T[NW] = {0};
int To[NW] = {0};
int Ts[NW] = {0};
double cor[MaxNT] = {0.0};
FILE *fp;

if (argc != 4){
        fprintf(stderr,"Usage: %s input_file time_step(ps) time_delay(ps)\n", argv[0]);
        exit(1); 
}
ts = atof(argv[2]);
Td = atof(argv[3])/ts; 

sprintf(inFile,"%s",argv[1]);
fp = fopen(inFile,"r");
fprintf(stderr,"open file %s\n",inFile);
l = 0;/*previous time index initialzed*/
while (fgets(buf,180,fp) != NULL){/*read until end of file*/
	sscanf(buf,"%lf %d %s %d %s %lf",&t,&nw,s1,&nO,s2,&r);
	if (nw > NW -1){
		fprintf(stderr,"nw = %d > %d. Increase NW\n",nw, NW-1);
		exit(1);
	}
	i = (t+0.1*ts)/ts;/*convert time to an integr index*/
	if (i == l){/*still in the same time step*/
		vec[nw] = 1;/*this water is inside the cavity*/
		T[nw] += 1;/*clock increament for this water molecule*/
	}
	if ( i == l+1){/*first molecule in the next time step*/
	    l = i;/*reset*/
	    if (i > 0){/*time correlation starting from 2nd time step*/
		for (j = 0; j< NW; j++){
			if (vec[j] == 1 && pvec[j] == 1)
				cor[T[j]] += 1.0;/* j remains inside*/
			if (vec[j] == 0 && pvec[j] == 1){/* j just exited*/
				Ts[j] = T[j];/*save for delay option*/
				T[j] = 0;/*j exit, clock reset*/
				To[j] += 1;/*for the delay option*/
			}
			if (vec[j] == 1 && pvec[j] == 0){/* j entered*/
				if (To[j] < Td){
					T[j] = Ts[j]+To[j];/*restore*/
					for (k = 0; k<To[j]; k++)
						cor[Ts[j]+k+1] +=1.0;
				}
			}
		}
	   }
		for (j = 0; j< NW; j++){/*update residence arrays*/
			pvec[j] = vec[j];/*save previous time vector*/
			vec[j] = 0;/*initialize for next time step*/
		}
		vec[nw] = 1;/*store current molecule*/
		T[nw] +=1;
	}
	if (i < l || i > l+1) {
		fprintf(stderr,"i = %d l= %d\n",i,l);
		exit(1);
	}
}/*end of file*/
avT = 0;
for (i = 2; i < MaxNT; i++){
	printf("%f %f\n",(i-2)*ts,cor[i]/cor[2]);
	avT += cor[i]/cor[2];
}
fprintf(stderr,"delayT = %f avT = %f\n",Td*ts, avT*ts);
}
