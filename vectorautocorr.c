/************Orinetational Correlation Function******************/

#include	<stdio.h>
#include	<stdlib.h>
#include	<math.h>
#define MaxX 4.49011  /*X-Dim of simulation box*/
#define MaxY 4.51261 /*Y-Dim of simulation box*/


#define SQR(x) ((x)*(x))
#define MaxTheta 180

/*Calculates Radial Distribution from Gromacs Data*/
/* Input data from file */
int main(int argc, char *argv[]) {
	
	if (argc != 6){
        printf("%s","Usage: %s InFile(.gro) SimLength(ps) TimeStep(ps) MaxTau(ps) OUTFILE(.txt)\n");
        exit(1);
    }
	
	else {
		
        int i,j,tau,MaxTau,nDiv,AtomID;
        double nT,TimeStep,Time,Theta,sum,MAG;
        double O1X,O1Y,O1Z,O2X,O2Y,O2Z,O3X,O3Y,O3Z;
        double x1,y1,z1,x2,y2,z2;
	    char infile[40],outfile[40],type[2],SolID[8],buf[180];
		FILE *fp;
        FILE *OUT;
        TimeStep=atof(argv[3]);
        nT = atoi(argv[2])/TimeStep;
        MaxTau = atoi(argv[4])/TimeStep;
        nDiv = nT/MaxTau;
        sprintf(infile,"%s",argv[1]);
        sprintf(outfile,"%s",argv[5]);
        OUT = fopen(outfile,"w");
		printf("Good artists copy\n");
        double *AC=(double*)malloc(MaxTau*sizeof(double));
		double *TC=(double*)malloc(MaxTau*sizeof(double));
        //Correlation Function Array//
        for(i=0;i<MaxTau;i++){
            TC[i]=AC[i]=0.0;
        }
		fp = fopen(infile,"r");
		//Calculate Correlation Function//
        for(i=0;i<nDiv;i++){
            double *ihat=(double*)malloc(MaxTau*sizeof(double));
            double *jhat=(double*)malloc(MaxTau*sizeof(double));
            double *khat=(double*)malloc(MaxTau*sizeof(double));
                for(j=0;j<MaxTau;j++){
                    //Get first two lines in Gro File//
                    fgets(buf,180,fp);
                    fgets(buf,180,fp);
                    //Get Gro File Data//
                    fgets(buf,180,fp);
                    sscanf(buf,"%s %2s %i %lf %lf %lf",&SolID,&type,&AtomID,&O1X,&O1Y,&O1Z);
                    fgets(buf,180,fp);
                    sscanf(buf,"%s %2s %i %lf %lf %lf",&SolID,&type,&AtomID,&O2X,&O2Y,&O2Z);
                    //Calculate Vector1//
                    x1 = O2X-O1X;
                    y1 = O2Y-O1Y;
                    z1 = O2Z-O1Z;
             
		/***Account for periodic boundary conditions in x-y plane***/
		/**if(x1 > MaxX/2) x1-=MaxX;

                if(x1 < -MaxX/2) x1+=MaxX;
                
                if(y1 > MaxY/2) y1-=MaxY;

                if(y1 < -MaxX/2) y1+=MaxY;

		***/

                    //Calculate Norm//
                    ihat[j] = x1;
                    jhat[j] = y1;
                    khat[j] = z1;
                    MAG = sqrt(ihat[j]*ihat[j]+jhat[j]*jhat[j]+khat[j]*khat[j]);
                    ihat[j]/=MAG;
                    jhat[j]/=MAG;
                    khat[j]/=MAG;
                    //Get Last line In Frame of Grofile//
                    fgets(buf,180,fp);
                }
                for(tau=0;tau<MaxTau;tau++){
					AC[tau] = 0.0;
                            for(j=0;j<MaxTau-tau;j++){
                                //Calculate Vector Magnitudes//
                                AC[tau] += ihat[j]*ihat[j+tau]+jhat[j]*jhat[j+tau]+khat[j]*khat[j+tau];
                            }
                }
			
			for(tau=0;tau<MaxTau;tau++){
				TC[tau]+=AC[tau]/(MaxTau-tau);	
			}
            free(ihat);
            free(jhat);
            free(khat);
		}
        fclose(fp);
        
        
		for(j=0;j<MaxTau;j++){
			fprintf(OUT,"%lf %lf\n",j*TimeStep,TC[j]/nDiv);
			printf("%lf %lf\n",j*TimeStep,TC[j]/nDiv);
		}
        free(AC);
    }
	
}
