/**********CALCULATES DIPOLE CORRELATION FUNCTION WITH DELAY******************/


#include	<stdio.h> /*an a1s2d3f$-output library*/
#include	<stdlib.h> /*other functions*/
#include	<math.h> /*math library*/
#include    <string.h> /*String Function Library*/


double X,Y,Z; //Dipole Vectors//
void dipole(double OX,double OY,double OZ,double H1X,double H1Y,double H1Z,double H2X,double H2Y,double H2Z)
{
    X = H1X+H2X-2*OX;
    Y = H1Y+H2Y-2*OY;
    Z = H1Z+H2Z-2*OZ;
}

//Correlation Dot Protduct Function//
double Corr;
void AutoCorr(double CX,double CY,double CZ)
{
     Corr=0.0;
     Corr=CX+CY+CZ;
}

//Vector Magnitude Function//
double MAG;
void MAGNITUDE(double X, double Y, double Z)
{
    MAG=0.0;
	MAG = sqrt(X*X+Y*Y+Z*Z);
}


int main(int argc,char *argv[]){
    // This little piece will stop execution if wrong input//
    if (argc != 10){
        printf("%s","Usage: %s (1)InputFile(.gro) (2)DistanceFile(.txt) (3)SimLength(ps) (4)TimeStep(ps) (5)MaxTau(ps) (6)Delay(ps) (7)Nsol (8)PoreRadius (9)OUTFILE\n");
        exit(1);
                 }
	else{
        /****Control Variables****/
        int i,j,k,tau; //Loop indices//
        /****Data Structure Info******/
        int nT,Nsol,Nelements,MaxTau,Ndiv,Delay; /***Number of solvent molecules,Number of Timesteps,Number of elements per stimestep,Max Lag Time***/
        double TimeStep,PoreRadius;
        /*Time Structure*/
        TimeStep = atof(argv[4]); /*Total simulation time in ps*/
        nT = atof(argv[3])/TimeStep; /*plus 1 for 0 timestep*/
        MaxTau = atoi(argv[5])/TimeStep; /*Max Lag, TotalNumber of data elements per slice*/
        Delay = atoi(argv[6])/TimeStep; /*Time Delay*/
        Ndiv = nT/MaxTau; /*Number of Time Slices to Analyze*/
        /*Solvent Structure*/
        Nsol = atoi(argv[7]);
		PoreRadius= atof(argv[8]);
        Nelements = Nsol*3;
        /*File Input*/
        char ID[2]="OW";
        char inFile[80],index[80],buf[80];
        FILE *fp,*fID;
		sprintf(inFile,"%s",argv[1]);
		sprintf(index,"%s",argv[2]);
		fID = fopen(index,"r"); /*Index file identifier*/
		fp = fopen(inFile,"r"); /*Gro data file identifier*/
        char OUTFILE[40];
		FILE *OUT;
		sprintf(OUTFILE,"%s",argv[9]);
		OUT = fopen(OUTFILE,"w"); /*Index file identifier*/
        
		//Initialize Correlation array//
		double *AC=(double*)malloc((MaxTau+1)*sizeof(double));/*Pore Autocorrelation*/
        double *BC=(double*)malloc((MaxTau+1)*sizeof(double));/*Pore Autocorrelation*/
		double *Dip_corrB=(double*)malloc((MaxTau+1)*sizeof(double));/*Total Bulk Dip Correlation*/
		double *Dip_corrI=(double*)malloc((MaxTau+1)*sizeof(double));/*Total Bulk Dip Correlation*/
		int *num1=(int*)malloc((MaxTau)*sizeof(int));/*Number of Bulk Water*/
	    int *num2=(int*)malloc((MaxTau)*sizeof(int));/*number of pore waters*/
        for(i=0;i<=MaxTau;i++) {
            Dip_corrI[i] = AC[i] = 0.0;
        }
		//Dipole Data Arrays//
		double **Dip_Delay;
		double **Dip;
		double **DipB;
        Dip_Delay = malloc((MaxTau+Delay)*sizeof(double));
		DipB=malloc((MaxTau+Delay)*sizeof(double));
        Dip = malloc((MaxTau+Delay)*sizeof(double));
		for(i=0;i<Delay+MaxTau;i++){
		Dip_Delay[i]=malloc((3*Nsol)*sizeof(double));
		Dip[i]=malloc((3*Nsol)*sizeof(double));
		DipB[i]=malloc((3*Nsol)*sizeof(double));
        }


		int AtomID,flag;
        double OX,OY,OZ,H1X,H1Y,H1Z,H2X,H2Y,H2Z,OWDist,Dist,CX,CY,CZ; //Dipole and Correlation variables//
        double dipx,dipy,dipz;
        char SolID[8],type[3],OWflag[2];
        //Get Distance File Headers//
        fgets(buf,180,fID);
        fgets(buf,180,fID);
		
        double time;
        for(i=0;i<Ndiv;i++){
			//Initial MultiD Arrays to 0//
                for(j=0;j<(MaxTau+Delay);j++){
                    for(k=0;k<3*Nsol;k++){
                        Dip[j][k] = 0.0;
                }
                }

                for(j=0;j<MaxTau+Delay;j++){
                    //Get first two lines in Gro File//
                    fgets(buf,180,fp);
                    fgets(buf,180,fp);
                    num1[j]=0;
				for(k=0;k<Nsol;k++){
					//Get Gro File Data//
                    fgets(buf,180,fp);
                    sscanf(buf,"%s %2s %i %lf %lf %lf",&SolID,&type,&AtomID,&OX,&OY,&OZ);
					fgets(buf,180,fp);
                    sscanf(buf,"%s %2s %i %lf %lf %lf",&SolID,&type,&AtomID,&H1X,&H1Y,&H1Z);
					fgets(buf,180,fp);
                    sscanf(buf,"%s %2s %i %lf %lf %lf",&SolID,&type,&AtomID,&H2X,&H2Y,&H2Z);
                        //Calculate Dipole Vectors//
                    //Distance Index//
                    fgets(buf,180,fID);
                    sscanf(buf,"%lf %i %lf\n",&time,&flag,&OWDist);
                    fgets(buf,180,fID);
                    sscanf(buf,"%lf %i %lf\n",&time,&flag,&Dist);
                    fgets(buf,180,fID);
                    sscanf(buf,"%lf %i %lf\n",&time,&flag,&Dist);
                    if(OWDist < Dist){
						dipole(OX,OY,OZ,H1X,H1Y,H1Z,H2X,H2Y,H2Z);
                        MAGNITUDE(X,Y,Z);
						Dip[j][3*(k+1)-3]=X/MAG; /*Make Unit Vector*/
						Dip[j][3*(k+1)-2]=Y/MAG;
						Dip[j][3*(k+1)-1]=Z/MAG;
                        num1[j]+=1; /*Number of water is pore at this timestep*/
                    }
                }
				fgets(buf,180,fp); //Get Last line In Frame of Grofile//
                }
			
            
			//Correlation Function//
            double sum,norm,clock,icor;
			for(tau=0;tau<MaxTau;tau++){
                AC[tau]=0.0;
                for(k=0;k<Nsol;k++){
                    icor=0.0;
                    sum=0.0;
                    clock=0.0;
                    norm=0.0;
                	for(j=0;j<MaxTau-tau;j++){
                                //Pore Residence Clock//
                                if(Dip[j][3*(k+1)-1]!=0.0){
                                    clock+=1.0;
                                    norm+=1.0;
                                }
            
                                //Pore Correlation//
                                if(clock > 0.0 && Dip[j+tau][3*(k+1)-1]!=0.0 ) {
									CX = Dip[j][3*(k+1)-3]*Dip[j+tau][3*(k+1)-3];
									CY = Dip[j][3*(k+1)-2]*Dip[j+tau][3*(k+1)-2];
									CZ = Dip[j][3*(k+1)-1]*Dip[j+tau][3*(k+1)-1];
									AutoCorr(CX,CY,CZ);
                                    sum+=Corr;
                            }
                        
                        
                                //Reset Clock//
                                if(clock > 0.0 && Dip[j][3*(k+1)-1]==0.0){
                                    clock=0.0;
                                }
                                
                                
                    }
                
            
                    icor += sum/norm;
                    
                    if(norm!=0.0) AC[tau]+=icor;

                    
                    
                }
            
			}
			
            
            //Collect in Corr array//
            for(tau=0;tau<MaxTau;tau++){
				Dip_corrI[tau] += AC[tau]/AC[0];
            }
		
		
		}
		
		/******End of Correlation Loop*****/
        fclose(fp);
		fclose(fID);
	
        /*Print Out Results*/
		for(i=0;i<MaxTau;i++){
			fprintf(OUT,"%lf %lf\n",i*TimeStep,Dip_corrI[i]/Ndiv);
            printf("%lf %lf\n",i*TimeStep,Dip_corrI[i]/Ndiv);
		}
        fclose(OUT);
	
		 for(i=0;i<MaxTau+Delay;i++){
			 free(Dip[i]);
			 free(Dip_Delay[i]);
		 }
		 free(Dip);
		 free(Dip_Delay);
		
    }
    
}
