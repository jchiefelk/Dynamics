/*Calculate when nanocavity contains water during the simulation*/ 


#include	<stdio.h> /*an input-output library*/
#include	<stdlib.h> /*other functions*/
#include	<math.h> /*math library*/
#include    <string.h> /*String Function Library*/
#define pi 3.1415926


int main(int argc,char *argv[]){
    // This little piece will stop execution if wrong input//
    if (argc != 7){
        printf("%s","Usage: %s distancefile(.txt) poreradius(nm) Simlength(ps) Timestep(ps) Nsol outfile(.txt)\n");
        exit(1);
                 }
	else{

		/*Control Variables*/
		int i,j,sum,ID;
		int nT = atof(argv[3])/atof(argv[4]); /*number of timesteps*/
		double timestep=atof(argv[4]);
		int Nsol = atoi(argv[5]);
		double PR = atof(argv[2]);
		char buf[180],infile[40],outfile[40];
		FILE *fp;
		FILE *out;
		sprintf(infile,"%s",argv[1]);
		sprintf(outfile,"%s",argv[6]);
		fp=fopen(infile,"r");
		out=fopen(outfile,"w");
		/*Occupancy List*/
		int *occ=(int*)malloc(nT*sizeof(int));
		double Dist,time; /*distance from cavity COM,simtime*/
		for(i=0;i<nT;i++){
			occ[i]=0;
			sum=0;
			for(j=0;j<Nsol;j++){
				//Distance Index//
				fgets(buf,180,fp);
				sscanf(buf,"%lf %i %lf\n",&time,&ID,&Dist);
			
				//Conditional Operation on Data//
				if(Dist < PR) sum+=1;
	
			}			
			occ[i]=sum;
			fprintf(out,"%f %i\n",timestep*i,occ[i]);
		}
		fclose(fp);
		free(occ);
		
	}
}




