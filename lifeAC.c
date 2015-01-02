
/****Calculates Lifetime Correlation Function***/
/****Create Matrix whose jk entry is 1 if a solvent molecule is less than the user
 specified distance*****/


#include	<stdio.h> /*an input-output library*/
#include	<stdlib.h> /*other functions*/
#include	<math.h> /*math library*/
#include    <string.h>


int main(int argc,char *argv[]){
    // This little piece will stop execution if wrong input//
    if (argc != 8){
        printf("%s","Usage: %s (1)Positions(.txt) (2)SimLength(ps) (3)TimeStep(ps) (4)MaxTau(ps) (5)PoreRadius(nm) (6)Nsol (7)OUT_file\n");
        exit(1);
    }
	else{
        
        
        int i,j,k,tau,maxtau,nT,Nsol,Ndiv; /*Loop Indices*/
        char infile[80],outfile[80],buf[80];
        double simlength,timestep,poreradius;
        simlength=atof(argv[2]); /*Length of Simulation in Picoseconds*/
        timestep=atof(argv[3]); /*Simulation Resolution*/
        nT = simlength/timestep; /* Total Number of Data Points */
        maxtau = atoi(argv[4])/timestep; /*Max Lag Time*/
        poreradius=atof(argv[5]);
        Nsol = atoi(argv[6]);
        Ndiv = nT/maxtau;

        
        //Input//
        FILE *fp;
        sprintf(infile,"%s",argv[1]);
        fp = fopen(infile,"r"); /*Index file identifier*/
        
        //Output//
        char OUTFILE[40];
        FILE *OUT;
        sprintf(OUTFILE,"%s",argv[7]);
        OUT = fopen(OUTFILE,"w"); /*Index file identifier*/
        
        //Lifetime Array//
        double *Life=(double*)malloc((maxtau)*sizeof(double));/*Pore Autocorrelation*/
        
        //AutoCorrelation//
        double *AC=(double*)malloc((maxtau)*sizeof(double));/*Pore Autocorrelation*/
        
        //Data Arrays//
        //Occupancy Stored as Binary//
        double **Data;
        Data = malloc((maxtau)*sizeof(double));
        for(i=0;i<maxtau;i++){
            Data[i]=malloc((Nsol)*sizeof(double));
        }
        
        //Number of Pore Solvent At time tau//
        double *porenum=(double*)malloc((maxtau)*sizeof(double));
        for(i=0;i<=maxtau;i++) {
           AC[i] = porenum[i] = 0.0;
        }
        
        
        //File Structure//
        int SolID;
        double Time,dist;
        
        /******************CORRELATION LOOP****************************/
        //Get Distance File Headers//
        fgets(buf,180,fp);
        fgets(buf,180,fp);
        for(i=0;i<Ndiv;i++){
    
            for(j=0;j<maxtau;j++){
                for(k=0;k<Nsol;k++){
                    Data[j][k] = 0;
                }
            }
            
            
            //Collect Data in Matrix//
            for(j=0;j<maxtau;j++){
                porenum[j]=0;
                    for(k=0;k<Nsol;k++){
                    fgets(buf,180,fp);
                    sscanf(buf,"%lf %i %lf\n",&Time,&SolID,&dist);
                        if(dist <= poreradius) {
                            Data[j][k]=1;
                            porenum[j]+=1;
                        }
                        
                    }

            }
            
            //Use data:Calculate Pore Water Correlation Function//
            double Corr,sum,clock,norm;
            for(tau=0;tau<=maxtau;tau++){
                AC[tau]=0.0;
                    for(k=0;k<Nsol;k++){
                        Corr=0.0;
                        sum=0.0;
                        clock=0.0;
                          for(j=0;j<maxtau-tau;j++){
                              
                                //Pore residence clock//
                                if(Data[j][k]==1.0){
                                    clock+=1.0; /*Clock Increment*/
                                    norm+=1.0; /*Calculates Total Number of Residencies */
                                }
                              
                                //Correlation//
                                if(clock > 0.0 && Data[j+tau][k]==1.0){
                                    sum+=1.0; /*Correlation*/
                                }
                              
                              
                                //Reset Clock//
                                if(clock > 0.0 && Data[j][k]==0.0){
                                    clock=0.0;
                                }
                              
                       
                          }
                        
                        sum/=norm;
                        
                        Corr += sum;
               
                        if(norm!=0.0) AC[tau]+=Corr;
                      
                    }
                
        
            }
         
            //Collect in Lifetime Array//
            for(tau=0;tau<maxtau;tau++){
                Life[tau] += AC[tau]/AC[0];
            }
            

        }
        fclose(fp);
        /*****************END OF CORRELATION LOOP******************/

        /*Print Out Results*/
        for(i=0;i<maxtau;i++){
            fprintf(OUT,"%lf %lf\n",i*timestep,Life[i]/Ndiv);
            printf("%lf %lf\n",i*timestep,Life[i]/Ndiv); /*For Test Purposes*/
        }
        free(Life);
        free(AC);
        fclose(OUT);
    

    }
    
	
}





