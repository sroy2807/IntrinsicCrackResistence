#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mt19937ar.h"   // suporting file to generate random number
#include "mt19937ar.c"   // suporting file to generate random number

#define Ns 100000    // size of the bundle
#define filename "Critical_Stress_Beta0.1_gamma4.0.dat"   // name of the output file
#define Number_Max Ns    // maximum number upto which the model evolves
#define gamma 4.0    // stress release range


FILE *fp;
double stress[Ns],Dif[Ns],R1,D;
double stresstmp[Ns],X,p;
double threshold[Ns],Min,Deviding_Factor,S,delta,beta,alpha,gamma1;
int broken[Ns],broken_prev[Ns],BF,NNR,NNL,D_NNR,D_NNL,D_TOT,P1,crack;
int nnl[Ns], nnr[Ns],Pmax,NN,R;
int Patch[Number_Max],P[Ns],L[Ns],Done[Ns],Q[Number_Max],Patch_Prev[Number_Max];
double Patch_Tot[Number_Max],External_Stress[Number_Max],Stress_Tot[Number_Max],Random_Number;

double threshold_distribution_powerlaw(double beta){    // threshold distribution (Eq.1 of the manuscript)
	R1=(2*genrand_real2()-1);
	return pow(10,(R1*beta));
	}

double break_minimum_stress(int *index){    // finiding the weakest link
    	int i=0;
    	double min=1000;
    	for(i=0;i<Ns;i++){Dif[i]=threshold[i]-stress[i];if((Dif[i]<min) && (broken[i]!=1)){min=Dif[i];*index=i;}}
    	return(Dif[*index]);
	}

void take_input(){   // random number generator
    	int stime;
    	long ltime;
    	ltime=time(NULL);
    	stime=(unsigned)ltime/2;
    	long seedval=stime;
    	init_genrand(seedval);
	}

int main(){
	
    	int i,j,T_Max,Max=1,T_Fail,DeltaTau,x;
	double Tmax_Tot=0,Tfail_Tot=0,DeltaTau_Tot=0,Y[Ns],C,Sum,C_ini,Sum_ini;
    	double stress_increment,absolute_stress,absolute_stress_tot,PF[Ns],Nsurvivors_tot,Number_Tot; 
    	int breaking_fiber,sample,flag;
    	int Nsurvivors,Nsurvivors_prev=0,N1,Nsurvivors_prev_stress,Number,Number1;
    	take_input();
	beta=0.1;   // disorder strength 
	Pmax=(3*((int)(pow(Ns,0.6666))));   // span of stress release range
	Sum=0;
	for(x=1;x<=Pmax;x++){Sum+=pow(x,-gamma);}    // normalization (eq.5 in the manuscript)
	C=2*Sum;   
	S=(1.0/(2*1.0*Pmax));
	for(crack=0;crack<10000;crack+=2){    // loop on initial crack length
		p=(double)(crack/2);
		Number=0;
		absolute_stress_tot=0;
		Nsurvivors_tot=0;
		for(sample=0;sample<10000;sample++){    // configuarations
			T_Fail=0;
			T_Max=0;
			for(i=0;i<Ns;i++){Done[i]=0;}
			Number=0;
			for(i=0;i<Ns;i++){
		 	       nnl[i]=i-1;
		        	nnr[i]=i+1;
		    		}
			nnl[0]=Ns-1;
			nnr[Ns-1]=0;
			R=(double)(Ns/2);
			for(i=0;i<Ns;i++){     // initialization
				threshold[i]=threshold_distribution_powerlaw(beta);
                		broken[i]=0;
				stress[i]=0;
				Y[i]=1.0;
            			}
			if(crack>0){        // initialization
				for(i=R-(double)(crack/2);i<R+(double)(crack/2);i++){
					threshold[i]=0.0;
					broken[i]=1;
					Done[i]=1;
					stress[i]=0;
					Y[i]=0.0;
					}
				}
			Sum_ini=0;
			for(x=1;x<=(Ns-crack);x++){Sum_ini+=pow(x,-gamma);}    // initial stress redistribution (eq.6 in the manuscript)
			C_ini=2*Sum_ini;
			flag=0; 
			for(i<R-(double)(crack/2);i>0;i--){
				flag++;
				stress[i]+=(stress[i]*(pow(flag,gamma)))/(C_ini);
				}
			flag=0; 
			for(i<R+(double)(crack/2);i<Ns;i++){
				flag++;
				stress[i]+=(stress[i]*(pow(flag,gamma)))/(C_ini);
				}
			absolute_stress=0;
			Nsurvivors=Ns;
			while(Nsurvivors>crack){     // stress increment
				Number1=Number;
				Number++;
				stress_increment=0.001;
				absolute_stress+=stress_increment;
				N1=Nsurvivors;   
				for(i=0;i<Ns;i++){stress[i]+=(stress_increment*Y[i]);}
				for(i=0;i<Ns;i++){
					if((threshold[i]<=stress[i]) && (broken[i]!=1)){
						broken[i]=1;
						Nsurvivors--;
						if(Nsurvivors==0){break;}
						nnl[nnr[i]]=nnl[i];
						nnr[nnl[i]]=nnr[i];
						}
					}
				while(Nsurvivors_prev!=Nsurvivors && Nsurvivors>crack){     // redistributions
					Nsurvivors_prev=Nsurvivors;
					for(i=0;i<Ns;i++){
						if(broken[i]==1 && Done[i]!=1){
							j=i;
							for(P1=1;P1<=Pmax;P1++){
								while(broken[j]==1){j=nnl[j];}
								stress[j]+=(stress[i]*(pow(P1,gamma)))/(C);
								j=nnl[j];
								}
							j=i;
							for(P1=1;P1<=Pmax;P1++){
								while(broken[j]==1){j=nnr[j];}
								stress[j]+=(stress[i]*(pow(P1,gamma)))/(C);
								j=nnr[j];
								}
							Done[i]=1;
							}
						}
					for(i=0;i<Ns;i++){
                        			if((threshold[i]<=stress[i]) && (broken[i]!=1)){
                            				broken[i]=1;
                            				Nsurvivors--;
                            				if(Nsurvivors==0){break;}
							nnl[nnr[i]]=nnl[i];
							nnr[nnl[i]]=nnr[i];
							}
						}
					}
				for(i=0;i<Ns;i++){Y[i]=stress[i]/absolute_stress;}
				}
			absolute_stress_tot=absolute_stress_tot+absolute_stress;
			Nsurvivors_tot=Nsurvivors_tot+N1;
			Number_Tot+=Number;
			}
		absolute_stress_tot=(double)(absolute_stress_tot/(sample));
		Nsurvivors_tot=(double)Nsurvivors_tot/(sample*Ns);
		Number_Tot=(double)Number_Tot/sample;
		fp=fopen(filename,"a");
      		fprintf(fp,"%d\t%f\n",crack,absolute_stress_tot);
       	fclose(fp); 
		}   
    	return 0;
	}
