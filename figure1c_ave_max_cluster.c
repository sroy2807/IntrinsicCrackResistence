#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mt19937ar.h"   // suporting file to generate random number
#include "mt19937ar.c"   // suporting file to generate random number

#define Ns 100000    // size of the bundle
#define filename "Cluster_Ave_Max_beta3.0_L100000.dat"   // name of the output file
#define Number_Max Ns    // maximum number upto which the model evolves

FILE *fp;
double stress[Ns],Dif[Ns],R1;
double stresstmp[Ns],X;
double threshold[Ns],Min,Deviding_Factor,S,delta,beta,p;
int broken[Ns],broken_prev[Ns],BF,NNR,NNL,D_NNR,D_NNL,D_TOT,P1,crack,Count1;
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

void take_input(){    // random number generator
    	int stime;
    	long ltime;
    	ltime=time(NULL);
    	stime=(unsigned)ltime/2;
    	long seedval=stime;
    	init_genrand(seedval);
	}

int main(){
	
    	int i,j,T_Max,T_Fail,DeltaTau,l0,i_last,patch;
	double Tmax_Tot=0,Tfail_Tot=0,DeltaTau_Tot=0,Y[Ns],Prob,Ave_Patch,Max_Patch,Ratio_Tot,Ratio,Max_Crack,Ini_Crack,Max,Max_Tot;
    	double stress_increment,absolute_stress,absolute_stress_tot,P[Ns],Nsurvivors_tot,Number_Tot,Crack_Ave
    	double Crack_Ave1,Crack_Ave_Tot,Crack_Ave1_Tot; 
    	int breaking_fiber,sample,flag,Count;
    	int Nsurvivors,Nsurvivors_prev=0,N1,Nsurvivors_prev_stress,Number,Number1;
    	take_input();
	Pmax=1;             // stress release range for LLS scheme 
	beta=0.7;           // strength of disorder
	S=(1.0/(2*1.0*Pmax));
	for(crack=1;crack<50;crack++){    // continuously varying intial length of the crack
		p=(0.5*crack);
		absolute_stress_tot=0;
		Nsurvivors_tot=0;
		for(Count=0;Count<Ns;Count++){P[Count]=0;}
		Crack_Ave_Tot=0.0;
		Max_Tot=0.0;
		for(sample=0;sample<10000;sample++){   // configuarations
			for(i=0;i<Ns;i++){Done[i]=0;}
			for(i=0;i<Ns;i++){    //neighbour list
				nnl[i]=i-1;
				nnr[i]=i+1;
			   	}
			nnl[0]=Ns-1;
			nnr[Ns-1]=0;
			R=(double)(Ns/2);
			if(crack==0){
        	    		for(i=0;i<Ns;i++){    // initialization
					threshold[i]=threshold_distribution_powerlaw(beta);
        	        		broken[i]=0;
					stress[i]=0;
					Y[i]=1.0;
        	    			}
				}
			else if(crack>0){    // initialization
				for(i=R;i<R+crack;i++){
					threshold[i]=0.0;
					broken[i]=1;
					Done[i]=1;
					stress[i]=0;
					Y[i]=0.0;
					}
				for(i=0;i<R;i++){
					threshold[i]=threshold_distribution_powerlaw(beta);
					broken[i]=0;
					stress[i]=0;
					Y[i]=1.0;
					}
				for(i=R+crack;i<Ns;i++){
					threshold[i]=threshold_distribution_powerlaw(beta);
					broken[i]=0;
					stress[i]=0;
					Y[i]=1.0;
					}
				}
			Y[R+crack]+=p;
			Y[nnl[R]]+=p;
			nnr[nnl[R]]=R+crack;
			nnl[R+crack]=nnl[R];
			absolute_stress=0;
			Nsurvivors=Ns;
			while(Nsurvivors>crack){  // stress increment
				for(i=0;i<Ns;i++){broken_prev[i]=broken[i];}
				stress_increment=0.0001;
				absolute_stress+=stress_increment;
				N1=Nsurvivors;   
				for(i=0;i<Ns;i++){stress[i]+=(stress_increment*Y[i]);}
				for(i=0;i<Ns;i++){
					if((threshold[i]<=stress[i]) && (broken[i]!=1)){
						broken[i]=1;
						i_last=i;
						Nsurvivors--;
						if(Nsurvivors==0){break;}
						nnl[nnr[i]]=nnl[i];
						nnr[nnl[i]]=nnr[i];
						}
					}
				while(Nsurvivors_prev!=Nsurvivors && Nsurvivors>crack){    // redistributions
					Nsurvivors_prev=Nsurvivors;
					for(i=0;i<Ns;i++){
						if(broken[i]==1 && Done[i]!=1){
							j=i;
							for(P1=0;P1<Pmax;P1++){while(broken[j]==1){j=nnl[j];}}
							NNL=j;
							if(NNL<BF){D_NNL=abs(BF-NNL);}
							else if(NNL>BF){D_NNL=(Ns-abs(BF-NNL));}
							j=i;
							for(P1=0;P1<Pmax;P1++){while(broken[j]==1){j=nnr[j];}}
							NNR=j;
							if(NNR>BF){D_NNR=abs(NNR-BF);}
							else if(NNR<BF){D_NNR=(Ns-abs(NNR-BF));}
							D_TOT=D_NNL+D_NNR;
							stress[NNL]+=(D_NNR/(1.0*D_TOT))*stress[i];
							stress[NNR]+=(D_NNL/(1.0*D_TOT))*stress[i];
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
			for(i=0;i<Ns;i++){Done[i]=0;}
			Crack_Ave=0.0;
			Max=0.0;
			patch=0;		
			for(i=0;i<Ns;i++){
				if(broken_prev[i]==1 && Done[i]==0){
					patch++;
					Done[i]=1;
					Count=1;
					for(j=i+1;j<Ns;j++){
						if(broken_prev[j]==1 && Done[j]==0){
							Done[j]=1;
							Count++;
							}
						else if(broken_prev[j]==0 || Done[j]==1){break;}
						}
					if(Count>Max){Max=Count;}
					Crack_Ave+=Count;
					}
				}
			i=R;
			Count1=1;
			for(j=i+1;j<Ns;j++){
				if(broken_prev[j]==1){Count1++;}
				else if(broken_prev[j]==0){break;}
				}
			i=R;
			for(j=i-1;j>=0;j--){
				if(broken_prev[j]==1){Count1++;}
				else if(broken_prev[j]==0){break;}
				}
			Crack_Ave1=Crack_Ave-Count1;
			Crack_Ave=((1.0*Crack_Ave)/(1.0*patch));
			if(patch>1){Crack_Ave1=((1.0*Crack_Ave1)/(1.0*(patch-1)));}
			else if(patch==1){Crack_Ave1=0;}
			Max_Tot+=Max;
			Crack_Ave_Tot+=Crack_Ave;
			Crack_Ave1_Tot+=Crack_Ave1;
			}
		Crack_Ave_Tot=((1.0*Crack_Ave_Tot)/(1.0*sample));
		Crack_Ave1_Tot=((1.0*Crack_Ave1_Tot)/(1.0*sample));
		Max_Tot=((1.0*Max_Tot)/(1.0*sample));
		fp=fopen(filename,"a");
		fprintf(fp,"%d\t%f\t%f\n",crack,Crack_Ave_Tot,Max_Tot);
		fclose(fp); 
		}
    	return 0;
	}
