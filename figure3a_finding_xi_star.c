#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mt19937ar.h"   // suporting file to generate random number
#include "mt19937ar.c"   // suporting file to generate random number

#define Ns 100000    // size of the bundle
#define filename "Xi_Star_L100000.dat"   // name of the output file
#define Number_Max Ns    // maximum number upto which the model evolves


FILE *fp;
double stress[Ns],Dif[Ns],R1;
double stresstmp[Ns],X,p;
double threshold[Ns],Min,Deviding_Factor,S,delta,beta;
int broken[Ns],broken_prev[Ns],BF,NNR,NNL,D_NNR,D_NNL,D_TOT,P1,crack;
int nnl[Ns], nnr[Ns],Pmax,NN,R;
int Broken[Number_Max],Patch[Number_Max],P[Ns],L[Ns],Done[Ns],Q[Number_Max],Patch_Prev[Number_Max];
double Broken_Tot[Number_Max],Patch_Tot[Number_Max],External_Stress[Number_Max],Stress_Tot[Number_Max],Random_Number;

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
	
    	int i,j,T_Max,T_Fail,DeltaTau;
	double Tmax_Tot=0,Tfail_Tot=0,DeltaTau_Tot=0,Y[Ns],Max,Max_tot;
    	double stress_increment,absolute_stress,absolute_stress_tot,PF[Ns],Nsurvivors_tot,Number_Tot; 
    	int breaking_fiber,sample,flag;
    	int Nsurvivors,Nsurvivors_prev=0,N1,Nsurvivors_prev_stress,Number,Number1;
    	take_input();
	for(beta=0.1;beta<2.0;beta+=0.05){    // loop for strength of disorder
		Pmax=1;
		S=(1.0/(2*1.0*Pmax));
		for(crack=0;crack<Ns;crack+=2){  // loop for intial crack-length
			p=(double)(crack/2);
			for(Number=0;Number<Number_Max;Number++){Q[Number]=0;}
			for(Number=0;Number<Number_Max;Number++){
				Patch_Tot[Number]=0;
				Broken_Tot[Number]=0;
				}
			Max_tot=0.0;
			for(sample=0;sample<10000;sample++){  // configuarations
				Max=0.0;
				for(Number=0;Number<Number_Max;Number++){Patch[Number]=0;}
				for(Number=0;Number<Number_Max;Number++){Broken[Number]=0;}
				Number=0;
				T_Fail=0;
				T_Max=0;
				for(i=0;i<Ns;i++){Done[i]=0;}
				Number=0;
				for(i=0;i<Ns;i++){    // neighbour list
					nnl[i]=i-1;
		        		nnr[i]=i+1;
		    			}
				nnl[0]=Ns-1;
				nnr[Ns-1]=0;
				R=(double)(Ns/2);
            			for(i=0;i<Ns;i++){   // initialization
					threshold[i]=threshold_distribution_powerlaw(beta);
                			broken[i]=0;
					stress[i]=0;
					Y[i]=1.0;
            				}
				if(crack>0){    // initialization 
					for(i=R-(double)(crack/2);i<R+(double)(crack/2);i++){
						threshold[i]=0.0;
						broken[i]=1;
						Done[i]=1;
						stress[i]=0;
						Y[i]=0.0;
						}
					}
				for(i=0;i<Ns;i++){if(broken[nnr[i]]==1 && broken[nnl[i]]==0){Y[i]+=p;}}
				for(i=0;i<Ns;i++){if(broken[nnr[i]]==0 && broken[nnl[i]]==1){Y[i]+=p;}}
				absolute_stress=0;
				Nsurvivors=Ns;
				while(Nsurvivors>crack){     // stress increment
					Number1=Number;
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
						Number++;
						Q[Number]++;
						Nsurvivors_prev=Nsurvivors;
						for(i=0;i<Ns;i++){
							if(broken[i]==1 && Done[i]!=1){
								BF=i;
								j=i;
								for(P1=0;P1<Pmax;P1++){
									while(broken[j]==1){j=nnl[j];}
									}
								NNL=j;
								D_NNL=abs(BF-NNL);
								j=i;
								for(P1=0;P1<Pmax;P1++){
									while(broken[j]==1){j=nnr[j];}
									}
								NNR=j;
								D_NNR=abs(NNR-BF);
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
						Broken[Number]=Ns-Nsurvivors;
						for(i=0;i<Ns;i++){if(broken[i]==0 && broken[i+1]==1){Patch[Number]++;}}
						}
					for(i=0;i<Ns;i++){Y[i]=stress[i]/absolute_stress;}
					if(Patch[Number]>Max){Max=Patch[Number];}
					}
				Max_tot+=Max;
				}
			Max_tot=((1.0*Max_tot)/(1.0*sample));
			if(Max_tot==1){break;}
			} 
		fp=fopen(filename,"a");
      		fprintf(fp,"%f\t%d\n",beta,crack);
       	fclose(fp); 
		}  
    	return 0;
	}
