/***************** Pricer in a recombining tree for CDSwaptions *************************/


/* ASSUMPTIONS: 
   1) Proportional model assume elasticity parameter \gamma = 1 
   2) Default intensity is defined \lambda(t) = \bar{r}(t)-r(t)
   3) Initial term structure of the forward credit spread is assumed to be flat: s(0,t) = lambda_0
   4) The volatility and the mean reversion term represented, respectively, by: vol and k_s are constants. 
   
    INPUT: Model parameters
	
	double vol;  volatitily (sigma)
	double k_s;  exogenous factor or mean reversion term
	
	
	int m;  number of elements in the partition of phi
	int N;  number of steps in the tree
	double strike; strike price of the CDSwaption
	double R; recovery rate
	
	string dateT; maturity date of the CDSwaption
	
	char typeoption; call 'C' or put option 'P', respectively, payer or receiver CDSwaption
	
	double r0; initial interest rate assumed to be constant over the time
	double lambda0; initial default intensity
	double phi0; initial accumulated variance for the forward credit spread
	
    
*/

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include<fstream>
#include<algorithm>
#include<assert.h>

using namespace std;


/* Node Declaration
 */
struct Nodes 
{
	int visit; // 0 - if the node wasn't visit , 1 - if the node was visit
    double lambda; // default intensity 
	double phi[2];// Array with two elements containing the minimum and maximum of phi, 
    double *optionPrice; // pointer to an array with the option price
    double *probUp; // pointer to an array with the up probabilities
};


int signal (double Z);

double calcJ(double mean, double sqrtt);

void buildLRStree( double vol,double k_s,int m, int N, double strike, char typeoption, double r0, double lambda0, 
                  double phi0,double R, const string& dateT, const vector<string>& dates, const vector<int>& numdate);


/*
 *******************Main Contains Menu***************************************************************
 */
 
int main()
{
	// declare and initialize the input parameters
    double vol=0.58, k_s=-0.1, strike=0.0625, r0=0.01, phi0=0.0, lambda0=0.06, R=0.4;
    
	string dateT= "20-12-2008";
	
    // number of partitions of phi and tree steps
    int m=25, N =64;
    
    // declare and initialize type of the option (call C ou put P)
    char typeoption='C';
    
    // extract the values of each column of the text file
    
    vector<string> dates;
    vector<int> daycount;
        
    ifstream theFile("dates.txt");
    
    string reading1;
	int reading2;  
        
    
    while(theFile>> reading1>> reading2){
	   	dates.push_back(reading1);
    	daycount.push_back(reading2);	
	}
   
    buildLRStree( vol, k_s, m, N, strike, typeoption,  r0, lambda0, phi0, R, dateT, dates, daycount);
    
    return 0;
     
}

void buildLRStree( double vol,double k_s, int m, int N, double strike, char typeoption, double r0, 
                   double lambda0, double phi0, double R,const string& dateT, const vector<string>& dates, 
				   const vector<int>& numdate){
	
		
	int posT = -1; // position on the txt file of the maturity date of the CDSwaption
	 
    for(unsigned int i=0; i < dates.size();i++){
    	if(dates[i]==dateT)
	 	posT=i;
    }
    
    cout << "Position of the expiry date of CDSoption: " << posT <<endl;
    assert(posT>-1); // verify if in fact exists the maturity date on the file ---> debug
    
    /*Local variables
	dt = time interval
	sqrt_dt = square root of the time interval
	*/
	
	double dt=0.0;
	double sqrt_dt=0.0;
	
	dt=((numdate[posT]-numdate[0])/365.0)/N; //time steps
	cout<<"tau "<<(numdate[posT]-numdate[0])/365.0 <<endl;
	cout<< "dt " << dt << endl;
	
	sqrt_dt= sqrt(dt);
	
    //Dynamically allocating memory in each node
	Nodes **tree = new Nodes *[N+1];
	
	for(int i=0; i<N+1 ; i++) // dimension of each line
	tree[i] = new Nodes [i+1];
			
	//Inicializing the root of the tree
	tree[0][0].lambda = lambda0;
    tree[0][0].phi[0]= phi0;
    tree[0][0].phi[1]= phi0;
    
    for (int i = 0; i < N+1 ; i++){
    	for (int j = 0; j <= i; j++){
    		tree[i][j].visit=0;
	 		tree[i][j].optionPrice= new double [m];
	 		tree[i][j].probUp = new double [m];
        }
    }
 
    
 /*** Calculating the values of lambda, phi and probUp for each node of the tree, except the probUp on the terminal node (i=N) ***/
    
for(int i=0; i<N ;i++){ 
    	for(int j=0; j<i+1 ;j++){
    		 
    		if(tree[i][j].phi[0]==tree[i][j].phi[1]){ //edges of the lattice
    			
				double mean1=0.0;
				int J1=0;   			
    			
				mean1= (k_s*(lambda0-tree[i][j].lambda) + tree[i][j].phi[0])/(vol * tree[i][j].lambda) - (vol/2);
			    J1=calcJ(mean1, sqrt_dt);
		        
				if(J1!=0){
					cout<< i << " " << j << " " << "J1 change is value for " << J1 <<endl;
					exit(0);
				}
				     
				for(int k=0; k<m ;k++)
	            tree[i][j].probUp[k]=(mean1*dt + (1-J1)*sqrt_dt)/(2*sqrt_dt); 
	        
	            // Insert the next values of lambda, if until the moment the next nodes aren't visit 
	            
    		    if(tree[i+1][j].visit==0)
			    tree[i+1][j].lambda= exp(log(tree[i][j].lambda) + vol*(J1+1) * sqrt_dt); 
			
			    if(tree[i+1][j+1].visit==0)	
			    tree[i+1][j+1].lambda= exp(log(tree[i][j].lambda) + vol*(J1-1) * sqrt_dt); 
			 		
               // Calculate the next value of phi
                double phi_next=0.0;
    		    phi_next= tree[i][j].phi[0] + (pow(vol,2) * pow(tree[i][j].lambda,2)- 2*k_s*tree[i][j].phi[0])*dt;
    		
    		    if(tree[i+1][j].visit==0){
    		    	tree[i+1][j].phi[0]= phi_next;
    		    	tree[i+1][j].phi[1]= phi_next;
    		    	tree[i+1][j].visit=1;
		    	}else{
		    		tree[i+1][j].phi[0] = min(phi_next, tree[i+1][j].phi[0] );
		    		tree[i+1][j].phi[1] = max(phi_next,tree[i+1][j].phi[1]);
			    }
			
			    if(tree[i+1][j+1].visit==0){
			    	tree[i+1][j+1].phi[0]= phi_next;
			    	tree[i+1][j+1].phi[1]= phi_next;
			    	tree[i+1][j+1].visit=1;
		    	}else{
		    		tree[i+1][j+1].phi[0]= min(phi_next, tree[i+1][j+1].phi[0] );
		    		tree[i+1][j+1].phi[1]= max(phi_next,tree[i+1][j+1].phi[1]);
	    		}

    		  			  	  
			} else{ //middle nodes
				
				double step=0.0; // step of the partition of phi
				step= (tree[i][j].phi[1]-tree[i][j].phi[0])/(m-1);
    		    
   		        vector<double> partPhi; // vector with m-partition of phi 
				partPhi.push_back(tree[i][j].phi[0]);
				  
    		    for(int k = 1; k < m-1 ; k++)
    		    partPhi.push_back(partPhi[k-1] + step);
    		    partPhi.push_back(tree[i][j].phi[1]);
    		                	
    		    vector<double> mean; // Now, the mean has m elements because we hve m values of phi 
    		    for(int k=0; k < m ; k++)
    		    mean.push_back((k_s * (lambda0 - tree[i][j].lambda) + partPhi[k] ) / (vol*tree[i][j].lambda) - (vol/2));
			    
			    vector<int> J;
			    int valueJ=0;
			    int check=0;
			    
			    for( int k=0; k<m ; k++){
			    	J.push_back(calcJ(mean[k], sqrt_dt));
			    	
					if(J[k]!=J[0]){
			    		check=1;
                        cout<< i << " " << j << " " << "Vector J change is value in" << k << " for " << J[k]<<endl;
					}		
			   	}
			   	
		        if(check==0)
				valueJ=J[0];
				else if(check==1)
				exit(0);
				
			   	    
    		    for(int k=0; k<m ;k++)
	            tree[i][j].probUp[k]=(mean[k]*dt + (1-J[k])*sqrt_dt)/(2*sqrt_dt); 
	        
                // Insert the next values of lambda, if until the moment the next nodes aren't visit 
    		    if( tree[i+1][j].visit==0)
    		    tree[i+1][j].lambda= exp(log(tree[i][j].lambda) + vol*(valueJ+1) * sqrt_dt); 
			
			    if(tree[i+1][j+1].visit==0)	
			    tree[i+1][j+1].lambda= exp(log(tree[i][j].lambda) + vol*(valueJ-1)* sqrt_dt); 
				
				vector <double> sucPhi; // vector with the successor values of partPhi (vector with dim=m)
    		    for(int k=0;k<m;k++)
   		    	sucPhi.push_back(partPhi[k]+ (pow(vol,2) * pow(tree[i][j].lambda,2)- 2 * k_s * partPhi[k] )*dt);
    		    	
			    /*Update the successor values of phi */
    		    if( tree[i+1][j].visit==1){ // By construction, we already visit this node 
    		    	
    		    	double currentmin=tree[i+1][j].phi[0];
			        double currentmax=tree[i+1][j].phi[1];
			            
			        for(int k=0;k<m;k++){
			        	if(sucPhi[k]<currentmin) //minimum value of sucPhi
			        	currentmin=sucPhi[k];
			        	
                        if(sucPhi[k]>currentmax) //maximum value of sucPhi
		    	        currentmax=sucPhi[k];
				    }
				   		    	
    		    	tree[i+1][j].phi[0]= currentmin;
			    	tree[i+1][j].phi[1]= currentmax;
			    	
			    }else{
			    	cout<< "You are doing the lattice in a wrong way-update(i+1,j)"<< endl;
			    }
			    
			    double minsucPhi=sucPhi[0]; //first value of sucPhi
 	            double maxsucPhi=sucPhi[m-1]; //last value of sucPhi
 	            
				 /*Calculation of the minimum and maximum values of the vector sucPhi*/
 	            
         		for(int k=0;k<m;k++){
					if(sucPhi[k]<minsucPhi) //minimum value of sucPhi
    		    	minsucPhi=sucPhi[k];
    		    	if(sucPhi[k]>maxsucPhi) //maximum value of sucPhi
    		    	maxsucPhi=sucPhi[k];
				}
			
			    if(tree[i+1][j+1].visit==0){ // By construction, we have't already visit this node
			    tree[i+1][j+1].phi[0]= minsucPhi;
			    tree[i+1][j+1].phi[1]= maxsucPhi;
			    tree[i+1][j+1].visit=1;
			    }else{
		    	cout<< "You are doing the lattice in a wrong way-update(i+1,j+1)"<< endl;
         	    }
				
			}
    		
        } //close j
	
	} //close i
	
	
/****** Terminal nodes *****/

/*Remember that posT  is the position of the maturity date T_k of CDSwaption on the textfile 
which coincides with the starting time of the forward CDS contract
*/	
    vector<double> delta; // interval of time between payments posT=T_k,...,T_n
		
    for(unsigned int k=0; k < numdate.size()-1;k++){
   	delta.push_back((numdate[k+1]-numdate[k])/365.0);
    
	}
	cout << "delta. size: " << delta.size() <<endl;
    	
    vector<double> survival; //survival factor S(0,T_i), for i=k,..,n-1,n
		
    for(unsigned int k=0; k < numdate.size();k++){
    survival.push_back(exp(-lambda0*((numdate[k]-numdate[0])/365.0)));
    //cout<< "survival0 " << survival[k] <<endl;
	}
    
    //Calculate the forward default swap rate at the valuation date:
    vector<double> bondPrice0; //default-free discount factor P(0,T_i+1), for i=k,..,n-1
		
    for(unsigned int k=0; k < numdate.size();k++){
      bondPrice0.push_back(exp(-r0*((numdate[k]-numdate[0])/365.0)));
    //cout<< "P0 " << bondPrice0[k]<<endl;
	}
                           	     	
    vector<double> defaultBondPrice0; // default zero coupon bond bar(P)(posT,T_i+1), i=k,..,n-1
    
	for( unsigned int k=0;k < numdate.size();k++){
		defaultBondPrice0.push_back(survival[k]*bondPrice0[k]);	
	//	cout<< "barP0: " << defaultBondPrice0[k]<<endl;
	}
	
	double num = 0.0;
	double den = 0.0;
  	
	for(unsigned int k=posT; k < numdate.size()-1;k++)
	num+= survival[k]*bondPrice0[k+1]-defaultBondPrice0[k+1];
	int count=posT+1;
	
	defaultBondPrice0.erase(defaultBondPrice0.begin(), defaultBondPrice0.begin()+count);
	
	int countDelta=posT;
	delta.erase(delta.begin(),delta.begin()+countDelta);
	for(unsigned int k=0; k < delta.size();k++)
	den+= delta[k]*defaultBondPrice0[k];
	//cout<<"discount0: "<<den<<endl;
	
    double fswaprate0=0.0;
	fswaprate0=(1-R)*(num/den);
    //cout<< "fswaprate0: " << fswaprate0 <<endl;
    
    vector<double> bondPrice; //default-free discount factor P(posT,T_i+1), for i=k,..,n-1
	    
	for(unsigned int k=0; k < numdate.size()-1;k++){
    	bondPrice.push_back(exp(-r0*((numdate[k+1]-numdate[posT])/365.0)));
    //	cout<<"P(t,T) "<<bondPrice[k]<<endl;
	}
	
	int count1=0;
	for( unsigned int k=0;k < bondPrice.size();k++){
		if(bondPrice[k]>=1)
		count1+=1;
	}
	bondPrice.erase(bondPrice.begin(),bondPrice.begin()+count1);
	cout << "P. size1: " << bondPrice.size() <<endl;
	
				
    
			 
/*
***---------Loop on the terminal nodes-------***

Goal: Calculate the values of option Price in each terminal node
*/
         
	for(int j=0; j < N+1 ;j++){
				
		if(tree[N][j].phi[0]==tree[N][j].phi[1]){ //edges of the lattice
		
	        vector<double> defaultRisk; // S(t,T) forward survival probability
	                          	     	
		    for(unsigned int k=0; k < numdate.size();k++){
	    	defaultRisk.push_back((survival[k]/survival[posT])*exp( ((1-exp(- k_s*((numdate[k]-numdate[posT])/365.0))) /k_s)* 
            (lambda0-tree[N][j].lambda)- 0.5*pow((1-exp(-k_s*((numdate[k]-numdate[posT])/365.0))),2)/pow(k_s,2)* tree[N][j].phi[0] ));      
			}
		    
	    	int count3=0;
						
			for( unsigned int k=0;k < defaultRisk.size();k++){
				if(defaultRisk[k]>1)
				count3+=1;
			}
			
			defaultRisk.erase(defaultRisk.begin(),defaultRisk.begin()+count3);

			for(unsigned int k=0; k < numdate.size();k++){     
	    	 cout<<"surv(t,T): "<<defaultRisk[k]<<endl;
			}
		    
			vector<double> defaultBondPrice; // default zero coupon bond bar(P)(posT,T_i+1), i=k,..,n-1
			
			for( unsigned int k=0;k < bondPrice.size();k++){
				defaultBondPrice.push_back(defaultRisk[k+1]*bondPrice[k]);
				cout<<"barP(t,T): "<<defaultBondPrice[k]<<endl;	
			}
//			cout << "P. size: " << bondPrice.size() <<endl;
//			cout << "DB. size: " << defaultBondPrice.size() <<endl;
//			cout << "Surv. size: " << defaultRisk.size() <<endl;
			
			
			/* Calculation of the forward swap rate \xi_{k,n} (fswaprate):
			
			auxfswap = numerator of the forward swap rate formula
			discount= denominator of the forward swap rate formula and also 
			          the discount factor of the payoff at T_k of the CDSoption
			          
			fswaprate = formula of the forward swap rate on remark 5.2
			*/
			
			double auxfswap = 0.0;
			double discount = 0.0;
			
			for(unsigned int k=0; k < bondPrice.size();k++){
				auxfswap+= defaultRisk[k]*bondPrice[k]-defaultBondPrice[k];
				discount+= delta[k]*defaultBondPrice[k];	
			}
		
			double fswaprate=0.0; 		
			fswaprate=(1-R)*(auxfswap/discount);
						
	        
			if(typeoption=='C'){ // if is a payer CDSwaption
				
				for(unsigned int k=0; k<m ;k++)
				tree[N][j].optionPrice[k]=discount* max( fswaprate - strike, 0.0);

	        }else if(typeoption=='P'){ // if is a receiver CDSwaption
	        	        	
				for(unsigned int k=0;k<m;k++)
				tree[N][j].optionPrice[k]=discount* max( strike - fswaprate , 0.0);
		    }
		 			
	    } else{ //middle nodes
	    	
			double step=0.0;		       		
            step=(tree[N][j].phi[1]-tree[N][j].phi[0])/(m-1);
            
            vector<double> partPhi;
			partPhi.push_back(tree[N][j].phi[0]);
				  
		    for(int k = 1; k < m-1 ; k++)
		    partPhi.push_back(partPhi[k-1] + step);
		    partPhi.push_back(tree[N][j].phi[1]);
         		          
            for(int k=0; k< m ;k++){ // for each value of the partition of phi
            	
				vector<double> defaultRisk;
				     	        
				for(int u=0 ; u < numdate.size();u++)
				defaultRisk.push_back((survival[u]/survival[posT])*exp( ((1-exp(- k_s*((numdate[u]-numdate[posT])/365.0))) /k_s)*
                (lambda0-tree[N][j].lambda)- 0.5*pow((1-exp(-k_s*((numdate[u]-numdate[posT])/365.0)))/k_s,2)* partPhi[k] ));          
                
                int count3=0;
				
				for( unsigned int u=0;u < numdate.size();u++){
					if(defaultRisk[u]>1)
					count3+=1;
				}
				
				defaultRisk.erase(defaultRisk.begin(),defaultRisk.begin()+count3);
			               	
			    vector<double> defaultBondPrice;
			    
				for(unsigned int u=0;u<bondPrice.size();u++)
	            defaultBondPrice.push_back(defaultRisk[u+1]*bondPrice[u]); 
// 	  
// 	            cout << "P. size: " << bondPrice.size() <<endl;
//			    cout << "DB. size: " << defaultBondPrice.size() <<endl;
//			    cout << "Surv. size: " << defaultRisk.size() <<endl;
			    
	            //Calculation of the forward swap rate:
	            double auxfswap = 0.0;
	            double discount = 0.0;
	           	
				for(int u=0; u< bondPrice.size() ;u++){
					auxfswap+= defaultRisk[u]*bondPrice[u]-defaultBondPrice[u];
					discount+=delta[u]*defaultBondPrice[u];
	            }
	            
	            double fswaprate=0.0;
	            fswaprate=(1-R)*(auxfswap/discount);
	            
	            
                if(typeoption=='C')
                tree[N][j].optionPrice[k]=discount* max( fswaprate - strike, 0.0);
				else if (typeoption=='P')
				tree[N][j].optionPrice[k]=discount* max( strike - fswaprate , 0.0);
				
		    }	
    	  }    			
	}//end of loop on the terminal nodes
	
 /*****------------Backward recursion----------------*****/
  
	for( int i=N-1; i>=0;i--){
	    for( int j=0; j<=i; j++)
       {
  	 	      
       	if(tree[i][j].phi[0]==tree[i][j].phi[1]){ //edge nodes of the lattice
       	
			// Calculate the next value of phi
            double phi_next=0.0;
            phi_next= tree[i][j].phi[0] + (pow(vol,2) * pow(tree[i][j].lambda,2)- 2*k_s*tree[i][j].phi[0])*dt;
    		
			if(phi_next==tree[i+1][j].phi[0] && phi_next==tree[i+1][j+1].phi[1]){
				for(int k=0;k<m;k++)
				tree[i][j].optionPrice[k]= (tree[i][j].probUp[k]*tree[i+1][j].optionPrice[0]+(1-tree[i][j].probUp[k])*
                                           tree[i+1][j+1].optionPrice[m-1]) *exp(-(tree[i][j].lambda+r0)*dt);    
			}else{
				cout<< " Doing in a wrong way the Linear interpolation: edge nodes of the lattice" <<endl;
				
			}		
      	   	   
		  
		} else{ //middle nodes
			
			/*At node (i,j) we will compare the values of next values of phi
			 with the values of phi in nodes (i+1,j) and (i+1,j+1)*/
						
        	double step=0.0;
        	step= (tree[i][j].phi[1]-tree[i][j].phi[0])/(m-1);
        	  	
            vector<double> partPhi;              
			partPhi.push_back(tree[i][j].phi[0]);
				  
		    for(int k = 1; k < m-1 ; k++)
		    partPhi.push_back(partPhi[k-1] + step);
		    partPhi.push_back(tree[i][j].phi[1]);
    		
    		vector <double> sucPhi; //next values of phi on node (i,j)
    		for(int k=0;k<m;k++)
			sucPhi.push_back(partPhi[k]+ (pow(vol,2) * pow(tree[i][j].lambda,2)- 2 * k_s * partPhi[k] )*dt);
					   	    		
    		double step2=0.0;          
			step2= (tree[i+1][j].phi[1]-tree[i+1][j].phi[0])/(m-1);
			
			vector<double> partPhi2; //values of the partition of phi in (i+1,j)
			partPhi2.push_back(tree[i+1][j].phi[0]);
				  
		    for(int k = 1; k < m-1 ; k++)
		    partPhi2.push_back(partPhi2[k-1] + step2);
		    partPhi2.push_back(tree[i+1][j].phi[1]);	
		    
		    //In theory the min and max pf partPhi2 are at positions 0 and m-1 --> just to check
		    double minpartPhi2=partPhi2[0];
		    double maxpartPhi2=partPhi2[m-1];
			int minpos2=0;
			int maxpos2=m-1;
			          
			for(int k=0;k<m;k++){
				if(partPhi2[k]<minpartPhi2){
					minpartPhi2=partPhi2[k];
					minpos2=k;
				} //minimum value of partPhi2
                
			        	
                if(partPhi2[k]>maxpartPhi2){
                	maxpartPhi2=partPhi2[k];
                	maxpos2=k;
				} //maximum value of partPhi2
                
		    }
		  		
			
    		double step3=0.0;
    		step3= (tree[i+1][j+1].phi[1]-tree[i+1][j+1].phi[0])/(m-1);

			vector<double> partPhi3; //values of the partition of phi in (i+1,j+1)
			partPhi3.push_back(tree[i+1][j+1].phi[0]);
			
			
    		for( int k = 1; k < m -1 ; k++)
    		partPhi3.push_back(partPhi3[k-1] + step3);
    		partPhi3.push_back(tree[i+1][j+1].phi[1]);
    				
    		//In theory, when J=0, the min and max of partPhi3 are at positions 0 and m-1 --> just to check
   		    double minpartPhi3=partPhi3[0];
		    double maxpartPhi3=partPhi3[m-1];
			int minpos3=0;
			int maxpos3=m-1;  
			    
			for(int k=0;k<m;k++){
				if(partPhi3[k]<minpartPhi3){
					minpartPhi3=partPhi3[k];
					minpos3=k;
				} //minimum value of partPhi3
                
			        	
                if(partPhi3[k]>maxpartPhi3){
                	maxpartPhi3=partPhi3[k];
                	maxpos3=k;
				} //maximum value of partPhi3
                
		    }
    				
    		/** Doing Backward recursion and Interpolation **/	
						
    		for(int k=0;k<m;k++){ 
			                   
				double gup=0.0;
                double gdown=0.0;
    			
    			if( (sucPhi[k] != minpartPhi2) && (sucPhi[k] != maxpartPhi3)){
    				
    				if(sucPhi[k] == partPhi2[k])
			        cout<<"error: sucphi(k)=partphi2(k)"<<endl;
			        if(sucPhi[k] == partPhi3[k])
			        cout<<"error: sucphi(k)=partphi3(k)"<<endl;
    				
    				int plus2 =0;
    				int minus2 =0;
    				for(int u=0; u<m; u++){
    					if((sucPhi[k] > partPhi2[u]) && (sucPhi[k] < partPhi2[u+1])){
    						minus2=u;
    						plus2=minus2+1;
   						}	
				    }
	    			
					int plus3 =0;
					int minus3 =0;
					for(int u=0; u<m; u++){
						if((sucPhi[k] > partPhi3[u]) && (sucPhi[k] < partPhi3[u+1])){
							minus3=u;
							plus3=minus3+1;		
					    }	
				    }
				    
			        gup=tree[i+1][j].optionPrice[minus2] + (sucPhi[k]-partPhi2[minus2])/(partPhi2[plus2]-partPhi2[minus2] ) *
				     (tree[i+1][j].optionPrice[plus2]-tree[i+1][j].optionPrice[minus2]);
			     
			        gdown=tree[i+1][j+1].optionPrice[minus3] + (sucPhi[k]-partPhi3[minus3])/(partPhi3[plus3]-partPhi3[minus3] ) * 
				     (tree[i+1][j+1].optionPrice[plus3]-tree[i+1][j+1].optionPrice[minus3]);
				     
					tree[i][j].optionPrice[k]=(tree[i][j].probUp[k] * gup + (1-tree[i][j].probUp[k])*gdown) * exp(-(tree[i][j].lambda+r0)*dt);
		        
			    }
			    
			          
		       if(sucPhi[k] == minpartPhi2 && sucPhi[k] != maxpartPhi3 && sucPhi[k] != minpartPhi3){ 
			   //coincide on the first midle node at i=N-1
			    if(sucPhi[k] == partPhi3[k])
			    cout<<"error: sucphi(k)=partphi3(k)"<<endl;
			    
					int plus3 =0;
					int minus3 =0;
					for(int u=0; u<m; u++){
						if((sucPhi[k] > partPhi3[u]) && (sucPhi[k] < partPhi3[u+1])){
							minus3=u;
							plus3=minus3+1;		
					    }	
				    }
		       	 
                gdown=tree[i+1][j+1].optionPrice[minus3] + (sucPhi[k]-partPhi3[minus3])/(partPhi3[plus3]-partPhi3[minus3] ) *
				 (tree[i+1][j+1].optionPrice[plus3]-tree[i+1][j+1].optionPrice[minus3]);

     			 tree[i][j].optionPrice[k]=(tree[i][j].probUp[k]*tree[i+1][j].optionPrice[minpos2]+
				  (1-tree[i][j].probUp[k])*gdown)*exp(-(tree[i][j].lambda+r0)*dt);

			    }
			    
			    if(sucPhi[k] == maxpartPhi3 && sucPhi[k] != minpartPhi2 && sucPhi[k] != maxpartPhi2){ 
				//coincide for the last edge node at i=N-1
				if(sucPhi[k] == partPhi2[k])
			    cout<<"error: sucphi(k)=partphi2(k)"<<endl;
			    
			    int plus2 =0;
     			int minus2 =0;
   				for(int u=0; u<m; u++){
   					if((sucPhi[k] > partPhi2[u]) && (sucPhi[k] < partPhi2[u+1])){
   						minus2=u;
   						plus2=minus2+1;
                    }	
			    }
			    			    
			    
			    gup=tree[i+1][j].optionPrice[minus2] + (sucPhi[k]-partPhi2[minus2])/(partPhi2[plus2]-partPhi2[minus2] ) * 
				(tree[i+1][j].optionPrice[plus2]-tree[i+1][j].optionPrice[minus2]);
				
     	        tree[i][j].optionPrice[k]=(tree[i][j].probUp[k]*gup + 
				 (1-tree[i][j].probUp[k])*tree[i+1][j+1].optionPrice[maxpos3] ) * exp(-(tree[i][j].lambda+r0)*dt);

			    }
			    
			    if( (sucPhi[k] == minpartPhi2) && (sucPhi[k] == maxpartPhi3)){
			    	cout<< "erro1"<<endl;
				}
				 if( sucPhi[k] == maxpartPhi2){
			    	cout<< "erro2"<<endl;
				}
				if( sucPhi[k] == minpartPhi3){
			    	cout<< "erro3"<<endl;
				}
			
	        } 		
          }
        }  	
  	}
  	
  	
  	for( int i=0 ; i<=N; i++){
  		for( int j=0;j<=i;j++){
//  			 cout<< i <<" "<< j << " "<<" dspread "<< tree[i][j].lambda<< endl;
//             cout <<i<<" "<< j << " "<<" phi0 "<<  tree[i][j].phi[0]<< endl;
//            cout <<i<<" "<< j << " "<<" phi1 "<< tree[i][j].phi[1]<< endl;
//           
  	       for(int k=0;k<m;k++){
    	    //   cout <<i<<" "<< j << " prob "<< tree[i][j].probUp[k]<< endl;
  	       	//cout <<i<<" "<< j << " "<< "optionPrice "<<tree[i][j].optionPrice[k]<< endl;
		  }
  	       
        }
    }  
    
       for(int k=0;k<m;k++){
       	//cout <<" "<<  " prob1 "<< tree[99][0].probUp[k]<< endl;
       	//cout <<" "<< " prob2 "<< tree[99][1].probUp[k]<< endl;
      	//cout<< " " << "probability " << tree[2][1].probUp[k];
      	//cout << "(2,1) "<< "optionPrice "<<tree[2][1].optionPrice[k]<< endl;
      }
      for(int k=0;k<m;k++ ){
      	//cout << "(3,2) "<< "optionPrice "<<tree[3][2].optionPrice[k]<< endl;
        cout << " "<< "optionPrice "<<tree[0][0].optionPrice[k]<< endl;
       	
	   }
       
       
  //return tree[0][0].optionPrice[0];
}

/* Auxilary Functions: Z and calcJ -----> for calculating J (LRS paper)*/
int signal ( double Z){
	
	int result;
	
	if(Z>0){
	result=1;
	} else if(Z<0)
	{
		result=-1;
	} else if( Z==0)
	{
		result=0;
	}
	
	return result;
	
	}
	


double calcJ(double mean, double sqrtt)
{
	int Z;
	double J;
	
	Z= floor(mean*sqrtt);
	
	if ( Z % 2 == 0)
	J=signal(Z)*Z;
	
	else
	J=signal(Z)*abs(Z) +1;
	
	
	return J; 
		
}


  
