/* INPUT:  Model parameters
	
	double vol; //volatitily
	double k_s; //exogenous factor
	(Proportional model assume elasticity parameter =1 )
	
	int m; //number of elements in the partition of phi
	int N; //number of steps in the tree
	double strike; // strike price of the CDSwaption
	double R; //recovery rate
	
	string dateT; maturity date of the CDSwaption
	
	char typeoption;// call 'C' or put option 'P' respectively payer or receiver cdswaption
	
	double r0; //initial interest rate assumed to be constant over the time
	double s0; //initial credit spread we assume s(0,t) flat initially
	double phi0; //accumulated variance
	
	int J1; // assume J (lrs paper) as constant

    
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
	int visit; //0- if the node wasn't already visit , 1- if the node was visit
    double s; //default intensity or instantaneous credit spread
	double phi[2];// Array with two elements containing the minimum and maximum of phi, which we calculate using the precedent generation
    double *optionPrice; // pointer to an array with the option price
    double *probUp; // pointer to an array with the up probabilities
};


int signal (double Z);

double calcJ(double mean, double sqrtt);

void buildLRStree( double vol,double k_s,int m, int N, double strike, char typeoption, double r0, double s0, 
                  double phi0,double R, const string& dateT, int J1, const vector<string>& dates, const vector<int>& numdate,
				  const vector<double>& discount, const vector<double>& survival);


/*
 *******************Main Contains Menu***************************************************************
 */
 
int main()
{
	// declare and initialize the input parameters
    double vol=0.94, k_s=-0.10, strike= 0.01 , r0=0.01, phi0=0.0, s0=0.02, R=0.4;
    
	string dateT= "20-12-2008";
    
    int m=3, N = 3, J1 = 0;
    
    //declare and initialize type of the option (put or call)
    char typeoption='C';
    
    //extract the values of each column of the text file
    
    vector<string> dates;
    vector<int> daycount;
    vector<double> zerobond;
    vector<double> zerodefault;
    
    ifstream theFile("dataIGdiff.txt");
    
    string reading1;
	int reading2;  
    double reading3;
    double reading4;
    
    
    while(theFile>> reading1>> reading2 >> reading3 >> reading4){
	
    	dates.push_back(reading1);
    	daycount.push_back(reading2);
    	zerobond.push_back(reading3);
    	zerodefault.push_back(reading4);
	
	}
   
    buildLRStree( vol, k_s, m, N, strike, typeoption,  r0, s0, phi0, R, dateT,J1, dates, daycount,zerobond, zerodefault);
    
    return 0;
     
}

void buildLRStree( double vol,double k_s, int m, int N, double strike, char typeoption, double r0, 
                   double s0, double phi0, double R,const string& dateT, int J1, const vector<string>& dates, 
				   const vector<int>& numdate, const vector<double>& discount, const vector<double>& survival){
	
		
	int posT = -1; // posT = position on the txt file of the exercise date of the CDSwaption
	 
    for(unsigned int i=0; i < dates.size();i++){
    	if(dates[i]==dateT)
	 	posT=i;
    }
    
    assert(posT>-1); //verify if in fact exist ---> debug
    
    /*Local variables
	dt = time interval
	sqrt_dt = square root of the time interval
	*/
	double dt=0.0;
	double sqrt_dt=0.0;
	
	dt=((numdate[posT]-numdate[0])/365.0)/N; //time steps
	
	sqrt_dt= sqrt(dt);
	
    //Dynamically allocating memory in each node
	Nodes **tree = new Nodes *[N+1];
	
	for(int i=0; i<=N; i++)
	tree[i] = new Nodes [i+1];
			
	//Inicializing the root of the tree
	tree[0][0].s = s0;
    tree[0][0].phi[0]= phi0;
    tree[0][0].phi[1]= phi0;
    
    for (int i = 0; i <= N ; i++){
    	for (int j = 0; j <= i; j++){
    		tree[i][j].visit=0;
	 		tree[i][j].optionPrice= new double [m];
	 		tree[i][j].probUp = new double [m];
        }
    }
 
    
 // ----> Calculating the values of s, phi and probUp for each node of the tree. Except the probUp on the terminal nodes(i=N)
    
	for(int i=0; i<=N-1;i++){
    	for(int j=0; j<=i;j++){
    		
   		    double mean1=0.0;
   		    double phi_next=0.0;
    		 
    		if(tree[i][j].phi[0]==tree[i][j].phi[1]){
    			    			
    			mean1= (k_s*(s0-tree[i][j].s) + tree[i][j].phi[0])/(vol * tree[i][j].s) - (vol/2);
			    //J1=calcJ(mean1, sqrt_dt);
//		        cout<< i << " " << j << " "<< "mean " << mean1 <<endl; 
//              cout<< i << " " << j << " "<< "meanJ " << calcJ(mean1,sqrt_dt) <<endl; 
//	            
				for(int k=0;k<m;k++)
	            tree[i][j].probUp[k]=(mean1*dt + (1-J1)*sqrt_dt)/(2*sqrt_dt); 
	        
	            // Insert the next values of s, if the nodes aren't visit yet
    		    if(tree[i+1][j].visit==0)
			    tree[i+1][j].s= exp(log(tree[i][j].s) + vol*(J1+1) * sqrt_dt); 
			
			    if(tree[i+1][j+1].visit==0)	
			    tree[i+1][j+1].s= exp(log(tree[i][j].s) + vol*(J1-1) * sqrt_dt); 
			 		
               //Calculate the next value of phi
    		    phi_next= tree[i][j].phi[0] + (pow(vol,2) * pow(tree[i][j].s,2)- 2*k_s*tree[i][j].phi[0])*dt;
    		
    		    if(tree[i+1][j].visit==0){
    		    	tree[i+1][j].phi[0]= phi_next;
    		    	tree[i+1][j].phi[1]=phi_next;
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

    		  			  	  
			} else{
				
				double step=0.0; // step of the partition of phi
				step= (tree[i][j].phi[1]-tree[i][j].phi[0])/(m-1);
    		    
   		        vector<double> partPhi; // partPhi= interval phi with m elemnts
				partPhi.push_back(tree[i][j].phi[0]);
				  
    		    for(int k = 1; k < m ; k++)
    		    partPhi.push_back(partPhi[k-1] + step);
    		                	
    		    vector<double> mean; 
    		    for(int k=0; k < m ; k++)
    		    mean.push_back((k_s * (s0 - tree[i][j].s) + partPhi[k] ) / (vol*tree[i][j].s) - (vol/2));
			    
//			    for( k=0; k < m ; k++)
//			    { J[k]= calcJ(mean[k], sqrt_dt);
//			     // cout<< i << j<< "arrayJ " << J[k]<<endl;
//    		 	}   
//

                 
    		    
    		    for(int k = 0; k< m ;k++)
	            tree[i][j].probUp[k]=(mean[k]*dt + (1-J1)*sqrt_dt)/(2*sqrt_dt); 
	        
				// Insert the next values of s, if the nodes aren't visit yet
				//ATENTION : different values in array J then many values of s ????
    		    if( tree[i+1][j].visit==0)
    		    tree[i+1][j].s= exp(log(tree[i][j].s) + vol*(J1+1) * sqrt_dt); 
			
			    if(tree[i+1][j+1].visit==0)	
			    tree[i+1][j+1].s= exp(log(tree[i][j].s) + vol*(J1-1)* sqrt_dt); 
				
				double sucPhi0=0.0;
                double sucPhim=0.0;    

    		    sucPhi0= partPhi[0]+ (pow(vol,2) * pow(tree[i][j].s,2) - 2 * k_s * partPhi[0] )*dt;
    		    sucPhim= partPhi[m-1]+ (pow(vol,2)* pow(tree[i][j].s,2)- 2 * k_s * partPhi[m-1] )*dt;
    		    
    		    if( tree[i+1][j].visit==0){
    		    	tree[i+1][j].phi[0]= sucPhi0;
    		    	tree[i+1][j].phi[1]= sucPhim;
    		    	tree[i+1][j].visit=1;
			    
			    }else{
			    	tree[i+1][j].phi[0]= min(sucPhi0, tree[i+1][j].phi[0] );
			    	tree[i+1][j].phi[1]= max(sucPhim, tree[i+1][j].phi[1]);
			    }
			
			    if(tree[i+1][j+1].visit==0){
			    tree[i+1][j+1].phi[0]= sucPhi0;
			    tree[i+1][j+1].phi[1]= sucPhim;
			    tree[i+1][j+1].visit=1;
			    }else{
				tree[i+1][j+1].phi[0]= min(sucPhi0, tree[i+1][j+1].phi[0] );
			    tree[i+1][j+1].phi[1]= max(sucPhim,tree[i+1][j+1].phi[1]);
			    }
				
			}
    		
        } //close j
	
	} //close i
	
	
    vector<double> bondPrice; //B(t,T) default-free discount factor
		
		for(unsigned int k=posT; k < discount.size()-1;k++)
		bondPrice.push_back(exp(-r0*((numdate[k+1]-numdate[posT])/365.0)));
				
		vector<double> delta; // delta = interval of time between T_k,...,T_n
		
		for(unsigned int k=posT; k < numdate.size()-1;k++)
		delta.push_back((numdate[k+1]-numdate[k])/365.0);
		

	 
/*---------------Loop on the terminal nodes----------------------------
On the terminal nodes: we calculate the values of optionPrice in each node*/
         
	for(int j=0; j<=N;j++){
				
		if(tree[N][j].phi[0]==tree[N][j].phi[1]){
	    		    		    	
	        vector<double> defaultRisk; //D(t,T) default discount factor
	                          	     	
		    for(unsigned int k=posT; k < discount.size();k++)
		    defaultRisk.push_back((survival[k]/survival[posT])*exp( ((1-exp(- k_s*((numdate[k]-numdate[posT])/365.0))) /k_s)* 
            (s0-tree[N][j].s)- 0.5*pow((1-exp(-k_s*((numdate[k]-numdate[posT])/365.0)))/k_s,2)* tree[N][j].phi[0] ));
                    
            
			vector<double> defaultBondPrice; // default zero coupon bond bar(B)(t,T)
			
			for( unsigned int k=0;k < bondPrice.size();k++)
			defaultBondPrice.push_back(defaultRisk[k+1]*bondPrice[k]);
			
 		    			
			/* Calculate the forward swap rate (fswaprate) 
			num = numerator of the forward swap rate formula
			den= denominator of the forward swap rate formula
			fswaprate = formula of the forward swap rate
			*/
			
			double num = 0.0;
			double den = 0.0;
			
  	
			for(unsigned int k=0; k < bondPrice.size();k++){
				num+= defaultRisk[k]*bondPrice[k]-defaultBondPrice[k];
	  	        den+= delta[k]*defaultBondPrice[k];	
			}
		
			
			double fswaprate=0.0; 		
			fswaprate=(1-R)*(num/den);
			
			
			double sumdefaultBondPrice=0.0;
	        for(unsigned int k=0; k< bondPrice.size();k++)
			sumdefaultBondPrice +=defaultBondPrice[k];
	        
			if(typeoption=='C'){ // if is a payer CDSwaption
				
				for(unsigned int k=0; k<m ;k++)
				tree[N][j].optionPrice[k]=sumdefaultBondPrice* max( fswaprate - strike, 0.0);

	        }else{ // if is a receiver CDSwaption
	        	        	
				for(unsigned int k=0;k<m;k++)
				tree[N][j].optionPrice[k]=sumdefaultBondPrice* max( strike - fswaprate , 0.0);
		    }
		 			
	    } else{
	    	
			double step=0.0;		       		
            step=(tree[N][j].phi[1]-tree[N][j].phi[0])/(m-1);
            
            vector<double> partPhi;
			partPhi.push_back(tree[N][j].phi[0]);
				  
		    for(int k = 1; k < m ; k++)
		    partPhi.push_back(partPhi[k-1] + step);
		    
         		          
            for(int k=0; k< m ;k++){ // for each value of the partition of phi
            	vector<double> defaultRisk;
     	        
            	
				for(int u=posT ; u < discount.size();u++)
				defaultRisk.push_back((survival[u]/survival[posT])*exp( ((1-exp(- k_s*((numdate[u]-numdate[posT])/365.0))) /k_s)*
                (s0-tree[N][j].s)- 0.5*pow((1-exp(-k_s*((numdate[u]-numdate[posT])/365.0)))/k_s,2)* partPhi[k] ));
                        
            	
			    vector<double> defaultBondPrice;
			    
				for(unsigned int u=0;u<bondPrice.size();u++)
	            defaultBondPrice.push_back(defaultRisk[u+1]*bondPrice[u]); 
	  	        	  
	            //---forward swap rate--
	            double num = 0.0;
	            double den = 0.0;
	           	
				for(int u=0; u< bondPrice.size() ;u++){
					num+= defaultRisk[u]*bondPrice[u]-defaultBondPrice[u];
					den+=delta[u]*defaultBondPrice[u];
	            }
	            
	            double fswaprate=0.0;
	            fswaprate=(1-R)*(num/den);
	            
	            
	            double sumdefaultBondPrice=0.0;
	            
	            for(unsigned int u=0; u< defaultBondPrice.size();u++)
			    sumdefaultBondPrice +=defaultBondPrice[u];
	            
                if(typeoption=='C')
                tree[N][j].optionPrice[k]=sumdefaultBondPrice* max( fswaprate - strike, 0.0);
				else
				tree[N][j].optionPrice[k]=sumdefaultBondPrice* max( strike - fswaprate , 0.0);
				
		    }	
			
    			
    	  }    			
			
			
	}//end of loop on the terminal nodes
	
 //------------------Backward recursion----------------------
 
// 
    for( int i=N-1; i>=0;i--){
	    for( int j=0; j<=i; j++)
       {
  	 	         
       	if(tree[i][j].phi[0]==tree[i][j].phi[1]){
       			     		
      		for(int k=0;k<m;k++)
      		tree[i][j].optionPrice[k]= (tree[i][j].probUp[k]*tree[i+1][j].optionPrice[0]+(1-tree[i][j].probUp[k])*
            tree[i+1][j+1].optionPrice[0]) *exp(-tree[i][j].s*dt);       	   
		  
		} else{
						
        	double step=0.0;
        	step= (tree[i][j].phi[1]-tree[i][j].phi[0])/(m-1);
        	   	
            vector<double> partPhi;
			partPhi.push_back(tree[i][j].phi[0]);
				  
		    for(int k = 1; k < m ; k++)
		    partPhi.push_back(partPhi[k-1] + step);
    		
    		vector <double> sucPhi;
    		for(int k=0;k<m;k++)
			sucPhi.push_back(partPhi[k]+ (pow(vol,2) * pow(tree[i][j].s,2)- 2 * k_s * partPhi[k] )*dt);
			
			   	    		
    		double step2=0.0;          
			step2= (tree[i+1][j].phi[1]-tree[i+1][j].phi[0])/(m-1);
			
			vector<double> partPhi2;
			partPhi2.push_back(tree[i+1][j].phi[0]);
				  
		    for(int k = 1; k < m ; k++)
		    partPhi2.push_back(partPhi2[k-1] + step2);
			
			
    		double step3=0.0;
    		step3= (tree[i+1][j+1].phi[1]-tree[i+1][j+1].phi[0])/(m-1);

			vector<double> partPhi3;
			partPhi3.push_back(tree[i+1][j+1].phi[0]);
			
    		for( int k = 1; k < m ; k++)
    		partPhi3.push_back(partPhi3[k-1] + step3);
			
    						
    		for(int k=0;k<m;k++){ //Interpolation.................................
			                   
    			
				double gup=0.0;
                double gdown=0.0;
    			
    			if( (sucPhi[k] < partPhi2[k]) && (sucPhi[k] > partPhi3[k])){
    				
    				
				gup=tree[i+1][j].optionPrice[k-1] + (sucPhi[k]-partPhi2[k-1])/(partPhi2[k]-partPhi2[k-1] ) *
				 (tree[i+1][j].optionPrice[k]-tree[i+1][j].optionPrice[k-1]);
				 
                gdown=tree[i+1][j+1].optionPrice[k] + (sucPhi[k]-partPhi3[k])/(partPhi3[k+1]-partPhi3[k] ) * 
				(tree[i+1][j+1].optionPrice[k+1]-tree[i+1][j+1].optionPrice[k]);
				
                tree[i][j].optionPrice[k]=(tree[i][j].probUp[k] * gup + (1-tree[i][j].probUp[0])*gdown) * exp(-tree[i][j].s*dt);
		  
			    }
			    
			          
		       if((sucPhi[k] > partPhi3[k]) && !(sucPhi[k] < partPhi2[k])){
		       	 
                gdown=tree[i+1][j+1].optionPrice[k] + (sucPhi[k]-partPhi3[k])/(partPhi3[k+1]-partPhi3[k] ) *
				 (tree[i+1][j+1].optionPrice[k+1]-tree[i+1][j+1].optionPrice[k]);

     			 tree[i][j].optionPrice[k]=(tree[i][j].probUp[k]*tree[i+1][j].optionPrice[k]+
				  (1-tree[i][j].probUp[k])*gdown)*exp(-tree[i][j].s*dt);

			    }
			    
			    if((sucPhi[k] < partPhi2[k])&& !(sucPhi[k] > partPhi3[k])){ 
			    
			    gup=tree[i+1][j].optionPrice[k-1] + (sucPhi[k]-partPhi2[k-1])/(partPhi2[k]-partPhi2[k-1] ) * 
				(tree[i+1][j].optionPrice[k]-tree[i+1][j].optionPrice[k-1]);
				
     	        tree[i][j].optionPrice[k]=(tree[i][j].probUp[k]*gup + 
				 (1-tree[i][j].probUp[k])*tree[i+1][j+1].optionPrice[k] ) * exp(-tree[i][j].s*dt);

			    }
			
	        } 		
          }
        }  	
  	}
  	
  	for( int i=0 ; i<=N ; i++){
  		for( int j=0;j<=i;j++){
  			 cout<< i <<" "<< j << " "<<" dspread "<< tree[i][j].s<< endl;
             cout <<i<<" "<< j << " "<<" phi0 "<<  tree[i][j].phi[0]<< endl;
             cout <<i<<" "<< j << " "<<" phi1 "<< tree[i][j].phi[1]<< endl;
                 
  	        for(int k=0;k<m;k++){
    	        cout <<i<<" "<< j << " prob "<< tree[i][j].probUp[k]<< endl;
  	        	cout <<i<<" "<< j << " "<< "optionPrice "<<tree[i][j].optionPrice[k]<< endl;
			  }
  	       
        }
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


  
