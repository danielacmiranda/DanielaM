/* INPUT: Parameter models
	
	double vol; //volatitily
	double k_s; //exogenous factor
	(double lambda; elasticity assume lambda=1 ---> for now this is not an input)
	
	int m; //number of elements in the partition of phi
	double T; //time to Maturity of the CDSwaption
	int N; //number of steps in the tree
	double strike; // strike price of the CDSwaption
	int cdsT; //maturity of the underlying forward CDS 
	double R; //recovery rate
	
	string dateT; maturity date of the CDSwaption
	
	char typeoption;// call 'C' or put option 'P' respectively payer or receiver cdswaption
	
	double s0; //credit spread we assume s(0,t) flat initially
	double phi0; //accumulated variance
	
	int J1;

    
*/

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include<fstream>
#include<algorithm>
#include<assert.h>

using namespace std;

const int m = 3;

/* Node Declaration
 */
struct Nodes 
{
	int visit; //0- if i don't already visit the node , 1- if i already visit the node
    double s; //default intensity or credit spread
	double phi[2];	 // Array with two elements containing the minimum anda maximum of phi, which we calculate using the precedent generation
    double *optionPrice;
    double *probUp;
};

/* Auxilary Functions: Z and calcJ -----> for calculating J (LRS paper)*/
int signal ( double Z){
	
	int result;
	
	if(Z>0)
	{
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
	
	if ( Z % 2 == 0){
		J=signal(Z)*Z;
	} else{
		J=signal(Z)*abs(Z) +1;
	}
	
	return J; 
		
}

/*
 *******************Main Contains Menu***************************************************************
 */
 
void buildLRStree( double vol,double k_s, int N, double strike, char typeoption, double r0, double s0, 
                  double phi0,double R, const string& dateT, int J1, const vector<string>& dates, const vector<int>& numdate,
				  const vector<double>& discount, const vector<double>& survival);

int main()
{
	// declare and initialise LRS parameters
    double vol=0.2, k_s=0.02, strike= 0.01 , r0=0.04, phi0=0.0, s0=0.04, R=0.4;
    string dateT= "20-09-2009";
    // declare and initialise tree paramaters (steps in tree)
    int N = 3, J1 = 0;
    //declare and initialize type of the option (put or call)
    char typeoption='C';
    
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
   
    buildLRStree( vol, k_s, N, strike, typeoption,  r0, s0, phi0, R, dateT,J1, dates, daycount,zerobond, zerodefault);
    
    return 0;
     
}

void buildLRStree( double vol,double k_s, int N, double strike, char typeoption, double r0, 
                   double s0, double phi0, double R,const string& dateT, int J1, const vector<string>& dates, 
				   const vector<int>& numdate, const vector<double>& discount, const vector<double>& survival){
	
	// Declaration of local variables
	double dt;
	double sqrt_dt;
	
	//Inicialization of local variables
	dt=((numdate[numdate.size()-1]-numdate[0])/365.0)/N; //time steps
	sqrt_dt= sqrt(dt);
	cout<< "int dt " << dt <<endl;

	//Dynamically allocating memory in each node
	Nodes **tree = new Nodes *[N];
	
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
 
    /********Local variables ************/
    //to calculate the succesor of phi
     double step2;
     double step3;
     //all the m values of phi(i,j)
     double partPhi2[m];
     double partPhi3[m];
	 
     //double J[m];
    
     double gplus;
     double gless;
     double sucPhi[m]; //allthe m successor values of phi(i+1,j) or phi(i+1,j+1)
     
 
 /******************/
    
    
     //double J1;
     
      
 
    for(int i=0; i<=N-1;i++){
    	for(int j=0; j<=i;j++){
    		
   		    double mean1=0.0;
   		    double phi_next=0.0;
    		 
    		if(tree[i][j].phi[0]==tree[i][j].phi[1]){
    			    			
    			mean1= (k_s*(s0-tree[i][j].s) + tree[i][j].phi[0])/(vol * tree[i][j].s) - (vol/2);
//			    J1=calcJ(mean1, sqrt_dt);
//			    cout<< i << " " << j << " "<< "mean " << mean1 <<endl; 

	            for(int k=0;k<m;k++)
	            tree[i][j].probUp[k]=(mean1*dt + (1-J1)*sqrt_dt)/(2*sqrt_dt); 
	        
	            // Insert the next values of r, if the nodes aren't visit yet
    		    if(tree[i+1][j].visit==0)
			    tree[i+1][j].s= exp(log(tree[i][j].s) + vol * sqrt_dt); 
			
			    if(tree[i+1][j+1].visit==0)	
			    tree[i+1][j+1].s= exp(log(tree[i][j].s) - vol * sqrt_dt); 
			 		
               //Calculate the successor of phi
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
				
				double step=0.0;
				step= (tree[i][j].phi[1]-tree[i][j].phi[0])/(m-1);
    		    
   		        vector<double> partPhi;
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
    		    
    		    for(int k = 0; k< m ;k++)
	            tree[i][j].probUp[k]=(mean[k]*dt + (1-J1)*sqrt_dt)/(2*sqrt_dt); 
	        
				// Insert the next values of r, if the nodes aren't visit yet
				//ATENTION : different values in array J then many values of r ???? many insert an if
    		    if( tree[i+1][j].visit==0)
    		    tree[i+1][j].s= exp(log(tree[i][j].s) + vol * sqrt_dt); 
			
			    if(tree[i+1][j+1].visit==0)	
			    tree[i+1][j+1].s= exp(log(tree[i][j].s) - vol* sqrt_dt); 
				
				double sucPhi0=0.0;
                double sucPhim=0.0;    

    		    sucPhi0= partPhi[0]+ (pow(vol,2) * pow(tree[i][j].s,2) - 2 * k_s * partPhi[0] );
    		    sucPhim= partPhi[m-1]+ (pow(vol,2)* pow(tree[i][j].s,2)- 2 * k_s * partPhi[m-1] );
    		    
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
	
	} //closei
	
	double delta;
 
	delta=(numdate[5]-numdate[4])/365.0;
    cout<< numdate[5]<<endl;
	cout<<"delta "<< delta<<endl;
	
	int posT = -1;
	 
    for(unsigned int i=0; i < dates.size();i++){
    	if(dates[i]==dateT)
	 	posT=i;
    }
    
    assert(posT>-1); //verify if in fact exist ---> debug
	 
    cout<< "position exercise date "<<posT<<endl;
    cout<< "position CDS maturity date " << dates.size()<<endl;
	 
//-----------------Loop on the terminal nodes----------------------------
         
	for(int j=0; j<=N;j++){
				
		vector<double> bondPrice;
		
		for(unsigned int k=posT; k < discount.size()-1;k++)
		bondPrice.push_back(exp(-r0*((numdate[k+1]-numdate[posT])/365.0)));
		cout<<j << "sizeB " << bondPrice.size()<<endl;	
		
	    if(tree[N][j].phi[0]==tree[N][j].phi[1]){
	    		    		    	
	        vector<double> defaultRisk;
	                          	     	
		    for(unsigned int k=posT; k < discount.size();k++)
		    defaultRisk.push_back((survival[k]/survival[posT])*exp( ((1-exp(- k_s*((numdate[k]-numdate[posT])/365.0))) /k_s)* 
            (s0-tree[N][j].s)- 0.5*pow((1-exp(-k_s*((numdate[k]-numdate[posT])/365.0)))/k_s,2)* tree[N][j].phi[0] ));
            
//            for(unsigned int k=posT; k < bondPrice.size();k++)
//            cout<<j << "D " << defaultRisk[k]<<endl;
//            cout<<j << "sizeD " << defaultRisk.size()<<endl;
            
			vector<double> defaultBondPrice;
			
			for( unsigned int k=0;k < bondPrice.size();k++)
			defaultBondPrice.push_back(defaultRisk[k+1]*bondPrice[k]);
			cout<<j << "sizebarB " << defaultBondPrice.size()<<endl;
			//---forward swap rate--
			double num = 0.0;
			double den = 0.0;
			
			
			for(unsigned int k=0; k < bondPrice.size();k++){
				num+= defaultRisk[k]*bondPrice[k]-defaultBondPrice[k];
	  	        den+= delta*defaultBondPrice[k];	
			}
			
			double fswaprate=0.0;			
			fswaprate=(1-R)*(num/den);
			
			double payoffpayer=0.0;
	        double payoffreceiver=0.0;
	        
			if(typeoption=='C'){
				for(unsigned int k=0; k< bondPrice.size();k++)
				payoffpayer +=defaultBondPrice[k]* max( fswaprate - strike, 0.0);
				
				for(unsigned int k=0; k<m ;k++)
				tree[N][j].optionPrice[k]=payoffpayer;
			
	        }else{
	        	for(unsigned int k=0; k< bondPrice.size() ;k++)
	        	payoffreceiver+= defaultBondPrice[k]* max( strike - fswaprate , 0.0);
	        	
				for(unsigned int k=0;k<m;k++)
				tree[N][j].optionPrice[k]=payoffreceiver;
		    }
		 			
	    } else{
	    	
			double step=0.0;		       		
            step=(tree[N][j].phi[1]-tree[N][j].phi[0])/(m-1);
            
            vector<double> partPhi;
			partPhi.push_back(tree[N][j].phi[0]);
				  
		    for(int k = 1; k < m ; k++)
		    partPhi.push_back(partPhi[k-1] + step);
               		          
            for(int k=0; k< m ;k++){
            	vector<double> defaultRisk;
				         	
            	for(int u=posT;u < discount.size()-1;u++)
            	defaultRisk.push_back((survival[u]/survival[posT])*exp( ((1-exp(- k_s*((numdate[u]-numdate[posT])/365.0))) /k_s)*
                (s0-tree[N][j].s)- 0.5*pow((1-exp(-k_s*((numdate[u+1]-numdate[posT])/365.0)))/k_s,2)* partPhi[k] ));
            
			    vector<double> defaultBondPrice;
			    
				for(int u=0;u<bondPrice.size();u++)
	            defaultBondPrice.push_back(defaultRisk[u+1]*bondPrice[u]); 
	  	        	  
	            //---forward swap rate--
	            double num = 0.0;
	            double den = 0.0;
	            
				for(int u=0; u< bondPrice.size() ;u++){
					num+= defaultRisk[u]*bondPrice[u]-defaultBondPrice[u];
					den+= delta*defaultBondPrice[u];
	            }
	            
	            double fswaprate=0.0;
	            fswaprate=(1-R)*(num/den);
	            
	            double payoffpayer=0.0;
	            double payoffreceiver=0.0;
	            for(int u=0; u< bondPrice.size() ;u++){
	            	payoffpayer += defaultBondPrice[u]* max( fswaprate - strike, 0.0);
	            	payoffreceiver+= defaultBondPrice[u]* max( strike - fswaprate , 0.0);
	            }
	            
                if(typeoption=='C')
                tree[N][j].optionPrice[k]=payoffpayer;
				else
				tree[N][j].optionPrice[k]=payoffreceiver;
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
          tree[i+1][j+1].optionPrice[m-1]) *exp(-tree[i][j].s*dt);       	   
		  
		} else{
						
        	double step=0.0;
        	step= (tree[i][j].phi[1]-tree[i][j].phi[0])/(m-1);
        	   	
            vector<double> partPhi;
			partPhi.push_back(tree[i][j].phi[0]);
				  
		    for(int k = 1; k < m ; k++)
		    partPhi.push_back(partPhi[k-1] + step);
    		
    		vector <double> sucPhi;
    		for(int k=0;k<m;k++)
    		sucPhi.push_back(partPhi[k]+ (pow(vol,2) * pow(tree[i][j].s,2)- 2 * k_s * partPhi[k] ));
    		
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
    		partPhi3.push_back(partPhi3[k-1] + step2);
    							
    		for(int k=0;k<m;k++){
    			
    			if( (sucPhi[k] < partPhi2[k]) && (sucPhi[k] > partPhi3[k])){
				//Interpolation.................................
			    double gup=0.0;
                double gdown=0.0;
                
				gup=tree[i+1][j].optionPrice[k-1] + (sucPhi[k]-partPhi2[k-1])/(partPhi2[k]-partPhi2[k-1] ) *
				 (tree[i+1][j].optionPrice[k]-tree[i+1][j].optionPrice[k-1]);
				 
                gdown=tree[i+1][j+1].optionPrice[k] + (sucPhi[k]-partPhi3[k])/(partPhi3[k+1]-partPhi3[k] ) * 
				(tree[i+1][j+1].optionPrice[k+1]-tree[i+1][j+1].optionPrice[k]);
				
                tree[i][j].optionPrice[k]=(tree[i][j].probUp[k] * gup + (1-tree[i][j].probUp[0])*gdown) * exp(-tree[i][j].s*dt);
		  
			    }
			    
			          
		       if((sucPhi[k] > partPhi3[k]) && !(sucPhi[k] < partPhi2[k])){
		       	
		       	double gdown=0.0;
                gdown=tree[i+1][j+1].optionPrice[k] + (sucPhi[k]-partPhi3[k])/(partPhi3[k+1]-partPhi3[k] ) *
				 (tree[i+1][j+1].optionPrice[k+1]-tree[i+1][j+1].optionPrice[k]);

     			 tree[i][j].optionPrice[k]=(tree[i][j].probUp[k]*tree[i+1][j].optionPrice[k]+
				  (1-tree[i][j].probUp[k])*gdown)*exp(-tree[i][j].s*dt);

			    }
			    
			    if((sucPhi[k] < partPhi2[k])&& !(sucPhi[k] > partPhi3[k])){ 
			    
			    double gup;
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
  			 cout<< i <<" "<< j << " "<<" rate "<< tree[i][j].s<< endl;
             cout <<i<<" "<< j << " "<<" phi0 "<<  tree[i][j].phi[0]<< endl;
             cout <<i<<" "<< j << " "<<" phi1 "<< tree[i][j].phi[1]<< endl;
                 
  	        for(int k=0;k<m;k++){
    	        //cout <<i<<" "<< j << " "<< tree[i][j].probUp[k]<< endl;
  	        	cout <<i<<" "<< j << " "<< tree[i][j].optionPrice[k]<< endl;
			  }
  	       
        }
    }  
  
  //return tree[0][0].optionPrice[0];
}


  
