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

using namespace std;

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
 *******************Main Contains Menu***************************************************************+
 */
 
void buildLRStree( double vol,double k_s, int m, double T, int N, double strike, double cdsT, char typeoption, double r0, double s0, double phi0,double R,string dateT, int J1, vector<string>& read1, vector<double>& read2, vector<double>& read3);

int main()
{
  // declare and initialise LRS parameters
  
  double vol=0.2, k_s=0.02, strike= 0.01 , r0=0.04, phi0=0.0, s0=0.04, T=3.0, cdsT=5.0, R=0.4;
  string dateT= "20-09-2009";
  
  // declare and initialise tree paramaters (steps in tree)
  int m=3,N=3,J1=0;
  
  //declare and initialize type of the option (put or call)
  char typeoption='C';
  //double result;
  
  
   vector<string> vecread1;
   vector<double> vecread2;
   vector<double> vecread3;
   
   ifstream theFile("dataIG.txt");//ifstream allows to read data from a file
	
   string reading1;
   double reading2;
   double reading3;
   
   
   while(theFile>> reading1>> reading2 >> reading3){
	
	vecread1.push_back(reading1);
	vecread2.push_back(reading2);
	vecread3.push_back(reading3);
	//cout<<reading1<< "," <<reading2<< "," << reading3 << endl;
	}
   
   buildLRStree( vol, k_s, m,  T, N, strike,  cdsT,  typeoption,  r0, s0, phi0, R, dateT,J1, vecread1, vecread2, vecread3);
   system("pause");
   return 0;
     
}

void buildLRStree( double vol,double k_s, int m, double T, int N, double strike, double cdsT, char typeoption, double r0, double s0, double phi0,double R,string dateT, int J1, vector<string>& read1, vector<double>& read2, vector<double>& read3){
	
	// Declaration of local variables
	double dt;
	double sqrt_dt;
	
	//Inicialization of local variables
	dt=T/N; //time steps
	cout<< "int dt " << dt <<endl;
	sqrt_dt= sqrt(dt);
	
	//....forwardprice of th bond
	
	//indices
	int i,j,k;
	
	
	//Dynamically allocating memory in each node of the tree or ( in each entrie of the lower triangular matrix)
	Nodes **tree = new Nodes *[N];
	
	for(i=0; i<=N; i++)
	 tree[i] = new Nodes [i+1];
	 
		
	//Inicializing the root of the tree
	 //tree[0][0].r = r0;
	 tree[0][0].s = s0;
	 tree[0][0].phi[0]= phi0;
	 tree[0][0].phi[1]= phi0;
	 
	 
	 for (int i = 0; i <= N ; i++)
    {
        for (int j = 0; j <= i; j++)
        { 
         tree[i][j].visit=0;
         tree[i][j].optionPrice= new double [m];
         tree[i][j].probUp = new double [m];
            
        }
    }
 
 /********Local variables ************/
 //to calculate the succesor of phi
 
     double mean1;
     //double J1;
     double phi_next;
     double step;
     double step2;
     double step3;
     double partPhi[m]; //all the m values of phi(i,j)
     double partPhi2[m];
     double partPhi3[m];
	 double mean[m];
     //double J[m];
     double sucPhi0;
     double sucPhim;
     double gplus;
     double gless;
     double sucPhi[m]; //allthe m successor values of phi(i+1,j) or phi(i+1,j+1)
     
 
 /******************/
  
 
    for(i=0; i<=N-1;i++){ //for each ith time step except the last one
    	
    	for(j=0; j<=i;j++){ //for each jth node in ith time step
    		
    		
    		if(tree[i][j].phi[0]==tree[i][j].phi[1])
			{
			 mean1= (k_s*(s0-tree[i][j].s) + tree[i][j].phi[0])/(vol * tree[i][j].s) - (vol/2);
//			 J1=calcJ(mean1, sqrt_dt);
//			 cout<< i << " " << j << " "<< "mean " << mean1 <<endl; 
		     //cout<< i << " " << j << "J1 " << J1 <<endl; 
	 
	        for(k=0;k<m;k++)
	        tree[i][j].probUp[k]=(mean1*dt + (1-J1)*sqrt_dt)/(2*sqrt_dt); 
	        
	        // Insert the next values of r, if the nodes aren't visit yet
    		if( tree[i+1][j].visit==0)
			tree[i+1][j].s= exp(log(tree[i][j].s) + vol * sqrt_dt); 
			
			if(tree[i+1][j+1].visit==0)	
			tree[i+1][j+1].s= exp(log(tree[i][j].s) - vol * sqrt_dt); 

			 		
            //Calculate the successor of phi
    		phi_next= tree[i][j].phi[0] + (pow(vol,2) * pow(tree[i][j].s,2)- 2*k_s*tree[i][j].phi[0])*dt;
    		
    		if( tree[i+1][j].visit==0){
    			tree[i+1][j].phi[0]= phi_next;
			    tree[i+1][j].phi[1]=phi_next;
			    tree[i+1][j].visit=1;
			}else{
				tree[i+1][j].phi[0]= min(phi_next, tree[i+1][j].phi[0] );
			    tree[i+1][j].phi[1]=max(phi_next,tree[i+1][j].phi[1]);
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
				
				step= (tree[i][j].phi[1]-tree[i][j].phi[0])/(m-1);
    		    
				partPhi[0]= tree[i][j].phi[0];
    		    for(k = 1; k < m ; k++)
    		    partPhi[k]= partPhi[k-1] + step;   
    		    
    		    for( k=0; k < m ; k++)
    		      mean[k]= (k_s * (s0 - tree[i][j].s) + partPhi[k] ) / (vol*tree[i][j].s) - (vol/2);
			    
//			    for( k=0; k < m ; k++)
//			    { J[k]= calcJ(mean[k], sqrt_dt);
//			     // cout<< i << j<< "arrayJ " << J[k]<<endl;
//    		 	}   
    		    
    		    for(k = 0; k< m ;k++)
	            tree[i][j].probUp[k]=(mean[k]*dt + (1-J1)*sqrt_dt)/(2*sqrt_dt); 
	        
				// Insert the next values of r, if the nodes aren't visit yet
				//ATENTION : different values in array J then many values of r ???? many insert an if
    		    if( tree[i+1][j].visit==0)
			        tree[i+1][j].s= exp(log(tree[i][j].s) + vol * sqrt_dt); 
			
			    if(tree[i+1][j+1].visit==0)	
			       tree[i+1][j+1].s= exp(log(tree[i][j].s) - vol* sqrt_dt); 

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
    		
        } //fecha j
	
	} //fecha i
	
	//Declaration:
     int posT;
     int pos;
     int u;
   vector<double> bondPrice;
   vector<double> defaultRisk;
   vector<double> defaultBondPrice;
     double num;
	 double den;
	 double delta;
	 double fswaprate;
	 double payoffpayer=0.0;
	 double payoffreceiver=0.0;
	 
	 delta=0.25;
//	  vector<double> delta; //quaterly
//	  
//	  delta.push_back();
    for(i=0; i < read1.size();i++){
    if(read1[i]==dateT)
      posT=i;
	 }
	 
	 cout<<posT<<endl;
	 
//-----------------Loop on the terminal nodes----------------------------
         
	for(j=0; j<=N;j++)
	
	{  
	    		
	  if(tree[N][j].phi[0]==tree[N][j].phi[1]){
		
		  for(k=posT; k < read2.size();k++)
         { 
          bondPrice.push_back(exp(-r0*(cdsT-T)/365));
          defaultRisk.push_back((read3[k]/read3[posT])*exp( ((1-exp(- k_s*(cdsT-T))) /k_s)* (s0-tree[N][j].s)- 0.5*pow((1-exp(-k_s*(cdsT-T)))/k_s,2)* tree[N][j].phi[0] ));
	     }
	     
//	     for(k=0; k < bondPrice.size();k++){
//	     	cout<<bondPrice[k]<<endl;
//	        cout<<defaultRisk[k]<<endl;
//		 }
//	     
	  
	      for(k=0;k<bondPrice.size();k++)
	     {
	  	  defaultBondPrice.push_back(defaultRisk[k+1]*bondPrice[k]); 
	  	// D(T,T_k) and B(T,T_k+1) for obtain defaultB(T, T_k+1)
	     }
	  
	  //---forward swap rate--
	  
	      for(k=0; k < bondPrice.size() ;k++)
	     {
	  	  num+= defaultRisk[k]*bondPrice[k]-defaultBondPrice[k];
	  	  den+= delta*defaultBondPrice[k];
	     }
	  
//	     cout<<"num "<<num<<endl;
//	     cout<<"den "<<den<<endl;
	     
	     fswaprate=(1-R)*(num/den);
//	     cout<< "fswaprate"<< fswaprate <<endl;
	  	     
	     if(typeoption=='C'){
			
			for(k=0; k< bondPrice.size() ;k++)
	  	    payoffpayer +=defaultBondPrice[k]* max( fswaprate - strike, 0.0);
	    
			for(k=0;k<m;k++)
			tree[N][j].optionPrice[k]=payoffpayer;
			
			
	      }else{
			for(k=0; k< bondPrice.size() ;k++)
	        payoffreceiver+= defaultBondPrice[k]* max( strike - fswaprate , 0.0);
	   
			for(k=0;k<m;k++)
			tree[N][j].optionPrice[k]=payoffreceiver;
		 }	
			
		} else{
			
			step = (tree[N][j].phi[1]-tree[N][j].phi[0])/(m-1);
			partPhi[0]= tree[N][j].phi[0];
    		for(k = 1; k < m ; k++)
    		partPhi[k]= partPhi[k-1] + step;
    		
    		for(k=0; k< m ;k++)
			{
    		   for(u=posT;u <read2.size();u++)
               { 
               bondPrice.push_back(exp(-r0*(cdsT-T)/365));
               defaultRisk.push_back((read3[u]/read3[posT])*exp( ((1-exp(- k_s*(cdsT-T))) /k_s)* (s0-tree[N][j].s)- 0.5*pow((1-exp(-k_s*(cdsT-T)))/k_s,2)* partPhi[k] ));
	           }
	  
	           for(u=0;u<bondPrice.size();u++)
	           {
	  	       defaultBondPrice.push_back(defaultRisk[u+1]*bondPrice[u]); 
	  	       // D(T,T_k) and B(T,T_k+1) for obtain defaultB(T, T_k+1)
	           }
	  
	         //---forward swap rate--
	  
	          for(u=0; u< bondPrice.size() ;u++)
	          {
	        	num+= defaultRisk[u]*bondPrice[u]-defaultBondPrice[u];
	        	den+= delta*defaultBondPrice[u];
	          }
	  
	         fswaprate=(1-R)*(num/den);
	  
	         for(u=0; u< bondPrice.size() ;u++)
	         {
	  	     payoffpayer += defaultBondPrice[u]* max( fswaprate - strike, 0.0);
	         payoffreceiver+= defaultBondPrice[u]* max( strike - fswaprate , 0.0);
	         }
	    
		    if(typeoption=='C'){
			
			tree[N][j].optionPrice[k]=payoffpayer;
			
		    }else{
			
			tree[N][j].optionPrice[k]=payoffreceiver;
			
		    }	
			
    			
    	  }    			
			
		}
			
	}//end of loop on the terminal nodes
	
 //	  //Backward recursion
	for(i=N-1;i>=0;i--)
     {
       for( j=0;j<=i;j++)
       {
       	
       	if(tree[i][j].phi[0]==tree[i][j].phi[1])
       	{
       	  for(k=0;k<m;k++ )
	      tree[i][j].optionPrice[k]= (tree[i][j].probUp[k]*tree[i+1][j].optionPrice[0]+(1-tree[i][j].probUp[k])*tree[i+1][j+1].optionPrice[m-1]) *exp(-tree[i][j].s*dt);       	   
		  
		} else{
        	
        	step= (tree[i][j].phi[1]-tree[i][j].phi[0])/(m-1);
        	partPhi[0]= tree[i][j].phi[0];
            
			for(k=1;k<m;k++)
    		partPhi[k]= partPhi[k-1] + step;
    		
    		for(k=0;k<m;k++)
    		{
    		  sucPhi[k]= partPhi[k]+ (pow(vol,2) * pow(tree[i][j].s,2)- 2 * k_s * partPhi[k] );
//    		  cout << i << " " << j << " sucphi " << sucPhi[k] << endl;
              }
    		
            
            
			step2= (tree[i+1][j].phi[1]-tree[i+1][j].phi[0])/(m-1);
			partPhi2[0]= tree[i+1][j].phi[0];
    		for( k = 1; k < m ; k++)
    		partPhi2[k]= partPhi2[k-1] + step2;
    		
    		step3= (tree[i+1][j+1].phi[1]-tree[i+1][j+1].phi[0])/(m-1);

			partPhi3[0]= tree[i+1][j+1].phi[0];
    		for( k = 1; k < m ; k++)
    		{
    			partPhi3[k]= partPhi3[k-1] + step3;
    			
			}
			
    		
    		    		
    		for(k=0;k<m;k++)
    		{ 
    		  

    		  if( (sucPhi[k] < partPhi2[k]) && (sucPhi[k] > partPhi3[k]))
			  { 
			    //Interpolation.................................
			   gplus=tree[i+1][j].optionPrice[k-1] + (sucPhi[k]-partPhi2[k-1])/(partPhi2[k]-partPhi2[k-1] ) * (tree[i+1][j].optionPrice[k]-tree[i+1][j].optionPrice[k-1]);
			   gless=tree[i+1][j+1].optionPrice[k] + (sucPhi[k]-partPhi3[k])/(partPhi3[k+1]-partPhi3[k] ) * (tree[i+1][j+1].optionPrice[k+1]-tree[i+1][j+1].optionPrice[k]);
               tree[i][j].optionPrice[k]=(tree[i][j].probUp[k] * gplus + (1-tree[i][j].probUp[0])*gless) * exp(-tree[i][j].s*dt);
		  
			  }
		      
		      if((sucPhi[k] > partPhi3[k]) && !(sucPhi[k] < partPhi2[k]))
		      { 
		        gless=tree[i+1][j+1].optionPrice[k] + (sucPhi[k]-partPhi3[k])/(partPhi3[k+1]-partPhi3[k] ) * (tree[i+1][j+1].optionPrice[k+1]-tree[i+1][j+1].optionPrice[k]);
//     	        cout << i << " " << j << " gless(21)" << gless << endl;
     			 tree[i][j].optionPrice[k]=(tree[i][j].probUp[k]*tree[i+1][j].optionPrice[k]+(1-tree[i][j].probUp[k])*gless) * exp(-tree[i][j].s*dt);
//                cout << i << " " << j << "optpriv  " << tree[i][j].optionPrice[k] << endl;

			  }
			  
			  if((sucPhi[k] < partPhi2[k])&& !(sucPhi[k] > partPhi3[k]))
		      { 
			    gplus=tree[i+1][j].optionPrice[k-1] + (sucPhi[k]-partPhi2[k-1])/(partPhi2[k]-partPhi2[k-1] ) * (tree[i+1][j].optionPrice[k]-tree[i+1][j].optionPrice[k-1]);
     	        tree[i][j].optionPrice[k]=(tree[i][j].probUp[k]*gplus + (1-tree[i][j].probUp[k])*tree[i+1][j+1].optionPrice[k] ) * exp(-tree[i][j].s*dt);

			  }
			
            
			  
			  } 		
			}
		    	            
    
          
    
		}
       	
       	
       }
	
 for(i=0;i<=N;i++){
 	
 	 for(j=0;j<=i;j++){
 	 cout<< i <<" "<< j << " "<<" rate "<< tree[i][j].s<< endl;
  	 for(k=0;k<m;k++)
      cout <<i<<" "<< j << " "<< tree[i][j].probUp[k]<< endl;
//     //cout <<i<<" "<< j << " "<< tree[i][j].optionPrice[k]<< endl;

   }
 }  
  
  //return tree[0][0].s;
}


  
