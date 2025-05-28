#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>     
#include <time.h>
#include<vector>
#include "Gauge_header.h"
#include <fstream>
#include <filesystem>
#include <cstdlib>

/***************************************************************************************************************************/

// THIS FUNCTION INIZIALIZE THE LATTICE IN THE FROZEN CONFIGURATION, I.E. WITH ALL THE LINKS +1
void lattice4_creation(vector<site4> &lattice, int n, int nt){ 
  
  site4 element; 

  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
        for(int k=0; k<n; k++){
            for(int l=0; l<nt; l++){
                element.position[0]=i;
                element.position[1]=j;
                element.position[2]=k;
                element.position[3]=l;

                element.value_direction[0]=1; // This is the link along the direction 0 from the site of coordinate (i,j,k,l)
                element.value_direction[1]=1; // This is the link along the direction 1 from the site of coordinate (i,j,k,l)
                element.value_direction[2]=1; // This is the link along the direction 2 from the site of coordinate (i,j,k,l)
                element.value_direction[3]=1; // This is the link along the direction 3 from the site of coordinate (i,j,k,l)

                lattice.push_back(element);

            }
        }
    }
  }


}

/***************************************************************************************************************************/

//THIS FUNCTION INIZIALIZE THE LATTICE IN A RANDOM CONFIGURATION
void lattice4Random_creation(vector<site4> &lattice, int n, int nt){
    site4 element; 
  srand (time(NULL));
  int spin[2]={1,-1}, random;

  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
        for(int k=0; k<n; k++){
            for(int l=0; l<nt; l++){
                element.position[0]=i;
                element.position[1]=j;
                element.position[2]=k;
                element.position[3]=l;

                random=rand() % 2;
                element.value_direction[0]=spin[random]; // This is the link along the direction 0 from the site of coordinate (i,j,k,l)
                random=rand() % 2;
                element.value_direction[1]=spin[random]; // This is the link along the direction 1 from the site of coordinate (i,j,k,l)
                random=rand() % 2;
                element.value_direction[2]=spin[random]; // This is the link along the direction 2 from the site of coordinate (i,j,k,l)
                random=rand() % 2;
                element.value_direction[3]=spin[random]; // This is the link along the direction 3 from the site of coordinate (i,j,k,l)

                lattice.push_back(element);

            }
        }
    }
  }

}
/***************************************************************************************************************************/

// THIS FUNCTION INIZIALIZE THE PLAQUETTES IN THE FROZEN CONFIGURATION, I.E. ALL THE PLAQUETTES ARE +1
void plaquette4_creation(vector<plaquette4> &plaquettes, vector<site4> &lattice, int lattice_size, int n, int nt){
  
  plaquettes.clear();
  
  plaquette4 plaq;

  for(int i=0;i<lattice_size; i++){
    
    //Here we save in the variable plaquette the coordinates of the site in bottom left, could be useful for check
    plaq.position[0]=lattice[i].position[0];
    plaq.position[1]=lattice[i].position[1];
    plaq.position[2]=lattice[i].position[2];
    plaq.position[3]=lattice[i].position[3];
    
    plaq.value[0]=1; //Paquette on the plane 0-1
    plaq.value[1]=1; //Paquette on the plane 0-2
    plaq.value[2]=1; //Paquette on the plane 0-3
    plaq.value[3]=1; //Paquette on the plane 1-2
    plaq.value[4]=1; //Paquette on the plane 1-3
    plaq.value[5]=1; //Paquette on the plane 2-3

    plaquettes.push_back(plaq);
  }

}

/***************************************************************************************************************************/

// THIS FUNCTION INIZIALIZE THE PLAQUETTES FROM A RANDOM CONFIGURATION OF THE LATTICE
void plaquette4Random_creation(vector<plaquette4> &plaquettes, vector<site4> &lattice, int lattice_size, int n, int nt){
  
  plaquettes.clear();
  
  plaquette4 plaq;

  for(int i=0;i<lattice_size; i++){
    
    //Here we save in the variable plaquette the coordinates of the site in bottom left, could be useful for check
    plaq.position[0]=lattice[i].position[0];
    plaq.position[1]=lattice[i].position[1];
    plaq.position[2]=lattice[i].position[2];
    plaq.position[3]=lattice[i].position[3];

    int x = lattice[i].position[0];
    int y = lattice[i].position[1];
    int z = lattice[i].position[2];
    int t = lattice[i].position[3];
    
    //Periodic boundary conditions implemented using the modulo operator
    int xplus = (x + 1 ) % n;
    int yplus = (y + 1 ) % n;
    int zplus = (z + 1 ) % n;
    int tplus = (t + 1 ) % nt;
    
    //The links are saved in the vector lattice. 
    //If we are interested on the links attached to the site (x,y,z,t), 
    // we have to take the (x * n^2 * nt + y * n * nt + z * nt + t)-th element of the vector Lattice .
    // Here we take the positions of the link attached to the sites (x+1, y, z, t), (x, y+1, z, t), ... 
    xplus = ((xplus * n + y) * n + z) * nt + t; 
    yplus = ((x * n + yplus) * n + z) * nt + t; 
    zplus = ((x * n + y) * n + zplus) * nt + t; 
    tplus = ((x * n + y) * n + z) * nt + tplus;     

    //Plaquette on the plane 0-1
    plaq.value[0]= -lattice[i].value_direction[0]*lattice[i].value_direction[1] * lattice[xplus].value_direction[1]*lattice[yplus].value_direction[0];
    //Plaquette on the plane 0-2
    plaq.value[1]= -lattice[i].value_direction[0]*lattice[i].value_direction[2] * lattice[xplus].value_direction[2]*lattice[zplus].value_direction[0];
    //Plaquette on the plane 0-3
    plaq.value[2]= -lattice[i].value_direction[0]*lattice[i].value_direction[3] * lattice[xplus].value_direction[3]*lattice[tplus].value_direction[0];
    //Plaquette on the plane 1-2
    plaq.value[3]= -lattice[i].value_direction[1]*lattice[i].value_direction[2] * lattice[yplus].value_direction[2]*lattice[zplus].value_direction[1];
    //Plaquette on the plane 1-3
    plaq.value[4]= -lattice[i].value_direction[1]*lattice[i].value_direction[3] * lattice[yplus].value_direction[3]*lattice[tplus].value_direction[1];
    //Plaquette on the plane 2-3
    plaq.value[5]= -lattice[i].value_direction[2]*lattice[i].value_direction[3] * lattice[zplus].value_direction[3]*lattice[tplus].value_direction[2];
    
    //The plaquettes are stored on the vector "plaquettes", in this way, the plaquettes have the same rule to take the position 
    //of the links in the vector Lattice, if we consider the upper-left site.
    plaquettes.push_back(plaq);
  }

}

/***************************************************************************************************************************/

// THIS FUNCTION SUM ALL THE PLAQUETTES
double sum4_all_plaquettes (vector<plaquette4> &plaquettes, int plaquettes_size){
  
  double sum=0;
  
  for(int i=0;i<plaquettes_size;i++){
    sum+=(double)(plaquettes[i].value[0]+plaquettes[i].value[1]+plaquettes[i].value[2]+plaquettes[i].value[3]+plaquettes[i].value[4]+plaquettes[i].value[5]);
  }
  return sum/(double)(6.0*plaquettes_size);

}

/***************************************************************************************************************************/

// THIS FUNCTION IS THE MARKOV-CHAIN
void minimize4_energy(vector<plaquette4> &plaq, int plaquettes_size ,vector<site4> &lattice, int lattice_size, int n, int nt, double beta, int iterations, string fileE, string fileP, int &numberWilson, int Wilson_lenght){

  ofstream open_file;
  ofstream temporary_file;
  
  //Here is where the files where the data-points of the energy are open 
  temporary_file.open(fileE+".txt" );
  temporary_file<<"#--------------------Beta "<<beta<<"--------------------"<<endl;

  open_file.open ("thermalization.txt",ios::app);
  open_file<<endl;
  open_file<<"#--------------------Beta "<<beta<<"--------------------"<<endl;
  open_file<<endl;
  

  int spatial_distance, time_distance; //For Wilson loop
  string file;
  int xminus, yminus, zminus, tminus;
  double deltaE, Boltzmann, sum, sum1;
  srand(time(NULL));

  for(int k=1;k<=iterations;k++){
    
    sum=sum4_all_plaquettes (plaq, plaquettes_size);
    Polyakov(nt, lattice, lattice_size, n, nt, beta, fileP);

    temporary_file<<k<<" "<<sum<<" "<<beta<<endl;
    open_file<<k<<" "<<sum<<" "<<endl;
    
    for(int j=0;j<lattice_size;j++){
      int x=lattice[j].position[0];
      int y=lattice[j].position[1];
      int z=lattice[j].position[2];
      int t=lattice[j].position[3];

      //The links are saved in the vector lattice. If we are interested on the links attached to the site (x,y,z,t), 
      // we have to take the (x * n^2 * nt + y * n * nt + z * nt + t)-th element of the vector Lattice .
      // Here we take the positions of the links attached to the sites (x-1, y, z, t), (x, y-1, z, t), ... 
      // And the same for the plaquettes, where (x, y, z, t) is the coordinate of the upper left site of the plaquette.
      // The factor +n or +nt is to avoid negative numbers.
      xminus=(x-1+n)%n;
      xminus= ((xminus * n + y) * n + z) * nt + t;
      yminus=(y-1+n)%n;
      yminus= ((x * n + yminus) * n + z) * nt + t;
      zminus=(z-1+n)%n;
      zminus= ((x * n + y) * n + zminus) * nt + t;
      tminus=(t-1+nt)%nt;
      tminus= ((x * n + y) * n + z) * nt + tminus;


      /****************************************/
      //Here we take only the plaquettes attached on the link along the direction 0 from the site (x, y, z, t)
      deltaE=2.0*(plaq[j].value[0]+plaq[j].value[1]+plaq[j].value[2]+plaq[yminus].value[0]+plaq[zminus].value[1]+plaq[tminus].value[2]);
      
      
      //Here we extract the pseudo-casual number
      Boltzmann = drand48(); 
      Boltzmann = (double)rand() / RAND_MAX;
      if(Boltzmann < exp(-beta * deltaE)){
        lattice[j].value_direction[0]*=-1; 
        plaq[j].value[0]*=-1;
        plaq[j].value[1]*=-1;
        plaq[j].value[2]*=-1;
        plaq[yminus].value[0]*=-1;
        plaq[zminus].value[1]*=-1;
        plaq[tminus].value[2]*=-1;
      }

      /****************************************/

      deltaE=2.0*(plaq[j].value[0]+plaq[j].value[3]+plaq[j].value[4]+plaq[xminus].value[0]+plaq[zminus].value[3]+plaq[tminus].value[4]);
      
      Boltzmann = drand48();
      if(Boltzmann < exp(-beta * deltaE)){
        lattice[j].value_direction[1]*=-1; 
        plaq[j].value[0]*=-1;
        plaq[j].value[3]*=-1;
        plaq[j].value[4]*=-1;
        plaq[xminus].value[0]*=-1;
        plaq[zminus].value[3]*=-1;
        plaq[tminus].value[4]*=-1;
      }
      /****************************************/

      deltaE=2.0*(plaq[j].value[1]+plaq[j].value[3]+plaq[j].value[5]+plaq[xminus].value[1]+plaq[yminus].value[3]+plaq[tminus].value[5]);

      Boltzmann = drand48();
      if(Boltzmann < exp(-beta * deltaE)){ 
        lattice[j].value_direction[2]*=-1; 
        plaq[j].value[1]*=-1;
        plaq[j].value[3]*=-1;
        plaq[j].value[5]*=-1;
        plaq[xminus].value[1]*=-1;
        plaq[yminus].value[3]*=-1;
        plaq[tminus].value[5]*=-1;
      }

      /****************************************/
      deltaE=2.0*(plaq[j].value[2]+plaq[j].value[4]+plaq[j].value[5]+plaq[xminus].value[2]+plaq[yminus].value[4]+plaq[zminus].value[5]);
      
      Boltzmann = drand48();
      if(Boltzmann < exp(-beta * deltaE)){
        lattice[j].value_direction[3]*=-1; 
        plaq[j].value[2]*=-1;
        plaq[j].value[4]*=-1;
        plaq[j].value[5]*=-1;
        plaq[xminus].value[2]*=-1;
        plaq[yminus].value[4]*=-1;
        plaq[zminus].value[5]*=-1;
      }

    }

     for(spatial_distance=3 ; spatial_distance<5 ; spatial_distance++){ //Wilson loop
      
       for(time_distance=1; time_distance<Wilson_lenght; time_distance++){
           file="Wilson_data/wilson_data"+to_string(numberWilson)+to_string(time_distance)+to_string(spatial_distance)+".txt";
           Wilson(time_distance, spatial_distance, lattice, n, nt, beta, file);
         }

     }
  }
  numberWilson++;
  temporary_file.close();
  open_file.close();

}

/***************************************************************************************************************************/

void Wilson(int T, int R, vector<site4> &lattice, int n, int nt, double beta, string file) {
    
    int lattice_size = lattice.size();
    int num_loop=0;
    double sum_loop=0;
    ofstream temporary_file;
    temporary_file.open(file,ios::app);

    for (int j = 0; j < lattice_size; j++) {

        int x = lattice[j].position[0];
        int y = lattice[j].position[1];
        int z = lattice[j].position[2];
        int t = lattice[j].position[3];
        int xfin=x;
        int yfin=y;
        int zfin=z;
        int tfin=t;
        int site_index;
        double loopx_value = 1.0;
        double loopy_value = 1.0;
        double loopz_value = 1.0;

        for (int i = 0; i < R; i++) {

            site_index = xfin * n * n * nt + y * n * nt + z * nt + t;
            loopx_value *= lattice[site_index].value_direction[0];  

            site_index = x * n * n * nt + yfin * n * nt + z * nt + t;
            loopy_value *= lattice[site_index].value_direction[1];  

            site_index = x * n * n * nt + y * n * nt + zfin * nt + t;
            loopz_value *= lattice[site_index].value_direction[2];  

            xfin = (xfin + 1 + n) % n;  
            yfin = (yfin + 1 + n) % n; 
            zfin = (zfin + 1 + n) % n;  
        }

        for (int i = 0; i < T; i++) {
            
            site_index = x * n * n * nt + y * n * nt + z * nt + tfin;
            loopx_value *= lattice[site_index].value_direction[3];  
            loopy_value *= lattice[site_index].value_direction[3];  
            loopz_value *= lattice[site_index].value_direction[3];  

            site_index = xfin * n * n * nt + y * n * nt + z * nt + tfin;
            loopx_value *= lattice[site_index].value_direction[3]; 

            site_index = x * n * n * nt + yfin * n * nt + z * nt + tfin; 
            loopy_value *= lattice[site_index].value_direction[3];  

            site_index = x * n * n * nt + y * n * nt + zfin * nt + tfin;
            loopz_value *= lattice[site_index].value_direction[3]; 

            tfin = (tfin + 1 + nt) % nt;  
        }
        
        xfin = (xfin - 1 + n) % n; 
        yfin = (yfin - 1 + n) % n; 
        zfin = (zfin - 1 + n) % n;

        for (int i = 0; i < R; i++) {

            site_index = xfin * n * n * nt + y * n * nt + z * nt + tfin;
            loopx_value *= lattice[site_index].value_direction[0];

            site_index = x * n * n * nt + yfin * n * nt + z * nt + tfin;
            loopy_value *= lattice[site_index].value_direction[1];

            site_index = x * n * n * nt + y * n * nt + zfin * nt + tfin;
            loopz_value *= lattice[site_index].value_direction[2];

            xfin = (xfin - 1 + n) % n; 
            yfin = (yfin - 1 + n) % n; 
            zfin = (zfin - 1 + n) % n;

        }
        
        sum_loop+=(double)(loopx_value+loopy_value+loopz_value)/(double)(3);
        num_loop++;
        
    }
    temporary_file<<sum_loop/(double)(num_loop)<<" "<<T<<" "<<R<<" "<<beta<<endl;
    
    temporary_file.close();
}

/***************************************************************************************************************************/

void Polyakov(int T, vector<site4> &lattice, int lattice_size, int n, int nt, double beta, string file){
 
  double polyakov_loop=1.0, sum_polyakov=0.0;
  double polyakov1=1.0, polyakov2=1.0, polyakov3=1.0;
  int normalization=0;
  ofstream Polyakov_file;
  Polyakov_file.open(file,ios::app);

   for(int x=0; x<n; x++){
     for(int y=0; y<n; y++){
       for(int z=0; z<n; z++){
        
         for(int t=0; t<T; t++){
           int site_index= x * n * n * nt + y * n * nt + z * nt + t;
           polyakov_loop*=lattice[site_index].value_direction[3];
         }
        
         normalization++;
         sum_polyakov+=polyakov_loop;
         polyakov_loop=1.0;
        
       }
     }
   }

  Polyakov_file<<sum_polyakov/(double)(normalization)<<" "<<beta<<endl;
  Polyakov_file.close();

}


