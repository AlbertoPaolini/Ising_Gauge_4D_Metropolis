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

int Indexposition(int x, int y, int z, int t, int n, int nt){
  return ( (x+n) %n ) * n * n * nt + ( (y+n) %n ) * n * nt + ( (z+n) %n ) * nt + ( (t+nt) %nt );
}
/***************************************************************************************************************************/

// THIS FUNCTION INIZIALIZE THE LATTICE IN THE FROZEN CONFIGURATION, I.E. WITH ALL THE LINKS +1
void lattice4_creation(vector<site4> &lattice, int n, int nt, vector< vector<int> > &map_nn){ 
  
  site4 element; 
  int nn;

  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
        for(int k=0; k<n; k++){
            for(int l=0; l<nt; l++){
                vector<int> nn_position;
                element.position[0]=i;
                element.position[1]=j;
                element.position[2]=k;
                element.position[3]=l;

                element.value_direction[0]=1; // This is the link along the direction 0 from the site of coordinate (i,j,k,l)
                element.value_direction[1]=1; // This is the link along the direction 1 from the site of coordinate (i,j,k,l)
                element.value_direction[2]=1; // This is the link along the direction 2 from the site of coordinate (i,j,k,l)
                element.value_direction[3]=1; // This is the link along the direction 3 from the site of coordinate (i,j,k,l)

                lattice.push_back(element);

                // Here we initialize a matrix of vector that contain the position of the near-neighbours of each site
                nn = Indexposition( i +1, j, k, l, n, nt );
                nn_position.push_back(nn);
                nn = Indexposition( i, j +1, k, l, n, nt );
                nn_position.push_back(nn);
                nn = Indexposition( i, j, k +1, l, n, nt );
                nn_position.push_back(nn);
                nn = Indexposition( i, j, k, l +1, n, nt ); 
                nn_position.push_back(nn);

                nn = Indexposition( i -1, j, k, l, n, nt );
                nn_position.push_back(nn);
                nn = Indexposition( i, j -1, k, l, n, nt );
                nn_position.push_back(nn);
                nn = Indexposition( i, j, k -1, l, n, nt );
                nn_position.push_back(nn);
                nn = Indexposition( i, j, k, l -1, n, nt );
                nn_position.push_back(nn);

                map_nn.push_back(nn_position);

            }
        }
    }
  }


}

/***************************************************************************************************************************/

//THIS FUNCTION INIZIALIZE THE LATTICE IN A RANDOM CONFIGURATION
void lattice4Random_creation(vector<site4> &lattice, int n, int nt, vector< vector<int> > &map_nn){
    site4 element; 
    int nn;
  srand (time(NULL));
  int spin[2]={1,-1}, random;

  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
        for(int k=0; k<n; k++){
            for(int l=0; l<nt; l++){
                vector<int> nn_position;
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

                // Here we initialize a matrix of vector that contain the position of the near-neighbours of each site
                nn = Indexposition( i +1, j, k, l, n, nt );
                nn_position.push_back(nn);
                nn = Indexposition( i, j +1, k, l, n, nt );
                nn_position.push_back(nn);
                nn = Indexposition( i, j, k +1, l, n, nt );
                nn_position.push_back(nn);
                nn = Indexposition( i, j, k, l +1, n, nt );
                nn_position.push_back(nn);

                nn = Indexposition( i -1, j, k, l, n, nt );
                nn_position.push_back(nn);
                nn = Indexposition( i, j -1, k, l, n, nt );
                nn_position.push_back(nn);
                nn = Indexposition( i, j, k -1, l, n, nt );
                nn_position.push_back(nn);
                nn = Indexposition( i, j, k, l -1, n, nt );
                nn_position.push_back(nn);

                map_nn.push_back(nn_position);

            }
        }
    }
  }

}

/***************************************************************************************************************************/

// THIS FUNCTION INIZIALIZE THE PLAQUETTES OF THE LATTICE
void plaquette4_creation(vector<int> &plaquettes, vector<site4> &lattice, int lattice_size, vector< vector<int> > &map_nn){
  
  plaquettes.clear();
  int position_plaquette = 0;
  int plaq;

  for(int i=0;i<lattice_size; i++){
    int xplus = map_nn[i][0];
    int yplus = map_nn[i][1];
    int zplus = map_nn[i][2];
    int tplus = map_nn[i][3];

    //Plaquette on the plane 0-1
    plaq = lattice[i].value_direction[0]*lattice[i].value_direction[1] * lattice[xplus].value_direction[1]*lattice[yplus].value_direction[0];
    lattice[i].plaquette_attached[0].push_back(position_plaquette); // Vector lattice[i].plaquette_attached[0] will contain the position in the vector "plaquette" of the plaquettes connected at the i-th link in the direction 0
    lattice[yplus].plaquette_attached[0].push_back(position_plaquette);
    lattice[i].plaquette_attached[1].push_back(position_plaquette);
    lattice[xplus].plaquette_attached[1].push_back(position_plaquette);
    plaquettes.push_back(plaq);
    position_plaquette ++;

    //Plaquette on the plane 0-2
    plaq = lattice[i].value_direction[0]*lattice[i].value_direction[2] * lattice[xplus].value_direction[2]*lattice[zplus].value_direction[0];
    lattice[i].plaquette_attached[0].push_back(position_plaquette);
    lattice[zplus].plaquette_attached[0].push_back(position_plaquette);
    lattice[i].plaquette_attached[2].push_back(position_plaquette);
    lattice[xplus].plaquette_attached[2].push_back(position_plaquette);
    plaquettes.push_back(plaq);
    position_plaquette ++;

    //Plaquette on the plane 0-3
    plaq = lattice[i].value_direction[0]*lattice[i].value_direction[3] * lattice[xplus].value_direction[3]*lattice[tplus].value_direction[0];
    lattice[i].plaquette_attached[0].push_back(position_plaquette);
    lattice[tplus].plaquette_attached[0].push_back(position_plaquette);
    lattice[i].plaquette_attached[3].push_back(position_plaquette);
    lattice[xplus].plaquette_attached[3].push_back(position_plaquette);
    plaquettes.push_back(plaq);
    position_plaquette ++;

    //Plaquette on the plane 1-2
    plaq = lattice[i].value_direction[1]*lattice[i].value_direction[2] * lattice[yplus].value_direction[2]*lattice[zplus].value_direction[1];
    lattice[i].plaquette_attached[1].push_back(position_plaquette);
    lattice[zplus].plaquette_attached[1].push_back(position_plaquette);
    lattice[i].plaquette_attached[2].push_back(position_plaquette);
    lattice[yplus].plaquette_attached[2].push_back(position_plaquette);
    plaquettes.push_back(plaq);
    position_plaquette ++;

    //Plaquette on the plane 1-3
    plaq = lattice[i].value_direction[1]*lattice[i].value_direction[3] * lattice[yplus].value_direction[3]*lattice[tplus].value_direction[1];
    lattice[i].plaquette_attached[1].push_back(position_plaquette);
    lattice[tplus].plaquette_attached[1].push_back(position_plaquette);
    lattice[i].plaquette_attached[3].push_back(position_plaquette);
    lattice[yplus].plaquette_attached[3].push_back(position_plaquette);
    plaquettes.push_back(plaq);
    position_plaquette ++;

    //Plaquette on the plane 2-3
    plaq = lattice[i].value_direction[2]*lattice[i].value_direction[3] * lattice[zplus].value_direction[3]*lattice[tplus].value_direction[2];
    lattice[i].plaquette_attached[2].push_back(position_plaquette);
    lattice[tplus].plaquette_attached[2].push_back(position_plaquette);
    lattice[i].plaquette_attached[3].push_back(position_plaquette);
    lattice[zplus].plaquette_attached[3].push_back(position_plaquette);
    plaquettes.push_back(plaq); //Ciclo for per plachette
    position_plaquette ++;
    
  }

}

/***************************************************************************************************************************/

// THIS FUNCTION GENERATE A FILE TO DEMONSTRATE THAT map_nn CONTAIN THE CORRECT NEAR-NEIGHBOUR
void test_map(vector<site4> &lattice, vector< vector<int> > &map_nn, int n, int nt){
  ofstream geometry_test;
  geometry_test.open("Geometry_Test.txt");
  for(int i=0; i<lattice.size(); i++){
    geometry_test<<"----------------------------------------------------------------------------------------------------------------------------------"<<endl;
    geometry_test<<"Lattice "<<n<<"^3 x "<<nt<<endl;
    geometry_test<<"("<<lattice[i].position[0]<<", "<<lattice[i].position[1]<<", "<<lattice[i].position[2]<<", "<<lattice[i].position[3]<<")-> Along x -> ";
    geometry_test<<"("<<lattice[map_nn[i][0]].position[0]<<", "<<lattice[map_nn[i][0]].position[1]<<", "<<lattice[map_nn[i][0]].position[2]<<", "<<lattice[map_nn[i][0]].position[3]<<") ";
    geometry_test<<"("<<lattice[map_nn[i][4]].position[0]<<", "<<lattice[map_nn[i][4]].position[1]<<", "<<lattice[map_nn[i][4]].position[2]<<", "<<lattice[map_nn[i][4]].position[3]<<") "<<endl;

    geometry_test<<"("<<lattice[i].position[0]<<", "<<lattice[i].position[1]<<", "<<lattice[i].position[2]<<", "<<lattice[i].position[3]<<")-> Along y -> ";
    geometry_test<<"("<<lattice[map_nn[i][1]].position[0]<<", "<<lattice[map_nn[i][1]].position[1]<<", "<<lattice[map_nn[i][1]].position[2]<<", "<<lattice[map_nn[i][1]].position[3]<<") ";
    geometry_test<<"("<<lattice[map_nn[i][5]].position[0]<<", "<<lattice[map_nn[i][5]].position[1]<<", "<<lattice[map_nn[i][5]].position[2]<<", "<<lattice[map_nn[i][5]].position[3]<<") "<<endl;

    geometry_test<<"("<<lattice[i].position[0]<<", "<<lattice[i].position[1]<<", "<<lattice[i].position[2]<<", "<<lattice[i].position[3]<<")-> Along z -> ";
    geometry_test<<"("<<lattice[map_nn[i][2]].position[0]<<", "<<lattice[map_nn[i][2]].position[1]<<", "<<lattice[map_nn[i][2]].position[2]<<", "<<lattice[map_nn[i][2]].position[3]<<") ";
    geometry_test<<"("<<lattice[map_nn[i][6]].position[0]<<", "<<lattice[map_nn[i][6]].position[1]<<", "<<lattice[map_nn[i][6]].position[2]<<", "<<lattice[map_nn[i][6]].position[3]<<") "<<endl;

    geometry_test<<"("<<lattice[i].position[0]<<", "<<lattice[i].position[1]<<", "<<lattice[i].position[2]<<", "<<lattice[i].position[3]<<")-> Along t -> ";
    geometry_test<<"("<<lattice[map_nn[i][3]].position[0]<<", "<<lattice[map_nn[i][3]].position[1]<<", "<<lattice[map_nn[i][3]].position[2]<<", "<<lattice[map_nn[i][3]].position[3]<<") ";
    geometry_test<<"("<<lattice[map_nn[i][7]].position[0]<<", "<<lattice[map_nn[i][7]].position[1]<<", "<<lattice[map_nn[i][7]].position[2]<<", "<<lattice[map_nn[i][7]].position[3]<<") "<<endl;
    geometry_test<<endl;
  }
  geometry_test.close();
}

/***************************************************************************************************************************/

// THIS FUNCTION SUM ALL THE PLAQUETTES
double sum4_all_plaquettes (vector<int> &plaquettes, int plaquettes_size){
  
  double sum=0;
  
  for(int i=0;i<plaquettes_size;i++){
    sum += (double)( plaquettes[i] );
  }
  return sum/(double)(plaquettes_size);

}


/***************************************************************************************************************************/

// THIS FUNCTION IS THE MARKOV-CHAIN
void minimize4_energy(vector< vector<int> > &map_nn, vector<int> &plaq, int plaquettes_size ,vector<site4> &lattice, int lattice_size, int n, int nt, double beta, int iterations, string fileE, string fileP, int &numberWilson, int Wilson_lenght){
  ofstream open_file;
  ofstream temporary_file;
  
  //Here is where the files where the data-points of the energy are open 
  temporary_file.open(fileE+".txt" );
  temporary_file<<"#--------------------Beta "<<beta<<"--------------------"<<endl;

  open_file.open ("thermalization.txt",ios::app);
  open_file<<endl;
  open_file<<"#--------------------Beta "<<beta<<"--------------------"<<endl;
  open_file<<endl;
  
  double deltaE, Boltzmann, sum;
  srand(time(NULL));

  for(int k=1;k<=iterations;k++){
    
    sum=sum4_all_plaquettes (plaq, plaquettes_size);
    Polyakov(lattice, lattice_size, n, nt, beta, fileP);

    temporary_file<<k<<" "<<sum<<" "<<beta<<endl;
    open_file<<k<<" "<<sum<<" "<<endl;
    
    for(int j=0;j<lattice_size;j++){

      for(int direction = 0; direction < 4; direction++){
        deltaE = 0.0;

        for(int posPlaq = 0; posPlaq < lattice[j].plaquette_attached[direction].size(); posPlaq++){
          deltaE += plaq[lattice[j].plaquette_attached[direction][posPlaq]];
        }
        deltaE *= 2.0;
        Boltzmann = drand48(); 
        if(Boltzmann < exp(-beta * deltaE)){
          lattice[j].value_direction[direction]*=-1;
          for(int posPlaq = 0; posPlaq < lattice[j].plaquette_attached[direction].size(); posPlaq++){
            plaq[lattice[j].plaquette_attached[direction][posPlaq]] *= -1;
          }
        }
      }

      // for(int R=3; R<4; R++){
      //   for(int T=2; T< Wilson_lenght; T++){
      //     Wilson( T, R, lattice, map_nn, beta, numberWilson );
      //   }
      // }

     }
  }
  numberWilson++;
  temporary_file.close();
  open_file.close();

}

/***************************************************************************************************************************/

void Wilson(int T, int R, vector<site4> &lattice, vector< vector<int> > &map_nn, double beta, int number) {
    
    int lattice_size = lattice.size();
    int num_loop=0;
    double sum_loop=0;
    ofstream temporary_file;
    temporary_file.open("Wilson_data/Wilson_"+to_string(number)+".txt",ios::app);

    for (int j = 0; j < lattice_size; j++) {

        int x_link_position = j;
        int y_link_position = j;
        int z_link_position = j;
        int t_link_position = j;
        double loopx_value = 1.0;
        double loopy_value = 1.0;
        double loopz_value = 1.0;

        for (int i = 0; i < R-1; i++) {

            loopx_value *= lattice[x_link_position].value_direction[0];  
            loopy_value *= lattice[y_link_position].value_direction[1];  
            loopz_value *= lattice[z_link_position].value_direction[2];  

            x_link_position = map_nn[x_link_position][0];
            y_link_position = map_nn[y_link_position][1];
            z_link_position = map_nn[z_link_position][2];  
        }

        for (int i = 0; i < T-1; i++) {
            
            loopx_value *= lattice[t_link_position].value_direction[3];  
            loopy_value *= lattice[t_link_position].value_direction[3];  
            loopz_value *= lattice[t_link_position].value_direction[3];  

            loopx_value *= lattice[x_link_position].value_direction[3]; 
            loopy_value *= lattice[y_link_position].value_direction[3];  
            loopz_value *= lattice[z_link_position].value_direction[3]; 

            t_link_position = map_nn[t_link_position][3]; 
            x_link_position = map_nn[x_link_position][3];
            y_link_position = map_nn[y_link_position][3];
            z_link_position = map_nn[z_link_position][3]; 
        }

        for (int i = 0; i < R-1; i++) {

            loopx_value *= lattice[x_link_position].value_direction[0];
            loopy_value *= lattice[y_link_position].value_direction[1];
            loopz_value *= lattice[z_link_position].value_direction[2];

            x_link_position = map_nn[x_link_position][4]; 
            y_link_position = map_nn[y_link_position][5]; 
            z_link_position = map_nn[z_link_position][6];

        }
        
        sum_loop+=(double)(loopx_value+loopy_value+loopz_value)/(double)(3);
        num_loop++;
        
    }
    temporary_file<<sum_loop/(double)(num_loop)<<" "<<T<<" "<<R<<" "<<beta<<endl;
    
    temporary_file.close();
}

/***************************************************************************************************************************/

void Polyakov(vector<site4> &lattice, int lattice_size, int n, int nt, double beta, string file){
 
  double polyakov_loop=1.0, sum_polyakov=0.0;
  int normalization=0;
  ofstream Polyakov_file;
  Polyakov_file.open(file,ios::app);

  for(int x=0; x<n; x++){
    for(int y=0; y<n; y++){
      for(int z=0; z<n; z++){
        
        for(int t=0; t<nt-1; t++){
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


