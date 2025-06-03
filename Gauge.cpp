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
using namespace std;

void MonteCarlo(vector< vector<int> > &map_nn, vector<int> &plaq, vector<site4> &lattice, int n, int nt, double &beta, double dbeta, int iterations, int Wilson_lenght, int number_file, int &numberWilson){
  
  cout<<beta<<endl;
  
  int plaquettes_size=plaq.size();
  int lattice_size=lattice.size();
  double energy, initial_energy;
  string fileE="Energy_data/energy_data"+to_string(number_file);
  string fileP="Polyakov_data/polyakov_data"+to_string(number_file)+".txt";

  initial_energy=sum4_all_plaquettes(plaq, plaquettes_size);
  minimize4_energy(map_nn, plaq, plaquettes_size,lattice, lattice_size, n, nt, beta, iterations, fileE, fileP, numberWilson, Wilson_lenght);
  energy=sum4_all_plaquettes(plaq, plaquettes_size);

  
  
  beta+=dbeta;


}


int main(){
  int removeFile=remove("thermalization.txt");
  int n, nt, StepsLow, StepsHigh, Wilson_lenght;
  
  cout<<"Lattice spatial size"<<endl;
  cin>>n;
  
  cout<<"Lattice time size"<<endl;
  cin>>nt;
  
  if(nt%2==0)
    Wilson_lenght=nt/2;
  else
    Wilson_lenght=(nt-1)/2;

  double beta;
  double beta_max, beta1, beta2;
  
  cout<<"Initial beta"<<endl;
  cin>>beta;

  cout<<"Final beta"<<endl;
  cin>>beta_max;

  cout<<"Range high density"<<endl;
  cin>>beta1;
  cin>>beta2;

  cout<<"How many steps in low density"<<endl;
  cin>>StepsLow;

  cout<<"How many steps in high density"<<endl;
  cin>>StepsHigh;

  double dbeta1, dbeta2, dbeta3;
  
  dbeta1=(beta1-beta)/(double)StepsLow;
  dbeta2=(beta2-beta1)/(double)StepsHigh;
  dbeta3=(beta_max-beta2)/(double)StepsLow;
  cout<<dbeta1<<" "<<dbeta2<<" "<<dbeta3<<endl;
  
  int iterations;
  cout<<"How many iterations in low density range?"<<endl;
  cin>>iterations;
  
  int iterationshigh;
  cout<<"How many iterations in high density range?"<<endl;
  cin>>iterationshigh;
  
  string choose;
  cout<<"Would you like a random initial lattice (R) o the frozen one (F) ?"<<endl;
  cin>>choose;
  
  vector<site4> lattice;
  vector< vector<int> > map_nn;
  vector<int> plaq;
  int lattice_size;

  if(choose=="R"){
      lattice4Random_creation(lattice, n, nt, map_nn);
      lattice_size=lattice.size();
      plaquette4_creation(plaq, lattice, lattice_size, map_nn);
  }

  else{
      lattice4_creation(lattice, n, nt, map_nn);
      lattice_size=lattice.size();
      plaquette4_creation(plaq, lattice, lattice_size, map_nn);
  }

  test_map(lattice, map_nn, n, nt);
  
  std::filesystem::create_directory("Energy_data");
  std::filesystem::create_directory("Wilson_data");
  std::filesystem::create_directory("Polyakov_data");

  int number_file=1;
  int numberWilson=1;

  while(beta<beta1){
    MonteCarlo(map_nn, plaq, lattice, n, nt, beta, dbeta1, iterations, Wilson_lenght, number_file, numberWilson);
    number_file++;
  }
    
  beta=beta1;

  while(beta<beta2){
    MonteCarlo(map_nn, plaq, lattice, n, nt, beta, dbeta2, iterationshigh, Wilson_lenght, number_file, numberWilson);
    number_file++;
  }

  beta=beta2;

  while(beta<=beta_max){
    MonteCarlo(map_nn, plaq, lattice, n, nt, beta, dbeta3, iterations, Wilson_lenght, number_file, numberWilson);
    number_file++;
  }

 
  
}
