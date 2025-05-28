# Monte Carlo simulation of the 4D Ising gauge model using Metropolis algorithm
Monte Carlo simulation based on Metropolis algorithm for the 4D Ising gauge model.

# 4D Ising gauge model
The Ising model is based on a set of discrete variables that take value ±1. The
most common Ising model is the Ising spin model, used to describe the the
magnetism in matter, where the variables are the spin values associated to the
sites of a lattice.
The Ising gauge model is a particular Ising model where the variables are de-
fined on the links between two sites.
In the 1970s, the Ising gauge model in 4D was proposed, due to its simplicity
and self-duality, as a tool to better understand lattice regularization of QFTs.
Unfortunately, the Ising gauge model exhibits a first-order phase transition,
that can cause problems of discontinuity in the continuum limit; thus it was
considered unsuitable for this purpose.
With advancements in technology, it is now possible to perform highly sophisticated Monte Carlo simulations, which are an essential tool in lattice calculations
where analytical solutions can’t be derived. This code has been developed to study the phase-transitions of the highly asymmetric 4D Ising gauge model, so
when one of the sides of the lattice is more shorter than the others.


## Metropolis algorithm
The Metropolis algorithm is one of the most known Monte Carlo method. 
The Metropolis algorithm is as follows:

1. Initialize the lattice in an arbitrary configuration with energy E1.
2. Pick one link and change its state, thus the new energy of the system is E2.
3. The new state is accepted with probability exp(-β(E2-E1)).
4. Steps 2-3 are repeated for all the links of the lattice. An iteration is when the code has attempted to switch all the links of the lattice.
5. Steps 2-4 are repeated until the set number of iterations is reached.

## Structure
- site4: This variable contain one site of the lattice with coordinates (x, y, z, t) stored in the vector position. In 4D, every site has 4 links 'forward', i.e. links that connect the site with the next first-neighbour ( for example (x, y, z, t)->(x+1 , y, z, t) along x-direction). The links are stored in the vector value_direction.
- plaquette4: This variable contain the plaquettes, which are the fundamental unity of the energy in the 4D Ising gauge model. The vector position contains the coordinates
  of the "upper left" site. Every site of the lattice, in 4D, has 6 plaquettes attached. Let's call the directions (x, y, z, t), the vector value contain the values of the plaquettes, respectively, along the planes xy, xz, xt, yz, yt, zt.
- lattice4_creation: This function initialize the lattice in the frozen-configuration, i.e. with all the links-value = +1.
- lattice4Random_creation: This function initialize the lattice in a random-configuration.
- plaquette4_creation: Initialize the plaquettes in the initial frozen-configuration.
- plaquette4Random_creation: Initialize the plaquettes in the initial random-configuration
- sum_all_plaquettes: This function sum all the plaquettes of the lattice.
- minimize4_energy: Perform the Metropolis algorithm.
- Wilson: Compute the rectangular Wilson loops. These loops have area RxT, with R = [3 , 4] and T = [1, 2, ..., Wilson_lenght]. In a lattice L^3 x Lt, where Lt is the small side, Wilson_lenght is L/2 if L is even, and (L-1)/2 if L is odd.
- Polyakov: Compute the Polyakov loop.

## Output
Three folders are generated from the simulation, one for the values of the density of interna energy (i.e. the average value of the plaquettes), one for the average value of the Polyakov loops and one for the Wilson loops. All the data-points are collected in the Metropolis process, thus must be uncorrelated using the blocking method. The python codes BlockingE.py, BlockingP.py and BlockingW.py perform the blocking method for, respectively, the energy, the Polyakov loop and the Wilson loops.

## How to launch the simulation
Open the terminal in folder where there is the code, so use the command


g++ -std=c++17 -stdlib=libc++ -o Gauge Gauge.cpp Gauge4_functions.cpp

Then, to start the simulation, use then command

./Gauge

The code will ask the following input:
- The size of the large sides of the lattice, identified as the space-like directions.
- The size of the smaller side of the lattice, identified as the time-like direction.
- The initial β-value.
- The final β-value.
- The range β-values where you are interested in taking more measurements. We will to this range as "high-density range".
- How many measurements you are interested in taking out of the high-density range. The suggestion is to take at least 2 measurements, in order to have the lattice properly thermalized when the range of interest is reached.
- How many measurements you are interested in taking in the high-density range.
- How many iterations are you interested to perform in the measurements out of the high-density range.
- How many iterations are you interested to perform in the high-density range.
- The initial configuration for the lattice.
