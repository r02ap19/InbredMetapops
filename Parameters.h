#pragma once

#include <fstream>

using namespace std;

class Parameters
{
public:
	Parameters();
	~Parameters();
	//Simulation
	int SimNr;
	int rep; //replicates
	int gen; //generations
	int out_int; //generations interval for population and trait outputs 
	int out_start;
	int ind_interval; //generations interval for individual outputs 
	int indRand_interval; //generation interval for random mating and output for ID
	int PopMut_interval;
	bool out_mutations;

	//Landscape
	bool env_stoch;
	int resol; //resolution (m)
	int x_max;
	int y_max;
	double prop_suitable; //proportion of suitable cells 
	double K; //carrying capacity (individuals/ha)

	int start_dispEvol;

	//Initialisation
	int min_seedX;
	int min_seedY;
	int max_seedX;
	int max_seedY;

	//Traits
	bool nearest; //TRUE = neartest neighbour dispersal 
	int L; //number of loci for each trait
	int dispEvol; //0 = no dispersal evolution; 1 = emigration p. evolving; 2 = disp. distance evolving; 3 = both evolving
	double* genot_mean; // initial genotypic mean for each trait (delta_tau, ep, dist, post_ep, alt)
	double* genot_std; // initial genotypic standard deviation for each trait
	double* genot_std2; //gentotypic standard deviation for when traits start evolving -AFFECTS MUTATION SD
	double mu; // haploid per allele mutation rate for each trait
	double mu_neutral;
	double* mu_std; // standard deviation for mutational effects for each trait (delta_tau, ep, dist)
	double rec_probability; //recombination probability for adaptive loci

	int nL; //number of neutral loci
	//Deleterious mutations
	int loadEffect; //0 = affects offsping survival; 1 = affects fecundity and fertilization or mating probability in males
	int initial_nMut; //initial no. of deleterious mutations
	double R; //genome map length (see Roze & Rousset 2009, JEB) - corresponds to recombination rate
	double Ud; //mutation rate for mildly deleterious mutations / diploid genome / generation (see Spigler et al. 2016, Evolution)
	double Ul; //mutation rate for lethal mutations / diploid genome / generation (see Spigler et al. 2016, Evolution)
	double Ub; //mutation rate for back mutations / diploid genome / generation 
	double mean_sd; //mean selection coefficient for mildly deleterious mutations
	double mean_hd; //mean dominance coefficient for mildly deleterious mutations
	double sl; //selection coefficient for lethal mutations (fix - not sampled from a distribution)	
	double hl; //dominance coefficient for lethal mutations
	double sb; //selection coefficient for beneficial mutations

	//DFE simulations
	int sh_dist; // 0 = s,h sampled as in Spigler et al 2016; 1 = s sampled from unireal (0.0,1.0), h sampled from unireal (0.0,0.); 2 = as 0, but with distributions for s and h switched. 3 = s sampled uniform, h = 0.5. 4 = s sampled from uniform, h = 0.2. 5 mean s=h=0 for all mutations.
	
	//Reproduction & survival 
	bool stoch_rep; //stochastic reproduction
	double fec; //mean fecundity

	//Dispersal
	double dispCost; //paid in terms of survival

	void outPara(string dir); //Parameter output
};