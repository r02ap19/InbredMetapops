#include "Parameters.h"

Parameters::Parameters() {
	//Simulation
	SimNr = 10080;
	rep = 1; //replicates
	gen = 200001; //generations
	out_int = 200; //generation interval for output (population and trait)
	out_start = 0; //output start generation
	ind_interval = 200; //generations interval for individual outputs 
	indRand_interval = 200;
	PopMut_interval = 200;
	out_mutations = 0;

	//Landscape
	env_stoch = false; //environmental stochasticity
	resol = 100; //resolution (m)
	x_max = 10;
	y_max = 20;
	prop_suitable = 1.0; //proportion of suitable cells 
	K = 50; //carrying capacity (individuals/ha)

	start_dispEvol = 0; //start dispersal evolution

	//Initialisation
	min_seedX = 0;
	min_seedY = 0;
	max_seedX = x_max;
	max_seedY = y_max;

	//Traits
	nearest = false; //nearest-neighbour dispersal
	L = 20; //number of loci for each trait
	dispEvol = 1; //0 = no dispersal evolution; 1 = emigration p. evolving; 2 = disp. distance evolving 3 = both evolving
	genot_mean = new double[2]{0.05, 200.0}; // initial genotypic mean for each trait (a, ep, dist, post-mating dispersal probability, probability of exhibiting altruistic behaviour)
	genot_std = new double[2]{0.1, 0.00000000001}; // initial genotypic standard deviation for each trait
	genot_std2 = new double[2]{0.1,  0.00000000001};
	mu = 0.001; // haploid per locus mutation rate for each trait
	mu_neutral = 0.001; // haploid per locus mutation probability for neutral loci
	mu_std = new double[2]{0.05, 0.0}; // standard deviation for mutational effects for each trait (a, ep, dist)

	rec_probability = 0.1; //recombination probability for adaptive loci

	nL = 500; //number of neutral loci

	//Deleterious mutations
	loadEffect = 0; //0 = affects offsping survival; 1 = affects fecundity and fertilization or mating probability in males; -9 = fix inbreeding depression hard-coded in PolyDisp.cpp
	initial_nMut = 0;
	R = 10.0; //genome (see Roze & Rousset 2009, JEB) - corresponds to recombination rate
	Ud = 0.1; //mutation rate for mildly deleterious mutations / diploid genome / generation (see Spigler et al. 2016, Evolution)
	Ul = 0.00000000000000001; //mutation rate for lethal mutations / diploid genome / generation (see Spigler et al. 2016, Evolution)
	Ub = 0.00000000000000001; //mutation rate for back mutations
	mean_sd = 0.05; //mean selection coefficient for mildly deleterious mutations
	mean_hd = 0.36; //mean dominance coefficient for mildly deleterious mutations
	sl = 1.0; //selection coefficient for lethal mutations (fix - not sampled from a distribution)
	hl = 0.02; //mean dominance coefficient for lethal mutations
	sb = -0.005;

	//Reproduction & survival 
	stoch_rep = true; //stochastic reproduction 
	fec = 12.0; //mean fecundit

	//Dispersal
	dispCost = 0.7; //paid in terms of survival prob. Can be subject to incrementation.

	//DFE Simulations
	sh_dist = 0;

}
//------------------------------
void Parameters::outPara(string name) {

	ofstream out;

	out.open(name.c_str());

	//Simulation
	out << "SimNr\t" << SimNr << endl;
	out << "rep\t" << rep << endl; //replicates
	out << "gen\t" << gen << endl; //generations
	out << "out_int\t" << out_int << endl;
	out << "out_start\t" << out_start << endl;
	out << "out_ind_interval\t" << ind_interval << endl;
	out << "indRand_interval\t" << indRand_interval << endl;
	out << "PopMut_interval\t" << PopMut_interval << endl;
	out << "out_mutations\t" << out_mutations << endl;

	//Landscape
	out << "envir_stochasticity\t" << env_stoch << endl;
	out << "resolution\t" << resol << endl; //resolution (m)
	out << "x_max\t" << x_max << endl;
	out << "y_max\t" << y_max << endl;
	out << "proportion_suitable_cells\t" << prop_suitable << endl; //proportion of suitable cells 	
	out << "K\t" << K << endl; //carrying capacity (individuals/ha)
	out << "dispersal_starts_evolving_at\t" << start_dispEvol << endl;

	//Initialisation
	out << "min_seedX\t" << min_seedX << endl;
	out << "min_seedY\t" << min_seedY << endl;
	out << "max_seedX\t" << max_seedX << endl;
	out << "max_seedY\t" << max_seedY << endl;

	//Traits
	out << "nearest-neigh\t" << nearest << endl;
	out << "nLoci_per_trait\t" << L << endl;
	out << "dispersal_evolving\t" << dispEvol << endl;
	out << "mean_emig_p\t" << genot_mean[0] << endl;
	out << "mean_disp_dist\t" << genot_mean[1] << endl;
	out << "sd_emig_p\t" << genot_std[0] << endl;
	out << "sd_disp_dist\t" << genot_std[1] << endl;
	out << "sd2_emig_p\t" << genot_std2[0] << endl;
	out << "mu\t" << mu << endl;
	out << "mu_neutral\t" << mu_neutral << endl;
	out << "mut_sd_emig_p\t" << mu_std[0] << endl;
	out << "mut_sd_disp_dist\t" << mu_std[1] << endl;
	out << "recombination\t" << rec_probability << endl;

	out << "nNeutral_loci\t" << nL << endl;
	//Deleterious mutations
	out << "load_effect\t" << loadEffect << endl;
	out << "initial_Nmut\t" << initial_nMut << endl;
	out << "genome_map_length\t" << R << endl;
	out << "Ud\t" << Ud << endl;
	out << "Ul\t" << Ul << endl;
	out << "Ub\t" << Ub << endl;
	out << "mean_selection_mild_mut\t" << mean_sd << endl;
	out << "mean_dominance_mild_mut\t" << mean_hd << endl;
	out << "selection_lethal_mut\t" << sl << endl;
	out << "dominance_lethal_mut\t" << hl << endl;
	out << "selection_beneficial\t" << sb << endl;
	out << "sh_dist\t" << sh_dist << endl;

	//Reproduction & survival 
	out << "fecundity\t" << fec << endl;
	out << "cost_dispersal\t" << dispCost << endl;
	out << "stoch_rep\t" << stoch_rep << endl;
	
	out.close();
}
//------------------------------
Parameters::~Parameters()
{
}