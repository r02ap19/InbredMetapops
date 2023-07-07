#include "PolyDisp.h"

//---------------------------------------------------------------------------

// Main function
#if LINUX
int main(int argc, char* argv[])
{
	// Get the current directory.
	char* buffer = getcwd(NULL, 0);
	dir = buffer;
	free(buffer);
	dir = dir + "/"; //Current directory path
	dirOut = dir + "Outputs/"; //Outpus folder path

	
	//para.SimNr = std::atoi(argv[1]);


	RunModel();

	cout << "Simulation completed" << endl;

	return 0;
}
#else
int _tmain(int argc, _TCHAR* argv[])
{
	
	// Get the current directory.
	char* buffer = _getcwd(NULL, 0);
	dir = buffer;
	free(buffer);
	dir = dir + "\\"; //Current directory path
	dirOut = dir + "Outputs\\"; //Outpus folder path

	extime = clock();

	RunModel();

	std::cout << "Simulation completed" << endl;

	return 0;
}
#endif

//---------------------------------------------------------------------------
const string Int2Str(const int x)
{
	ostringstream o;
	if (!(o << x)) return "ERROR";
	return o.str();
}
//---------------------------------------------------------------------------
const string Float2Str(const double x)
{
	ostringstream o;
	if (!(o << x)) return "ERROR";
	return o.str();
}
//---------------------------------------------------------------------------
void RunModel(void) {
	std::cout << "Simulation nr. " << para.SimNr << endl;

	para.dispCost = fluc(rdgen);

	para.mu_std[0] = para.genot_std[0] / sqrt(2.0 * (double)para.L);
	para.mu_std[1] = para.genot_std[1] / sqrt(2.0 * (double)para.L);


	string name = dirOut + "Sim" + Int2Str(para.SimNr) + "_Para.txt";
	para.outPara(name);


	if (para.out_int > 0) {
		outPop_header();
	}
	if (para.ind_interval > 0) outInds_header();


	if (para.indRand_interval > 0) outRanInds_header();
	if (para.PopMut_interval > 0) outPopMut_header();

	//Loop through REPLICATES --------------------------------------------
	for (r = 0; r < para.rep; r++) {
		std::cout << "rep = " << r << "==================" << endl;

		if (para.out_mutations) outMut_header();

		//create a new landscape
		landscape(para.x_max, para.y_max, para.prop_suitable);

		//Initialisation
		initialisation(para.x_max, para.y_max, para.min_seedX, para.min_seedY, para.max_seedX, para.max_seedY);
		//cout << "initialisation ok" << endl;

		mat_trials = false;

		//Loop through GENERATIONS ------------------------------------
		for (g = 0; g < para.gen; g++) {

			if (g % 50 == 0) {
				std::cout << "gen = " << g << endl;
				extime = clock() - extime;
				std::cout << "time = " << (float)extime / CLOCKS_PER_SEC << " sec" << endl;
				extime = clock();
			}

			//update traits parameters
			if (g == para.start_dispEvol) {
				if (para.dispEvol == 1 || para.dispEvol == 3) para.mu_std[0] = para.genot_std2[0] / sqrt(2.0 * (double)para.L);
				if (para.dispEvol > 1) para.mu_std[1] = para.genot_std2[1] / sqrt(2.0 * (double)para.L);
			}

			//Reproduction
			reproduction();

			//Dispersal 
			dispersal();

			//Density-dependent survival
			if (para.stoch_rep) survival();

			if (g > 0 && (g % para.indRand_interval == 0)) {
				betweenPop_mating();
				mat_trials = true;
			}
			if (Ntot < 2) break;
		}

		//Delete populations and the landscape
		for (int x = 0; x < para.x_max; x++) {
			for (int y = 0; y < para.y_max; y++) {
				if (pop[x][y] != NULL) {
					pop[x][y]->deleteAdults();
					delete pop[x][y]; pop[x][y] = NULL;
				}
				delete land[x][y]; land[x][y] = NULL;
			}
			delete[] pop[x]; pop[x] = NULL;
			delete[] land[x]; land[x] = NULL;
		}
		delete[] pop; pop = NULL;
		delete[] land; land = NULL;

		if (muts.is_open()) muts.close();
	}

	if (pops.is_open()) pops.close();
	if (trait.is_open()) trait.close();
	if (inds.is_open()) inds.close();
	if (ranInd.is_open()) ranInd.close();
	if (popmut.is_open()) popmut.close();
	if (kinship.is_open()) kinship.close();
}
//---------------------------------------------------------------------------
void landscape(int xmax, int ymax, double prop) {

	int x, y, i, ncells;
	std::uniform_int_distribution<> ranx(0, xmax - 1);
	std::uniform_int_distribution<> rany(0, ymax - 1);

	land = new Landscape * *[xmax];
	for (int j = 0; j < xmax; j++) {
		land[j] = new Landscape * [ymax];
		for (int jj = 0; jj < ymax; jj++) {
			land[j][jj] = new Landscape();
		}
	}

	pop = new Population * *[xmax];
	for (int j = 0; j < xmax; j++) {
		pop[j] = new Population * [ymax];
		//do not intialise the single pop. - do it only when there is a population
		for (int jj = 0; jj < ymax; jj++) pop[j][jj] = NULL;
	}

	ncells = (int)(xmax * ymax * prop); //nr. of suitable cells 

	//Distribute suitable cells randomly
	i = 0;
	do {
		do {
			x = ranx(rdgen);
			y = rany(rdgen);
		} while (land[x][y]->suitable == 1);
		land[x][y]->suitable = 1;
		land[x][y]->local_K = para.K;
		i++;
	} while (i < ncells);

}
//---------------------------------------------------------------------------
void initialisation(int xmax, int ymax, int minSx, int minSy, int maxSx, int maxSy) {

	//Normal distributions for traits initialisation
	std::normal_distribution<> distrEmig(para.genot_mean[0] / (2.0 * (double)para.L), para.genot_std[0] / sqrt(2.0 * (double)para.L));
	std::normal_distribution<> distrDist(para.genot_mean[1] / (2.0 * (double)para.L), para.genot_std[1] / sqrt(2.0 * (double)para.L));

	//Inintialise populations in all the suitable cells in the defined area
	for (int x = minSx; x < maxSx; x++) {
		for (int y = minSy; y < maxSy; y++) {
			if (land[x][y]->suitable) {
				pop[x][y] = new Population(x, y);
				pop[x][y]->initialise_pop(land[x][y]->local_K, k, para, distrEmig, distrDist, s_mild, neutral, position);
			}
		}
	}

	cy_max = para.max_seedY;
	cy_min = 0;
}

void reproduction(void) {
	int n_offs;
	int n_mut;
	int m, dad;
	int a_males; //number of available mates
	int n_offs_b; //number of offspring before any selection to calculate frequency of deleterious mutations
	double p_surv; //survival probability 

	Individuals* ind;
	vector<int> a_mates; //vector of indeces of available males
	vector<Individuals>::iterator iter;

	//Normal distributions for trait mutations
	std::normal_distribution<> MutEm(0.0, para.mu_std[0]);
	std::normal_distribution<> MutDist(0.0, para.mu_std[1]);

	std::poisson_distribution<> n_mildmut(para.Ud);
	std::poisson_distribution<> n_lethal(para.Ul);
	std::poisson_distribution<> n_backmut(para.Ub);

	//Loop through the landscape
	for (int x = 0; x < para.x_max; x++) {
		for (int y = cy_min; y < cy_max; y++) {

			if (pop[x][y] != NULL) {

				pop[x][y]->Noffs = 0;
				pop[x][y]->Foffs = 0;
				pop[x][y]->Moffs = 0;

				//clear Map of population's deleterious mutations
				if (!pop[x][y]->popMuts.empty()) pop[x][y]->popMuts.clear();
				n_offs_b = 0;

				if (pop[x][y]->Nf > 0 && pop[x][y]->Nm > 0) {

					//Loop through females for mating and survival
					for (iter = pop[x][y]->females.begin(); iter != pop[x][y]->females.end(); iter++) {
						p_surv = 1.0; //survival probability

						if (iter->reproduce) {
							//MATING-------------------------------------------------------------								
							if (pop[x][y]->Nm > 1) {
								//available males
								a_males = pop[x][y]->Nm;
								for (int i = 0; i < pop[x][y]->Nm; i++) {
									if (pop[x][y]->males[i].reproduce) a_mates.push_back(i);
									else a_males--;
								}

								if (a_males > 0) {
									//sample the first male
									std::uniform_int_distribution<> mat(0, a_males - 1);
									m = mat(rdgen);
									iter->mates.push_back(a_mates[m]);
									iter->n_mates++;
									a_males--;
									a_mates.erase(a_mates.begin() + m);

								}
							}
							else {
								if (pop[x][y]->males[0].reproduce) {
									//mate with the only male available
									iter->mates.push_back(0);
									iter->n_mates++;
								}
							}

							if (iter->n_mates > 0) {
								if (iter->alive) {

									//produce 10 offs from random mating within-population for ID and heterosis calculations
									if (mat_trials) {
										//father sampling distribution
										std::uniform_int_distribution<> sdad(0, pop[x][y]->Nm - 1);

										for (int i = 0; i < 10; i++) {
												//cout << "doing extra-mating " << endl;
												ind = new Individuals(para, Bern(rdgen), x, y);
												dad = sdad(rdgen);

												//inheritance
												inheritance(ind, *iter, pop[x][y]->males[dad]);
												//output individual
												ind->outRanInd(0, r, g, &ranInd);
												ind->deleteInd();

												delete ind;
										}
									}


									//OFFSPRING (if female survives)---------------------------------------------------
									std::bernoulli_distribution surv(p_surv);
									if (surv(rdgen)) {
										if (para.stoch_rep) n_offs = pois(rdgen);
										else n_offs = (int)para.fec; //deterministic reproduction

										pop[x][y]->Noffs += n_offs;
										//father sampling distribution
										std::uniform_int_distribution<> sdad(0, iter->n_mates - 1);

										if (para.stoch_rep) {
											for (int i = 0; i < n_offs; i++) {
												ind = new Individuals(para, Bern(rdgen), x, y);
												if (iter->n_mates > 1) dad = iter->mates[sdad(rdgen)];
												else dad = iter->mates[0];


												//inheritance
												inheritance(ind, *iter, pop[x][y]->males[dad]);

												//neutral mutation
												n_mut = n_neutmut(rdgen);
												if (n_mut > 0) ind->neutral_mutation(n_mut, neutr_position, neutral);

												//back mutations
												n_mut = n_backmut(rdgen);
												if (n_mut > 0 && ind->chromo.nMut > n_mut) ind->back_mutation(n_mut);

												//mildly deleterius mutations
												n_mut = n_mildmut(rdgen);
												if (n_mut > 0) ind->delet_mutation(n_mut, k, position, s_mild, uniform, para);

												//lethal mutations
												n_mut = n_lethal(rdgen);
												if (n_mut > 0) ind->delet_mutation(n_mut, position, para.sl, para.hl);

												//mutations to adaptive alleles
												if (g > para.start_dispEvol - 1) adaptiveMutations(ind, MutEm, MutDist);

												n_offs_b++;

												if (ind->w > 0.0) {

													if (para.loadEffect == 0) {
														if (ind->w > 1.0) p_surv = 1.0;
														else p_surv = ind->w;
														std::bernoulli_distribution survOffs(p_surv);
														ind->alive = survOffs(rdgen);
													}
													if (para.loadEffect == -9) {
														p_surv = std::exp(-1.0 * (ind->h / (double)para.nL));
														std::bernoulli_distribution survOffs(p_surv);
														ind->alive = survOffs(rdgen);
													}
													if (ind->alive) {
														if (ind->sex) {
															pop[x][y]->Jfemales.push_back(*ind);
															pop[x][y]->Foffs++;
														}
														else {
															pop[x][y]->Jmales.push_back(*ind);
															pop[x][y]->Moffs++;
														}
													}
													else {
														ind->deleteInd();
														pop[x][y]->Noffs--;
													}
												}
												else {
													ind->deleteInd();
													pop[x][y]->Noffs--;
												}

												delete ind;
											}
										}
									}
								}
							}
							if (!a_mates.empty()) a_mates.clear();
						}
					}
				}
				pop[x][y]->N = 0;
				pop[x][y]->Nf = 0;
				pop[x][y]->Nm = 0;
				pop[x][y]->deleteAdults();
			}
		}
	}
	mat_trials = false;
}

//---------------------------------------------------------------------------
//#if TEST
void inheritance(Individuals* pup, Individuals mom, Individuals dad) {
	int rdn, rdn2, pos, pos2;
	int hom;
	int n_crossovers;
	double cross;
	std::map<double, mutation>::iterator iter, iter2;
	std::set <double> recomSites;
	std::set<double>::iterator itercross;

	//Neutral loci and homozygosity
	rdn = Bern(rdgen);
	rdn2 = Bern(rdgen);

	pos = 0; pos2 = 0;
	for (int i = 0; i < para.nL; i++) {
		if (pos == i) {
			if (rdn) rdn = 0;
			else rdn = 1;
			pos += 1 + geom(rdgen);
		}
		if (pos2 == i) {
			if (rdn2) rdn2 = 0;
			else rdn2 = 1;
			pos2 += 1 + geom(rdgen);
		}
		pup->markers.push_back(mom.markers[i * 2 + rdn]);
		pup->markers.push_back(dad.markers[i * 2 + rdn2]);
		if (mom.markers[i * 2 + rdn] == dad.markers[i * 2 + rdn2]) pup->h++;
	}

	//Recombination of deleterious mutations (see Roze & Rousset 2009, JEB - Appendix 2)
	//Inherit homologue 1 from the mother----------------------------------------------

	if (mom.chromo.nMut > 0) {
		//sample no. of crossovers
		n_crossovers = crossn(rdgen);
		//sample crossover positions
		for (int i = 0; i < n_crossovers; i++) {
			cross = position(rdgen);
			recomSites.insert(cross);
		}
		itercross = recomSites.begin(); //iterator through crossover positions

		//sample starting homologue
		hom = Bern(rdgen);
		iter = mom.chromo.mutations.begin();

		//no mutations before cross-overs positions
		while (n_crossovers > 0 && *itercross < iter->first) {
			if (hom == 0) hom++;
			else hom--;
			itercross++;
			n_crossovers--;
		}
		for (iter = mom.chromo.mutations.begin(); iter != mom.chromo.mutations.end(); iter++) {
			//cross-overs
			while (n_crossovers > 0 && *itercross < iter->first) {
				if (hom == 0) hom++;
				else hom--;
				itercross++;
				n_crossovers--;
			}

			//if mutation is on the right homologue inherit it, otherwise ignore it
			if (iter->second.homol == hom || iter->second.homol == 2) {
				pup->chromo.mutations[iter->first] = iter->second;
				pup->chromo.mutations[iter->first].homol = 0; //inherit first homologue from mon

				pup->chromo.nMut++;
				//calculate fitness considering the mutation as heterozygote
				pup->w *= (1.0 - iter->second.h * iter->second.s);
			}

		}
		if (!recomSites.empty()) recomSites.clear();
	}
	else {
		hom = Bern(rdgen);
	}

	if (dad.chromo.nMut > 0) {
		//sample no. of crossovers
		n_crossovers = crossn(rdgen);
		//sample crossover positions
		for (int i = 0; i < n_crossovers; i++) {
			cross = position(rdgen);
			recomSites.insert(cross);
		}
		itercross = recomSites.begin(); //iterator through crossover positions
										//sample starting homologue
		hom = Bern(rdgen);
		iter = dad.chromo.mutations.begin();

		//no mutations before cross-overs positions
		while (n_crossovers > 0 && *itercross < iter->first) {

			if (hom == 0) hom++;
			else hom--;
			itercross++;
			n_crossovers--;
		}
		for (iter = dad.chromo.mutations.begin(); iter != dad.chromo.mutations.end(); iter++) {
			//crossovers
			while (n_crossovers > 0 && *itercross < iter->first) {
				if (hom == 0) hom++;
				else hom--;
				itercross++;
				n_crossovers--;
			}

			//if mutation is on the right homologue inherit it, otherwise ignore it
			if (iter->second.homol == hom || iter->second.homol == 2) {

				iter2 = pup->chromo.mutations.find(iter->first);
				//if mutation is already present --> it is homozygous
				if (iter2 != pup->chromo.mutations.end()) {
					iter2->second.homol = 2; //mutation is homozygote
					pup->chromo.Nho++;
					//change fitness effect
					pup->w /= (1.0 - iter2->second.h * iter2->second.s);
					pup->w *= (1.0 - iter2->second.s);
				}
				else { //mutation is heterozygote
					pup->chromo.mutations[iter->first] = iter->second;
					pup->chromo.mutations[iter->first].homol = 1;
					pup->chromo.nMut++;
					//fitness effect
					pup->w *= (1.0 - iter->second.h * iter->second.s);
				}
			}
		}
		if (!recomSites.empty()) recomSites.clear();
	}
	//Recombination of the adaptive traits
	if (para.dispEvol) {
		if (para.dispEvol > 0) {
			switch (para.dispEvol) {
			case 1:
				rdn = Bern(rdgen);
				rdn2 = Bern(rdgen);
				for (int i = 0; i < para.L; i++) {
					if (rdn) rdn -= recomb(rdgen);
					else rdn += recomb(rdgen);
					if (rdn2) rdn2 -= recomb(rdgen);
					else rdn2 += recomb(rdgen);
					pup->traits.ep_f[i * 2] = mom.traits.ep_f[i * 2 + rdn];
					pup->traits.ep_f[i * 2 + 1] = dad.traits.ep_f[i * 2 + rdn2];
					*pup->traits.g_ep_f += pup->traits.ep_f[i * 2] + pup->traits.ep_f[i * 2 + 1];
				}
				//phenotype
				*pup->traits.p_ep_f = *pup->traits.g_ep_f;
				if (*pup->traits.p_ep_f < 0.0) *pup->traits.p_ep_f = 0.0;
				if (*pup->traits.p_ep_f > 1.0) *pup->traits.p_ep_f = 1.0;
				break;
			case 2:
				rdn = Bern(rdgen);
				rdn2 = Bern(rdgen);
				for (int i = 0; i < para.L; i++) {
					if (rdn) rdn -= recomb(rdgen);
					else rdn += recomb(rdgen);
					if (rdn2) rdn2 -= recomb(rdgen);
					else rdn2 += recomb(rdgen);
					pup->traits.dist_f[i * 2] = mom.traits.dist_f[i * 2 + rdn];
					pup->traits.dist_f[i * 2 + 1] = dad.traits.dist_f[i * 2 + rdn2];
					*pup->traits.g_dist_f = pup->traits.dist_f[i * 2] + pup->traits.dist_f[i * 2 + 1];
				}
				//phenotype
				*pup->traits.p_dist_f = *pup->traits.g_dist_f;
				if (*pup->traits.p_dist_f < 0.0) *pup->traits.p_dist_f = 0.0;
				break;
			case 3:
				rdn = Bern(rdgen);
				rdn2 = Bern(rdgen);
				for (int i = 0; i < para.L; i++) {
					if (rdn) rdn -= recomb(rdgen);
					else rdn += recomb(rdgen);
					if (rdn2) rdn2 -= recomb(rdgen);
					else rdn2 += recomb(rdgen);
					pup->traits.ep_f[i * 2] = mom.traits.ep_f[i * 2 + rdn];
					pup->traits.ep_f[i * 2 + 1] = dad.traits.ep_f[i * 2 + rdn2];
					*pup->traits.g_ep_f += pup->traits.ep_f[i * 2] + pup->traits.ep_f[i * 2 + 1];
				}
				//phenotype
				*pup->traits.p_ep_f = *pup->traits.g_ep_f;
				if (*pup->traits.p_ep_f < 0.0) *pup->traits.p_ep_f = 0.0;
				if (*pup->traits.p_ep_f > 1.0) *pup->traits.p_ep_f = 1.0;

				rdn = Bern(rdgen);
				rdn2 = Bern(rdgen);
				for (int i = 0; i < para.L; i++) {
					if (rdn) rdn -= recomb(rdgen);
					else rdn += recomb(rdgen);
					if (rdn2) rdn2 -= recomb(rdgen);
					else rdn2 += recomb(rdgen);
					pup->traits.dist_f[i * 2] = mom.traits.dist_f[i * 2 + rdn];
					pup->traits.dist_f[i * 2 + 1] = dad.traits.dist_f[i * 2 + rdn2];
					*pup->traits.g_dist_f = pup->traits.dist_f[i * 2] + pup->traits.dist_f[i * 2 + 1];
				}
				//phenotype
				*pup->traits.p_dist_f = *pup->traits.g_dist_f;
				if (*pup->traits.p_dist_f < 0.0) *pup->traits.p_dist_f = 0.0;

				break;
			}
		}
	}
}

void adaptiveMutations(Individuals* pup, std::normal_distribution<> MutEm, std::normal_distribution<> MutDist) {
	int nmut;
	int allele;

	if (para.dispEvol) {
		if (para.dispEvol > 0) {
			switch (para.dispEvol) {
			case 1:
				nmut = n_amut(rdgen);
				if (nmut > 0) {
					for (int i = 0; i < nmut; i++) {
						allele = uni_loci(rdgen);
						pup->traits_mutation(1, allele, MutEm);
					}
				}
				break;
			case 2:
				nmut = n_amut(rdgen);
				if (nmut > 0) {
					for (int i = 0; i < nmut; i++) {
						allele = uni_loci(rdgen);
						pup->traits_mutation(2, allele, MutDist);
					}
				}
				break;
			case 3:
				nmut = n_amut(rdgen);
				if (nmut > 0) {
					for (int i = 0; i < nmut; i++) {
						allele = uni_loci(rdgen);
						pup->traits_mutation(1, allele, MutEm);
					}
				}
				nmut = n_amut(rdgen);
				if (nmut > 0) {
					for (int i = 0; i < nmut; i++) {
						allele = uni_loci(rdgen);
						pup->traits_mutation(2, allele, MutDist);
					}
				}
				break;
			}
		}
	}
}


//--------------------------------------------------------------------------
void dispersal(void) {
	int max_y_disp;
	int maxy;
	int new_x, new_y;
	double x_rand, y_rand;
	double R1, dist, rndAngle;
	double emp; //emigration probabilty 
	double mdist;
	vector<Individuals>::iterator iter;

	std::uniform_real_distribution<> unireal_disp(0.0, 0.999);
	std::uniform_real_distribution<> unireal_dispB(0.0000001, 1.0);

	max_y_disp = para.y_max;

	maxy = cy_max;
	for (int x = 0; x < para.x_max; x++) {
		for (int y = cy_min; y < maxy; y++) {
			if (land[x][y]->suitable && pop[x][y] != NULL) {
				if (pop[x][y]->Noffs > 0) {
					//females
					for (iter = pop[x][y]->Jfemales.begin(); iter != pop[x][y]->Jfemales.end(); iter++) {
						if (iter->alive) {
							emp = 0.0;
							if (para.dispEvol == 1 || para.dispEvol == 3) emp = *iter->traits.p_ep_f;
							else emp = para.genot_mean[0];
							//Disperses?
							if ( unireal(rdgen) < emp) {
								if (unireal(rdgen) > para.dispCost) { //survives?
									if (para.nearest) { //nearest-neighbour
										std::uniform_int_distribution<> sample_xy(-1, 1);
										do {
											do {
												new_x = x + sample_xy(rdgen);
												new_y = y + sample_xy(rdgen);
											} while (new_x == x && new_y == y);
										} while (new_x < 0 || new_x >(para.x_max - 1) || new_y < cy_min || new_y >(max_y_disp - 1));

									}
									else {
										//dispersal distance
										if (para.dispEvol == 2 || para.dispEvol == 3) mdist = *iter->traits.p_dist_f;
										else mdist = para.genot_mean[1];
										//sample new location
										x_rand = unireal_disp(rdgen);
										y_rand = unireal_disp(rdgen);
										do {
											do {
												R1 = unireal_dispB(rdgen);
												dist = (-1.0 * mdist) * std::log(R1);
												rndAngle = unireal(rdgen) * 2.0 * PI;
												new_x = (int)(dist * cos(rndAngle) / (double)para.resol + x_rand + x);
												new_y = (int)(dist * sin(rndAngle) / (double)para.resol + y_rand + y);
											} while (new_x == x && new_y == y);
										} while (new_x < 0 || new_x >(para.x_max - 1) || new_y < cy_min || new_y >(max_y_disp - 1));
									}
									//settle or die
									pop[x][y]->Foffs--;
									pop[x][y]->Noffs--;
									if (land[new_x][new_y]->suitable) {
										iter->dispersed = true;
										if (pop[new_x][new_y] == NULL) pop[new_x][new_y] = new Population(new_x, new_y);
										pop[new_x][new_y]->tmp_females.push_back(*iter);
										pop[new_x][new_y]->Foffs++;
										pop[new_x][new_y]->Noffs++;

										if (new_y > cy_max) cy_max = new_y + 1; //update maximum y
									}
									else iter->deleteInd();
								}
								else {
									pop[x][y]->Foffs--;
									pop[x][y]->Noffs--;
									iter->deleteInd();
								}
							}
							else { //resident
								pop[x][y]->tmp_females.push_back(*iter);
							}
						}
						else iter->deleteInd();
					}
					//males
					for (iter = pop[x][y]->Jmales.begin(); iter != pop[x][y]->Jmales.end(); iter++) {

						if (iter->alive) {
							emp = 0.0;
							if (para.dispEvol == 1 || para.dispEvol == 3) {
								emp = *iter->traits.p_ep_f;
							}
							else emp = para.genot_mean[0];
							//Disperses?
							if (unireal(rdgen) < emp) {
								if (unireal(rdgen) > para.dispCost) { //survives
								//iter->w -= para.dispCost;
									if (para.nearest) { //nearest-neighbour
										std::uniform_int_distribution<> sample_xy(-1, 1);
										do {
											do {
												new_x = x + sample_xy(rdgen);
												new_y = y + sample_xy(rdgen);
											} while (new_x == x && new_y == y);
										} while (new_x < 0 || new_x >(para.x_max - 1) || new_y < cy_min || new_y >(max_y_disp - 1));

									}
									else {
										//dispersal distance
										if (para.dispEvol == 2 || para.dispEvol == 3) {
											mdist = *iter->traits.p_dist_f;
										}
										else mdist = para.genot_mean[1];
										//sample new location
										x_rand = unireal_disp(rdgen);
										y_rand = unireal_disp(rdgen);
										do {
											do {
												R1 = unireal_dispB(rdgen);
												dist = (-1.0 * mdist) * std::log(R1);
												rndAngle = unireal(rdgen) * 2.0 * PI;
												new_x = (int)(dist * cos(rndAngle) / (double)para.resol + x_rand + x);
												new_y = (int)(dist * sin(rndAngle) / (double)para.resol + y_rand + y);
											} while (new_x == x && new_y == y);
										} while (new_x < 0 || new_x >(para.x_max - 1) || new_y < cy_min || new_y >(max_y_disp - 1));
									}
									//settle or die
									pop[x][y]->Moffs--;
									pop[x][y]->Noffs--;
									if (land[new_x][new_y]->suitable) {
										iter->dispersed = true;
										if (pop[new_x][new_y] == NULL) pop[new_x][new_y] = new Population(new_x, new_y);
										pop[new_x][new_y]->tmp_males.push_back(*iter);
										pop[new_x][new_y]->Moffs++;
										pop[new_x][new_y]->Noffs++;

										if (new_y > cy_max) cy_max = new_y + 1; //update maximum y
									}
									else iter->deleteInd();
								}
								else {
									pop[x][y]->Moffs--;
									pop[x][y]->Noffs--;
									iter->deleteInd();
								}
							}
							else { //resident
								pop[x][y]->tmp_males.push_back(*iter);
							}
						}
						else iter->deleteInd();
					}
				}
				pop[x][y]->Jfemales.clear();
				pop[x][y]->Jmales.clear();
			}
		}
	}
}
//--------------------------------------------------------------------------
//Density-dependent survival
void survival(void) {
	int maxy;
	double ps;

	vector<Individuals>::iterator iter;
	std::map<double, mutation>::iterator iter2;

	Ntot = 0;

	maxy = cy_max;

	for (int x = 0; x < para.x_max; x++) {
		for (int y = cy_min; y < maxy; y++) {
			if (land[x][y]->suitable && pop[x][y] != NULL) {
				if (pop[x][y]->Noffs > 0) {

					pop[x][y]->set2zero(); //set population stats to zero
					ps = std::fmin(land[x][y]->local_K / (double)pop[x][y]->Noffs, 1.0);

					//ps = std::fmin(land[x][y]->local_K / (double)pop[x][y]->Noffs, 1.0);
					std::bernoulli_distribution survive(ps);
					//cout << pop[x][y]->alt_K << endl;
					for (iter = pop[x][y]->tmp_females.begin(); iter != pop[x][y]->tmp_females.end(); iter++) {
						
						if (iter->alive && survive(rdgen)) {

							pop[x][y]->females.push_back(*iter);
							pop[x][y]->Nf++;
							pop[x][y]->N++;
							pop[x][y]->computeSums(para, *iter);

							//sum deleterious mutation to calculate population frequencies
							if (g % para.PopMut_interval == 0) {
								for (iter2 = iter->chromo.mutations.begin(); iter2 != iter->chromo.mutations.end(); iter2++) {
									pop[x][y]->addMutation(iter2->first, iter2->second.s, iter2->second.h);
									if (iter2->second.homol == 2) pop[x][y]->addMutation(iter2->first, iter2->second.s, iter2->second.h);
								}
							}
							//output mutations
							if (para.out_mutations)(iter->outMut(g, &muts));

							//output individuals
							if (para.ind_interval > 0 && g % para.ind_interval == 0)
								iter->outInd(para, r, g, &inds);
						}
						else iter->deleteInd();
					}
					for (iter = pop[x][y]->tmp_males.begin(); iter != pop[x][y]->tmp_males.end(); iter++) {
						
						if (iter->alive && survive(rdgen)) {

							pop[x][y]->males.push_back(*iter);
							pop[x][y]->Nm++;
							pop[x][y]->N++;
							pop[x][y]->computeSums(para, *iter);

							//sum deleterious mutation to calculate population frequencies
							if (g % para.PopMut_interval == 0) {
								for (iter2 = iter->chromo.mutations.begin(); iter2 != iter->chromo.mutations.end(); iter2++) {
									pop[x][y]->addMutation(iter2->first, iter2->second.s, iter2->second.h);
									if (iter2->second.homol == 2) pop[x][y]->addMutation(iter2->first, iter2->second.s, iter2->second.h);
								}
							}
							//output mutations
							if (para.out_mutations)(iter->outMut(g, &muts));

							//output individuals
							if (para.ind_interval > 0 && g % para.ind_interval == 0)
								iter->outInd(para, r, g, &inds);
						}
						else iter->deleteInd();
					}

					if (pop[x][y]->N > 0) {

						pop[x][y]->computeStats(para);
						//output populations and traits
						if (g > para.out_start - 1 && (g % para.out_int == 0)) {
							pop[x][y]->outPop(r, g, &pops);
							if ((para.dispEvol) && g > para.start_dispEvol - 1)
								pop[x][y]->outTrait(para, r, g, &trait);
						}
						//output frequency of deleterious mutations
						if (g > para.out_start - 1 && g > 0 && g % para.PopMut_interval == 0) pop[x][y]->outMutations(pop[x][y]->N, r, g, &popmut);

						pop[x][y]->Noffs = 0;
						pop[x][y]->Foffs = 0;
						pop[x][y]->Moffs = 0;
						pop[x][y]->tmp_females.clear();
						pop[x][y]->tmp_males.clear();

						Ntot += pop[x][y]->N;
					}
					else {
						delete pop[x][y];
						pop[x][y] = NULL;
					}

				}
				else {
					delete pop[x][y];
					pop[x][y] = NULL;
				}
			}
		}
	}
}

void betweenPop_mating(void) {
	int maxy = cy_max;
	int xx, yy;
	int dad;
	int count_N; //counts the number of non-empty cells.

	Individuals* ind;

	vector<Individuals>::iterator iter;

	std::uniform_int_distribution<> sample_x(0, para.x_max - 1);
	std::uniform_int_distribution<> sample_y(cy_min, cy_max - 1);
	
	count_N = 0;
	for (int x = 0; x < para.x_max; x++) {
		for (int y = cy_min; y < maxy; y++) {
			if (pop[x][y] != NULL)(count_N += 1);
		}
	}
	
	if (count_N >= 4) {
		//Loop through the landscape
		for (int x = 0; x < para.x_max; x++) {
			for (int y = cy_min; y < maxy; y++) {
				if (pop[x][y] != NULL) {
					if (pop[x][y]->Nf > 0) {
						//Loop through females and generate 10 offspring each with males from other random populations
						for (iter = pop[x][y]->females.begin(); iter != pop[x][y]->females.end(); iter++) {

							for (int i = 0; i < 10; i++) {
								do {
									xx = sample_x(rdgen);
									yy = sample_y(rdgen);
								} while (pop[xx][yy] == NULL || pop[xx][yy]->Nm < 1 || (xx == x && yy == y));
								std::uniform_int_distribution<> sample_mate(0, pop[xx][yy]->Nm - 1);
								dad = sample_mate(rdgen);

								ind = new Individuals(para, Bern(rdgen), x, y);
								//inheritance
								inheritance(ind, *iter, pop[xx][yy]->males[dad]);

								//output individual
								ind->outRanInd(1, r, g + 1, &ranInd);

								ind->deleteInd();

								delete ind;
							}
						}
					}
				}
			}
		}
	}
}


void outPop_header(void) {
	string name;

	name = dirOut + "Sim" + Int2Str(para.SimNr) + "_Pops.txt";
	pops.open(name.c_str());
	pops << "rep\tgen\tx\ty\tNf\tNm";
	pops << "\tW\tWstd\tHomoz\tHstd\tDelMut\tDelMutSd";
	pops << endl;
}
//----------------------------------------------------------------------------------------
void outTrait_header(void) {
	string name;

	name = dirOut + "Sim" + Int2Str(para.SimNr) + "_Traits.txt";
	trait.open(name.c_str());

	trait << "rep\tgen\tx\ty";
	if (para.dispEvol > 0) {
		switch (para.dispEvol) {
		case 1:
			trait << "\tfep\tfepStd";
			break;
		case 2:
			trait << "\tfdist\tfdistStd";
			break;
		case 3:
			trait << "\tfep\tfepStd";
			trait << "\tfdist\tfdistStd";
			break;
		}
	}
	trait << endl;
}
//----------------------------------------------------------------------------------------
void outMut_header(void) {
	string name;

	name = dirOut + "Sim" + Int2Str(para.SimNr) + "_Muts_rep" + Int2Str(r) + ".txt";
	muts.open(name.c_str());
	muts << "gen\ts\th\thomo" << endl;
}
//----------------------------------------------------------------------------------------
void outInds_header(void) {
	string name;

	name = dirOut + "Sim" + Int2Str(para.SimNr) + "_Inds.txt";
	inds.open(name.c_str());
	inds << "rep\tgen\tx\ty\tdispersed\tW\thomoz";
	inds << "\tNmut\tNho";
	if (para.dispEvol > 0) {
		switch (para.dispEvol) {
		case 1:
			inds << "\tpEp";
			break;
		case 2:
			inds << "\tpDist";
			break;
		case 3:
			inds << "\tgEp\tpEp\tgDist\tpDist";
			break;
		}
	}
	inds << endl;
}
//----------------------------------------------------------------------------------------
void outRanInds_header(void)
{
	string name;

	name = dirOut + "Sim" + Int2Str(para.SimNr) + "_RanInds.txt";
	ranInd.open(name.c_str());
	//"out" indicates whether it is a "between population" off (1) or a "within population" off (0)
	ranInd << "rep\tgen\tx\ty\tout\tW\thomoz\tNmut\tNho";
	ranInd << endl;
}
//----------------------------------------------------------------------------------------
void outPopMut_header(void)
{
	string name;
	name = dirOut + "Sim" + Int2Str(para.SimNr) + "_PopMut.txt";
	popmut.open(name.c_str());

	popmut << "rep\tgen\tx\ty\ts\th\tfreq" << endl;
}



