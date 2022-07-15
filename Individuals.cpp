#include "Individuals.h"

Parameters para2;

std::random_device rd2;
std::mt19937 gen(rd2());

std::bernoulli_distribution findhom(0.5); //for sampling the homologue
std::uniform_int_distribution<> loci(0, (para2.L - 1));

Individuals::Individuals(Parameters para, bool sx, int xx, int yy) {
	dispersed = false;
	alive = true;
	reproduce = true;
	sex = sx;
	x = xx; y = yy;
	xnatal = xx; ynatal = yy;
	noffs = 0;
	n_mates = 0;
	w = 1.0;
	h = 0;
	traits.initialise(para);

}
//--------------------------------------------------------------------------
Individuals::~Individuals() {

}
//--------------------------------------------------------------------------
void Individuals::deleteInd(void) {
	traits.deleteTraits();
	chromo.deleteChromo();
	if (!mates.empty()) mates.clear();
	if (!markers.empty()) markers.clear();
}
//--------------------------------------------------------------------------

void Individuals::initialise(double kk, Parameters para, std::normal_distribution<> NormEm, std::normal_distribution<> NormDist, std::gamma_distribution<> finds,
	std::uniform_real_distribution<> neutr, std::uniform_real_distribution<> findpos) {

	int hom;
	double pos;
	double ss;
	double hh;

	//initialise neutral markers
	for (int i = 0; i < (para.nL * 2); i++) {
		markers.push_back(neutr(gen));
	}
	if (para.initial_nMut > 0) {

		for (int i = 0; i < para.initial_nMut; i++) {
			//sample homologue
			hom = findhom(gen);
			//sample selection coefficient & dominance coefficient
			if (para.sh_dist == 0) {
				ss = finds(gen);
				std::uniform_real_distribution<> findh(0.0, std::exp(-kk * ss));
				hh = findh(gen);
			}
		//sample position
			pos = findpos(gen);
			//add mutation
			w *= chromo.addDelMutation(hom, pos, ss, hh);

		}
	}

	if (para.dispEvol > 0) {
		//if dispersal is not sex limited use only the female trait
		switch (para.dispEvol) {
		case 1:
			for (int i = 0; i < (2 * para.L); i++) {
				traits.ep_f[i] = NormEm(gen);
				*traits.g_ep_f += traits.ep_f[i];
			}
			*traits.p_ep_f = *traits.g_ep_f;
			if (*traits.p_ep_f < 0.0) *traits.p_ep_f = 0.0;
			if (*traits.p_ep_f > 1.0) *traits.p_ep_f = 1.0;
			break;
		case 2:
			for (int i = 0; i < (2 * para.L); i++) {
				traits.dist_f[i] = NormDist(gen);
				*traits.g_dist_f += traits.dist_f[i];
			}
			*traits.p_dist_f = *traits.g_dist_f;
			if (*traits.p_dist_f < 0.0) *traits.p_dist_f = 0.0;
			break;
		case 3:
			for (int i = 0; i < (2 * para.L); i++) {
				traits.ep_f[i] = NormEm(gen);
				*traits.g_ep_f += traits.ep_f[i];
				traits.dist_f[i] = NormDist(gen);
				*traits.g_dist_f += traits.dist_f[i];
			}
			*traits.p_ep_f = *traits.g_ep_f;
			if (*traits.p_ep_f < 0.0) *traits.p_ep_f = 0.0;
			if (*traits.p_ep_f > 1.0) *traits.p_ep_f = 1.0;
			*traits.p_dist_f = *traits.g_dist_f;
			if (*traits.p_dist_f < 0.0) *traits.p_dist_f = 0.0;
			break;
		}
	}
}
//--------------------------------------------------------------------------
void Individuals::neutral_mutation(int nmut, std::uniform_int_distribution<> findpos, std::uniform_real_distribution<> mutate) {
	int allele;

	//std::random_device rd;
	//std::mt19937 gen(rd());

	for (int i = 0; i < nmut; i++) {
		allele = findpos(gen);
		//if locus was homozygote, discount for that
		if (allele % 2 == 0.0) { //1st allele
			if (markers[allele] == markers[allele + 1]) h--;
		}
		else { //2nd allele
			if (markers[allele - 1] == markers[allele]) h--;
		}
		//assume that you won't get an homozygote from continuous mutation 
		markers[allele] = mutate(gen);
	}
}
//--------------------------------------------------------------------------
void Individuals::benef_mutation(int nmut, double s, std::uniform_real_distribution<> findpos) {
	int hom;
	double pos;
	double ss;
	std::bernoulli_distribution findhom(0.5); //for sampling the homologue

	for (int i = 0; i < nmut; i++) {
		//sample homologue
		hom = findhom(gen);
		//sample selection coefficient
		ss = s;
		//sample position
		pos = findpos(gen);
		//add mutation
		w *= chromo.addDelMutation(hom, pos, ss, 1.0); //assume complete dominance of beneficial mutations
	}
}
//--------------------------------------------------------------------------
void Individuals::back_mutation(int nmut)
{
	int rdn;
	std::map<double, mutation, std::less<double>>::iterator it;

	std::bernoulli_distribution hom(0.5);

	for (int i = 0; i < nmut; i++) {
		std::uniform_int_distribution<> samp(0, chromo.mutations.size() - 1);
		rdn = samp(gen);
		it = chromo.mutations.begin();
		std::advance(it, rdn);

		if (it->second.homol == 2) { //mutation is homozygote
			it->second.homol = hom(gen);
			chromo.Nho--;
			//update fitness
			w /= (1.0 - it->second.s);
			w *= (1.0 - it->second.h * it->second.s);
		}
		else { //mutation is heterozygote
			chromo.mutations.erase(it);
			chromo.nMut--;
			//update fitness
			w /= (1.0 - it->second.h * it->second.s);
		}
	}
}
//--------------------------------------------------------------------------
void Individuals::delet_mutation(int nmut, double kk, std::uniform_real_distribution<> findpos,
	std::gamma_distribution<> finds, std::uniform_real_distribution<> uniform, Parameters para) {
	int hom;
	double pos;
	double ss;
	double hh;

	for (int i = 0; i < nmut; i++) {
		//sample homologue
		hom = findhom(gen);
		if (para.sh_dist == 0) {
			ss = finds(gen);
			std::uniform_real_distribution<> findh(0.0, std::exp(-kk * ss));
			hh = findh(gen);
		}
		//sample position
		pos = findpos(gen);
		//add mutation
		w *= chromo.addDelMutation(hom, pos, ss, hh);
	}
}
//--------------------------------------------------------------------------
void Individuals::delet_mutation(int nmut, std::uniform_real_distribution<> findpos,
	double ss, double hh) {
	int hom;
	double pos;
	for (int i = 0; i < nmut; i++) {
		hom = findhom(gen);
		pos = findpos(gen);
		w *= chromo.addDelMutation(hom, pos, ss, hh);
	}
}
//--------------------------------------------------------------------------
void Individuals::traits_mutation(int trait, int allele, std::normal_distribution<> Norm) {

	switch (trait) {
	case 1: //female emigration probability
		*traits.g_ep_f -= traits.ep_f[allele];
		traits.ep_f[allele] += Norm(gen);
		*traits.g_ep_f += traits.ep_f[allele];
		*traits.p_ep_f = *traits.g_ep_f;
		if (*traits.p_ep_f < 0.0) *traits.p_ep_f = 0.0;
		if (*traits.p_ep_f > 1.0) *traits.p_ep_f = 1.0;
		break;
	case 2: //female dispersal distance
		*traits.g_dist_f -= traits.dist_f[allele];
		traits.dist_f[allele] += Norm(gen);
		*traits.g_dist_f += traits.dist_f[allele];
		*traits.p_dist_f = *traits.g_dist_f;
		if (*traits.p_dist_f < 0.0) *traits.p_dist_f = 0.0;
		break;
	case 3: //male emigration probability
		*traits.g_ep_m -= traits.ep_m[allele];
		traits.ep_m[allele] += Norm(gen);
		*traits.g_ep_m += traits.ep_m[allele];
		*traits.p_ep_m = *traits.g_ep_m;
		if (*traits.p_ep_m < 0.0) *traits.p_ep_m = 0.0;
		if (*traits.p_ep_m > 1.0) *traits.p_ep_m = 1.0;
		break;
	case 4: //male dispersal distance
		*traits.g_dist_m -= traits.dist_m[allele];
		traits.dist_m[allele] += Norm(gen);
		*traits.g_dist_m += traits.dist_m[allele];
		*traits.p_dist_m = *traits.g_dist_m;
		if (*traits.p_dist_m < 0.0) *traits.p_dist_m = 0.0;
		break;
	}
}



//--------------------------------------------------------------------------
void Individuals::outMut(int g, std::ofstream* out) {

	map<double, mutation>::iterator iter;

	for (iter = chromo.mutations.begin(); iter != chromo.mutations.end(); iter++) {
		*out << g << "\t" << iter->second.s << "\t" << iter->second.h << "\t";
		if (iter->second.homol < 2) *out << 0;
		else *out << 1;
		*out << endl;
	}
}
//--------------------------------------------------------------------------
void Individuals::outRanInd(int outb, int r, int g, std::ofstream* out)
{
	* out << r << "\t" << g << "\t" << x << "\t" << y << "\t" << outb << "\t" << w << "\t" << h << "\t" << chromo.nMut << "\t" << chromo.Nho << endl;
}

//--------------------------------------------------------------------------
void Individuals::outInd(Parameters para, int r, int g, std::ofstream* out) {
	*out << r << "\t" << g << "\t" << x << "\t" << y << "\t" << dispersed;
	*out << "\t" << w << "\t" << h;
	* out << "\t" << chromo.nMut << "\t" << chromo.Nho;
	if (para.dispEvol > 0) {
		switch (para.dispEvol) {
		case 1:
			*out << "\t" << *traits.p_ep_f;
			break;
		case 2:
			*out << "\t" << *traits.g_dist_f << "\t" << *traits.p_dist_f;
			break;
		case 3:
			*out << "\t" << *traits.g_ep_f << "\t" << *traits.p_ep_f << "\t" << *traits.g_dist_f << "\t" << *traits.p_dist_f;
			break;
		}
	}
	*out << endl;
}

void Individuals::outKinship(int g, int gg, std::ofstream* out)
{
	map<double, mutation>::iterator iter;

	for (iter = chromo.mutations.begin(); iter != chromo.mutations.end(); iter++) {
		*out << g << "\t" << iter->second.s << "\t" << iter->second.h << "\t" << iter->first << "\t" << gg << "\t" << iter->second.homol;
		*out << endl;
	}
}

