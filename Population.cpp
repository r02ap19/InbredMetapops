#include "Population.h"
#include "Individuals.h"


std::random_device rd3;
std::mt19937 gene(rd3());

Population::Population(int xx, int yy) {
	x = xx;
	y = yy;
	N = 0;
	Noffs = 0;
	Nf = 0;
	Nm = 0;
	Foffs = 0;
	Moffs = 0;
	DelMut = 0;
	meanDelMut = 0.0;
	DelMutSd = 0.0;
	W = 0.0;
	Wsd = 0.0;
	Homoz = 0.0; Hsd = 0.0;

	mFep = 0.0; mMep = 0.0; mFdist = 0.0; mMdist = 0.0;
	sFep = 0.0; sMep = 0.0; sFdist = 0.0; sMdist = 0.0;

	if (!females.empty()) females.clear();
	if (!males.empty()) males.clear();
	if (!Jfemales.empty()) Jfemales.clear();
	if (!Jmales.empty()) Jmales.clear();
	if (!tmp_females.empty()) tmp_females.clear();
	if (!tmp_males.empty()) tmp_males.clear();
	if (!popMuts.empty()) popMuts.clear();
}
//----------------------------------------
Population::~Population() {
	if (!females.empty()) females.clear();
	if (!males.empty()) males.clear();
	if (!Jfemales.empty()) Jfemales.clear();
	if (!Jmales.empty()) Jmales.clear();
	if (!tmp_females.empty()) tmp_females.clear();
	if (!tmp_males.empty()) tmp_males.clear();
}
//----------------------------------------
void Population::initialise_pop(double k, double kk, Parameters para, std::normal_distribution<> NormEm, std::normal_distribution<> NormDist,
	std::gamma_distribution<> finds, std::uniform_real_distribution<> neutr, std::uniform_real_distribution<> findpos) {

	N = (int)k;

	std::uniform_real_distribution<> unireal(0.0, 1.0);

	//initialise females
	for (int i = 0; i < (int)(N / 2); i++) {
		//for (int i = 0; i < 1; i++) {
		females.push_back(Individuals(para, true, x, y));
		females[i].initialise(kk, para, NormEm,NormDist, finds, neutr, findpos);
		Nf++;
		W += 1.0; //at initialisation every individual has viability = -1

	}
	//initialise males
	for (int i = 0; i < (int)N / 2; i++) {
		//for (int i = 0; i < 1; i++) {
		males.push_back(Individuals(para, false, x, y));
		males[i].initialise(kk, para, NormEm,NormDist, finds, neutr, findpos);
		Nm++;
		W += 1.0; //at initialisation every individual has viability = -1

	}

}
//----------------------------------------
void Population::computeSums(Parameters para, Individuals ind) {
	DelMut += ind.chromo.nMut;
	DelMutSd += (double)ind.chromo.nMut * (double)ind.chromo.nMut; //sum of squares
	W += ind.w;
	Wsd += ind.w * ind.w; //sum of squares
	//neutral homozygosity
	Homoz += (double)ind.h / (double)para.nL;
	Hsd += ((double)ind.h / (double)para.nL) * ((double)ind.h / (double)para.nL);//sum of squares
	
	if (para.dispEvol > 0) {
		switch (para.dispEvol) {
		case 1:
			mFep += *ind.traits.p_ep_f; sFep += (*ind.traits.p_ep_f) * (*ind.traits.p_ep_f);
			break;
		case 2:
			mFdist += *ind.traits.p_dist_f; sFdist += (*ind.traits.p_dist_f) * (*ind.traits.p_dist_f);
			break;
		case 3:
			mFep += *ind.traits.p_ep_f; sFep += (*ind.traits.p_ep_f) * (*ind.traits.p_ep_f);
			mFdist += *ind.traits.p_dist_f; sFdist += (*ind.traits.p_dist_f) * (*ind.traits.p_dist_f);
			break;
		}
	}
}
//----------------------------------------
void Population::computeStats(Parameters para) {

	DelMutSd = std::sqrt((DelMutSd - ((double)DelMut * (double)DelMut) / (double)N) / (double)N);
	meanDelMut = (double)DelMut / (double)N;

	Wsd = std::sqrt((Wsd - (W * W) / (double)N) / (double)N);
	W /= (double)N;

	Hsd = std::sqrt((Hsd - (Homoz * Homoz) / (double)N) / (double)N);
	Homoz /= (double)N;

	if (para.dispEvol > 0) {
		switch (para.dispEvol) {
		case 1:
			sFep = std::sqrt((sFep - (mFep * mFep) / (double)N) / (double)N);
			mFep /= (double)N;
			break;
		case 2:
			sFdist = std::sqrt((sFdist - (mFdist * mFdist) / (double)N) / (double)N);
			mFdist /= (double)N;
			break;
		case 3:
			sFep = std::sqrt((sFep - (mFep * mFep) / (double)N) / (double)N);
			mFep /= (double)N;
			sFdist = std::sqrt((sFdist - (mFdist * mFdist) / (double)N) / (double)N);
			mFdist /= (double)N;
			break;
		}
	}
}
//----------------------------------------
void Population::set2zero(void) {
	DelMut = 0;
	DelMutSd = 0.0;
	W = 0.0;
	Wsd = 0.0;
	Homoz = 0.0; Hsd = 0.0;
	mFep = 0.0; mMep = 0.0; mFdist = 0.0; mMdist = 0.0;
	sFep = 0.0; sMep = 0.0; sFdist = 0.0; sMdist = 0.0;
}
//----------------------------------------
void Population::deleteAdults(void) {
	vector<Individuals>::iterator iter;

	if (!females.empty()) {
		for (iter = females.begin(); iter != females.end(); iter++) iter->deleteInd();
		females.clear();
	}
	if (!males.empty()) {
		for (iter = males.begin(); iter != males.end(); iter++) iter->deleteInd();
		males.clear();
	}
}
//----------------------------------------
void Population::outPop(int r, int g, std::ofstream* out) {
	*out << r << "\t" << g << "\t" << x << "\t" << y << "\t" << Nf << "\t" << Nm;
	* out << "\t" << W << "\t" << Wsd << "\t";
	*out << Homoz << "\t" << Hsd << "\t" << meanDelMut << "\t" << DelMutSd;
	* out << endl;
}
//----------------------------------------
void Population::outTrait(Parameters para, int r, int g, std::ofstream* out) {
	*out << r << "\t" << g << "\t" << x << "\t" << y;
	//*out << "\t" << mMates << "\t" << sMates;
	if (para.dispEvol > 0) {
		switch (para.dispEvol) {
		case 1:
			*out << "\t" << mFep << "\t" << sFep;
			break;
		case 2:
			*out << "\t" << mFdist << "\t" << sFdist;
			break;
		case 3:
			*out << "\t" << mFep << "\t" << sFep;
			*out << "\t" << mFdist << "\t" << sFdist;
			break;
		}
	}
	*out << endl;
}

void Population::addMutation(double pos, double ss, double hh)
{
	pop_muts muts;

	muts.h = hh;
	muts.s = ss;

	if (popMuts.find(pos) == popMuts.end()) {
		muts.count = 1;
		popMuts[pos] = muts;
	}
	else popMuts[pos].count++;

}

void Population::outMutations(int n, int r, int g, std::ofstream* out)
{
	map<double, pop_muts>::iterator iter;

	for (iter = popMuts.begin(); iter != popMuts.end(); iter++) {
		*out << r << "\t" << g << "\t" << x << "\t" << y << "\t";
		*out << iter->second.s << "\t" << iter->second.h << "\t";
		*out << (double)iter->second.count / (2.0 * (double)n) << endl;
	}
}









