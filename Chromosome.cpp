#include  "Chromosome.h"


Chromosome::Chromosome() {
	nMut = 0;
	Nho = 0;
}
//-------------------------------
Chromosome::~Chromosome() {
}
//-------------------------------
void Chromosome::deleteChromo(void) {
	if (!mutations.empty()) mutations.clear();
}

double Chromosome::addDelMutation(int hom, double pos, double ss, double hh) {
	double v;
	mutation mut;
	mut.homol = hom;
	mut.s = ss;
	mut.h = hh;

	if (mutations.find(pos) == mutations.end()) {
		nMut++;
		mutations[pos] = mut;
		v = (1.0 - hh * ss);
	}
	else v = 1.0;

	return v;
}
