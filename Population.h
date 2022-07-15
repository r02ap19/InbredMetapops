#pragma once

#include <stdio.h>
#include <vector>
#include <iostream>

#include "Individuals.h"
#include "Parameters.h"

using namespace std;

struct pop_muts {

	int count; //mutation occurrence to calculate frequency when output
	double s; //selection coefficient
	double h; //dominance coefficient
};

//map containing all mutations in the population: map<position, pop_muts> 
typedef std::map<double, pop_muts, std::less<double>> MapPopMuts;

class Population {
public:
	Population(int, int); //pass coordinates
	~Population();
	int x, y;
	int N, Noffs;
	int Nf, Nm; //total females and males
	int Foffs, Moffs; //nr. of female and male offs
	int DelMut; //nr. of deleterious mutations
	double meanDelMut;
	double DelMutSd; //standard deviation
	double W, Wsd; //mean population (calculated after survival) and standard deviation
	double Homoz, Hsd; //mean neutral homozygosity and standard deviation

	double mFep, mMep, mFdist, mMdist; //mean trait phenotypic values
	double sFep, sMep, sFdist, sMdist; //standard deviations

	vector<Individuals> females, males, Jfemales, Jmales;
	vector<Individuals> tmp_females, tmp_males;

	MapPopMuts popMuts;

	void initialise_pop(double, double, Parameters, std::normal_distribution<>, std::normal_distribution<>, 
		std::gamma_distribution<>, std::uniform_real_distribution<>, std::uniform_real_distribution<>);
	void deleteAdults(void);
	void computeSums(Parameters, Individuals);
	void computeStats(Parameters);
	void set2zero(void);
	void outPop(int, int, std::ofstream*);
	void outTrait(Parameters, int, int, std::ofstream*);
	void addMutation(
		double, //position
		double, //s
		double //h
	);
	void outMutations(int, int, int, std::ofstream*);
};
