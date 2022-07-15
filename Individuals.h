#pragma once

#include <stdio.h>
#include <fstream>
#include <math.h>
#include <random>
#include <vector>
#include <set>
#include <iostream>

#include "Chromosome.h"
#include "Traits.h"

using namespace std;


class Individuals {
public:
	Individuals(Parameters, bool, int, int);
	~Individuals();
	bool dispersed; //1 = it has dispersed
	bool alive;
	bool sex; //1 = female
	bool reproduce;
	int x, y;    //cell coordinates
	int xnatal, ynatal; //coordinates of natal cell
	int noffs; //number of offspring
	int n_mates; //number of mates
	int h; //number of homozygote neutral loci

	double w; //individual viability

	Traits traits;

	Chromosome chromo;

	vector<int> mates; //mates' indeces
	vector<double> markers; //continuous neutral markers to determine inbreeding

	void initialise(double kk, Parameters para, std::normal_distribution<>, std::normal_distribution<>, 
		std::gamma_distribution<> finds, std::uniform_real_distribution<>, std::uniform_real_distribution<>);
	void traits_mutation(int, int, std::normal_distribution<>);
	void neutral_mutation(int, std::uniform_int_distribution<>, std::uniform_real_distribution<>);
	void benef_mutation(int, double, std::uniform_real_distribution<>);
	void back_mutation(int);
	void delet_mutation(int, double, std::uniform_real_distribution<>, std::gamma_distribution<>, std::uniform_real_distribution<>, Parameters para);
	void delet_mutation(int, std::uniform_real_distribution<>, double, double);

	void deleteInd(void);
	void outInd(Parameters para, int r, int g, std::ofstream* out); //individual output
	void outMut(int g, std::ofstream* out);
	void outRanInd(int outb, int r, int g, std::ofstream* out); //individual output for calculating ID
	void outKinship(int g, int gg, std::ofstream* out);
private:
};

