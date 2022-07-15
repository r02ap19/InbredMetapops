#pragma once

#include <map>
#include <set>
#include <stdio.h>
#include <fstream>
#include <iostream>

#include "Parameters.h"


struct mutation {
	int homol; //0 = only on homologue 1; 1 = only on homologue 2; 2 = on both homologues (hence homozygote)
	double s; //selection coefficient
	double h; //dominance coefficient
};

//map containing all mutations: map<position, mutation> 
typedef std::map<double, mutation, std::less<double>> MapMuts;

typedef std::map<double, int, std::less<double>> MapPos;

//map containing continuous allele coding for a trait: map<position, allele>
typedef std::map<double, double, std::less<double>> MapTrait;


class Chromosome {
public:
	Chromosome();
	~Chromosome();
	int nMut; //no. of deleterious mutations
	int Nho; //number of homozygote mutations

	MapMuts mutations; //list of all deleterious mutations, their coefficients and whether they are homo or heterozygote 

	double addDelMutation(
		int, //homologue
		double, //position
		double, //s
		double //h
	);

	void deleteChromo();
private:
};