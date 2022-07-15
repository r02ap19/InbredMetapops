#pragma once

#include <vector>
#include <stdio.h>

#include "Parameters.h"

class Traits
{
public:
	Traits();
	~Traits();


	double* ep_f, * ep_m; //female and male emigration probability 
	double* dist_f, * dist_m; //female and male mean dispersal distance

	//genotypic values
	double* g_ep_f, * g_ep_m;
	double* g_dist_f, * g_dist_m;

	//phenotypic values
	double* p_ep_f, * p_ep_m;
	double* p_dist_f, * p_dist_m;

	void initialise(Parameters);
	void deleteTraits(void);

private:
};


