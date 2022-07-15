#include "Traits.h"


Traits::Traits() {
	ep_f = NULL;
	ep_m = NULL;
	dist_f = NULL;
	dist_m = NULL;
	g_ep_f = NULL;
	g_ep_m = NULL;
	g_dist_f = NULL;
	g_dist_m = NULL;
	p_ep_f = NULL;
	p_ep_m = NULL;
	p_dist_f = NULL;
	p_dist_m = NULL;
}
//--------------------------------------------------------------------------------------------
void Traits::initialise(Parameters para) {

	if (para.dispEvol > 0) {
		switch (para.dispEvol) {
		case 1:
			ep_f = new double[2 * para.L];
			break;
		case 2:
			dist_f = new double[2 * para.L];
			break;
		case 3:
			ep_f = new double[2 * para.L];
			dist_f = new double[2 * para.L];
			break;
		}
	}
	//genotypic and phenotypic values
	if (para.dispEvol > 0) {
		switch (para.dispEvol) {
		case 1:
			g_ep_f = new double(0.0);
			p_ep_f = new double(0.0);
			break;
		case 2:
			g_dist_f = new double(0.0);
			p_dist_f = new double(0.0);
			break;
		case 3:
			g_ep_f = new double(0.0);
			p_ep_f = new double(0.0);
			g_dist_f = new double(0.0);
			p_dist_f = new double(0.0);
			break;
		}
	}
}
//--------------------------------------------------------------------------------------------
Traits::~Traits() {

}
//--------------------------------------------------------------------------------------------
void Traits::deleteTraits(void) {
	if (ep_f != NULL) { delete[] ep_f; ep_f = NULL; }
	if (ep_m != NULL) { delete[] ep_m; ep_m = NULL; }
	if (dist_f != NULL) { delete[] dist_f; dist_f = NULL; }
	if (dist_m != NULL) { delete[] dist_m; dist_m = NULL; }

	if (g_ep_f != NULL) { delete g_ep_f; g_ep_f = NULL; }
	if (g_ep_m != NULL) { delete g_ep_m; g_ep_m = NULL; }
	if (g_dist_f != NULL) { delete g_dist_f; g_dist_f = NULL; }
	if (g_dist_m != NULL) { delete g_dist_m; g_dist_m = NULL; }

	if (p_ep_f != NULL) { delete p_ep_f; p_ep_f = NULL; }
	if (p_ep_m != NULL) { delete p_ep_m; p_ep_m = NULL; }
	if (p_dist_f != NULL) { delete p_dist_f; p_dist_f = NULL; }
	if (p_dist_m != NULL) { delete p_dist_m; p_dist_m = NULL; }
}