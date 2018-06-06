/* Subgraph Matching Kernel (Supplementary Functions)
 * Author: Ioannis Siglidis <y.siglidis@gmail.com>
 * License: BSD 3 clause"
 * Code taken from: http://www.partow.net/programming/hashfunctions/#APHashFunction
 */
#include "../include/functions.hpp"
#include <list>
#include <cmath>
#include <stdlib.h>

using namespace std;

double *cv;
double **ce;
double *totalValue;
unsigned int k;

void sm_core(double value, list<int> c, list<int> p, int* d, int lBound, int uBound) {
	
	while (!p.empty()) {
		int i = p.front();
		p.pop_front();

		double nValue = value * cv[i];
		double* iEdgeValue = ce[i];
		for (list<int>::const_iterator it = c.begin(); it != c.end(); it++) {
			nValue *= abs(iEdgeValue[*it]);
		}

		totalValue[c.size()] += nValue;

		if (c.size()+1 < k) {
			c.push_back(i);
			
			// prepare candidate set for recursive call
			list<int> newP;
			for (list<int>::const_iterator it = p.begin(); it != p.end(); it++) {
				int v = *it;
				if (iEdgeValue[v] != 0)
					newP.push_back(v);
			}
			
			int newUBound = uBound;
			int newLBound = lBound;
			if (lBound <= uBound) {
				int tmp;
				while (iEdgeValue[d[newUBound]] == 0 && --newUBound > lBound);
				while (iEdgeValue[d[newLBound]] == 0 && ++newLBound < newUBound);
				int nm = newLBound-1;
				while (++nm <= newUBound) {
					if (iEdgeValue[d[nm]] < 0) continue;
					if (iEdgeValue[d[nm]] > 0)
						newP.push_back(d[nm]);
					// swap
					tmp = d[newLBound];
					d[newLBound] = d[nm];
					d[nm] = tmp;
					newLBound++;
				}
			}
			
			sm_core(nValue, c, newP, d, newLBound, newUBound);
			c.pop_back();
		}
	}
}
	

void sm_core_init(double value, int* d, int nv, int kappa, double *cost_vertices, double **cost_edges, double *total_value) {
	
	cv = cost_vertices;
	ce = cost_edges;
	totalValue = total_value;
	k = (unsigned int) kappa;

	list<int> c;
	int lBound = 0;
	int uBound = nv-1;

	for (int it=lBound; it<=uBound; it++) {
		int i = d[it];
		double nValue = value * cv[i];
		
		totalValue[0] += nValue;

		if (k > 1) {
			c.push_back(i);
			
			// prepare candidate set for recursive call
			list<int> p;
			int tmp;
			int newUBound = uBound;
			double *iEdgeValue = ce[i];
			while (iEdgeValue[d[newUBound]] == 0 && --newUBound > it);
			int newLBound = it;
			while (++newLBound <= newUBound && iEdgeValue[d[newLBound]] == 0);
			int nm = newLBound-1;
			while (++nm <= newUBound) {
				if (iEdgeValue[d[nm]] < 0) continue;
				if (iEdgeValue[d[nm]] > 0)
					p.push_back(d[nm]);
				// swap
				tmp = d[newLBound];
				d[newLBound] = d[nm];
				d[nm] = tmp;
				newLBound++;
			}
			
			sm_core(nValue, c, p, d, newLBound, newUBound);
			c.pop_back();
		}
	}
}
