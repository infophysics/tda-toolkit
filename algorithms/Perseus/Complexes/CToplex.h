/*
 * CToplex.h
 *
 * Cubical Toplex Class Definitions
 */

#ifndef CTOPLEX_H_
#define CTOPLEX_H_

#define CUBE_RGMAP map<Point<PS>*, Cell<C,BT>*, ptcomplex<PS> >
#define CUBE_MAP map<vector<num>*, CUBE_RGMAP*, addcomp >
#define PB_MAP map<Point<PS>*, BT, ptcomplex<PS> >
#define BM_MAP multimap<BT,Point<PS>*>

#include <map>
#include <utility>
#include <vector>
#include <set>

#include "../Global/Combinatorics.h"
#include "../Cells/Cell.hpp"
#include "Morse.h"

// vector lexico comparison for map-making
struct addcomp
{
	bool operator () (vector<num>* v1, vector<num>* v2)
	{
		//cout<<" here "; cin.get();
		if (v1 == v2) return false;
		if (v1 == NULL && v2 != NULL) return true;
		if (v2 == NULL) return false;

		if (v1->size() < v2->size()) return true;
		if (v1->size() > v2->size()) return false;
		// sizes are equal
		for (unsigned int i=0; i<v1->size(); i++)
		{
			if (v1->at(i) < v2->at(i)) return true;
			else if (v1->at(i) > v2->at(i)) return false;
			else continue;
		}
		return false; // equal
	}
};


// C is the ring of coefficients to make chains, PS is the space
// from which vertex points are selected.
template <typename C = int, typename PS = int, typename BT = int>
class CToplex
{
public:
	// Data

	// this is the entire structure: we map the addin vector (i.e.
	// set of dimensions in which to increment the anchor point) to
	// another map.
	// this "range map" for a given addin vector maps anchor points
	// to the cube with that anchor and the original addin

	// see definition of CUBE_MAP above!

	CUBE_MAP mymap;

	// Functions

	// Constructors
	CToplex()
	{
		mymap.clear();
	}

	CToplex(const num topdim, const vector<Point<PS>*>& anchors, bool face = false)
	{
		buildTopCubeMap(topdim, anchors);
		if (face) addCubeFaces();
	}

	~CToplex()
	{
		Destroy(); // kill the anchor points and addin vectors, but not the cubes
	}

	friend ostream& operator << (ostream& out, const CToplex<C,PS,BT>& toprint)
	{
		// iterate over map's addin vectors:
		typename CUBE_MAP::const_iterator aditer;
		CUBE_RGMAP* currg;
		// and over anchors:
		typename CUBE_RGMAP::const_iterator rgiter;
		typename vector<num>::const_iterator vit;

		for (aditer = toprint.mymap.begin(); aditer != toprint.mymap.end(); ++aditer)
		{
			currg = aditer->second;
			out<<"Addin: [";
			for (vit = aditer->first->begin(); vit != aditer->first->end(); ++vit) out<<*vit<<" ";
			out<<"]\n";
			for (rgiter = currg->begin(); rgiter != currg->end(); ++rgiter)
			{
				if (rgiter->second == NULL) continue;
				out<<"     Anc: "<<*(rgiter->first)<<": "<<*(rgiter->second)<<"\n";
			}
		}
		return out;
	}


	bool buildTopCubeMap(const num, const vector<Point<PS>*>&, const vector<BT>*); // populate map using anchors
	void addCubeFaces(); // add faces given starting cubes
	bool singleCubeFaces(typename CUBE_MAP::reverse_iterator, typename CUBE_RGMAP::iterator); // builds faces for a single cube
	void insertFacePairOf(Cell<C,BT>*, vector<num>*,Point<PS>*,Point<PS>*,const C); // updates map to add point
	bool makeCellComplex(const num, const vector<Point<PS>*>&, MComplex<C,BT>&); // makes mcomplex
	void writeComplex(MComplex<C,BT>&,bool); // writes from map to complex
	void writeComplexExcept(MComplex<C,BT>&,const BT&,bool); // writes all but cells of given birth
	void Destroy(); // wipes out the structure, and optionally allocated memory
	pair<num,bool> readFileInfo(ifstream&,num&,PB_MAP&);
	pair<num,bool> makeFromFile(ifstream&); // reads toplex info
	MComplex<C,BT>* readMovie(ifstream&); // reads movie toplex info
	void removeAllBut(const BT&);
};



#endif /* CTOPLEX_H_ */
