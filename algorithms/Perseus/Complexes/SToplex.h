/*
 * Toplex.h
 *
 * simplicial and cubical toplex definitions
 */

#ifndef TOPLEX_H_
#define TOPLEX_H_

# include <set>
# include <map>
# include <algorithm>

# include "../DebugV.h"
# include "../Cells/Simplex.hpp"

# include "PComplex.hpp"

# define SIMP_RG multimap<num, Simplex<C,PS,BT>*>
# define SIMP_MAP map<num, SIMP_RG*, decnum>


// C is the ring coefficients over which chains are built, PS
// is the space from which points are allocated
template <typename C=int, typename PS=double, typename BT=int>
class SToplex
{
public:
	// this is really the whole stoplex structure. each simplex is
	// assigned a hash, which sums the hashes of its vertices to a
	// number 'h'.
	// we populate a map sending dimension to a range
	// multimap. the key of this range map is a number, i.e. the hash h.
	// h then forms the key for all simplices of the given dimension
	// which have hash h.

	// see definition of SIMP_MAP and SIMP_RG above!

	SIMP_MAP mymap;
	num topdim; // highest dimension encountered
	num psize; // number of points

	PT_SET verlist; // vertex list


	// methods

	// constructor
	SToplex()
	{
		mymap.clear();
		topdim = psize = 0;
	}
	SToplex(const vector<Simplex<C,PS,BT>*>& topsimps)
	{
		buildTopSimpMap(topsimps);
	}

	virtual ~SToplex()
	{
		Destroy(false);
	}

	friend ostream& operator << (ostream& out, const SToplex<C,PS,BT>& toprint)
	{
		typename SIMP_MAP::const_iterator miter;
		typename SIMP_RG::const_iterator liter;
		out<<"\n++++++++S_TOP PRINT!+++++++++\n";
		for (miter = toprint.mymap.begin(); miter != toprint.mymap.end(); ++miter)
		{
			out<<"dim: ["<<miter->first<<"]:\n"; // print hash
			// iterate over range map
			for (liter = miter->second->begin(); liter != miter->second->end(); ++liter)
			{
				// liter->first would be the hash here
				out<<"       "<<*(liter->second)<<" born "<<liter->second->birth<<" @ hash: "<<liter->first<<'\n';
			}
		}
		return out;
	}

	void printSizeInfo (ostream& out = cout)
	{
		typename SIMP_MAP::const_iterator miter;
		typename SIMP_RG::const_iterator liter;
		out<<"\n++++++++S_TOP SZ PRINT!+++++++++\n";
		for (miter = mymap.begin(); miter != mymap.end(); ++miter)
		{
			out<<"dim: ["<<miter->first<<"]: size = "<<miter->second->size()<<"\n"; // print hash
		}
	}

	num getHashMod() const;
	pair<num,bool> makeFromFile(ifstream&, bool, bool, bool); // create from file
	Simplex<C,PS,BT>* isIn(const PT_SET&); // checks for containment
	// makes cell complex out of a vector of topdim simplices, stores in sec arg
	bool makeCellComplex(const vector<Simplex<C,PS,BT>*>&, MComplex<C,BT>&, bool);
	bool buildTopSimpMap (const vector<Simplex<C,PS,BT>*>&,bool); // populates map, and verlist optionally
	bool addSimpFaces(bool); // adds faces descending from given dimension
	void writeComplex(MComplex<C,BT>&) const; //makes simplicial cell complex
	void Destroy(bool); //cleans out the map structure, optionally deletes simplices!
	void getVerts(PT_SET&); // extracts all vertices;
	num populateMap(const vector<Simplex<C,PS,BT>*>&, bool); // builds map with faces etc.
	bool addCoface(Simplex<C,PS,BT>*); // adds a single simplex provided its boundary simplices already exist
	bool makeFacesOf(Simplex<C,PS,BT>*,bool); // makes faces of single simplex
	num capDimension(map<PT_SET,BT>&,num); // updates simplex vector by capping top dimension!

    bool barySub(SToplex<C,PS,BT>&); // stores barycentric subdivision of complex in input!
    bool subdivide (Simplex<C,PS,BT>*, map<Simplex<C,PS,BT>*,Simplex<C,PS,BT>*>&,set<PT_SET>&); // bary-subdivides a simplex
    bool insertVertex (Simplex<C,PS,BT>*&); // insert vertex: that is, zero dim simplex

	bool topWrapper(const vector<Simplex<C,PS,BT>*>&, string );

	pair<num,bool> makeRicci (ifstream&, ifstream&);
	bool readEdgeBirths(ifstream&, map<Point<PS>*, num, ptcomplex<PS> >&);
	bool inheritBirths(const num);
};

# endif // TOPLEX_H
