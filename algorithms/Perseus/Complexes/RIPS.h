/*
 * RIPS.h
 *
 *  Created on: Dec 2, 2010
 *      Author: Vidit
 */

#ifndef RIPS_H_
#define RIPS_H_

#include "SToplex.hpp"
// upper triangular matrix to store point-point neighbor relations
#define RROW vector<pair<num,BT> >
#define RIPMAT vector<RROW*>
// vector of vertices
#define VERTS vector<Point<PS>*>

// structure to store high dimensional (> 0) simplices
#define SLIST vector<Simplex<C,PS,BT>*>
#define INSRIP map<num, SLIST*>


template <typename C=int, typename PS=double, typename BT=int>
class RIPS : public SToplex<C,PS,BT>
{
public:
	// DATA
	RIPMAT isnbr; // contains upper-triangular dist matrix info
	VERTS vert; // contains vertices

	vector<double> rads; // contains radius for each point
	vector<BT> birth; // vector of birth times in case of

	double epsilon, stepsize;
	num numsteps;

	// FUNCTIONS

	RIPS()
	{
		InitRIPS();
	}

	void InitRIPS()
	{
		epsilon = 1; stepsize = 1;
		numsteps = 0;
		isnbr.clear();
		vert.clear();
		rads.clear();
		// in case we begin to allocate points...
		srand(time(NULL));
	}

	~RIPS()
	{
		clear();
	}

	void allocateNbrMatrix();
	num addVertex(const BT&);
	bool addEdge(const num, const num, const BT&);

	BT getEdgeBirth(num, num) const;
	void clear(); // clean rips structure!
	pair<num,bool> makeFromFile(ifstream&); // make rips complex where births are inherited from vertices
	pair<num,bool> makeFromFile_GrowBalls(ifstream&,bool,bool,bool); // make rips complex by growing epsilon balls around points
	pair<num,bool> makeFromDistMatrixFile(ifstream&,bool,bool); // make rips complex from matrix of pairwise distances
	pair<num,bool> makeFromTimeSeriesFile(ifstream&); // make rips complex from matrix of pairwise distances

    num witness(double); // prune vertex via witnesses. call before building nbrmatrix!
	void getLowerNbrs(num, vector<num>&) const; // get lower order neighbors!
	void getLowerNbrs_GrowBalls(num, vector<num>&, vector<BT>&, bool) const;
	void incRIPS(const num); // builds RIPS complex incrementally
	void incRIPS_GrowBalls(const num, const BT&);
	void makeRandom_EdgeShuffle(const num, bool);

	bool AddRIPSCofaces(Simplex<C,PS,BT>*, const vector<num>&, num, INSRIP&, num);
	bool AddRIPSCofaces_GrowBalls(Simplex<C,PS,BT>*, const vector<num>&, const vector<BT>&, num, INSRIP&);

	void makeNbrMatrix(); // populates upper triangular matrix of point neighbor relations
	void makeNbrMatrix_GrowBalls();

	void clearPoints(); // removes memory allocated to vertices

	void showPoints(const vector<num>&) const; // debug help for nbr lists!
	void storeCoface(const Simplex<C,PS,BT>*, INSRIP&) const;
	void insertToToplex(INSRIP&);
	void ComputePersistence(num,string,bool);
	void ComputePersistence_GrowBalls(num,string);

	void ComputePersistence(num ptdim, map<num, vector<pair<BT,BT> > >&);
	void showNbrMatrix(ostream& out) const;
};

#endif /* RIPS_H_ */
