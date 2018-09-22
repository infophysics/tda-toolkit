/*
 * PComplex.h
 *
 *  Persistent Complex
 *      Author: Vidit
 */

#ifndef PCOMPLEX_H_
#define PCOMPLEX_H_

# include "../Global/Global.h"
# include "../Cells/PCell.h"

// vector of persistent cells!
# define PCVEC vector<PCell<C,BT>*>
// store betti numbers for each frame
#define BETTI_STR map<BT, vector<BT>*>

// store (birth, death) intervals for each dimension
#define INTVEC vector<pair<BT,BT> >
#define INT_STR map<num, INTVEC*>

#define PCMAP map<Cell<C,BT>*, PCell<C,BT>*, CELLORD >

template <typename C = int, typename BT = int>
class PComplex
{
public:
	// complex data: persistent cells ordered by birth +  dimension
	PCVEC klist;
	INT_STR ints;
	BETTI_STR betti;

	// functions:
	PComplex()
	{
		ints.clear();
		betti.clear();
		klist.clear();
	}

	~PComplex()
	{
		destroyCells();
		destroyPersData();
	}


	friend ostream& operator << (ostream& out, const PComplex<C,BT>& toprint)
	{
		out<<"\n\nPCOMPLEX PRINTER...";
		for (num i = 0; i < (num) toprint.klist.size(); ++i)
		{
			out<<"\n   "<<*(toprint.klist.at(i))<<" bd "<<toprint.klist.at(i)->getBD();
		}
		return out;
	}

	num makeFromComplex(const Complex<C,BT>&, bool);
	bool makeBDChains (PCMAP&,bool);

	// Persistence Routines
	void showBetti(ostream&); // shows betti numbers for each frame to std out
	map<num, vector<pair<BT,BT> > > getInts() const; // returns persistence intervals
	void showInts(ostream&); // shows persistence intervals for all dimensions
	bool makeOutputFiles(const string&); // writes output to files
	void initPersData(const Complex<C,BT>&);
	num incrementBetti(const BT&, const BT&, const num&);
	void destroyCells();
	void destroyPersData();
	void REMOVE_PIVOT_ROW(const PCell<C,BT>*, PCCHAIN&, bool);
	void COMPUTE_INTERVALS(const Complex<C,BT>&, bool, bool);

};

#endif /* PCOMPLEX_H_ */
