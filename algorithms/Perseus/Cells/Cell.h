/*
 * PCell.h
 * Contains Class Definitions for Cells
 */


#ifndef CELL_H_
#define CELL_H_

// how to order cells!
#define CELLORD cellord<C,BT>


// define a chain to be a map from cell pointer to coefficient... so the chain 2A + 3B
// maps the cell pointer A to coefficient 2 and B to coefficient 3. This makes searching
// and removal and all chain operations logarithmic
#define CHAIN map<D,C>
#define CCHAIN Chain<C,Cell<C,BT>*>
#define CSTRUCT map<Cell<C,BT>*,C>

# include <vector>
# include <set>
# include <list>
# include <utility>
# include <iostream>
# include <map>

// contains pre-defined values for "default" birth time and "banned" birth time
#include "../Global/Pflags.h"
// defines a size_type-compatible typedef to index complex data structures
#include "../Global/Size.h"

using namespace std;

template <typename C, typename BT> class Chain; // forward declaration, see class Chain below

// this class describes a single cell in a filtered complex
// where C is a class of Ring coefficients and BT is the
// filtration index space
template <typename C, typename BT = int>
class Cell
{
public:
	/* BASIC CELL DATA */
	num dim; // dimension of cell
	num ind; // index into cell complex
	CCHAIN bd; // boundary of cell
	CCHAIN cb;	// coboundary of cell

	Cell<C,BT>* srcell; // source cell in cubical or simplicial data

	/* COMPLEX DATA */
	//typename list<Cell<C,BT>*>::iterator mypos;  // position in complex where inserted

	bool isIn; // is in complex?
	bool isFaced;

	BT birth;		// birth time
	num k_index;	// index into linear persistent complex
	//bool marked;

	/* DISCRETE MORSE THEORY DATA */
	bool isCrit; // is it critical?
	bool isQueued; // added to co-reduction queue?

	CCHAIN* gpath; // descending gradient path, chain with critical cells of co-index 1
	CCHAIN* kgen; // descending king chain to complete the critical cell into a cycle
	//CCHAIN* TLST; // persistence debugger
	/* BASIC FUNCTIONS */

	// silly constructor
	Cell()
	{
		isFaced = false;
		dim = 0; ind = 0;
		bd.clear() ; cb.clear();
		// initialize coreduction flags
		isCrit = isQueued = isIn = false;
		// initialize gradient paths
		gpath = NULL;
		kgen = NULL;
		birth = INITBT;

		srcell = NULL;
	}

	// useful constructor: d = dimension, i = index into cell complex
	Cell(const num d, const num i=0)
	{
		isFaced = false;
		dim = d; ind = i;
		bd.clear(); cb.clear();
    	// morse data
		isCrit = isQueued = false;
		gpath = NULL;
		kgen = NULL;
		birth = INITBT;

		srcell = NULL;
	}

	// copy constructor
	Cell(const Cell<C,BT>& other)
	{
		isFaced = other.isFaced;
		dim = other.dim; ind = other.ind;
		bd = other.bd;
		cb = other.cb;
		isIn = other.isIn;
		//mypos = other.mypos;
		// morse data
		isCrit = other.isCrit;
		isQueued = other.isQueued;
		// shallow copy?
		gpath = other.gpath;
		kgen = other.kgen;
		birth = other.birth;

		srcell = other.srcell;
	}

	// destructor
	virtual ~Cell()
	{
		clear();
	}

	// prints cell... using the "print" function is a standard hack to
	// support "virtual" printing... a class inheriting from Cell can
	// define its own print() function to override output
	friend ostream& operator << (ostream& out, const Cell<C,BT>& toPrint)
	{
		toPrint.print(out);
		return out;
	}

	// wipes out cell
	void clear();
	// erases dynamically allocated components from memory
	void Destroy();
	// print helper... add "virtual" to print cubes or simplices instead!
	virtual void print(ostream& out) const;

	// accessors for dimension and boundary/coboundary chains
	virtual num getDim() const { return (num)dim; }
	virtual num getInd() const {return (num)ind;}
	virtual const CCHAIN& getBD() const { return bd; }
	virtual const CCHAIN& getCB() const { return cb; }

	virtual bool isMarked() const {return false;}
	virtual void showGens(ostream&, bool, const string) const;

	// initialize/copy persistence data
	void InitP();
	void CopyP(const Cell<C,BT>&);

	// modifiers for boundary and coboundary chains
	virtual void setBD(const CCHAIN&, bool);
	virtual void setCB(const CCHAIN&, bool);

	// operators
	virtual bool operator == (const Cell<C,BT>&) const; // compares cells
	virtual Cell<C,BT>& operator = (const Cell<C,BT>&); // copy

	// add link to boundary or coboundary with option to update birth times
	// so that no cell is born _after_ a cell in its boundary
	virtual void addBDLink(Cell<C,BT>* toadd, const C& coeff, bool, bool);
	virtual void addCBLink(Cell<C,BT>* toadd, const C& coeff, bool, bool);
};

/******************* Cell orderings **********************/

// fast and dirty: pointer comparison
template <typename C, typename BT>
struct cellordptr
{
	bool operator () (const Cell<C,BT>* c1, const Cell<C,BT>* c2) const
	{
		return (c1 < c2) ? true : false;
	}
};

// nicer to read, but less efficient: comparison by filtration order
template <typename C, typename BT>
struct cellord
{
	bool operator () (const Cell<C,BT>* c1, const Cell<C,BT>* c2) const
	{
		if (c1 == NULL) return true;
		else if (c2 == NULL) return false;
		if (c1->birth < c2->birth) return true; // birth order
		else if (c1->birth > c2->birth) return false; // birth order
		// reaching here means same birth time
		if (c1->dim < c2->dim) return true; // dim order
		else if (c1->dim > c2->dim) return false; // dim order
		// reaching here means same birth + dim
		if (c1->ind < c2->ind) return true; // ind order
		else if (c1->ind > c2->ind) return false; // ind order
		// reaching here means all is same!
		return false;
	}
};

// comparison by number of unit incidences in the BOUNDARY.
// we want to PREFER cells with fewer unit-boundary incidences
template <typename C, typename BT>
struct bdszord
{
	bool operator () (const Cell<C,BT>* c1, const Cell<C,BT>* c2) const
	{
		if (c1 == NULL) return true;
		else if (c2 == NULL) return false;
		// if there is no difference between number of possible reduction pairs,...
		if (c1->getBD().cbsum() == c2->getBD().cbsum())
        {
            // greedily remove more non-units first as critical
            return (c1->getBD().nonUnitCount() >  c2->getBD().nonUnitCount());
        }
        // otherwise prefer to greedily optimize reduction pairs
        else return (c1->getBD().cbsum() < c2->getBD().cbsum());
	}
};

template <typename C, typename BT>
struct cbszord
{
	bool operator () (const Cell<C,BT>* c1, const Cell<C,BT>* c2) const
	{
		if (c1 == NULL) return true;
		else if (c2 == NULL) return false;
		// if there is no difference between possible coreduction pairs...
        if (c1->getCB().bdsum() == c2->getCB().bdsum())
        {
            return (c1->getCB().nonUnitCount() > c2->getCB().nonUnitCount());
        }
        // otherwise greedily optimize coreduction pairs
        else return (c1->getCB().bdsum() < c2->getCB().bdsum());
	}
};

// comparison by number of unit incidences in the BOUNDARY.
// we want to PREFER cells with fewer unit-boundary incidences
template <typename C, typename BT>
struct unitbdord
{
	bool operator () (const Cell<C,BT>* c1, const Cell<C,BT>* c2) const
	{
		if (c1 == NULL) return true;
		else if (c2 == NULL) return false;
		// if unit inc count in bd of c2 exceeds that of c2 return true
		if (c1->getBD().unitCount() <= c2->getBD().unitCount()) return true;
		// in the opposite case return false
		if (c1->getBD().unitCount() > c2->getBD().unitCount()) return false;
		// if they are equal, prefer the one with the larger sized boundary:
		return (c1->getBD().size() < c2->getBD().size());
	}
};

// comparison by number of unit incidences in the BOUNDARY.
// we want to PREFER cells with fewer unit-boundary incidences
template <typename C, typename BT>
struct unitcbord
{
	bool operator () (const Cell<C,BT>* c1, const Cell<C,BT>* c2) const
	{
		if (c1 == NULL) return true;
		else if (c2 == NULL) return false;
		// if unit inc count in bd of c2 exceeds that of c2 return true
		if (c1->getCB().unitCount() <= c2->getCB().unitCount()) return true;
		// in the opposite case return false
		if (c1->getCB().unitCount() > c2->getCB().unitCount()) return false;
		// if they are equal, prefer the one with the larger sized boundary:
		return (c1->getCB().size() < c2->getCB().size());
	}
};

// comparison by number of unit incidences in the BOUNDARY.
// we want to PREFER cells with fewer unit-boundary incidences
template <typename C, typename BT>
struct unitcbordold
{
	bool operator () (const Cell<C,BT>* c1, const Cell<C,BT>* c2) const
	{
		if (c1 == NULL) return true;
		else if (c2 == NULL) return false;
		// if sizes are equal...
		if (c1->getCB().size() == c2->getCB().size())
        {
                return (c1->getBD().unitCount() <= c2->getBD().unitCount());
        }
        // otherwise prefer the big one
		return (c1->getCB().size() > c2->getCB().size());
	}
};

// comparison by number of unit incidences in the CO-BOUNDARY
template <typename C, typename BT>
struct unitbdordold
{
	bool operator () (const Cell<C,BT>* c1, const Cell<C,BT>* c2) const
    {
        if (c1 == NULL) return true;
		else if (c2 == NULL) return false;
		//if sizes are equal...

		// if sizes are equal...
       if (c1->getBD().size() == c2->getBD().size())
        {
                return (c1->getCB().unitCount() < c2->getCB().unitCount());
        }
        // otherwise prefer the big one
		return (c1->getBD().size() > c2->getBD().size());
	}
};

// This class describes C-linear combinations of objects of
// type D, where D must be pointers to objects that have:
// 1. num k_index, for "getMaxIndex"
// 2. function getBD(), getCB() returning chains over D
// 3. bool marked
// 4. print(ostream& out) defined
template <typename C, typename D>
class Chain
{
public:
	/* BASIC DATA */

	// chain structure, see definition above. we map cell pointers to coefficients
	CHAIN cells;

	/* BASIC FUNCTIONS */

	// constructors
	Chain()
	{
		cells.clear();
	}
	Chain(const CHAIN& c)
	{
		cells = c;
	}

	Chain(const Chain<C,D>& other)
	{
		cells.clear();
		cells = other.cells;
	}
	virtual ~Chain()
	{
		clear();
	}

	// output function to print chains
	friend ostream& operator <<(ostream& out, const Chain<C,D>& toPrint)
	{
		typename CHAIN::const_iterator i;
		for(i = toPrint.begin(); i != toPrint.end(); ++i)
		{
			out<<"("<<i->second<<")"<<*(i->first); //<<" @ "<<i->first;
			out<<" + ";
		}
	return out;
	}

	/* CHAIN MATH */
	virtual bool addLink (const D&, const C&); // add a single cell with coefficient given
	virtual Chain<C,D>& operator = (const Chain<C,D>& other);  // assign
	virtual const Chain<C,D> operator + (const Chain<C,D>&) const;  // add
	virtual const Chain<C,D> operator - (const Chain<C,D>&) const; // subtract
	virtual Chain<C,D>& operator += (const Chain<C,D>&); // add and assign
	virtual Chain<C,D>& operator -= (const Chain<C,D>&); // subtract and assign
	virtual Chain<C,D> scale (const C&) const; // returns chain scaled by coefficient
	virtual void scaleMe(const C&); // scales chain directly
	virtual void makeZp(num);

	virtual Chain<C,D> getBD() const; // returns boundary of chain
	virtual Chain<C,D> getCB() const; // returns co-boundary of chain


	// other functions
	num size() const {return cells.size();}
	bool removeCell(const D&); // removes cell
	bool isIn(const D&, typename CHAIN::iterator&); // finds cell

	// iterators for chains!!
	typename CHAIN::const_iterator begin() const {return cells.begin();}
	typename CHAIN::const_iterator end() const {return cells.end();}

	// eradicate chain data
	void clear();

	// get coefficient of given cell... 0 if not found
	C getCoeff(const D&) const;

	// removes chains as well as pointers to incident cells
	void Destroy(bool);

	// removes unmarked cells for persistence algo
	// to use this, typename D must have isMarked() function!
	num removeUnmarked(bool);
	// returns cell of maximal index
	num getMaxIndex(D&) const;
	num getDim() const; // get dimension of chain???

	num unitCount() const; // get number of unit incidence coefficients in chain!
	num nonUnitCount() const; // count non unit coeffs
    num bdsum() const;
    num cbsum() const;
};

#endif /* CELL_H_ */
