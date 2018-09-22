/*
 * PCell.h
 *
 *  Created on: Jan 8, 2011
 *  Persistent cell class, to be coded
 *      Author: Vidit
 */

#ifndef PCELL_H_
#define PCELL_H_

#define PCCHAIN Chain<C,PCell<C,BT>*>
# include "Cell.hpp"

template <typename C, typename BT = int>
class PCell
{
private:
	bool marked;
public:
	num dim;
// Persistence Data!
	BT birth;
	PCCHAIN bdry;
	PCCHAIN* TLIST;
	num kindex;

	// more hackery....
	CCHAIN* kgen;
	//BT death;		// death time
	//num T_Index;	// uh...

// functions:
	// constructor!
	PCell()
	{
		dim = 0;
		marked = false;
		bdry.clear();
		birth = INITBT;
		TLIST = NULL;
		kindex = 0;
		//birth = INITBT;
		//death = INITBT;
		//T_Index = 0;
		kgen = NULL;
	}

	// copy constructor
	PCell(const PCell<C,BT>& other)
	{
		dim = other.dim;
		TLIST = other.TLIST;
		bdry = other.getBD();
		kindex = other.kindex;

		kgen = other.kgen;
	}

	~PCell()
	{
		bdry.clear();
		if (TLIST != NULL)
		{
			delete TLIST;
			TLIST = NULL;
		}
	}

	friend ostream& operator << (ostream& out, const PCell<C,BT>& toprint)
	{
		toprint.print(out);
		return out;
	}
//functions
	num getInd()const {return kindex;}
	void makeFromCell(const Cell<C,BT>*);
	void print(ostream&) const;
	num getDim() const;
	const PCCHAIN& getBD() const;
	const PCCHAIN& getCB() const; // DUMMY: DO NOT USE.
	bool addBdryLink(PCell<C,BT>*&, const C&);
	bool isMarked() const {return marked;}
	void mark(){marked = true;}
};

#endif /* PCELL_H_ */
