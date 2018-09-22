/*
 * PCell.hpp
 * Contains Function Definitions for PCell.h
 */

#ifndef CELL_HPP_
#define CELL_HPP_

# include "../DebugV.h"
# include "Cell.h"

/***********************************************
 *  Cell FUNCTIONS
 ***********************************************/

// basic

// ghetto unit check for integers... eventually this will
// be a property of the ring over which we are doing the
// algebra
template <typename C>
bool isIntUnit(const C& coeff)
{
	return (coeff==1 || coeff == -1);
}


// initialize persistence data DEPRECATED
template <typename C, typename BT>
void Cell<C,BT>::InitP()
{

}

// copy persistence data DEPRECATED
template <typename C, typename BT>
void Cell<C,BT>::CopyP(const Cell<C,BT>& other)
{

}


// adds a link (i.e. cell-coefficient pair given by toadd-coeff) to the boundary of the calling
// cell. if propdownbt is enabled then toadd's birth time is suitably decremented to be <= the
// calling cell's birth time. if mkcb is enabled, we also add the calling cell to the coboundary
// of toadd
template <typename C, typename BT>
void Cell<C,BT>::addBDLink(Cell<C,BT>* toadd, const C& coeff, bool mkcb = true, bool propdownbt = true)
{
	//cout<<"\n\n bdlinker: for calling cell: "<<*this;
	if (coeff != 0 && toadd != NULL)
	{
		//cout<<" to bd cell "<<*toadd<<"\n";
		bd.addLink(toadd, coeff);
		if(mkcb) toadd->cb.addLink(this, coeff);
		if (propdownbt && ((toadd->birth > birth) || toadd->birth == INITBT)) toadd->birth = birth;
		//cout<<"****bd is now: "<<this->getBD();
	}

}


// adds a link (i.e. cell-coefficient pair given by toadd-coeff) to the coboundary of the calling
// cell. if propupbt is enabled then toadd's birth time is suitably incremented to be >= the
// calling cell's birth time. if mkbd is enabled, we also add the calling cell to the boundary
// of toadd
template <typename C, typename BT>
void Cell<C,BT>::addCBLink(Cell<C,BT>* toadd, const C& coeff, bool mkbd = true, bool propupbt = true )
{
	if (coeff != 0 && toadd != NULL)
	{
		cb.addLink(toadd, coeff);
		if (mkbd) toadd->bd.addLink(this,coeff);
		if (propupbt && ((toadd->birth < birth) || toadd->birth == INITBT)) toadd->birth = birth;
	}
}

// wipes out cell
template <typename C, typename BT>
void Cell<C,BT>::clear()
{
    bd.clear();
	cb.clear();
	Destroy();
	//cout<<"\n    wiped cell clean!"; cin.get();
}

// erases dynamically allocated components (gradient path) from memory
template <typename C, typename BT>
void Cell<C,BT>::Destroy()
{
	if (gpath != NULL) {gpath->clear(); delete gpath;}
	if (kgen != NULL) {kgen->clear(); delete kgen;}
}

template <typename C, typename BT>
void Cell<C,BT>::showGens(ostream& out = cout, bool samebirth = false, const string buf = "------ ") const
{
}

// prints cell as a triple of dimension, index and birth time
template <typename C, typename BT>
void Cell<C,BT>::print(ostream& out) const
{
	out<<"["<<birth<<";"<<getDim()<<","<<ind<<"]";
	if (isCrit) out<< " (*)";
}

// sets calling cell's boundary to given chain toset. if mkcb is enabled, we also
// update the coboundary of each cell in the support of toset to include the calling
// cell
template <typename C, typename BT>
void Cell<C,BT>::setBD(const CCHAIN& toset, bool mkcb=true)
{
	// first set boundary to equal the argument chain...
	bd = toset;

	// and every cell in this boundary chain now gets its coboundary
	// changed to include the calling cell
	if (mkcb)
	{
		typename CSTRUCT::const_iterator i;
		for (i = getBD().begin(); i != getBD().end(); ++i)
		{
			i->first->cb.cells[this] = i->second;
		}
	}
}

// create co-boundary chain
template <typename C, typename BT>
void Cell<C,BT>::setCB(const CCHAIN& toset, bool mkbd=true)
{
	// first set coboundary equal to argument chain
	cb = toset;

	// and every cell in this coboundary chain now gets its boundary
	// changed to include the calling cell
	if (mkbd)
	{
		typename CSTRUCT::const_iterator i;
		for (i = cb.begin(); i != cb.end(); ++i)
		{
			i->first->bd.cells[this] = i->second;
		}
	}
}

// assignment operator
template <typename C, typename BT>
Cell<C,BT>& Cell<C,BT>::operator = (const Cell<C,BT>& other)
{
	dim = other.dim; ind = other.ind;
	bd = other.getBD();
	cb = other.getCB();
	// copy persistence junk as well?
	this->CopyP(other);
	return *this;
	// maybe also morse data?
}
// comparison operator
template <typename C, typename BT>
bool Cell<C,BT>::operator == (const Cell<C,BT>& other) const
{
	return (dim==other.dim && ind==other.ind && birth==other.birth);
}

// checks if face admits a free face collapse, optionally checks for
// birth time equality in case of persistence
template <typename C, typename BT>
bool isFreeFace(const Cell<C,BT>* face, Cell<C,BT>*& coface, bool checkbt=true)
{
	if (face == NULL) return false;
	if (face->getCB().size()==1) // otherwise not a free face!
	{
		// check for unit incidence!
		coface = face->getCB().begin()->first;
		if (isIntUnit(face->getCB().begin()->second))
		{
			// if the coface has other boundary elements, do not collapse!s
			if(checkbt)
			{
				if(face->birth != coface->birth) return false;
			}
			return true;
		}
	}
	return false;
}

// checks if face admits a free coface collapse, and modifies
// optionally checks for birth time equality in case of persistence
// checks if face admits a free face collapse, optionally checks for
// birth time equality in case of persistence
template <typename C, typename BT>
bool isFreeCoface(const Cell<C,BT>* coface, Cell<C,BT>*& face, bool checkbt=true)
{
	if (coface == NULL) return false;
	if (coface->getBD().size()==1) // otherwise not a free coface!
	{
		// check for unit incidence!
		face = coface->getBD().begin()->first;
		if (isIntUnit(coface->getBD().begin()->second))
		{
			// if the coface has other boundary elements, do not collapse!s
			if(checkbt)
			{
				if(face->birth != coface->birth) return false;
			}
			return true;
		}
	}
	return false;
}

template <typename C, typename BT>
bool isKQPair(const Cell<C,BT>* k, const Cell<C,BT>* q, bool checkbt=true)
{
	// nothing serious here now, but eventually:
	// check if the incidence coefficient is a unit!

	//cout<<"\n check kqp... ";
	if (k->getBD().size()==1) // if there is only one boundary element of k
	{
		//cout<<" ok size... ";
		if (q == k->getBD().begin()->first) // and if it equals q
		{
			// replace with UNIT INCIDENCE TEST!
			if (isIntUnit(k->getBD().begin()->second))
			{
				if(checkbt)
				{
					if(k->birth != q->birth) return false;
				}
				return true;
			}
		}
	}
	return false;
}

template <typename C, typename BT>
bool isQKPair(const Cell<C,BT>* q, const Cell<C,BT>* k, bool checkbt=true)
{
	// nothing serious here now, but eventually:
	// check if the incidence coefficient is a unit!

	//cout<<"\n check kqp... ";
	if (q->getCB().size()==1) // if there is only one coboundary element of q
	{
		//cout<<" ok size... ";
		if (k == q->getCB().begin()->first) // and if it equals k
		{
			// replace with UNIT INCIDENCE TEST!
			if (isIntUnit(q->getCB().begin()->second))
			{
				if(checkbt)
				{
					if(q->birth != k->birth) return false;
				}
				return true;
			}
		}
	}
	return false;
}



/********************************************
 *  Chain Functions
 ********************************************/

// adds a link consisting of a cell poiner and a coefficient to
// the calling chain
template <typename C, typename D>
bool Chain<C,D>::addLink (const D& ptr, const C& coeff) // add a single cell with coefficient given
{
	if (ptr == NULL) return false;
	if (coeff != 0)
	{
		typename CHAIN::iterator i;
		// if the cell is in the support of this chain...
		if(isIn(ptr,i))
		{
			// increment the coefficient of that cell
			i->second += coeff;
			if (i->second == 0) // but remove if zero!
			{
			    cells.erase(i);
			}
		}
		// otherwise, just add it!
		else cells.insert(make_pair( (D)ptr, (C) coeff ) );
	}

	return true;
}

// assignment operator
template <typename C, typename D>
Chain<C,D>& Chain<C,D>::operator = (const Chain<C,D>& other)
{
	if (this != &other) // check self-assignment
	{
		cells.clear();
		cells = other.cells;
	}
	return *this;
}

// addition operator
template <typename C, typename D>
const Chain<C,D> Chain<C,D>::operator + (const Chain<C,D>& other) const// add
{
	Chain<C,D> result(*this);
	result += other;
	return result;
}

// subtraction operator
template <typename C, typename D>
const Chain<C,D> Chain<C,D>::operator - (const Chain<C,D>& other) const // subtract
{
	Chain<C,D> result = *this;
	result -= other;
	return result;
}

// add and assign
template <typename C, typename D>
Chain<C,D>& Chain<C,D>::operator += (const Chain<C,D>& other) // add and assign
{
	typename CHAIN::const_iterator i;
	for(i=other.cells.begin(); i != other.cells.end();++i)
	{
		this->addLink(i->first,i->second);
	}
	return *this;
}

// subtract and assign
template <typename C, typename D>
Chain<C,D>& Chain<C,D>::operator -= (const Chain<C,D>& other) // subtract and assign
{
	typename CHAIN::const_iterator i;
	for(i=other.cells.begin(); i != other.cells.end();++i)
	{
		this->addLink(i->first,-i->second);
	}
	return *this;
}

// returns chain multiplied by a scalar coefficient
template <typename C, typename D>
Chain<C,D> Chain<C,D>::scale (const C& coeff) const // scale by coefficient
{
	Chain<C,D> toret;
	typename CHAIN::const_iterator i;
	for(i=cells.begin(); i != cells.end(); ++i)
	{
		toret.addLink(i->first,i->second * coeff);
	}
	return toret;
}

template <typename C, typename D>
void Chain<C,D>::makeZp(num p = 2) // make Z_p version: kills even coeff cells and makes odd coeff = 1
{
    vector<typename CHAIN::iterator> tokill; // store everything to kill

    typename CHAIN::iterator i;
	for(i=cells.begin(); i != cells.end(); ++i)
	{
		i->second = i->second % p; // assign parity
		if (i->second == 0) tokill.push_back(i); // we will kill this soon...
	}

	// loop over iterators to links that we want to remove
	typename vector<typename CHAIN::iterator>::iterator kit;
	for (kit = tokill.begin(); kit != tokill.end(); ++kit)
	{
	    cells.erase(*kit); // and remove them...
	}
}

// scales this chain directly
template <typename C, typename D>
void Chain<C,D>::scaleMe (const C& coeff) // scale by coefficient
{
    typename CHAIN::iterator i;
	for(i=cells.begin(); i != cells.end(); ++i)
	{
		i->second *= coeff;
	}
}

// removes the cell pointed to by tokill if found,
// otherwise returns false...
template <typename C, typename D>
bool Chain<C,D>::removeCell(const D& tokill)
{
	if (tokill == NULL) return false;
	return cells.erase((D)tokill);
}

template <typename C, typename D>
num Chain<C,D>::removeUnmarked(bool trace=false)
{

	// vector of iterators to ALL unmarked cells in the chain!
	vector<typename CHAIN::iterator> unmit;

	// loop over chain...
	typename CHAIN::iterator i; // to iterate over cells
	for (i=cells.begin(); i!=cells.end(); ++i)
	{
		if(trace) cout<<"\n Check: "<<*(i->first)<<" at "<<i->first;

		if (!(i->first->isMarked())) // if this cell is unmarked,...
		{
			if (trace) cout<<" UNmarked!!!";
			unmit.push_back(i); // store iterator in vector...
		}
		else
		{
			if (trace) cout<<" \n marked: "<<*(i->first)<<" at "<<i->first;
		}
	}
	// now loop over vector of iterators to remove and remove them!
	typename vector<typename CHAIN::iterator>::iterator vit;
	for (vit = unmit.begin(); vit!= unmit.end(); ++vit)
	{
		cells.erase(*vit);
	}
	return (num) unmit.size();
}

// for persistence, used by REMOVE_PIVOT_ROW, sadly linear in size
template <typename C, typename D>
num Chain<C,D>::getMaxIndex(D& maxcell) const
{
	num curmax = 0;
	typename CHAIN::const_iterator link;
	for (link = begin(); link != end(); ++link)
	{
		if (link->first->getInd() > curmax)
		{
			maxcell = link->first;
			curmax = maxcell->getInd();
		}
	}
	return curmax;
}

// returns coefficient of argument cell, 0 if not in support
template <typename C, typename D>
C Chain<C,D>::getCoeff(const D& tofind) const
{
	typename CHAIN::const_iterator locate = cells.find(tofind);
	if (locate == end()) return 0;
	return locate->second;
}

// returns true if the cell pointed to by the argument lies in the
// support of the calling chain, false otherwise. It modifies pos
// to the position of the cell if found
template <typename C, typename D>
bool Chain<C,D>::isIn(const D& tofind, typename CHAIN::iterator& pos)
{
	// can't find squat if it is null
	if (tofind == NULL) return false;
	// find by key!
	pos = cells.find((D)tofind);
	return (pos == cells.end()) ? false : true;
}

template <typename C, typename D>
void Chain<C,D>::clear()
{
	cells.clear();
}

template <typename C, typename D>
Chain<C,D> Chain<C,D>::getBD() const
{
	Chain<C,D> toret;
	typename CHAIN::const_iterator curcell;
	for (curcell = begin(); curcell != end(); ++curcell)
	{
		toret += curcell->first->getBD().scale(curcell->second);
	}
	return toret;
}

template <typename C, typename D>
Chain<C,D> Chain<C,D>::getCB() const
{
	Chain<C,D> toret;
	typename CHAIN::const_iterator curcell;
	for (curcell = begin(); curcell != end(); ++curcell)
	{
		toret += curcell->first->getCB().scale(curcell->second);
	}
	return toret;
}

template <typename C, typename D>
num Chain<C,D>::getDim() const
{
    if (size()==0) return 0; // should return -1, really.
    return begin()->first->getDim();
}

template <typename C, typename D>
void Chain<C,D>::Destroy(bool killcell = false)
{
	// i iterates over map structure
	if (killcell)
	{	typename CHAIN::iterator i;
		for (i=begin(); i!= end(); ++i)
		{
			// remove this cell!
			delete i->first;
		}
		// and clear map
	}
	clear();
}

// return number of non-unit-coefficients.
template <typename C, typename D>
num Chain<C,D>::nonUnitCount() const
{
    num toret = 0;
    typename CHAIN::const_iterator i;

    for (i = begin(); i != end(); ++i)
    {
        if (!isIntUnit(i->second)) toret++;
    }
    return toret;
}


// return number of unit-coefficients.
template <typename C, typename D>
num Chain<C,D>::unitCount() const
{
    num toret = 0;
    typename CHAIN::const_iterator i;

    for (i = begin(); i != end(); ++i)
    {
        if (isIntUnit(i->second)) toret++;
    }
    return toret;
}

// returns the number of cells in this chain
// whose coboundaries have size 2
template <typename C, typename D>
num Chain<C,D>::cbsum() const
{
    double toret = 0;
    typename CHAIN::const_iterator i;

    for (i = begin(); i != end(); ++i)
    {
        //if (i->first->getCB().size() == 2) toret++;
        toret += (i->first->getCB().size() - 1);
    }
    return toret;
}

// returns the number of cells in this chain
// whose boundaries have size 2
template <typename C, typename D>
num Chain<C,D>::bdsum() const
{
    double toret = 0;
    typename CHAIN::const_iterator i;

    for (i = begin(); i != end(); ++i)
    {
        //if (i->first->getBD().size() == 2) toret++;
        toret += (i->first->getBD().size() - 1);
    }
    return toret;
}

#endif /* CELL_HPP_ */
