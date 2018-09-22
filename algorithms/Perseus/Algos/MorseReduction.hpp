/*
 * Morse.h
 * Routines to perform Discrete Morse Theoretic reductions via free coface
 * collapses of a Cell Complex
 */

#ifndef MORSERED_HPP_
#define MORSERED_HPP_

#include "../DebugV.h"
#include "../Complexes/Morse.hpp"

// processes a given cell as critical during co-reduction
template <typename C, typename BT>
void MComplex<C,BT>::makeCritical_Cored(Cell<C,BT>* tocrit, bool makekgen = false)
{
	tocrit->isCrit = true;

	if (GPTALK)
	{
		debout<<" crit gpath: ";
		if (tocrit->gpath==NULL) debout<<" NULL";
		else debout<< *(tocrit->gpath);
	}

	// save this cell's address in all remaining (unit) co-boundaries' critbd
	updateGpathNearAce_Cored(tocrit);

	// if making generator... well, this information is almost already in kgen of the ace!
	// just need to add the ace itself with coefficient 1. Bam!


	if (makekgen)
	{
		if (tocrit->kgen == NULL) // maybe there are no descending gradient paths!
		{
			tocrit->kgen = new CCHAIN;
		}
		tocrit->kgen->addLink(tocrit,1);
		if (KGTALK) debout<<"\n setting kgen of ace "<<*tocrit<<" to "<<tocrit->kgen;
	}


	typename CLIST::iterator pos;
	// and now remove this cell from the complex as critical...
	//cout<<" removing... "; cin.get();
	removeCell(tocrit, pos, true);
	//cout<<"done!"; cin.get();
}

// processes a given cell as critical during reduction
template <typename C, typename BT>
void MComplex<C,BT>::makeCritical_Red(Cell<C,BT>* tocrit, bool makekgen = false)
{
	tocrit->isCrit = true;

	if (GPTALK)
	{
		debout<<" crit gpath: ";
		if (tocrit->gpath==NULL) debout<<" NULL";
		else debout<< *(tocrit->gpath);
	}

	// save this cell's address in all remaining (unit) co-boundaries' critbd
	updateGpathNearAce_Red(tocrit);

	// if making the generator chain... start populating with the current ace, coefficient 1.

	if (makekgen)
	{
		if (tocrit->kgen == NULL) // maybe there are no ascending gradient paths!
		{
			tocrit->kgen = new CCHAIN;
		}
		else tocrit->kgen->clear();
		tocrit->kgen->addLink(tocrit,1);
		if (KGTALK) debout<<"\n setting kgen of ace "<<*tocrit<<" to "<<tocrit->kgen;
	}


	typename CLIST::iterator pos;
	// and now remove this cell from the complex as critical...
	//cout<<" removing... "; cin.get();
	removeCell(tocrit, pos, true);
	//cout<<"done!"; cin.get();
}


// handles gradient path relaying to coboundaries of ace (i.e. critical cell)
template <typename C, typename BT>
void MComplex<C,BT>::updateGpathNearAce_Cored(Cell<C,BT>* ace)
{
	// save this cell's address in all remaining (unit) co-boundaries' critbd
	typename CSTRUCT::const_iterator cbc;
	for (cbc = ace->getCB().begin(); cbc != ace->getCB().end(); ++cbc)
	{
		//if (!isIntUnit(cbc->second)) continue;
		if (cbc->first->gpath == NULL)
		{
			cbc->first->gpath = new CCHAIN;
			if (GPMEMTALK) debout<<"\n      Assigning new critbd for: "<<*(cbc->first);
			//cin.get();
		}
		cbc->first->gpath->addLink(ace,cbc->second);

		if (GPTALK)
		{
			debout<<"\n    adding to critbd of: "<<*(cbc->first)<<" which is now "<<*(cbc->first->gpath);
			//cin.get();
		}
	}
}


// handles gradient path relaying to coboundaries of ace (i.e. critical cell)
template <typename C, typename BT>
void MComplex<C,BT>::updateGpathNearAce_Red(Cell<C,BT>* ace)
{
	// save this cell's address in all remaining (unit) co-boundaries' critbd
	typename CSTRUCT::const_iterator bdc;
	for (bdc = ace->getBD().begin(); bdc != ace->getBD().end(); ++bdc)
	{
		//if (!isIntUnit(cbc->second)) continue;
		if (bdc->first->gpath == NULL)
		{
			bdc->first->gpath = new CCHAIN;
			if (GPMEMTALK) debout<<"\n      Assigning new critbd for: "<<*(bdc->first);
			//cin.get();
		}
		bdc->first->gpath->addLink(ace,bdc->second);

		if (GPTALK)
		{
			debout<<"\n    adding to critcb of: "<<*(bdc->first)<<" which is now "<<*(bdc->first->gpath);
			//cin.get();
		}
	}
}


// performs ffc's: call AFTER markUncritical_Red!
template <typename C, typename BT>
num MComplex<C,BT>::freeCofaceCollapser()
{
    // find the frame:
    typename COMPLEX::iterator bit;

    FRAME* curf = bit->second;
    typename FRAME::iterator fiter;

    CLIST* curlist;
    typename CLIST::iterator liter;

    Cell<C,BT>* curcell;  // current cell
    Cell<C,BT>* freecoface; // potential free face?

    num dim;
    num toret = 0; // how many cell pairs removed?
    C coeff; // incidence coefficient


    vector<typename CLIST::iterator> tokill;
    typename CLIST::iterator killhere;

    typename vector<typename CLIST::iterator>::iterator removeus;
    BT birth;

    //cout<<"\n collapser"; cin.get();
    for (bit = this->clist.begin(); bit != this->clist.end(); ++bit)
    {
        birth = bit->first;
        curf = bit->second;
        // loop forwards over frame, meaning bottom dimension first.
        for (fiter = curf->begin(); fiter != curf->end(); ++fiter)
        {
            dim = fiter->first;
            if (dim == (this->getTopDim(birth)).first) break; // no cofaces of top dimension cells

            // loop over list of cells of this dim
            curlist = fiter->second;
            //cout<<"\n dim: "<<dim; cin.get();

            for (liter = curlist->begin(); liter != curlist->end(); ++liter)
            {
                curcell = *liter; // extract cell.
                if (curcell == NULL) continue;
                // if boundary has only one face...
                if (curcell->getCB().size()==1)
                {
                    // check unit incidence:
                    coeff = curcell->getCB().begin()->second;
                    if (isIntUnit(coeff))
                    {
                        freecoface = curcell->getCB().begin()->first;

                        if (freecoface->birth != curcell->birth) continue;

                        //cout<<"\n     pair: "<<*curcell<<" and "<<*freeface;
                        toret++;
                        removeCell(curcell,killhere,false);
                        tokill.push_back(killhere);
                        removeCell(freecoface,killhere,false);
                        tokill.push_back(killhere);
                    }
                }
            }
        }

        //cout<<"\n cleaning! "<<toret; cin.get();

        // now clean out list:
        for (removeus = tokill.begin(); removeus != tokill.end(); ++removeus)
        {
            delete **removeus;
            **removeus = NULL;
        }
        tokill.clear();

    }

    return toret;
}


// performs ffc's downwards from given frame: call AFTER markUncritical_Red!
template <typename C, typename BT>
num MComplex<C,BT>::freeFaceCollapser()
{
    // find the frame:
    typename COMPLEX::iterator bit;


    FRAME* curf;
    typename FRAME::reverse_iterator fiter;

    CLIST* curlist;
    typename CLIST::iterator liter;

    Cell<C,BT>* curcell;  // current cell
    Cell<C,BT>* freeface; // potential free face?

    num dim;
    num toret = 0; // how many cell pairs removed?
    C coeff; // incidence coefficient


    vector<typename CLIST::iterator> tokill;
    typename CLIST::iterator killhere;

    typename vector<typename CLIST::iterator>::iterator removeus;
    BT birth;

    // loop over complex
    for (bit = this->clist.begin(); bit != this->clist.end(); ++bit)
    {
        birth = bit->first;
        curf = bit->second;
        // loop backwards over frame, meaning top dimension first.
        for (fiter = curf->rbegin(); fiter != curf->rend(); ++fiter)
        {
            dim = fiter->first;
            if (dim == (this->getBotDim(birth)).first) break; // no faces of bottom dimension cells

            // loop over list of cells of this dim
            curlist = fiter->second;
            //cout<<"\n dim: "<<dim; cin.get();

            for (liter = curlist->begin(); liter != curlist->end(); ++liter)
            {
                curcell = *liter; // extract cell.
                if (curcell == NULL) continue; // ignore nulls
                if (!curcell->isIn) continue; // and removeds

                // if boundary has only one face...
                if (curcell->getBD().size()==1)
                {
                    // check unit incidence:
                    coeff = curcell->getBD().begin()->second;
                    if (isIntUnit(coeff))
                    {
                        freeface = curcell->getBD().begin()->first;

                        if (freeface->birth != curcell->birth) continue;
                        //cout<<"\n     pair: "<<*curcell<<" and "<<*freeface;
                        toret++;
                        removeCell(curcell,killhere,false);
                        tokill.push_back(killhere);
                        removeCell(freeface,killhere,false);
                        tokill.push_back(killhere);
                    }
                }
            }
        }
    }

    //cout<<"\n cleaning! "<<toret; cin.get();

    // now clean out list:
    for (removeus = tokill.begin(); removeus != tokill.end(); ++removeus)
    {
        delete **removeus;
        **removeus = NULL;
    }
    tokill.clear();

    return toret;
}

// enqueues all cells in a given chain to the queue structure:
template <typename C, typename BT>
void MComplex<C,BT>::enqueue(const CCHAIN& mychain, deque<Cell<C,BT>*>& Queue)
{
    // want to enqueue in general, BUT: if hyperq is on AND the Queue is larger
    // than some fraction of the entire complex, maybe we should just bail?
        CSTRUCT myst = mychain.cells;
        typename CSTRUCT::const_iterator chit;
        for (chit = myst.begin(); chit != myst.end(); ++chit)
        {
            if (!chit->first->isQueued)
            {
                Queue.push_back(chit->first);
                chit->first->isQueued = true;
                if (CORETALK) debout<<"\n      Enq: "<<*(chit->first);
            }
        }
}


// relays gradient path info from cell picked as queen to all remaining coboundaries.
// additionally enques all coboundary cells with isQueued = false to myQ!
template <typename C, typename BT>
void MComplex<C,BT>::updateGpathNearQueen_Cored(Cell<C,BT>* queen, deque<Cell<C,BT>*>& Queue, bool rightdim)
{
    rightdim = true;
    enqueue(queen->getCB(),Queue);

	// loop over coboundaries of queen....
	typename CSTRUCT::const_iterator cbs;
	for (cbs = queen->getCB().begin(); cbs != queen->getCB().end(); ++cbs)
	{
		// don't proceed unless the queen has the right dimension and a nontrivial gpath
		if (!rightdim || queen->gpath == NULL) continue;

		// inform coboundaries of the queen's gpath, allocating memory if needed
		if (cbs->first->gpath == NULL)
		{
			cbs->first->gpath = new CCHAIN;
			if (GPMEMTALK) debout<<"\n   Assigning new critbd of "<<*(cbs->first);
 		}
		*(cbs->first->gpath) += queen->gpath->scale(cbs->second);
		if (GPTALK) debout<<"\n   Updating gpath of "<<*(cbs->first)<<" to "<<*(cbs->first->gpath);
	}
}





// relays gradient path info from cell picked as queen to all remaining coboundaries.
// additionally enques all coboundary cells with isQueued = false to myQ!
template <typename C, typename BT>
void MComplex<C,BT>::updateGpathNearKing_Red(Cell<C,BT>* king, deque<Cell<C,BT>*>& Queue, bool rightdim)
{
    rightdim = true;
    // enqueue boundary of king:
    enqueue(king->getBD(),Queue);

	// loop over boundaries of king....
	typename CSTRUCT::const_iterator bds;
	for (bds = king->getBD().begin(); bds != king->getBD().end(); ++bds)
	{
		// don't proceed unless the king has the right dimension and a nontrivial gpath
		if (!rightdim || king->gpath == NULL) continue;

		// inform coboundaries of the queen's gpath, allocating memory if needed
		if (bds->first->gpath == NULL)
		{
			bds->first->gpath = new CCHAIN;
			if (GPMEMTALK) debout<<"\n   Assigning new critcb of "<<*(bds->first);
 		}
		*(bds->first->gpath) += king->gpath->scale(bds->second);
		if (GPTALK) debout<<"\n   Updating gpath of "<<*(bds->first)<<" to "<<*(bds->first->gpath);
	}
}

// relays king chain info from cell picked as queen to all remaining co-boundaries.
// !!! ONLY CALL THIS FUNCTION WHEN THE QUEEN HAS THE SAME DIM AS THE UNDERLYING ACE!!!
template <typename C, typename BT>
void MComplex<C,BT>::updateKgenNearQueen_Cored(Cell<C,BT>* king, Cell<C,BT>* queen, const C& coeff)
{
	// loop over coboundaries of queen....
	typename CSTRUCT::const_iterator cbs;
	for (cbs = queen->getCB().begin(); cbs != queen->getCB().end(); ++cbs)
	{
		// inform coboundaries of the queen: wake up, gens are coming
		if (cbs->first->kgen == NULL)
		{
			cbs->first->kgen = new CCHAIN;
			if (KGTALK) debout<<"\n   Assigning new kgen of "<<*(cbs->first);
 		}

		// if there is a descending king chain already,...
		if (king->kgen != NULL)
		{
			*(cbs->first->kgen) += *(king->kgen); // pass it upwards from the queen
		}
		cbs->first->kgen->addLink(king,-(cbs->second/coeff)); // and add this king to it.

		if (KGTALK) debout<<"\n   Updating kgen of "<<*(cbs->first)<<" to "<<*(cbs->first->kgen);
	}
	// remove king's kgen? not here, because we only call this function when queen has the right
	// dimension. the removal of the king's kgen should happen in all cases... in markPair();
}

// relays king chain info from cell picked as king to all remaining boundaries.
// !!! ONLY CALL THIS FUNCTION WHEN THE KING HAS THE SAME DIM AS THE UNDERLYING ACE!!!
template <typename C, typename BT>
void MComplex<C,BT>::updateKgenNearKing_Red(Cell<C,BT>* queen, Cell<C,BT>* king, const C& coeff)
{
	// loop over boundaries of king...
	typename CSTRUCT::const_iterator bds;
	for (bds = king->getBD().begin(); bds != king->getBD().end(); ++bds)
	{
		// inform boundaries of the king
		if (bds->first->kgen == NULL)
		{
			bds->first->kgen = new CCHAIN;
			if (KGTALK) debout<<"\n   Assigning new kgen of "<<*(bds->first);
 		}

		// if there is an ascending king chain already,...
		if (queen->kgen != NULL)
		{
			*(bds->first->kgen) += *(queen->kgen); // pass it onwnwards to other bdrys of the king
		}
		bds->first->kgen->addLink(king, -(bds->second/coeff)); // and add this king to it.

		if (KGTALK) debout<<"\n   Updating kgen of "<<*(bds->first)<<" to "<<*(bds->first->kgen);
	}
	// remove queen's kgen? not here, because we only call this function when queen has the right
	// dimension. the removal of the queen's kgen should happen in all cases... in markPair();
}



// inherit queen's gradient path from the king's
template <typename C, typename BT>
bool MComplex<C,BT>::inheritGpathFromKing_Cored(Cell<C,BT>* queen, Cell<C,BT>* king, bool rightdim)
{
    rightdim = true;
	// if there is no gradient path to inherit or we are in the wrong dimension,...
	if (king->gpath == NULL)
	{
		// then null out queen's gradient path!
		if (queen->gpath != NULL)
		{
			queen->gpath->clear();
			delete queen->gpath;
			queen->gpath = NULL;
		}
		// and exit
		return false;
	}

	// if we're here, then there is a nontrivial gradient path to inherit
	// and we are in the right dimension!
	C kqcoe = king->getBD().begin()->second; // coefficient of KQ incidence

	// it is possible the queen' gradient path might already have critical info, but we
	// can ignore it! the pre-stored gradient path for a queen of dimQ will be a
	// chain of dimQ-1, which may be ignored since the paths that count are
	// pairs of type [dimQ, dimQ+1]...
	if (queen->gpath != NULL)
	{
		// ignore pre-existing gpath of queen
		queen->gpath->clear();
	}
	else // otherwise allocate new memory
	{
		if (GPMEMTALK) debout<<"\n   Assigning NEW gpath of "<<*queen;
		queen->gpath = new CCHAIN;
	}
	// inherit king's gp if non null...
	typename CSTRUCT::const_iterator critdn; // loops over critical boundary of king
	for (critdn = king->gpath->begin(); critdn != king->gpath->end(); ++critdn)
	{
		queen->gpath->addLink(critdn->first, -critdn->second / kqcoe); // and adds those to queen's gpath
	}
	if (GPTALK) debout<<"\n   Queen's Gpath is "<<*(queen->gpath);

	// now we safely remove king's gpath
	king->gpath->clear();
	delete king->gpath;
	king->gpath = NULL;

	if (GPMEMTALK) debout<<"\n          deleting critbd of "<<*king;
	return true;
}


// inherit king's gradient path from the queen's
template <typename C, typename BT>
bool MComplex<C,BT>::inheritGpathFromQueen_Red(Cell<C,BT>* king, Cell<C,BT>* queen, bool rightdim)
{
    rightdim = true;
	// if there is no gradient path to inherit or we are in the wrong dimension,...
	if (queen->gpath == NULL)
	{
		// then null out king's gradient path!
		if (king->gpath != NULL)
		{
			king->gpath->clear();
			delete king->gpath;
			king->gpath = NULL;
		}
		// and exit
		return false;
	}

	// if we're here, then there is a nontrivial gradient path to inherit
	// and we are in the right dimension!
	C qkcoe = queen->getCB().begin()->second; // coefficient of KQ incidence

	// it is possible the king' gradient path might already have critical info, but we
	// can ignore it! the pre-stored gradient path for a king of dim K will be a
	// chain of dimK - 1, which may be ignored since the paths that count are
	// pairs of type [dimK+1, dimK]...
	if (king->gpath != NULL)
	{

        if (king->gpath->size() != 0)
        {
            if (GPTALK)
            {
                debout<<"\n clearing gpath of king "<<*king<<" from "<<*(king->gpath);
                debout<<"\n to inherit queen's gpath, which is "<<*(queen->gpath);
            }
            king->gpath->clear();
        }
		// ignore pre-existing gpath of queen
    }
	else // otherwise allocate new memory
	{
		if (GPMEMTALK) debout<<"\n   Assigning NEW gpath of "<<*king;
		king->gpath = new CCHAIN;
	}
	// inherit queen's gp if non null...
	typename CSTRUCT::const_iterator critup; // loops over critical coboundary of queen
	for (critup = queen->gpath->begin(); critup != queen->gpath->end(); ++critup)
	{
		king->gpath->addLink(critup->first, -critup->second / qkcoe); // and adds those to king's gpath
	}
	if (GPTALK) debout<<"\n   King's Gpath is "<<*(king->gpath);

	// now we safely remove king's gpath
	queen->gpath->clear();
	delete queen->gpath;
	queen->gpath = NULL;

	if (GPMEMTALK) debout<<"\n          deleting critcb of "<<*queen;
	return true;
}


// excises queen-king (reduction) pairs... "queen" is the Queen, "king" is
// its only coboundary element, Q is the reduction queue, rootdim is the
// dimension of the root for this instance of reduce. returns true in case
// gpaths were manipulated
template <typename C, typename BT>
bool MComplex<C,BT>::markPair_Red(Cell<C,BT>* queen, Cell<C,BT>* king, deque<Cell<C,BT>*>& Queue, num rootdim,
                    typename CLIST::iterator& qpos, typename CLIST::iterator& kpos, bool makekgens=false)
{
	if (GPTALK || CORETALK || AKQTALK)
	{
		debout<<"\n  Making pair: "<<*queen<<" "<<*king;
		//cin.get();
	}

	// In every case, eliminate king's kgen information for it is useless!
	if (king->kgen != NULL)
	{
		king->kgen->clear();
		delete king->kgen;
		king->kgen = NULL;
	}

	bool rightdim = (king->getDim() == rootdim) ? true : false;

	// store queen's critcb information in king' gradient path if the dimension matches
	inheritGpathFromQueen_Red(king, queen, rightdim);

    CCHAIN qbd = queen->getBD();
    removeCell(queen, qpos, false); // excise queen as uncritical!

	C qkcoe = queen->getCB().begin()->second; // coefficient of QK incidence

	// Call to keep track of generators as king chains
	if (makekgens)
	{
		if(rightdim) updateKgenNearKing_Red(queen, king, qkcoe);
		// now remove queen's kgen
		if (queen->kgen != NULL)
		{
			queen->kgen->clear();
			delete queen->kgen;
			queen->kgen = NULL;
		}
	}

	// update gradient paths and enqueue remaining boundaries of King:
	//cout<<"\ntry qgp update"; cin.get();
	updateGpathNearKing_Red(king, Queue, rightdim);


	// now we safely remove the king's gradient path
	// MAY WANT TO STORE IT ANYWAY
	if (king->gpath != NULL)
	{
		delete king->gpath;
		king->gpath = NULL;
		if (GPMEMTALK) debout<<"\n          deleting gpath of "<<*king;
	}

    //enqueue bd of queen:
    enqueue (qbd, Queue);

	// now remove the king itself from the complex as uncritical
	//cout<<"\n kill king! "<<*queen; cin.get();
	removeCell(king, kpos, false);
	return true;
}



// excises king-queen (coreduction) pairs... "king" is the King, "queen" is
// its only boundary element, Q is the coreduction queue, rootdim is the
// dimension of the root for this instance of coreduce. returns true in case
// gpaths were manipulated
template <typename C, typename BT>
bool MComplex<C,BT>::markPair_Cored(Cell<C,BT>* king, Cell<C,BT>* queen, deque<Cell<C,BT>*>& Queue, num rootdim,
                    typename CLIST::iterator& kpos, typename CLIST::iterator& qpos, bool makekgens=false)
{
	if (GPTALK || CORETALK || AKQTALK)
	{
		debout<<"\n  Making pair: "<<*king<<" "<<*queen;
		//cin.get();
	}

	// In every case, eliminate queens' kgen information for it is useless!
	if (queen->kgen != NULL)
	{
		queen->kgen->clear();
		delete queen->kgen;
		queen->kgen = NULL;
	}

	bool rightdim = (queen->getDim() == rootdim) ? true : false;

	// store king's critbd information in queens' gradient path if the dimension matches
	inheritGpathFromKing_Cored(queen, king, rightdim);
	//cout<<"\n kill king! "<<*king; cin.get();
	// excise the king as uncritical

//EXPERIMENT 1
   CCHAIN kcb = king->getCB();

	removeCell(king, kpos, false);
	C kqcoe = king->getBD().begin()->second; // coefficient of KQ incidence

	// Call to keep track of generators as king chains
	if (makekgens)
	{
		if(rightdim) updateKgenNearQueen_Cored(king, queen, kqcoe);
		// now remove king's kgen
		if (king->kgen != NULL)
		{
			king->kgen->clear();
			delete king->kgen;
			king->kgen = NULL;
		}
	}

	// update gradient paths and enqueue remaining coboundaries of Queen:
	//cout<<"try qgp update"; cin.get();
	updateGpathNearQueen_Cored(queen, Queue, rightdim);

	// now we safely remove the queen's gradient path
	// MAY WANT TO STORE IT ANYWAY
	if (queen->gpath != NULL)
	{
		delete queen->gpath;
		queen->gpath = NULL;
		if (GPMEMTALK) debout<<"\n          deleting gpath of "<<*queen;
	}

    // enqueue coboundary of king
    enqueue(kcb,Queue);


	// now remove the queen itself from the complex as uncritical
	//cout<<"\n kill queen! "<<*queen; cin.get();
	removeCell(queen, qpos, false);
	return true;
}



// removes a free face pair from the complex, returns an iterator to the
// next cell of lower dimension
template <typename C, typename BT>
void MComplex<C,BT>::removeFreeFace(Cell<C,BT>* torem, bool deep=true)
{
	//typename COMPLEX::iterator toret;
	Cell<C,BT>* coface = torem->cb.begin()->first;

	if (FFTALK)
	{
	    cout<<"\n\n Removing FFC: "<<*torem<<" "<<*coface;
	    cout<<" coface has boundary: "<<coface->getBD();
	}


	typename CLIST::iterator cellpos, coffpos;
	removeCell(torem, cellpos, false); // remove as uncritical
	removeCell(coface, coffpos, false); // and remove coface also as uncritical

	if (FFTALK)
	{
	    cout<<"\n\n Removed!!";
	    cout<<" now, coface has boundary: "<<coface->getBD();
	    cin.get();
	}



	// certainly, the gradient paths are useless independent of whether we want
	// a deep erase or not.
	if ((*cellpos)->gpath != NULL)
	{
		delete (*cellpos)->gpath;
		(*cellpos)->gpath = NULL;
	}
	if ((*coffpos)->gpath != NULL)
	{
		delete (*coffpos)->gpath;
		(*coffpos)->gpath = NULL;
	}

	if (deep) // erase cells altogether from memory if deep, set to null!
	{
		delete *cellpos; *cellpos = NULL;
		delete *coffpos; *coffpos = NULL;
	}
}

// removes a free coface pair from the complex, returns an iterator to the
// next cell of lower dimension
template <typename C, typename BT>
void MComplex<C,BT>::removeFreeCoface(Cell<C,BT>* torem, bool deep=true)
{
	//typename COMPLEX::iterator toret;
	Cell<C,BT>* face = torem->bd.begin()->first;
	typename CLIST::iterator cellpos, coffpos;
	removeCell(torem, cellpos, false); // remove torem as uncritical
	removeCell(face, coffpos, false); // and remove its face also as uncritical
	// certainly, the gradient paths are useless independent of whether we want
	// a deep erase or not.
	if ((*cellpos)->gpath != NULL)
	{
		delete (*cellpos)->gpath;
		(*cellpos)->gpath = NULL;
	}
	if ((*coffpos)->gpath != NULL)
	{
		delete (*coffpos)->gpath;
		(*coffpos)->gpath = NULL;
	}

	if (deep) // erase cells altogether from memory if deep, and set to null!
	{
		delete *cellpos; *cellpos = NULL;
		delete *coffpos; *coffpos = NULL;
	}
}

// finds next critical cell after performing free coface collapses
template <typename C, typename BT>
num MComplex<C,BT>::getNextCrit_Red (const BT& birth, Cell<C,BT>*& newace, bool deep = true)
{
 	num remsize = 0; // number of free face pairs removed!
	newace = NULL; // initialize to null.

	Cell<C,BT>* face = NULL;

	// extract dimension with next critical cell
	pair<num,bool> mdret = this->getMaxDim(birth);
	// check if there are no more cells left in this dimension
	if (mdret.second == false) return 0;

	// okay, we want to find the next critical cell for dimension mdret.first
	// at frame corresponding to birth... so,

	typename COMPLEX::const_iterator findb = this->clist.find(birth);
	FRAME* curf = findb->second;

	typename FRAME::const_iterator findd;
	CLIST* curl;
	typename CLIST::const_iterator lit;
	num finalind;
	//Cell<C,BT>* cof = NULL; //potential coface for free face collapsing!

	do // while there are uncritical cells...
	{
		mdret = getMaxDim(birth); // EXISTING MIN UNCRIT DIMENSION!!
		newace = NULL; // reset new ace to NULL...

		if (mdret.second == false) // bad mindim, so bad birth!
		{
			return 0; // exit at NULL newace
		}

		// okay, we now have a valid mindim with last cell uncritical...
		// this cell is at finalind, to be computed

		findd = curf->find(mdret.first); // locate critical dimension
		curl = findd->second; // the list containing the critical cell...
		// get the index of the highest-index unremoved cell in this list...
		finalind = (this->getSizeInfo(birth,mdret.first)).second;

		if (minindex > finalind)
		{
			debout << "min fin violation at birth = "<<birth<<" mindim = "<<mdret.first;
			debout << " m = "<<minindex<<" but f = "<<finalind;
			minindex = 0;
		}

		// now traverse list between minindex and finalind to get next crit cell
		for (lit = curl->begin() + minindex; lit <= curl->begin() + finalind; ++lit)
		{
			newace = *lit;
			if (newace == NULL) continue; // ignore nulls
			if (!(newace->isIn)) continue; // and removeds
			if (newace->isCrit) // flag critical cell order error
			{
				//cout<<"\n    disorderly crit cell: "<<newace->ind<<" but minindex = "<<minindex ;
				//cin.get();
				continue;
			}
			// if we reach here, the cell is _not_ critical and may be checked for ffc!
			if (isFreeCoface(newace, face))
			{
				remsize++; // one more pair removed!
				removeFreeCoface(newace,deep);

				// want to check if there are any more uncrit cells left at this dimension!
				// if not, just break out of the FOR loop and reevaluate min dimension... this
				// happens in the WHILE() at the end of the DO loop!
				if (!this->hasUncritCells(birth,mdret.first))
				{
				    // reset minindex since we are changing dimension
					minindex = 0;
					break;
				}
				// otherwise, just continue with the list...
				continue;
			}
			// if we get here, we have a critical cell candidate stored in newace!
			minindex = newace->ind;
			return remsize;
		}
		// if we get here, the min dim has changed. re-evaluate it!
	} while (hasUncritCells(birth,getMaxDim(birth).first));

	// if we get here, there is no hope of new ace
	newace = NULL;
	return remsize;
}

// finds next critical cell after performing free face collapses
template <typename C, typename BT>
num MComplex<C,BT>::getNextCrit_Cored (const BT& birth, Cell<C,BT>*& newace, bool deep = true)
{
	num remsize = 0; // number of free face pairs removed!
	newace = NULL; // initialize to null.

	Cell<C,BT>* coface = NULL;

	// extract dimension with next critical cell
	pair<num,bool> mdret = this->getMinDim(birth);
	// check if there are no more cells left in this dimension
	if (mdret.second == false) return 0;

	// okay, we want to find the next critical cell for dimension mdret.first
	// at frame corresponding to birth... so,

	typename COMPLEX::const_iterator findb = this->clist.find(birth);
	FRAME* curf = findb->second;

	typename FRAME::const_iterator findd;
	CLIST* curl;
	typename CLIST::const_iterator lit;
	num finalind;

	do // while there are uncritical cells...
	{
		mdret = getMinDim(birth); // EXISTING MIN UNCRIT DIMENSION!!
		newace = NULL; // reset new ace to NULL...

		if (mdret.second == false) // bad mindim, so bad birth!
		{
			return 0; // exit at NULL newace
		}

		// okay, we now have a valid mindim with last cell uncritical...
		// this cell is at finalind, to be computed

		findd = curf->find(mdret.first); // locate critical dimension
		curl = findd->second; // the list containing the critical cell...
		// get the index of the highest-index unremoved cell in this list...
		finalind = (this->getSizeInfo(birth,mdret.first)).second;

		if (minindex > finalind)
		{
			debout << "min fin violation at birth = "<<birth<<" mindim = "<<mdret.first;
			debout << " m = "<<minindex<<" but f = "<<finalind;
			minindex = 0;
		}

		// now traverse list between minindex and finalind to get next crit cell
		for (lit = curl->begin() + minindex; lit <= curl->begin()+ finalind; ++lit)
		{
			newace = *lit;
			if (newace == NULL) continue; // ignore nulls
			if (!(newace->isIn)) continue; // and removeds
			if (newace->isCrit) // flag critical cell order error
			{
				//cout<<"\n    disorderly crit cell: "<<*curcell<<" but minindex = "<<minindex ;
				//cin.get();
				continue;
			}
			// if we reach here, the cell is _not_ critical and may be checked for ffc!
			if (isFreeFace(newace, coface))
			{
				//cout<<"\n about to collapse::: "<<*newace;
				remsize++; // one more pair removed!
				removeFreeFace(newace,deep);

				// want to check if there are any more uncrit cells left at this dimension!
				// if not, just break out of the FOR loop and reevaluate min dimension... this
				// happens in the WHILE() at the end of the DO loop!
				if (!this->hasUncritCells(birth,mdret.first))
				{
					// reset minindex since we are changing dimension
					minindex = 0;
					break;
				}
				// otherwise, just continue with the list...
				continue;
			}
			// if we get here, we have a critical cell candidate stored in newace!
			minindex = newace->ind;
			return remsize;
		}
		// if we get here, the min dim has changed. re-evaluate it!
	} while (hasUncritCells(birth,getMinDim(birth).first));

	// if we get here, there is no hope of new ace
	newace = NULL;
	return remsize;
}


// reduces starting from first highest dimensional available uncritical cell
// returns number of QK pairs excised!

// this version returns a pair, the first element of which is the number of
// cell pairs excised (as free faces or cofaces) and the second element is
// a flag indicating:
// -1: the entire frame is excised or nonexistent.
// 1: no more collapses found to current root, but more could be possible
// -2: no more uncritical cells remain in current frame.

template <typename C, typename BT>
pair<num,int> MComplex<C,BT>::Reduce(const BT& birth = INITBT, const BT& ignore = BANBT,
                                     bool deep = true, bool tracer = false)
{
    isred = true;
	//if (birth == 1) tracer = true;

	if(this->fsize(birth) != 0)
	{
		if (CORETALK || CMFTALK || tracer)
		{
			debout<<"\n         REDUCTION ["<<birth<<"], with ";
			debout<<" bot "<<this->getBotDim(birth).first;
			debout<<" and max "<<getMaxDim(birth).first;

            //this->printFrame(1);
            cin.get();
		}
	}
	else return make_pair(0,-1); // frame does not exist

	num numrem = 0; //ffres.first; // number removed via free face collapse

	Cell<C,BT>* ace; // potential candidate to reduce
	numrem += getNextCrit_Red(birth,ace,deep);
	// we now select an ace of maximal dimension that does
	// not allow a simple free coface collapse to a lower
	// dimensional coboundary cell

	if (ace == NULL) return make_pair(numrem,-1);
	vector<typename CLIST::iterator> toKill;

	// now the routine

	deque<Cell<C,BT>*> Queue; // create empty queue of cell pointers
	Queue.clear();
	Queue.push_back(ace); // add root to queue!

	if (CORETALK || AKQTALK || tracer)
	{
		debout<<"\n trying ACE at index: "<<ace->ind;
		cin.get();
		debout<<*ace<<" bd size: "<<(ace->getBD()).size();
	}

	vector<Cell<C,BT>*> flushMe; // this list holds all uncritical cells that are
	// enqueued but NOT reduced... we will "flush" them by setting their queued
	// status to false once coreduction completes, so they may be enqueued in the
	// next coreduction cycle!

	// make root critical, and remove it from boundaries of higher cells
	//cout<<"critting!"; cin.get();
	makeCritical_Red(ace, !deep);
	//cout<<"\n+++++ Cored with ace: "<<*ace;

	// king is the next popped element from the queue, queen is its
	// potential coreduction pair partner of lower dimension
	Cell<C,BT>* king; Cell<C,BT>* queen;

	while(!Queue.empty())
	{
		//if (CORETALK) printQ(Queue);
		queen = Queue.at(0); // pop first element: queen
		Queue.pop_front();
		if (queen->birth == ignore) continue; // do not reduce from ignored birthtime!!
		if (CORETALK || tracer)
		{
			debout<<"\n Popping... "<<*queen;
			debout<<" with cb size: "<<(queen->cb).size()
                  <<" and bd size: "<<(queen->bd).size();
			cin.get();
		}

		if (queen->getCB().size() == 1) // if queen has nontrivial coboundary,...
		{
			king = queen->cb.begin()->first; // first cell in boundary
			// check for kq pair
			if (isQKPair(queen, king))
			{
				if (AKQTALK || tracer)
				{
					debout<<"\n  King: "<<*king<<" at "<<king;
					debout<<", and Queen: "<<*queen<<" at "<<queen;
					cin.get();
				}
				//cout<<"\nTrying to mark pair!"; cin.get();
				typename CLIST::iterator kpos, qpos;
				markPair_Red(queen, king, Queue, ace->getDim(),qpos, kpos, !deep); // make king queen pair. This requires us to know of the Queue and
				//cout<<"\nmarked!!!"; cin.get();

				numrem++; 	// increment count of cells removed!
				// add option to add king to generator list or something
				// if deep, delete cells pointed to by king and queen
				if (deep)
				{
					toKill.push_back(qpos);
					toKill.push_back(kpos);
				}
			}
			// we don't have unit incidence... what to do?
		} // end if queen->cb has size 1
		else if (this->hyperq || queen->getCB().size()==0)
		{
			if (CORETALK || tracer)
			{
				debout<<", No king; Enq bd...";
			}
			//enqueue the boundary of the current cell:
			enqueue(queen->getBD(),Queue);
			// if cell is critical, remove links to its boundary
			if (queen->isCrit) queen->bd.clear();
			flushMe.push_back(queen); // otherwise prepare to flush Q status
		} // else if king->bd empty
		//else flushMe.push_back(queen);
	}  // while Q nonempty
	if (MTALKBIG || CORETALK || tracer) debout<<"\n   red-remove #: "<<numrem;
	// flush the isqueued statuses here:
	flushQ(flushMe);

	// delete memory associated to QK pairs!
	typename vector<typename CLIST::iterator>::iterator killme;
	for (killme = toKill.begin(); killme != toKill.end(); ++killme)
	{
		delete **killme;
		**killme = NULL;
	}
	return make_pair(numrem,1); // number of collapses performed
} // coreduce

// coreduces starting from first lowest dimensional available uncritical cell
// returns number of KQ pairs excised! Originally by Mrozek!

// this version returns a pair, the first element of which is the number of
// cell pairs excised (as free faces or cofaces) and the second element is
// a flag indicating:
// -1: the entire frame is excised or nonexistent.
// 1: no more collapses found to current root, but more could be possible
// -2: no more uncritical cells remain in current frame.

template <typename C, typename BT>
pair<num,int> MComplex<C,BT>::CoReduce(const BT& birth = INITBT, const BT& ignore = BANBT,
                                       bool deep = true, bool tracer = false)
{
    isred = false;
	//typename COMPLEX::iterator cit = this->clist.find(birth);
	//cout<<"\n++++++++++NEW COREDUCE++++++++++++++++";

	if(this->fsize(birth) != 0)
	{
		if (CORETALK || CMFTALK || tracer)
		{
			debout<<"\n         COREDUCTION ["<<birth<<"]\n";
			debout<<" top "<<this->getTopDim(birth).first;
			debout<<" vs min "<<getMinDim(birth).first; cin.get();
		}
	}
	else return make_pair(0,-1); // frame does not exist

	num numrem = 0; //ffres.first; // number removed via free face collapse

	Cell<C,BT>* ace; // potential candidate to reduce
	numrem += getNextCrit_Cored(birth,ace,deep);
	// we now select an ace of minimal dimension that does
	// not allow a simple free face collapse from a higher
	// dimensional coboundary cell

	if (ace == NULL) return make_pair(numrem,-1);
	vector<typename CLIST::iterator> toKill;

	// now the routine

	deque<Cell<C,BT>*> Queue; // create empty queue of cell pointers
	Queue.clear();
	Queue.push_back(ace); // add root to queue!

	if (CORETALK || AKQTALK || tracer)
	{
		debout<<"\n trying ACE: ";
		cin.get();
		debout<<*ace<<" with cb size: "<<(ace->getCB().size());
	}

	vector<Cell<C,BT>*> flushMe; // this list holds all uncritical cells that are
	// enqueued but NOT reduced... we will "flush" them by setting their queued
	// status to false once coreduction completes, so they may be enqueued in the
	// next coreduction cycle!

	// make root critical, and remove it from boundaries of higher cells
	//cout<<"critting!"; cin.get();
	makeCritical_Cored(ace, !deep);
	//cout<<"\n+++++ Cored with ace: "<<*ace;

	// king is the next popped element from the queue, queen is its
	// potential coreduction pair partner of lower dimension
	Cell<C,BT>* king; Cell<C,BT>* queen;

	while(!Queue.empty())
	{
		//if (CORETALK) printQ(Queue);
		king = Queue.at(0); // pop first element: king
		Queue.pop_front();
		if (king->birth == ignore) continue; // do not reduce from ignored birthtime!!
		if (CORETALK || tracer)
		{
			debout<<"\n Popping... "<<*king;
			debout<<" with bd size: "<<(king->bd).size()
                  <<" and cb size: "<<(king->cb).size();
			cin.get();
		}

		if (king->getBD().size() == 1) // if king has nontrivial boundary,...
		{
			queen = king->bd.begin()->first; // first cell in boundary
			// check for kq pair
			if (isKQPair(king, queen))
			{
				if (AKQTALK || tracer)
				{
					debout<<"\n  King: "<<*king<<" at "<<king;
					debout<<", and Queen: "<<*queen<<" at "<<queen;
					cin.get();
				}
				//cout<<"\nTrying to mark pair!"; cin.get();
				typename CLIST::iterator kpos, qpos;
				markPair_Cored(king, queen, Queue, ace->getDim(),kpos, qpos, !deep); // make king queen pair. This requires us to know of the Queue and
				//cout<<"\nmarked!!!"; cin.get();

				numrem++; 	// increment count of cells removed!
				// add option to add king to generator list or something
				// if deep, delete cells pointed to by king and queen
				if (deep)
				{
					toKill.push_back(kpos);
					toKill.push_back(qpos);
				}
			}
			// we don't have unit incidence... what to do?
		} // end if king->bd has size 1
		else if (this->hyperq || king->getBD().size()==0)
		{
			if (CORETALK || tracer)
			{
				debout<<", No queen; Enq cb...";
			}
			//enqueue the coboundary of the queue element
			enqueue (king->getCB(),Queue);
			//king->isQueued = false;
			// if cell is critical, clear out its co-boundary
			if (king->isCrit) king->cb.clear();
			flushMe.push_back(king); // otherwise prepare to flush Q status
		} // else if king->bd empty
	}  // while Q nonempty
	if (MTALKBIG || CORETALK || tracer) debout<<"\n   cored-remove #: "<<numrem;
	// flush the isqueued statuses here:
	flushQ(flushMe);

	// delete memory associated to KQ pairs!
	typename vector<typename CLIST::iterator>::iterator killme;
	for (killme = toKill.begin(); killme != toKill.end(); ++killme)
	{
		delete **killme;
		**killme = NULL;
	}
	return make_pair(numrem,1); // number of collapses performed
} // coreduce


// makes morse complex using coreduce iteratively
// statuses: (the second integer returned value):
template <typename C, typename BT>
pair<num,bool> MComplex<C,BT>::RedMorsify(const BT& birth, bool deep = true,
                                          BT except = BANBT,bool tracer=false)
{
   //if (birth == 5) tracer = true;


    minindex = 0;
	typename COMPLEX::iterator cit;
	bool ct;
	if (CMFTALK || TIMEINFO || tracer)
	{
		debout<<"\n      redmorsify: "<<birth;
		//cin.get();
	}

	pair<num,int> coret; // stores return value of coreduce
	pair<num,bool> oldmdret, newmdret; // stores result of critical dimension computations
	int numrem = 0;
	num numiter = 0;
	do
	{
		// old critical dimension
		ct = false;
		oldmdret = this->getMaxDim(birth); // get max crit dimension BEFORE reduce is called

		coret = Reduce(birth, except, deep, tracer); // if deep, then we erase cell pairs from memory
        numrem += coret.first;

        // set traces here:
        //tracer = (++numiter % 5000 == 0) ? true : false;

		if (CMFTALK || tracer)
		{
			debout<<"\n        reduce returned: "<<coret.first<<" "<<coret.second;
			//this->printSizeInfo(true);
			//cin.get();
		}
        //tracer = false;

		switch(coret.second)
		{
		case -1: // frame does not exist...
			//cout<<"Frame Wiped"; cin.get();
			if (CMFTALK || tracer)
			{
				debout<<"\n               red wiped frame for birth time "<<birth;
				// reset min dimension for next frame, time to move on!
			}
			minindex = 0;
            //cin.get();
			return make_pair(numrem, false);
		case 1: // frame exists, we can coreduce some more
			// increment number of removed cells
			// extract minimum uncrit dimension...
			newmdret = getMaxDim(birth); // get max dimension again, AFTER coreduce
			if (newmdret.second == false)
			{
				// no more reductions, everything left is critical
				return make_pair(numrem,false);
			}
			// have we changed max dimension? if so, reset minindex to zero
			if (newmdret.first != oldmdret.first)
			{
				//cout<<"\n old minindex was "<<minindex;
				//cout<<"\n different mindim: "<<newmdret.first<<", not "<<oldmdret.first;
            minindex = 0;
			}

			if (CMFTALK || tracer)
			{
				debout<<" and removed "<<coret.first<<" pairs";
			}
			if (!hasUncritCells(birth))
			{
				if (CMFTALK || tracer)
				{
					debout<<"\n       returning without iteration, there are no uncrits left at this birth!! ";
					//cin.get();
				}
				if (TIMEINFO)
				{
					debout<<"\n****/redmorsify";
				}
				return make_pair(numrem,true);
			}
			else if (CMFTALK || tracer)
			{
				debout<<"\n          more uncrit cells remain, carrying on coreduce! \n";
				//cin.get();
			}
			ct = true;
			break;
		}
		if (ct) { /*cout<<"\n ctn"; cin.get();*/ continue;}
		else
		{
			break;
		}
	} while(1); // redundant condition!

	//cout<<"\n>>>>>>>> SHOULD NOT HAPPEN ";
	return make_pair(numrem, true);
}

// makes morse complex using coreduce iteratively
// statuses: (the second integer returned value):
template <typename C, typename BT>
pair<num,bool> MComplex<C,BT>::CoredMorsify(const BT& birth, bool deep = true,
                                            BT except = BANBT, bool tracer = false)
{
    minindex = 0;
	typename COMPLEX::iterator cit;
	bool ct;
	if (CMFTALK || TIMEINFO || tracer)
	{
		debout<<"\n****coredmorsify: "<<birth;
		cin.get();
	}

	pair<num,int> coret; // stores return value of coreduce
	pair<num,bool> oldmdret, newmdret; // stores result of critical dimension computations
	int numrem = 0;
	num numiter = 0;
	do
	{
		// old critical dimension
		ct = false;
		oldmdret = this->getMinDim(birth); // get min crit dimension BEFORE coreduce is called

        coret = CoReduce(birth, except, deep, tracer); // if deep, then we erase cell pairs from memory
        numrem += coret.first;

        // set traces here:
       //tracer = (++numiter % 5000 == 0) ? true : false;

		if (CMFTALK || tracer)
		{
			debout<<"\n        coreduce returned: "<<coret.first<<" "<<coret.second<<" in iteration "<<numiter;
			//this->printSizeStruct();
			//f;
			//cin.get();
		}
		//tracer = false;

		switch(coret.second)
		{
		case -1: // frame does not exist...
			if (CMFTALK || tracer)
			{
				debout<<"\n               cored wiped frame for birth time "<<birth;
				// reset min dimension for next frame, time to move on!
				minindex = 0;
				//cin.get();
			}
			return make_pair(numrem, false);
		case 1: // frame exists, we can coreduce some more
			newmdret = getMinDim(birth); // get min dimension again, AFTER coreduce

			if (newmdret.second == false)
			{
				// no more reductions, everything left is critical
				return make_pair(numrem,false);
			}

			// have we changed min dimension? if so, reset minindex to zero
			if (newmdret.first != oldmdret.first)
			{
				//cout<<"\n old minindex was "<<minindex;
				//cout<<"\n different mindim: "<<newmdret.first<<", not "<<oldmdret.first;
				minindex = 0;
			}

			if (CMFTALK || tracer)
			{
				debout<<" and removed "<<coret.first<<" pairs";
			}
			if (!hasUncritCells(birth))
			{
				if (CMFTALK || tracer)
				{
					debout<<"\n       returning without iteration, there are no uncrits left at this birth!! ";
					//cin.get();
				}
				if (TIMEINFO)
				{
					debout<<"\n****/coredmorsify";
				}
				return make_pair(numrem,true);
			}
			else if (CMFTALK)
			{
				debout<<"\n          more uncrit cells remain, carrying on coreduce! \n";
				//cin.get();
			}
			ct = true;
			break;
		}
		if (ct) { /*cout<<"\n ctn"; cin.get();*/ continue;}
		else
		{
			break;
		}
	} while(1); // redundant condition!

	//cout<<"\n>>>>>>>> SHOULD NOT HAPPEN ";
	return make_pair(numrem, true);
}

// repeatedly performs morse coreduction until no new cells are removed!
template <typename C, typename BT>
num MComplex<C,BT>::MorseCoreduce(bool savegens = false, bool tracer = false)
{
	if (TIMEINFO || tracer)
	{
		debout<<"\n   MorseCoreduce";
        cin.get();
	}

	typename COMPLEX::const_iterator cit;
	pair<num,bool> numrem;
	//vector<BT> origbirths;
	BT birth;
	num toret = 0;
	num oldsz;

	// store all encountered birth times in origbirths...
	for(cit = this->clist.begin(); cit != this->clist.end(); ++cit)
	{
		birth = cit->first;
		if (this->fsize(birth) == 0) continue; // if you have killed the frame...
		// new birth, reset minindex!
		minindex = 0;
		if (BIGMORSTALK || tracer)
		{
			debout<<"\n\n*****BIRTH: "<<birth;
			//cin.get();
		}
		oldsz = this->size();

		numrem = CoredMorsify(birth, !savegens, BANBT, tracer); // deep removal if generators are not needed
        markUncritical_Cored(birth); // if frame has crit cells left, update!


		// if half the number of cells have been removed and we don't want to save generators...
		// recopy the lists to eliminate nulls
		if (!savegens && 2*this->size() <= oldsz)
		{
			memorySaver(birth);
		}

		toret += numrem.first;

		if (BIGMORSTALK || tracer)
		{
			debout<<"\n***Cmfy returned: "<<numrem.first<<" "<<numrem.second;
			debout<<", pairs removed so far: "<<toret<<"*******\n"; //cin.get();
			//cout<<*this;
		}
	}
	if(TIMEINFO || tracer)
	{
		debout<<"*/morsecoreduce";
		//cin.get();
	}
	return toret;
}

// repeatedly performs morse reduction until no new cells are removed!
template <typename C, typename BT>
num MComplex<C,BT>::MorseReduce(bool savegens = false, bool tracer = false)
{
	if (TIMEINFO || BIGMORSTALK || tracer)
	{
		debout<<"\n   MorseReduce";
		cin.get();
		//debout<<*this;
	}

	typename COMPLEX::const_reverse_iterator cit;
	pair<num,bool> numrem;
	//vector<BT> origbirths;
	BT birth;
	num toret = 0;
	num oldsz;

	// store all encountered birth times in origbirths...
	for(cit = this->clist.rbegin(); cit != this->clist.rend(); ++cit)
	{
		birth = cit->first;
		if (this->fsize(birth) == 0) continue; // if you have killed the frame...
		// new birth, reset minindex!
		minindex = 0;
		if (BIGMORSTALK || tracer)
		{
			debout<<"\n\n*****BIRTH: "<<birth;
			//cin.get();
		}
		oldsz = this->size();
		numrem = RedMorsify(birth, !savegens, BANBT,tracer); // deep removal if generators are not needed!
        markUncritical_Red(birth); // if frame has crit cells left, update!


		// if half the number of cells have been removed and we don't want to save generators...
		// recopy the lists to eliminate nulls
		if (!savegens && 2*this->size() <= oldsz)
		{
		    if (BIGMORSTALK || tracer) debout<<"\n mem-saving for birth "<<birth;
			memorySaver(birth);
		}

		toret += numrem.first;

		if (BIGMORSTALK || tracer)
		{
			debout<<"\n***Rmfy returned: "<<numrem.first<<" "<<numrem.second;
			debout<<", pairs removed so far: "<<toret<<"*******\n"; //cin.get();
			//debout<<*this;
		}
	}
	if(TIMEINFO || tracer)
	{
		debout<<"*/morsereduce";
		//cin.get();
	}
	return toret;
}


template <typename C, typename BT>
num MComplex<C,BT>::TopAlternator(bool redfirst=true, double thresh = 0.2,
                                  bool savegens = false, bool tracer = false)
{
    if (tracer)
    {
        cout << "\nTOP-ALTERNATOR";
        cin.get();
    }
    num toret = 0;
    if (this->size() == 0) // nothing to do...?
    {
        thresh --;
        return 0;
    }

    if (redfirst)
    {
       toret += MorseReduce(savegens,tracer);
       //toret += MorseCoreduce(savegens);
    }
    else
    {
        toret += MorseCoreduce(savegens,tracer);
        //toret += MorseReduce(savegens,thresh);
    }
    //cin.get();
    return toret;
}


// alternates coreduction and reduction
template <typename C, typename BT>
num MComplex<C,BT>::Alternator(const BT& birth = INITBT, bool savegens = false)
{
	if (TIMEINFO)
	{
		debout<<"\n*AltReduce";
		//cin.get();
	}

	typename COMPLEX::iterator cit = this->clist.find(birth);
	if (cit == this->clist.end()) return 0;

	pair<num,bool> numrem;

	num toret = 0;
	num thisiter = 1, lastiter = 1;
    num oldsz = this->fsize(birth);

	if (this->fsize(birth) == 0) return 0; // if you have killed the frame...

	// new birth, reset minindex!
	minindex = 0;
	oldsz = this->size();
	num topmass = this->getTopDim(birth).first;
	num botmass = this->getBotDim(birth).first;

    bool redfirst = (topmass > botmass);

	while (thisiter != 0 && lastiter != 0)
	{
	    lastiter = thisiter;
	    thisiter = 0;

        if (redfirst)
        {
            numrem = RedMorsify(birth, !savegens);
            markUncritical_Red(birth);
            thisiter = numrem.first + freeCofaceCollapser(birth) + freeFaceCollapser(birth);
            numrem = CoredMorsify(birth, !savegens);
            markUncritical_Cored(birth);
            thisiter += numrem.first + freeFaceCollapser(birth) + freeCofaceCollapser(birth);
            debout<<"\n***** removed "<<thisiter<<" red -> cored pairs!";
        }
        else
        {
            numrem = CoredMorsify(birth, !savegens);
            markUncritical_Cored(birth);
            thisiter = numrem.first + freeFaceCollapser(birth) + freeCofaceCollapser(birth);
            numrem = RedMorsify(birth, !savegens);
            markUncritical_Red(birth);
            thisiter += numrem.first + freeCofaceCollapser(birth) + freeFaceCollapser(birth);
            debout<<"\n***** removed "<<thisiter<<" cored -> red pairs!";
        }
        toret += thisiter;
        redfirst = !redfirst;

        if (!savegens && 2*this->fsize(birth) <= oldsz)
		{
			memorySaver(birth);
		}

	}
	if (BIGMORSTALK)
	{
		debout<<"\n***Altrem returned: "<<toret;
        //cout<<*this;
	}

	if(TIMEINFO)
	{
		debout<<"*/altreduce";
		if (MAKEBPS) cin.get();
	}
	return toret;
}





// repeatedly performs morse coreduction until no new cells are removed,
// EXCEPT for cells of birthtime given by exception!
template <typename C, typename BT>
num MComplex<C,BT>::MorseCoreduceExcept(const BT& start = INITBT, const BT& exception = BANBT, bool savegens = false)
{
	if (MOVTALK)
	{
		debout<<"\n*MorseCoreduce-Except "<<exception;
		//cin.get();
	}
	// reset min critical cell index, everything starts from scratch now

	typename COMPLEX::const_iterator cit;
	pair<num,bool> numrem;
	//vector<BT> origbirths;
	BT birth;
	num toret = 0;

	// store all encountered birth times in origbirths...
	for(cit = this->clist.lower_bound(start); cit != this->clist.end(); ++cit)
	{
		// new birth: reset minindex
		minindex = 0;
		birth = cit->first;
		if (birth >= exception) break; // halt if birth exceeds exceptional birth!
		if (MOVTALK)
		{
			debout<<"\n\n*****BIRTH: "<<birth;
			//cin.get();
		}
		// now we morse-reduce
		numrem = CoredMorsify(birth, !savegens, exception); // deep removal if generators are not needed!
		// if this frame has been swallowed, move on to the next one!
//		if (this->fsize(birth)==0)
//		{
//			//cout<<"\nswallowed: "<<birth;
//			continue;
//		}
		markUncritical_Cored(birth); // if frame has crit cells left, update!
		//cout<<*this;
		//verify(birth);
		toret += numrem.first;

		if (MOVTALK)
		{
			debout<<"\n***Cmfy returned: "<<numrem.first<<" "<<numrem.second;
			debout<<", pairs removed so far: "<<toret<<"*******\n"; //cin.get();

		}
	}
	if(TIMEINFO)
	{
		debout<<"*/morsecoreduce";
		cin.get();
	}
	return toret;
}


// verifies status of gradient paths etc.
template <typename C, typename BT>
bool MComplex<C,BT>::verify(const BT& birth)
{
	// seek frame
	typename COMPLEX::iterator finder = this->clist.find(birth);
	if (finder == this->clist.end()) return false;

	FRAME* frm = finder->second; // extract frame
	typename FRAME::iterator fiter;
	// list variables
	CLIST* curl;
	typename CLIST::iterator liter;
	Cell<C,BT>* curcell;
	for (fiter = frm->begin(); fiter != frm->end(); ++fiter)
	{
		curl = fiter->second;
		for (liter = curl->begin(); liter != curl->end(); ++liter)
		{
			curcell = *liter;
			if (curcell == NULL) continue;
			if (!(curcell -> isIn)) continue;

			if (GPTALK && curcell->gpath != NULL)
			{
				debout<<"\nBad gpath on Cell: "<<*curcell;
				//cin.get();
				return false;
			}
		}
		//cout<<"\n gpath integrity verified!";
		return true;
	}
	return false;
}


// call this after redmorsify  to mark cells uncritical for next round
template <typename C, typename BT>
bool MComplex<C,BT>::markUncritical_Red(const BT& stop)
{
    bool tracer = false;
    //if (stop == 1) tracer = true;

	if (UCTALK || CMFTALK || tracer)
	{
		debout<<"\nUncriting for birth "<<stop<<" fsize "<<this->fsize(stop);
        if (MAKEBPS) cin.get();
		//cin.get();
	}
	typename COMPLEX::const_iterator cit = this->clist.find(stop);
	if (cit == this->clist.end()) return false; // frame not found!
	//if (this->fsize(stop) == 0) return false; // no cells in this frame!

	// otherwise, we have a frame whose cells we want to make uncritical!
	FRAME* touncrit = cit->second;
	typename FRAME::const_reverse_iterator fit;

	CLIST* curlist;
	typename CLIST::const_iterator lit;
	Cell<C,BT>* curcell;

	// iterate over frame...
	for (fit = touncrit->rbegin(); fit != touncrit->rend(); ++fit)
	{
		curlist = fit->second;
		for (lit = curlist->begin(); lit != curlist->end(); ++lit)
		{
			curcell = *lit;
			if (curcell ==  NULL) continue; // ignore nulls
			if (!(curcell->isIn))
			{
				//cout<<"\n------------------------ outcell???? "<<*curcell;
				continue; // and removeds
			}

			if (!curcell->isCrit && !hyperq)
			{
			    // this should never happen.
				cout<<"\n URGENT ERROR::: unprocessed cell??? "<<*curcell;
			}


			curcell->isCrit = false; // reset critical status
			// restore boundary from gradient path
			//curcell->bd.clear();
			if (curcell->gpath != NULL)
			{
				if (GPMEMTALK || tracer)
				{
				    debout<<"\n... restoring cb from gp = "
				                   <<*(curcell->gpath)<<" of "<<*curcell;
                    //cin.get();
				}
				//if (MAKEBPS) cin.get();
				curcell->setCB(*(curcell->gpath));
				delete curcell->gpath;
				curcell->gpath = NULL;
			}
		}
	}
	return true;
}

// call this after coredmorsify  to mark cells uncritical for next round
template <typename C, typename BT>
bool MComplex<C,BT>::markUncritical_Cored(const BT& stop)
{
	//num count = 0;
	if (CMFTALK)
	{
		debout<<"\n\n********Marking Uncritical for time: "<<stop;
		if (MAKEBPS) cin.get();
	}
	if (UCTALK)
	{
		debout<<"\nUncriting for birth "<<stop<<" fsize "<<this->fsize(stop);
		//cin.get();
	}
	typename COMPLEX::const_iterator cit = this->clist.find(stop);
	if (cit == this->clist.end()) return false; // frame not found!
	//if (this->fsize(stop) == 0) return false; // no cells in this frame!

	// otherwise, we have a frame whose cells we want to make uncritical!
	FRAME* touncrit = cit->second;
	typename FRAME::const_iterator fit;

	CLIST* curlist;
	typename CLIST::const_iterator lit;
	Cell<C,BT>* curcell;

	// iterate over frame...
	for (fit = touncrit->begin(); fit != touncrit->end(); ++fit)
	{
		curlist = fit->second;
		for (lit = curlist->begin(); lit != curlist->end(); ++lit)
		{
			curcell = *lit;
			if (curcell ==  NULL) continue; // ignore nulls
			if (!(curcell->isIn))
			{
				//cout<<"\n------------------------ outcell???? "<<*curcell;
				continue; // and removeds
			}

			if (!curcell->isCrit && !hyperq)
			{
				cout<<"\n URGENT ERROR: unprocessed cell??? "<<*curcell;
			}


			curcell->isCrit = false; // reset critical status
			// restore boundary from gradient path
			//curcell->bd.clear();
			if (curcell->gpath != NULL)
			{
				if (GPMEMTALK) debout<<"\n... restoring bd from gp = "<<*(curcell->gpath)<<" of "<<*curcell;
				//if (MAKEBPS) cin.get();
				curcell->setBD(*(curcell->gpath));
				delete curcell->gpath;
				curcell->gpath = NULL;
			}
		}
	}
	return true;
}

template <typename C, typename BT>
num MComplex<C,BT>::MorseWrapper_Cored(bool savegens = false, double thresh = 0.3)
{
	num premct, oldsz; // number of pairs removed in each iteration of the morse tower
	num tht = 1; // height of morse tower
	num totrem = 0;
	float fracrem = 0;
	do
	{
		//cout<<"in..."; cin.get();
		oldsz = this->size();
		premct = MorseCoreduce(savegens);
		//cout<<*this;
		totrem += premct;
		fracrem = (oldsz == 0) ? 1 : float(2*premct)/float(oldsz);
		if (FLOWTALK)
		{
			debout<<"\n******Performed "<<premct<<" Morse pair coreductions at height "<<tht<<", ";
			debout<<" fraction reduced = "<<fracrem;
			if (MAKEBPS) cin.get();
		}
		//if (fracrem < thresh) break;
		tht++;
		//cout<<" fracrem "<<fracrem<<" vs thresh "<<thresh; cin.get();
	} while (fracrem > thresh && premct != 0);
	return totrem;
}

template <typename C, typename BT>
num MComplex<C,BT>::MorseWrapper_Red( bool savegens = false, double thresh = 0.3)
{
	num premct, oldsz; // number of pairs removed in each iteration of the morse tower
	num tht = 1; // height of morse tower
	num totrem = 0;
	float fracrem = 0;
	do
	{
		//cout<<"in..."; cin.get();
		oldsz = this->size();
		premct = MorseReduce(savegens);
		totrem += premct;
		fracrem = (oldsz == 0) ? 1 : float(2*premct)/float(oldsz);
		//cout<<ccomp;
		//cin.get();
		if (FLOWTALK)
		{
			debout<<"\n******Performed "<<premct<<" Morse pair reductions at height "<<tht<<", ";
			debout<<" fraction reduced = "<<fracrem;
			if (MAKEBPS) cin.get();
			//ccomp.printList(3,2);
		}
		tht++;
		//cout<<" fracrem "<<fracrem<<" vs thresh "<<thresh; cin.get();
	} while (fracrem > thresh && premct != 0);
	return totrem;
}

template <typename C, typename BT>
void MComplex<C,BT>::MorseWrapperExcept_Cored(const BT& start, const BT& except, bool savegens = false,
                                                double thresh = 0.3)
{
	num premct, oldsz; // number of pairs removed in each iteration of the morse tower
	num tht = 1; // height of morse tower
	float fracrem;
	do
	{
		//cout<<"in..."; cin.get();
		oldsz = this->size();
		premct = MorseCoreduceExcept(start, except, savegens);
		fracrem = (oldsz == 0) ? 1 : float(2*premct)/float(oldsz);
		//cout<<ccomp;
		//cin.get();
		if (FLOWTALK)
		{
			debout<<"\n******Performed "<<premct<<" Morse pair coreductions at height "<<tht<<", ";
			debout<<" fraction reduced = "<<fracrem;
			if (fracrem < thresh) break;
			//ccomp.printList(3,2);
		}
		if (MAKEBPS) cin.get();
		tht++;
	} while (premct > 0);
}

template <typename C, typename BT>
void MComplex<C,BT>::MorseWrapperExcept_Red(const BT& start, const BT& except, bool savegens = false,
                                            double thresh = 0.3)
{
	num premct, oldsz; // number of pairs removed in each iteration of the morse tower
	num tht = 1; // height of morse tower
	float fracrem;
	do
	{
		//cout<<"in..."; cin.get();
		oldsz = this->size();
		premct = MorseReduceExcept(start, except, savegens);
		fracrem = (oldsz == 0) ? 1 : float(2*premct)/float(oldsz);
		//cout<<ccomp;
		//cin.get();
		if (FLOWTALK)
		{
			debout<<"\n******Performed "<<premct<<" Morse pair reductions at height "<<tht<<", ";
			debout<<" fraction reduced = "<<fracrem;
			if (fracrem < 0.3) break;
			//ccomp.printList(3,2);
            if (MAKEBPS) cin.get();
		}
		tht++;
	} while (premct > 0);
}



template <typename C, typename BT>
void CoreduceAndUpdate (MComplex<C,BT>*& mycomp, bool savegens = false,
                        double thresh = 0.5, bool superopt = false, ostream& out = cout)
{
	MComplex<C,BT>* rcomp = NULL; // reduced complex!
	num premct, oldsz; // number of pairs removed in each iteration of the morse tower
	num tht = 1; // height of morse tower
	float fracrem;
	do
	{
		//cout<<"in..."; cin.get();
		oldsz = mycomp->size();
		premct = mycomp->MorseCoreduce(savegens);
		fracrem = (oldsz == 0) ? 1 : float(2*premct)/float(oldsz);
		//cout<<ccomp;
		//cin.get();
		if (FLOWTALK)
		{
			out<<"\n+++coreductions: "<<oldsz<<" -> "<<oldsz - 2*premct
               <<", fraction removed = "<<fracrem<<" at height "<<tht;
			//cout<<*mycomp;
			//cin.get();
			//ccomp.printList(3,2);
			if (MAKEBPS) cin.get();

		}
		tht++;

		if (fracrem > thresh && !savegens) // we want to keep reducing!
		{
			if (fracrem >= 0.5 || superopt) // want to move over?
			{
			    // allocate space...
                rcomp = new MComplex<C,BT>;
                // move over to new reduced complex
                mycomp->moveOver_Sort(*rcomp,false); // optimize for coreduction
                //mycomp->moveOver(*rcomp);
                //cout<<" moved over, new size: "<<rcomp->size(); cin.get();
                // remove old complex
                delete mycomp;
                // set old pointer to new complex
                mycomp = rcomp;
			}
		}
		else break; // we are done, not enough cells were reduced
	} while (premct > 0);
}


template <typename C, typename BT>
void AlternateAndUpdate_makeStat (MComplex<C,BT>*& mycomp,  vector<num>& svec,
                                  double thresh = 0.0)
{
    svec.clear();
    double movewhen = 0.0;
    svec.push_back(mycomp->size()); // initialize svec with unreduced size...
    MComplex<C,BT>* rcomp = NULL; // reduced complex!
    num oldsz, newsz;
    bool redcycle = false; // if false, we coreduce first otherwise we reduce first
    float fracrem = 1, oldfrem = 1;
    num tht = 0; // tower height

    do
    {
        tht++;
        oldfrem = fracrem;
        oldsz = mycomp->size();
        mycomp->TopAlternator(redcycle,thresh);
        newsz = mycomp->size();

        fracrem = (oldsz == 0) ? 1 : float(oldsz-newsz)/float(oldsz);

        //out<<" \n+++"<<mtype<<" "<<oldsz<<" -> "<<newsz<<", fraction removed "<<fracrem<<" at height "<<tht;

        redcycle = !redcycle; // if reducing, next cycle is coreduce and vice versa

        if ( (fracrem > thresh) || (oldfrem > thresh) ) // we want to keep reducing!
		{

		    // non-trivial size reduction: note down size in svec!
		    if (oldsz != newsz) svec.push_back(newsz);

			if (fracrem > movewhen) // want to move over?
			{
			    // allocate space...
                rcomp = new MComplex<C,BT>;
                // move over to new reduced complex
                mycomp->moveOver(*rcomp);
                //cout<<" moved over, new size: "<<rcomp->size(); cin.get();
                //remove old complex
                delete mycomp;
                // set old pointer to new complex
                mycomp = rcomp;
			}
		}
		else break; // we are done, not enough cells were reduced
    } while (fracrem > thresh || oldfrem > thresh);
}



template <typename C, typename BT>
num AlternateAndUpdate (MComplex<C,BT>*& mycomp,  bool savegens = false, double thresh = 0.5,
                         bool redfirst = false, bool superopt = false, ostream& out = cout)
{
    double movewhen = (superopt) ? 0.0 : 0.5;
    MComplex<C,BT>* rcomp = NULL; // reduced complex!
    num oldsz, newsz;
    bool tracer = false; // trace off by default
    float fracrem = 1, oldfrem = 1;
    string mtype;
    num tht = 0;

    do
    {
        tht++;

        oldfrem = fracrem;
        mtype = (redfirst) ? "reductions:   " : "coreductions: ";
        oldsz = mycomp->size();

        mycomp->TopAlternator(redfirst,thresh,savegens,tracer);
        newsz = mycomp->size();
        fracrem = (oldsz == 0) ? 1 : float(oldsz-newsz)/float(oldsz);

        out<<" \n+++"<<mtype<<" "<<oldsz<<" --> "<<newsz
           <<", fraction removed "<<fracrem<<" at height "<<tht;

        if (fracrem <= thresh && oldfrem <= thresh) return tht; // we are done, not enough cells were reduced

		// otherwise, toggle cored <-> red
		redfirst = !redfirst;
		// optimize complex for chosen reduction type:
        if (superopt && !savegens)
		{
             rcomp = new MComplex<C,BT>;
            mycomp->moveOver(*rcomp);
            delete mycomp;
            mycomp = rcomp;
		}
    } while (1);
}

template <typename C, typename BT>
void ReduceAndUpdate (MComplex<C,BT>*& mycomp, bool savegens = false,
                      double thresh = 0.5, bool superopt = false, ostream& out = cout)
{
	MComplex<C,BT>* rcomp = NULL; // reduced complex!
	num premct, oldsz; // number of pairs removed in each iteration of the morse tower
	num tht = 1; // height of morse tower
	float fracrem;
	do
	{
		//cout<<"in..."; cin.get();
		oldsz = mycomp->size();
		premct = mycomp->MorseReduce(savegens);
		fracrem = (oldsz == 0) ? 1 : float(2*premct)/float(oldsz);
		//cout<<ccomp;
		//cin.get();
		if (FLOWTALK)
		{
			out<<"\n+++reductions:   "<<oldsz<<" -> "<<oldsz - 2*premct;
			out<<", fraction removed = "<<fracrem<<" at height = "<<tht;
			//cout<<*mycomp;
			//cin.get();
			//ccomp.printList(3,2);
			if (MAKEBPS) cin.get();

		}
		tht++;

		if (fracrem > thresh && !savegens) // we want to keep reducing!
		{
			if (fracrem >= 0.5 || superopt) // want to move over?
			{
			    // allocate space...
                rcomp = new MComplex<C,BT>;
                // move over to new reduced complex
                mycomp->moveOver_Sort(*rcomp, true); // optimize for reduction
                //cout<<" moved over, new size: "<<rcomp->size(); cin.get();
                // remove old complex
                delete mycomp;
                // set old pointer to new complex
                mycomp = rcomp;
			}
		}
		else break; // we are done, not enough cells were reduced
	} while (premct > 0);
}

template <typename C, typename BT>
void AlternateAndUpdate_Shuffle (MComplex<C,BT>*& mycomp,  num attempts=0, ostream& out = cout)
{
    MComplex<C,BT>* rcomp = NULL; // reduced complex!
    num oldsz, newsz;
    bool redfirst = false; // trace off by default
    float fracrem = 1, oldfrem = 1;
    string mtype;
    num tht = 0;

    for (num i=0; i < attempts; i++)
    {
        tht++;


        // first reduce!
        oldfrem = fracrem;
        mtype = (redfirst) ? "reductions:   " : "coreductions: ";
        oldsz = mycomp->size();

        mycomp->TopAlternator(redfirst,0.0,false,false);
        newsz = mycomp->size();
        fracrem = (oldsz == 0) ? 1 : float(oldsz-newsz)/float(oldsz);

        out<<" \n+++"<<mtype<<" "<<oldsz<<" --> "<<newsz
           <<", fraction removed "<<fracrem<<" at height "<<tht;

       redfirst = !redfirst;

       // then shuffle!

       rcomp = new MComplex<C,BT>;
       mycomp->moveOver_Shuffle(*rcomp);
       delete mycomp;
       mycomp = rcomp;
    }
		// optimize complex for chosen reduction type:
}


// flushes out queue statuses of cells in this list
template <typename C, typename BT>
void flushQ(const vector<Cell<C,BT>*>& toflush)
{
    //cout<<"\n flushQ size: "<<toflush.size();
   //cin.get();
	typename vector<Cell<C,BT>*>::const_iterator lit;
	for(lit = toflush.begin(); lit!=toflush.end(); ++lit)
	{
		(*lit)->isQueued = false;
	}
}

#endif /* MORSERED_HPP_ */
