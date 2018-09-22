/*
 * PComplex.hpp
 * Basic Complex Methods: Insertion, Removal, Check for containment
 */


#ifndef COMPLEX_HPP_
#define COMPLEX_HPP_

#include "Complex.h"

// returns size pair for given birth and dimension, i.e. a pair, the first
// number gives the total number of cells of this birth and dimension still
// in the complex, the second gives the largest indexed cell still present

// NO CHECKS ARE MADE HERE. you should make sure that there have been cells
// of the given birth and dimension via clist or sizeinfo lookup prior to
// calling this function
template <typename C, typename BT>
pair<num,num> Complex<C,BT>::getSizeInfo(const BT& birth, const num& dim) const
{
	pair<num,num> toret = make_pair(0,0);
	typename DMAP::const_iterator mit = sizeinfo.find(birth);
	if (mit == sizeinfo.end()) return toret;
	DSIZE* cursz = mit->second;
	typename DSIZE::const_iterator sit = cursz->find(dim);
	if (sit == cursz->end()) return toret;
	return sit->second;
}

// how many cells of this birth and dimension?
template <typename C, typename BT>
num Complex<C,BT>::lsize(const BT birth, const num dim) const
{
	// search for this birth time
	typename COMPLEX::const_iterator cit = clist.find(birth); // complex iterator
	if (cit == clist.end()) return (num) 0; // frame not found!

	// assign frame, birth time found!
	FRAME* curframe = cit->second;
	typename DMAP::const_iterator dit = sizeinfo.find(birth);

	// search for this dimension
	typename FRAME::const_iterator fit = curframe->find(dim); // frame iterator
	if (fit == curframe->end()) return (num) 0; // frame found, but dimension not found

	typename DSIZE::const_iterator dsit = dit->second->find(dim);

/*	// return list size for cells still considered "in"
	CLIST* curlist = fit->second;
	typename CLIST::const_iterator curcell;
	num toret = 0;

	for(curcell = curlist->begin(); curcell != curlist->end(); ++curcell)
	{
		if(*curcell == NULL) continue; // ignore null cells...
		if(!(*curcell)->isIn) continue; // or 'removed' cells...
		// otherwise, count up the size!
		toret++;
	}*/
	//both frame and dim found, so return from sizeinfo
	return dsit->second.first;
}

// how many cells in this frame?
template <typename C, typename BT>
num Complex<C,BT>::fsize(const BT birth) const
{
	// search for this birth time
	typename COMPLEX::const_iterator cit = clist.find(birth); // complex iterator
	if (cit == clist.end()) return (num) 0; // frame not found!

	// get frame
	FRAME* curframe = cit->second;
	num toret = 0;
	// iterate over frame...
	typename FRAME::const_iterator fit;
	for (fit = curframe->begin(); fit != curframe->end(); ++fit)
	{
		toret += lsize(birth,fit->first); // call list size for this birth with current dimension
	}
	return toret;
}

// reallocates a new vector for cells located at (birth,dim) in
// the complex structure that removes all the NULLs... call when
// getSizeInfo(birth,dim) is much less than the size() of the
// relevant vector. returns the size saving in terms of number
// of cells-pointers.
template <typename C, typename BT>
num Complex<C,BT>::pruneList(const BT& birth, const num dim)
{
	num toret = 0;
	typename COMPLEX::const_iterator compit = clist.find(birth);
	if (compit == clist.end()) return toret;

	FRAME* curf = compit->second;
	typename FRAME::iterator frit = curf->find(dim);
	if (frit == curf->end()) return toret;

	// if we reach here the cell vector to prune is frit->second
	CLIST* toprune = frit->second;
	typename CLIST::iterator liter;

	CLIST* newlist = new CLIST; // to copy over...

	for (liter = toprune->begin(); liter != toprune->end(); ++liter)
	{
		if (*liter == NULL) // if we want to ignore this cell...
		{
			toret++; // add up space saving...
			continue; // and go on without copying
		}
		// otherwise,... copy and reset index
		newlist->push_back(*liter);
		(*liter)->ind = (num)(newlist->size()-1);
	}
	// memory saving on new list...
	CLIST(*newlist).swap(*newlist);
	// and on the old
	toprune->clear();
	CLIST(*toprune).swap(*toprune);
	delete toprune;
	// assign the new pruned list to the complex
	frit->second = newlist;

	// update sizeinfo structure now... the size should be the same as before
	// but the index of the last element is now the index of the last element
	// in newlist
	((*(sizeinfo[birth]))[dim]).second = (num)(newlist->size()-1);
	return toret;
}


// how many cells in the complex?
template <typename C, typename BT>
num Complex<C,BT>::size() const
{
	num toret = 0;
	typename COMPLEX::const_iterator cit;
	for (cit = clist.begin(); cit != clist.end(); ++cit)
	{
		toret += fsize(cit->first);
	}
	return toret;
}


// size computations! how many cells of this dimension across all birth times
// strictly less than cap (optional, defaults to answering across all birth times)?
template <typename C, typename BT>
num Complex<C,BT>::dsize(const num dim, BT cap = BANBT) const
{
	typename COMPLEX::const_iterator cit; // complex iterator
	typename FRAME::const_iterator fit; // frame iterator
	FRAME* curframe;

	typename COMPLEX::const_iterator findcap = clist.find(cap);

	num toret = 0;

	for (cit = clist.begin(); cit != findcap; ++cit)
	{
		curframe = cit->second;
		fit = curframe->find(dim); // locate cells of this dimension
		if (fit == curframe->end()) continue; // if none, move on to next frame

		toret += lsize(cit->first,dim); // add size of list for this dimension

	}
	return toret;

}


template <typename C, typename BT>
void Complex<C,BT>::printSizeInfo(bool dimdetail = false, bool framedetail = true, ostream& out = cout) const
{
	out<<"\n**********Complex Size: "<<size()<<"*************";

	if (framedetail)
	{
		typename COMPLEX::const_iterator cit; // complex iterator
		typename FRAME::const_iterator fit;
		FRAME* curframe;
		for (cit = clist.begin(); cit != clist.end(); ++cit)
		{
			out<<"\n       Birth: ["<<cit->first<<"], Size: "<<fsize(cit->first);
			if (dimdetail)
			{
				curframe = cit->second;
				for (fit = curframe->begin(); fit != curframe->end(); ++fit)
				{
					out<<"\n                     Dim: <"<<fit->first<<">, Size: "<<lsize(cit->first, fit->first);
				}
			}
		}
	}
}




// check del-square = 0, complex property!
template <typename C, typename BT>
bool Complex<C,BT>::checkComplex(ostream& out = cout) const
{
	typename COMPLEX::const_iterator cit; // complex iterator
	typename FRAME::const_iterator fit; // frame iterator
	typename CLIST::const_iterator lit; // list iterator

	FRAME* curframe;
	CLIST* curlist;
	Cell<C,BT>* curcell;
	bool retval = true;
	typename CSTRUCT::const_iterator chit;

	for (cit = clist.begin(); cit != clist.end(); ++cit)
	{
		//out<<"\n BCHK: "<<cit->first;

		curframe = cit->second;
		for (fit = curframe->begin(); fit != curframe->end(); ++fit)
		{

			//out<<"\n     DCHK "<<fit->first;

			curlist = fit->second;
			for (lit = curlist->begin(); lit != curlist->end(); ++lit)
			{
				curcell = *lit;

				if(curcell == NULL) continue; // ignore null
				if(!curcell->isIn) continue; // and removeds

				if (curcell->getDim() == 0) break; // no need to check dim 0!
				// check boundary operator
				if (curcell->getBD().getBD().size() != 0)
				{
					out<<"\n\n Violation of Complex Structure!!!";
					out<<"\n @ cell "<<*curcell<<" with bd "<<curcell->getBD();
					out<<" \n and bd^2 is "<<curcell->getBD().getBD()<<"\n";
					out<<" note index = "<<curcell->ind;
					//cin.get();
					retval = false;
				}
				//check coboundary operator
				if (curcell->getCB().getCB().size() != 0)
				{
					out<<"\n\n Violation of Complex Structure!!!";
					out<<"\n @ cell "<<*curcell<<" with cb "<<curcell->getCB();
					out<<" \n and cb^2 is "<<curcell->getCB().getCB()<<"\n";
					out<<" note index = "<<curcell->ind;
					//cin.get();
					retval = false;
				}
				// check birth times at boundary:
				for (chit = curcell->bd.begin(); chit != curcell->bd.end(); ++chit)
				{
				    if (chit->first->birth > curcell->birth)
				    {
				        out<<"\n\n Violation of Filtered Structure!!!";
				        out<<"\n @ cell "<<*curcell<<" born "<<curcell->birth
				           <<" with bd "<<(curcell->getBD())<<" cells born at "
				           <<chit->first->birth<<"\n";
				        retval = false;
				        cin.get();
				    }
                }

			}
		}
	}
	//cout<<" Dunn ";
	return retval;
}


// quick insert... inserts an entire vector of given "birth" and "dimension".
template <typename C, typename BT>
void Complex<C,BT>::quicksert(const BT birth, const num dim, const vector<Cell<C,BT>*>& toins)
{
	//cout<<"\n\n********QUICK INSERT at ["<<birth<<", "<<dim<<"] ***************";
	FRAME* curf; // current frame to insert in, corresponding to birth
	CLIST* curl; // current list to insert in, corresponding to dimension

	// look for this birth time
	typename COMPLEX::iterator finder = clist.find(birth);
	typename FRAME::iterator dimloc;

	DSIZE* curdsz;

	if (finder == clist.end()) // if not found...
	{
		//cout<<"\nFrame for birth "<<birth<<" created!";
		curf = new FRAME; // make new frame
		curl = new CLIST; // and new list
		curf->insert(curf->begin(),make_pair(dim,curl)); // add list to frame for correct dimension
		clist.insert(clist.begin(),make_pair(birth,curf)); // and frame to complex at right birth

		// create size info!
		curdsz = new DSIZE;
		(*curdsz)[dim] = make_pair(0,0); // initialize to nothing for this new birth + dimension...
		sizeinfo[birth] = curdsz; // initialize to nothing for this new birth + dimension...
	}
	else // if found,...
	{
		//cout<<"\nFrame for birth "<<birth<<" already existed!";
		curf = finder->second; // extract frame!
		curdsz = sizeinfo[birth]; // read list of dims from sizeinfo

		// now search for list corresponding to this dimension
		dimloc = curf->find(dim);
		if (dimloc == curf->end()) // not found!
		{
			//cout<<"\nDim for birth "<<birth<<" already existed!";
			curl = new CLIST;
			curf->insert(curf->end(),make_pair(dim,curl));
			// size info!
			(*curdsz)[dim] = make_pair(0,0); // initialize to nothing for this new dimension...
		}
		else
		{
			curl = dimloc->second;
		}
	}


	//cout<<"\n>>>>>>vec ins... "; cin.get();
	// now we merely need to start inserting into curl!

	// first use assign to save memory (???)
	curl->reserve(toins.size());
	curl->assign(toins.begin(),toins.end());

	// now loop over to change indices etc.
	for (num i = 0; i < (num)toins.size(); ++i)
	{
		curl->at(i)->isIn = true;
		curl->at(i)->ind = i;
	}
	(*curdsz)[dim].first += toins.size(); // # cells added incremented
	(*curdsz)[dim].second += toins.size()-1; // highest index incremented
	//cout<<">>>>>done!!!"; cin.get();
	// memory saver?
	CLIST(*curl).swap(*curl);
}

// prints cells of dimension dim born at time birth.
template <typename C, typename BT>
void Complex<C,BT>::printList(const BT& birth, const num dim, ostream& out = cout) const
{
	typename COMPLEX::const_iterator cit; // complex iterator
	typename FRAME::const_iterator fit; // frame iterator
	typename CLIST::const_iterator lit; // list iterator

	FRAME* curframe;
	CLIST* curlist;

	Cell<C,BT>* curcell;

	cit = clist.find(birth);
	if (cit != clist.end())
	{
		curframe = cit->second;
		fit = curframe->find(dim);
		if (fit != curframe->end())
		{
			curlist = fit->second;
			out<<"\n*** Frame "<<birth<<", Dim: "<<fit->first<<"  ["<<lsize(cit->first,fit->first)<<"]\n";
			for (lit = curlist->begin(); lit != curlist->end(); ++lit)
			{
				curcell = *lit;
				if (curcell == NULL) continue;
				if (!(curcell->isIn)) continue;
				//cin.get();
				out<<"\n   #"<<curcell->ind<<". "<<*curcell  //<<" // qstat = "<<curcell->isQueued
				   <<"\n            bd: "<<curcell->getBD();
				   //<<"\n            cb: "<<curcell->getCB()<<"\n";
			}
		}
		else out<<"\n*** List for dim "<<dim<<" empty in frame "<<birth;
		out<<"\n*****************************\n";
	}
	else out <<"\n***Frame "<<birth<<" is empty!!";
}


// prints cells of every dimension born at time "birth"
template <typename C, typename BT>
void Complex<C,BT>::printFrame(const BT& birth, ostream& out = cout) const
{
	typename COMPLEX::const_iterator cit; // complex iterator
	typename FRAME::const_iterator fit; // frame iterator
	typename CLIST::const_iterator lit; // list iterator

	FRAME* curframe;

	cit = clist.find(birth);
	if (cit != clist.end())
	{
		out<<"\n     >>>***BIRTH "<<cit->first<<", #Dims: ["<<cit->second->size()<<"] ***<<<";
		curframe = cit->second;
		for (fit = curframe->begin(); fit != curframe->end(); ++fit)
		{
			printList(birth,fit->first,out);
		}
		out<<"---------------------------\n";
	}
	else out <<"\n***Frame "<<birth<<" is empty!!";
}

// inserts a cell pointed to by "toin"
template <typename C, typename BT>
bool Complex<C,BT>::insertCell(Cell<C,BT>* toin)
{
	if (toin == NULL) return false; // can't add null cell
	if (CTALK) cout<<"\n\n INSERT: "<<*toin<<" at "<<toin;
	BT born = toin->birth;
	num dim = toin->getDim();

	// the frame and list into which we will insert this cell
	FRAME* myframe; // frame of right birth time
	CLIST* mylist; // list of right dimension in frame

	// see if the frame for this birth time exists...
	typename COMPLEX::iterator findbt = clist.find(born);
	DSIZE* cursz;

	if (findbt == clist.end()) //if we don't have this frame,...
	{
		if (CTALK) cout<<"\n\n   New Frame for birth "<<born;
		myframe = new FRAME; // make the frame!

		if (CTALK) cout<<"\n       ... and new list for dim "<<dim;
		mylist = new CLIST; // and then add a list for this dimension
		// insert list for this dimension into the frame
		myframe->insert(make_pair(dim,mylist));
		// and this frame into the complex
		clist.insert(make_pair(born,myframe));

		// also create size information!
		cursz = new DSIZE;
		(*cursz)[dim] = make_pair(1,0); // 1 cell of this dim, with index 0
		sizeinfo[born] = cursz; // add to sizeinfo for this frame

	}
	else // already have this frame
	{
		myframe = findbt->second;
		// ok, but do we have a list for this dimension?
		typename FRAME::iterator finddim = myframe->find(dim);
		cursz = sizeinfo[born]; // this frame must already be in complex size info
		if (finddim == myframe->end()) // if this dim's list does not exist...
		{
			if (CTALK) cout<<"\n   New list for birth/dim "<<born<<"/"<<dim;
			mylist = new CLIST; // make the list!
			myframe->insert(make_pair(dim,mylist));
			//and size info...
			(*cursz)[dim] = make_pair(1,0); // 1 cell of this dim, with index 0
		}
		else // frame exists and list exists!
		{
			mylist = finddim->second;
			if (CTALK) cout<<"\n    insert birth/dim "<<born<<"/"<<dim;
			// size:
			(*cursz)[dim].first++; // update number of cells of this dimension
			(*cursz)[dim].second++; // and the index of the last one!
		}
	}

	mylist->push_back(toin); // insert this cell into list
	toin->isIn = true;
	toin->ind = mylist->size()-1; // set cell's index

	return true;
}

// if you inserted cells one by one, you may want to squeeze the vectors and save
// some memory!
template <typename C, typename BT>
void Complex<C,BT>::squeeze()
{
	typename COMPLEX::iterator compit;
	FRAME* curf;
	typename FRAME::iterator fiter;
	CLIST* curl;

	for (compit = clist.begin(); compit != clist.end(); ++compit)
	{
		curf = compit->second;
		for (fiter = curf->begin(); fiter != curf->end(); ++fiter)
		{
			curl = fiter->second;
			// swap trick
			CLIST(*curl).swap(*curl);
		}
	}

}


// checks for containment of given cell pointer... returns iterator to the cell
// pointer if found at the desired index.
template <typename C, typename BT>
bool Complex<C,BT>::isIn(const Cell<C,BT>* toCheck, typename CLIST::iterator& pos)
{
	if (toCheck == NULL) return false; // null, so ignore
	if (!(toCheck->isIn)) return false; // not in, who cares...

	BT birth = toCheck->birth;
	num dim = toCheck->getDim(); // dimension of cell to find

	typename COMPLEX::const_iterator cit = clist.find(birth);
	if (cit == clist.end()) return false;

	typename FRAME::const_iterator fit = cit->second->find(dim);
	if (fit == cit->second->end()) return false;

	CLIST* curlist = fit->second; // the cells for this dimension!
	if (toCheck->ind >= (num)curlist->size()) return false;

	// okay, now a simple check to see if the vector element at this index matches...
	pos = curlist->begin()+toCheck->ind;
	if(*pos == toCheck) return true;

	return false;
}


// removes a cell from the complex along with references to it from
// boundary and co-boundary cells...
// if remc is true, the cell is actually removed, from the complex, else not
// if fmbd is true, it is removed from its coboundaries' boundaries, else not
// if fmcb is true, it is removed from its boundaries' coboundaries, else not
// an iterator to the cell is returned (if found!) via cpos...
template <typename C, typename BT>
bool Complex<C,BT>::removeCell(Cell<C,BT>* tokill, typename CLIST::iterator& cpos, bool remc = true,
		                       bool fmbd = true, bool fmcb = true)
{
	if (tokill == NULL) return false; // can't kill the null-cell, it is already dead...
	//cout<<"\nRemoving: "<<*tokill; cin.get();


	// get iterator to the cell pointer in this complex
	if (!isIn(tokill,cpos)) return false;

	// now we have the cell pointer as *pos if it was actually found in the complex
	if(remc) // remove from complex
	{
		// just change isIn value to indicate this cell is dead!
		tokill->isIn = false;
		BT birth = tokill->birth;
		num dim = tokill->getDim();

		// we also must update size info... first the number of cells, which is easy
		num newtotsize = getSizeInfo(birth,dim).first - 1; //decrement #cells of this dim
		(*(sizeinfo[birth]))[dim].first = newtotsize;

		if (newtotsize == 0) (*(sizeinfo[birth]))[dim].second = 0;
		// and if we killed the max-indexed one while there are more cells left...
		else if(newtotsize > 0 && tokill->ind == getSizeInfo(birth,dim).second)
		{
			//cout<<"\n size Adjustment for: "<<*tokill;
			Cell<C,BT>* curcell;
			num newpos = 0;
			for (newpos = 1; newpos <= tokill->ind; ++newpos)
			{
				//cin.get();
				//cout<<"\n newpos "<<newpos;
				curcell = *(cpos - newpos);
				if (curcell == NULL)
				{
					continue; // ignore wiped
				}
				if (!(curcell -> isIn))
				{
					continue; // ignore deleted
				}
				break;
			}
			// finally, we know how much to shrink max index by!
			//cout<<"set sizeinfo... "; cin.get();
			(*(sizeinfo[birth]))[dim].second = tokill->ind - newpos;
		}
	}

	typename CSTRUCT::const_iterator i;
	// remove it from its co-boundaries' boundaries, if needed!
	if(fmbd)
	{
		for (i = tokill->cb.begin(); i != tokill->cb.end();++i)
		{
			i->first->bd.removeCell(tokill);
		}
	}
	// remove it from its boundaries' co-boundaries, if needed!
	if (fmcb)
	{
		for (i = tokill->bd.begin(); i != tokill->bd.end();++i)
		{
			i->first->cb.removeCell(tokill);
		}
	}
	return true;
}

// returns bottom dimension of calling frame!
template <typename C, typename BT>
pair<num,bool> Complex<C,BT>::getBotDim(const BT& born) const
{
    typename COMPLEX::const_iterator mit = clist.find(born);
	if (mit == clist.end())  // no such frame
	{
		return make_pair(0,false);
	}
	// frame located...
	FRAME* curf = mit->second;
	// we now want to find the lowest-dimensional cell list in this frame with nonzero size, so...
	// start at the lowest dim using a forward iterator,
	typename FRAME::iterator dimloc;
	for(dimloc = curf->begin(); dimloc != curf->end(); ++dimloc)
	{
		// and return the first dim you get with nonzero cells in it
		//cout<<"\n+++++++ list size of birth, dim: "<<born<<", "<<dimloc->first<<" is "<<lsize(born,dimloc->first); cin.get();
		if (lsize(born,dimloc->first) > 0) return make_pair(dimloc->first, true);
	}
	//if you are here, this frame has zero top dim
	return make_pair(0, false);
}

// returns top dimension of the calling frame!
template <typename C, typename BT>
pair<num,bool> Complex<C,BT>::getTopDim(const BT& born) const
{
	//cout<<"\n\n     returning topdim of frame: "<<born;
	typename COMPLEX::const_iterator mit = clist.find(born);
	if (mit == clist.end())  // no such frame
	{
		//cout<<"... not found! ";
		return make_pair(0,false);
	}
	// frame located...
	FRAME* curf = mit->second;

	// we now want to find the highest-dimensional cell list in this frame with nonzero size, so...
	// start at the highest dim using a reverse iterator,
	typename FRAME::reverse_iterator dimloc;
	for(dimloc = curf->rbegin(); dimloc != curf->rend(); ++dimloc)
	{
		// and return the first dim you get with nonzero cells in it
		//cout<<"\n+++++++ list size of birth, dim: "<<born<<", "<<dimloc->first<<" is "<<lsize(born,dimloc->first); cin.get();
		if (lsize(born,dimloc->first) > 0) return make_pair(dimloc->first, true);
	}
	//if you are here, this frame has zero top dim
	return make_pair(0, false);
}

// gets top dimension of the whole complex!
template <typename C, typename BT>
pair<num,bool> Complex<C,BT>::getTopDim() const
{
	//cout<<"\ndetermining topdim... ";
	typename COMPLEX::const_iterator mit;

	bool flag = false;
	num toret = 0;
	for(mit = clist.begin(); mit != clist.end(); ++mit)
	{
		flag = true;
		//cout<<"\n    comp "<<toret<<" vs "
		//	<<getTopDim(mit->first).first<<" at "<<mit->first
		//	<<" with flag: "<<getTopDim(mit->first).second;
		if (toret < getTopDim(mit->first).first)
		{
			toret = getTopDim(mit->first).first;
		}
	}
	return make_pair(toret,flag);
}

// gets bottom dimension of the whole complex!
template <typename C, typename BT>
pair<num,bool> Complex<C,BT>::getBotDim() const
{
	//cout<<"\ndetermining topdim... ";
	typename COMPLEX::const_iterator mit;

	bool flag = false;
	num toret = 0;
	for(mit = clist.begin(); mit != clist.end(); ++mit)
	{
		flag = true;
		//cout<<"\n    comp "<<toret<<" vs "
		//	<<getTopDim(mit->first).first<<" at "<<mit->first
		//	<<" with flag: "<<getTopDim(mit->first).second;
		if (toret < getBotDim(mit->first).first)
		{
			toret = getBotDim(mit->first).first;
		}
	}
	return make_pair(toret,flag);
}


// suspension by a point, dim indicates the "order" of suspension,
// for example: 1 makes edges from the "god point", two makes faces,
// and so on
/*template <typename C, typename BT>
bool Complex<C,BT>::suspend(const num dim, const BT gbirth = INITBT)
{
	if (dim < 0) return false;
	if (dim > topdim) dim = topdim; // make sure nothing is out of bounds
	// first create a zero cell.. the "god point" from which this complex will hang
	Cell<C,BT>* godpt = new Cell<C,BT>(0,0);
	godpt->birth = gbirth;

	if (dim >= 1) // otherwise there is nothing to connect!
	{
		typename COMPLEX::iterator cur; // iterates over frames

		// stores the to-insert structure, each 0-cell is mapped to
		// a list of its edge with the god-point
		map<Cell<C,BT>*,Cell<C,BT>*> connect;

		Cell<C,BT>* curedge;

		// now loop over zero cells of complex
		for(cur = begin(0); cur != end(0); ++cur)
		{
			// and to each zero cell, associate an edge from the
			// god point to that zero cell with correct boundary...
			curedge = new Cell<C,BT>(1,0);
			// add each zero cell from the complex to curedge's boundary
			curedge->addBDLink(*cur,1);
			(*cur)->addCBLink(curedge,1);
			// and also the god point goes to curedge's boundary
			curedge->addBDLink(godpt,-1);
			godpt->addCBLink(curedge,-1);
			// associate curedge to cur in the map
			connect[*cur] = curedge;
		}


		// at this point, we have the one skeleton of the suspension
		// from the god point built. Maybe we need higher faces, so
		// we build them!

		// iterate over boundary of current cell!
		typename CHAIN::const_iterator bditer;

		for (num i = 1; i < dim; i++)
		{
			for (cur = begin(i); cur != end(i); ++cur)
			{
				// want this to correspond to cur in the map.
				curedge = new Cell<C,BT>(i+1,0);
				// but we need to inherit its boundary from the
				// mapped images of cur's boundary! First, cur itself...
				curedge->addBDLink(*cur,1);
				(*cur)->addCBLink(curedge,1);

				// now iterate over cur's boundary
				for (bditer = (*cur)->bd.begin(); bditer != (*cur)->bd.end(); ++bditer)
				{
					// and add to curedge's boundary the images of cur's boundary's cells
					curedge->addBDLink(connect[bditer->first],bditer->second*-1);
					connect[bditer->first]->addCBLink(curedge,bditer->second*-1);
				}
				// finally, map cur to its connection with the god point
				connect[*cur] = curedge;
			}
		}
		// now use the map to insert all god point connected cells!
		typename map<Cell<C,BT>*,Cell<C,BT>*>::const_iterator miter;
		for (miter = connect.begin(); miter != connect.end(); ++miter)
		{
			insertCell(miter->second);
		}
	}
	insertCell(godpt);
	return true;
}*/

// clears complex without deleting cells!
template <typename C, typename BT>
void Complex<C,BT>::clear()
{
	Destroy(false);
}

template <typename C, typename BT>
void Complex<C,BT>::Destroy(bool killcells = true)
{
    //cout<<"\n complex destroyer... "; cin.get();
	typename COMPLEX::const_iterator cit; // complex iterator
	typename FRAME::const_iterator fit; // frame iterator
	typename CLIST::iterator lit; // list iterator

	FRAME* curf;
	CLIST* curl;

	for (cit = clist.begin(); cit != clist.end(); ++cit)  // loop over dimension...
	{
		curf = cit->second;
		for (fit = curf->begin(); fit != curf->end(); ++fit)
		{
			curl = fit->second;
			if (killcells)
			{
			    //cout<<"\n    killing cells... "; cin.get();
				for (lit = curl->begin(); lit != curl->end(); ++lit)
				{
					if (*lit != NULL)
					{
						delete *lit;
						//*lit = NULL;
					}
				}
			}
			curl->clear();
			CLIST(*curl).swap(*curl);
			delete curl;
		}
		curf->clear();
		delete curf;
	}
	clist.clear();

	// now remove size info
	typename DMAP::iterator dmit;
	for (dmit = sizeinfo.begin(); dmit != sizeinfo.end(); ++dmit)
	{
		dmit->second->clear(); // clear the range map
		delete dmit->second; // and delete the pointer!
	}
	sizeinfo.clear();
	//cout<<"\n killing done!!! "; cin.get();
}

// helper function to output the complex:
template <typename C, typename BT>
void Complex<C,BT>::printSizeStruct(ostream& out = cout) const
{
	typename DMAP::const_iterator dit;
	DSIZE* cursz;
	typename DSIZE::const_iterator sit;

	for (dit = sizeinfo.begin(); dit != sizeinfo.end(); ++dit)
	{
		out<<"\n Birth: "<< dit->first;
		cursz = dit->second;
		for (sit = cursz->begin(); sit != cursz->end(); ++sit)
		{
			out<<"\n     Dim: "<<sit->first;
			out<<" <"<<sit->second.first<<","<<sit->second.second<<">";
		}
	}
}



// helper function to output the complex:
template <typename C, typename BT>
void Complex<C,BT>::print (ostream& out = cout) const
{
	typename COMPLEX::const_iterator cit; // complex iterator

	out<<"\n\n ::::Complex:::: ["<<clist.size()<<"]\n";
	for (cit = clist.begin(); cit != clist.end(); ++cit)  // loop over dimension...
	{
		printFrame(cit->first,out);
		out<<"\n*****************************\n";
	}
}

template <typename C, typename BT>
void Complex<C,BT>::clearDim(num dim)
{
	typename COMPLEX::const_iterator cit; // complex iterator
	typename FRAME::const_iterator fit; // frame iterator

	typename DMAP::iterator dmit;
	DSIZE* cursz;
	for (cit = clist.begin(); cit != clist.end(); ++cit)  // loop over dimension...
	{
		fit = cit->second->find(dim);
		if(fit != cit->second->end()) //found!
		{
			fit->second->clear();
			// update dimension info
			dmit = sizeinfo.find(cit->first);
			cursz = dmit->second;
			(*cursz)[dim].first = 0;
			(*cursz)[dim].second = 0;
		}
	}
}


#endif /* COMPLEX_HPP_ */
