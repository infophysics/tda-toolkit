/*
 * Morse.hpp
 *
 * Basic methods for morse complexes
 */

#ifndef MORSE_HPP_
#define MORSE_HPP_


// copies over non-removed and non-null cells to
// "other" complex after clearing it out. if forred
// is true then the cells are re-ordered by the
// unit incidence count in the boundary, otherwise
// in the co-boundary.
template <typename C, typename BT>
void MComplex<C,BT>::moveOver_Sort(MComplex<C,BT>& other, bool forred)
{

    // set structures to store cells ordered by unit-incidence: boundary or
    // coboundary: if (forred) then boundary, else co-boundary.

    multiset<Cell<C,BT>*,unitbdord<C,BT> > bdsort;
    multiset<Cell<C,BT>*,unitcbord<C,BT> > cbsort;

    // and their iterators:
    typename multiset<Cell<C,BT>*, unitbdord<C,BT> >::iterator bditer; // for reduce
    typename multiset<Cell<C,BT>*, unitcbord<C,BT> >::iterator cbiter; // for coreduce

	other.Destroy(false); // remove all prior cells from other complex
	// iterate over complex...
	typename COMPLEX::const_iterator compit;
	FRAME* curf;
	typename FRAME::const_iterator frit;
	CLIST* curl;
	typename CLIST::const_iterator clit;
	Cell<C,BT>* curcell;
	// iterate over complex
	for (compit = this->clist.begin(); compit != this->clist.end(); ++compit)
	{
		// get frame
		curf = compit->second;

		// iterate over frame
		for (frit = curf->begin(); frit != curf->end(); ++frit)
		{
		    // extract list
			curl = frit->second;
			// iterate over list
			for (clit = curl->begin(); clit != curl->end(); ++clit)
			{
				// extract current cell:
				curcell = *clit;
				// if this cell is not null...
				if (curcell == NULL) continue;
				// and not removed
				if (!(curcell->isIn)) continue;
				// then add the cell to the new complex


				// depending on whether we are ordering for reduction or coreduction,
				// choose the appropriate set:

				if (forred) bdsort.insert(curcell);
				else cbsort.insert(curcell);
            }

            // insert from set into other complex:
            if (forred)
            {
                for (bditer = bdsort.begin(); bditer != bdsort.end(); ++bditer)
                {
                    other.insertCell(*bditer);
                    //other.printList(compit->first,frit->first);
                }
            }
            else
            {
                for (cbiter = cbsort.begin(); cbiter != cbsort.end(); ++cbiter)
                {
                    other.insertCell(*cbiter);
                    //other.printList(compit->first,frit->first);
                }
            }

            bdsort.clear();
            cbsort.clear();
			// clear list
			curl->clear();
			CLIST(*curl).swap(*curl);
			delete curl;
		}

		// clear frame
		curf->clear();
		delete curf;
	}

	// and of course... clear this complex!
	this->clist.clear();

}

// returns the boundary matrix with C coeffs for the subcomplex up to
// and including level "birth" (default: all) in the dimensions updim ---> updim-1.
// should call after moveOver() so that the cell-indices work out. Returns the
// density: number of nonzero entries / total number of entries
template <typename C, typename BT>
double MComplex<C,BT>::boundaryMatrix (const num& updim, vector<vector<C> >& bmat,
                                     BT birth = BANBT) const
{
    double nonzero_count = 0; // number of nonzero entries in bdry matrix
    num upsz = this->dsize(updim, birth); // number of n-cells
    num dnsz = this->dsize(updim-1, birth); // number of (n-1)-cells
    long double matsz = upsz * dnsz; // total number of entries in bdry matrix


    bmat.assign(upsz, vector<C>(dnsz,INITBT));

    typename COMPLEX::const_iterator compit;
	typename COMPLEX::const_iterator bcap;
    bcap = this->clist.find(birth);

	FRAME* curf;
	typename FRAME::const_iterator frit;
	CLIST* curl;
	typename CLIST::const_iterator clit;
	Cell<C,BT>* curcell;

	typename CSTRUCT::const_iterator bditer; // loops over critical coboundary of queen

	for (compit = this->clist.begin(); compit != bcap; ++compit)
    {
        curf = compit->second;
        if (curf == NULL) continue;

        frit = curf->find(updim); // find list for this dimension
        if (frit == curf->end()) continue;

        curl = frit->second;
        if (curl == NULL) continue;

        // now go over the list...
        for (clit = curl->begin(); clit != curl->end(); ++clit)
        {
            curcell = *clit; // get current cell...
            nonzero_count += curcell->getBD().size();
            // iterating over its boundary,...
            for (bditer = curcell->getBD().begin(); bditer != curcell->getBD().end(); ++bditer)
            {
                // populate the coefficient!
                bmat.at(curcell->getInd()).at(bditer->first->getInd()) = bditer->second;
            }
        }
    }
    return (matsz == 0) ? 1 : (nonzero_count / matsz);
}


// returns the density of the boundary matrix with C coeffs for the subcomplex up to
// and including level "birth" (default: all) in the dimensions updim ---> updim-1.

template <typename C, typename BT>
double MComplex<C,BT>::getBoundaryDensity (const num& updim, BT birth = BANBT) const
{
    double nonzero_count = 0; // number of nonzero entries in bdry matrix
    num upsz = this->dsize(updim, birth); // number of n-cells
    num dnsz = this->dsize(updim-1, birth); // number of (n-1)-cells
    long double matsz = upsz * dnsz; // total number of entries in bdry matrix

    if (matsz == 0) return 1;

    // iterators for the complex...
    typename COMPLEX::const_iterator compit;
	typename COMPLEX::const_iterator bcap;
    bcap = this->clist.find(birth);

	FRAME* curf;
	typename FRAME::const_iterator frit;
	CLIST* curl;
	typename CLIST::const_iterator clit;
	Cell<C,BT>* curcell;

	typename CSTRUCT::const_iterator bditer; // loops over critical coboundary of queen

    // start iteration
	for (compit = this->clist.begin(); compit != bcap; ++compit)
    {
        curf = compit->second;
        if (curf == NULL) continue;

        frit = curf->find(updim); // find list for this dimension
        if (frit == curf->end()) continue;

        curl = frit->second;
        if (curl == NULL) continue;

        // now add up nontrivial boundary components!
        for (clit = curl->begin(); clit != curl->end(); ++clit)
        {
            curcell = *clit; // get current cell...
            nonzero_count += curcell->getBD().size();
        }
    }
    return (nonzero_count / matsz);
}

// returns the fraction of nonzero entries in the boundary matrix with C coeffs for the subcomplex up to
// and including level "birth" (default: all) in the dimensions updim ---> updim-1.

template <typename C, typename BT>
double MComplex<C,BT>::getUnitFraction (const num& updim, BT birth = BANBT) const
{
    double nonzero_count = 0; // number of nonzero entries in bdry matrix
    double unit_count = 0;
    num upsz = this->dsize(updim, birth); // number of n-cells
    num dnsz = this->dsize(updim-1, birth); // number of (n-1)-cells
    long double matsz = upsz * dnsz; // total number of entries in bdry matrix

    if (matsz == 0) return 1;

    // iterators for the complex...
    typename COMPLEX::const_iterator compit;
	typename COMPLEX::const_iterator bcap;
    bcap = this->clist.find(birth);

	FRAME* curf;
	typename FRAME::const_iterator frit;
	CLIST* curl;
	typename CLIST::const_iterator clit;
	Cell<C,BT>* curcell;

	typename CSTRUCT::const_iterator bditer; // loops over critical coboundary of queen

    // start iteration
	for (compit = this->clist.begin(); compit != bcap; ++compit)
    {
        curf = compit->second;
        if (curf == NULL) continue;

        frit = curf->find(updim); // find list for this dimension
        if (frit == curf->end()) continue;

        curl = frit->second;
        if (curl == NULL) continue;

        // now add up nontrivial boundary components!
        for (clit = curl->begin(); clit != curl->end(); ++clit)
        {
            curcell = *clit; // get current cell...
            nonzero_count += curcell->getBD().size();
            unit_count += curcell->getBD().unitCount();
        }
    }
    return (nonzero_count == 0) ? 1: (unit_count / nonzero_count);
}



// copies over non-removed and non-null cells to
// "other" complex in order after clearing it out. if forred
// is true then the cells are re-ordered by the
// unit incidence count in the boundary, otherwise
// in the co-boundary.
template <typename C, typename BT>
void MComplex<C,BT>::moveOver(MComplex<C,BT>& other)
{

	other.Destroy(false); // remove all prior cells from other complex
	// iterate over complex...
	typename COMPLEX::const_iterator compit;
	FRAME* curf;
	typename FRAME::const_iterator frit;
	CLIST* curl;
	typename CLIST::const_iterator clit;
	Cell<C,BT>* curcell;
	// iterate over complex
	for (compit = this->clist.begin(); compit != this->clist.end(); ++compit)
	{
		// get frame
		curf = compit->second;

		// iterate over frame
		for (frit = curf->begin(); frit != curf->end(); ++frit)
		{
		    // extract list
			curl = frit->second;
			// iterate over list
			for (clit = curl->begin(); clit != curl->end(); ++clit)
			{
				// extract current cell:
				curcell = *clit;
				// if this cell is not null...
				if (curcell == NULL) continue;
				// and not removed
				if (!(curcell->isIn)) continue;
				// then add the cell to the new complex
                other.insertCell(curcell);
            }
			// clear list
			curl->clear();
			CLIST(*curl).swap(*curl);
			delete curl;
		}

		// clear frame
		curf->clear();
		delete curf;
	}

	// and of course... clear this complex!
	this->clist.clear();

}

// copies over non-removed and non-null cells to
// "other" complex in random order after clearing it out. if forred
// is true then the cells are re-ordered by the
// unit incidence count in the boundary, otherwise
// in the co-boundary.
template <typename C, typename BT>
void MComplex<C,BT>::moveOver_Shuffle(MComplex<C,BT>& other)
{

	other.Destroy(false); // remove all prior cells from other complex
	// iterate over complex...
	typename COMPLEX::const_iterator compit;
	FRAME* curf;
	typename FRAME::const_iterator frit;
	CLIST* curl;
	typename CLIST::const_iterator clit;
	Cell<C,BT>* curcell;

	srand(time(NULL)); // for random shuffle!

	// iterate over complex
	for (compit = this->clist.begin(); compit != this->clist.end(); ++compit)
	{
		// get frame
		curf = compit->second;

		// iterate over frame
		for (frit = curf->begin(); frit != curf->end(); ++frit)
		{
		    // extract list
			curl = frit->second;
			random_shuffle(curl->begin(),curl->end());
			// iterate over list
			for (clit = curl->begin(); clit != curl->end(); ++clit)
			{
				// extract current cell:
				curcell = *clit;
				// if this cell is not null...
				if (curcell == NULL) continue;
				// and not removed
				if (!(curcell->isIn)) continue;
				// then add the cell to the new complex
                other.insertCell(curcell);
            }
			// clear list
			curl->clear();
			CLIST(*curl).swap(*curl);
			delete curl;
		}

		// clear frame
		curf->clear();
		delete curf;
	}

	// and of course... clear this complex!
	this->clist.clear();

}

template <typename C, typename BT>
void MComplex<C,BT>::showGens (const BT& born, const num& dim, bool resbirth = false, ostream& out = cout) const
{
	out<<"\n\n+++++++++ GENERATORS FOR BIRTH: ["<<born<<"] DIMENSION <"<<dim<<"> ++++++++++++";

	typename COMPLEX::const_iterator findbirth = this->clist.find(born);

	num mpsz = 0;

	if (findbirth != this->clist.end())
	{
		FRAME* bframe = findbirth->second; // frame for this birth time
		typename FRAME::const_iterator finddim = bframe->find(dim); // locates cells of this dimension

		if(finddim != bframe->end())
		{
			CLIST* dlist = finddim->second;
			// okay, now loop over this list and get gen structures output
			typename CLIST::const_iterator curcell;
			for (curcell = dlist->begin(); curcell != dlist->end(); ++curcell)
			{
				if (*curcell == NULL) continue; // ignore nulls
				if (!(*curcell)->isIn) continue; // and removeds

				if ((*curcell)->gen != NULL) mpsz += (*curcell)->gen->size();

				out<<"\n++Root: "<<**curcell<<" ";
				(*curcell)->showGens(out, resbirth);
			}
		}
	}
	//cout<<"\n mpsz: "<<mpsz;
}


template <typename C, typename BT>
bool MComplex<C,BT>::insertCell(Cell<C,BT>* toins)
{
	// try insertion into complex using base class's insert
	bool toret = this->Complex<C,BT>::insertCell(toins);
	BT birth = toins->birth;

	if (toret) // if we successfully inserted,...
	{
		if (!isred)
		{
            // at the minimal uncritical dimension or less
            if (toins->getDim() <= getMinDim(birth).first)
            {
                // update min index if the new index is smaller
                if (toins->ind < minindex) minindex = toins->ind;
            }
		}
		else
		{
            // at the minimal uncritical dimension or less
            if (toins->getDim() >= getMaxDim(birth).first)
            {
                // update min index if the new index is smaller
                if (toins->ind < minindex) minindex = toins->ind;
            }
		}
	}
	return toret;
}

template <typename C, typename BT>
void MComplex<C,BT>::quicksert(const vector<Cell<C,BT>*>& toins)
{
	if (toins.size() > 0)
	{
		BT birth = toins.at(0)->birth;
		num dim = toins.at(0)->getDim();
		num index = toins.at(0)->ind;

		//cout<<"\n ins:...";
		this->Complex<C,BT>::quicksert(birth, dim, toins);
		if (!isred && dim <= getMinDim(birth).first)
		{
			if (minindex > index) minindex = index;
		}
		else if (isred && dim >= getMaxDim(birth).first)
		{
			if (minindex > index) minindex = index;
		}
	}
}

// returns the MAXIMAL dimension containing UNCRITICAL cells. For reduction!
template <typename C, typename BT>
pair<num,bool> MComplex<C,BT>::getMaxDim(const BT& btime)
{
    // plan: reverse_iterate over dimensions, skip those with size 0, and
	// consult critinfo structure to see when the (birth,dim) pair
	// has nonzero dimension... and then check if the last cell of
	// the current dim is uncritical, etc...
	if (this->sizeinfo.find(btime) == this->sizeinfo.end()) return make_pair(0,false);
	// no such birth time located!
	DSIZE* cursz = this->sizeinfo[btime];
    typename DSIZE::reverse_iterator dimit;
    // reverse-iterate over sizeinfo structure
    for (dimit = cursz->rbegin(); dimit != cursz->rend(); ++dimit)
    {
        if (this->lsize(btime,dimit->first)==0) continue; // ignore dimensions with size 0!
        // check for uncritical cells...
        if (hasUncritCells(btime,dimit->first)) // if there are uncrit cells in this dimension,...
        {
			return make_pair(dimit->first, true); // return this dimension + true!
		}
	}
	// if you reached here, the complex is exhausted, the max dim is the bot dim
	return make_pair(this->getBotDim(btime).first,false);

}

// returns the MINIMAL dimension containing UNCRITICAL cells. For coreduction!
template <typename C, typename BT>
pair<num,bool> MComplex<C,BT>::getMinDim(const BT& btime)
{
	// plan: iterate over dimensions, skip those with size 0, and
	// consult critinfo structure to see when the (birth,dim) pair
	// has nonzero dimension... and then check if the last cell of
	// the current dim is uncritical, etc...
	if (this->sizeinfo.find(btime) == this->sizeinfo.end()) return make_pair(0,false);
	// no such birth time located!
	DSIZE* cursz = this->sizeinfo[btime];

	typename DSIZE::iterator dimit;
	// loop over the sizeinfo structure
	for (dimit = cursz->begin(); dimit != cursz->end(); ++dimit)
	{
		if (this->lsize(btime,dimit->first)==0) continue; // ignore dimensions with size 0!

		// check if the critinfo index is less than the max index
		// the first term is the "next uncritical cell", the second is the "largest index" cell
		// still present in the complex.

		//cout<<" calling uncrit... for birth "<<btime<<" and dim "<<dimit->first; cin.get();
		if (hasUncritCells(btime,dimit->first)) // if there are uncrit cells in this dimension,...
		{
			return make_pair(dimit->first, true); // return this dimension + true!
		}
	}
	// if you reached here, the complex is exhausted, the min dim is the top dim
	return make_pair(this->getTopDim(btime).first,false);
}

// prunes list with < 50% cell occupancy
template <typename C, typename BT>
void MComplex<C,BT>::memorySaver(BT birth = INITBT)
{
	typename COMPLEX::const_iterator compit;
	FRAME* curf;
	typename FRAME::const_iterator frit;
	CLIST* curl;
	typename CLIST::const_iterator clit;

	BT curbirth;
	num curdim;
	bool samebirth;

	vector<pair<BT, num> > prunelists;

	for (compit = this->clist.begin(); compit != this->clist.end(); ++compit)
	{
		samebirth = false;
		curbirth = compit->first;

		if (curbirth == birth)
		{
			samebirth = true;
		}

		curf = compit->second; // get frame
		for (frit = curf->begin(); frit != curf->end(); ++frit)
		{
			curdim = frit->first;
			curl = frit->second; // get list
			// now... compare list raw size to lsize and if it is
			// more than twice as much then add it to the pruned
			// vector candidates
			if (2*this->lsize(curbirth,curdim) <= (num)curl->size() )
			{
				if (!isred && samebirth && getMinDim(curbirth).first == curdim)
				{
					// if we are pruning the min crit dimension for the current birth time,...
					minindex = 0;
				}
				else if (isred && samebirth && getMaxDim(curbirth).first == curdim)
				{
					// if we are pruning the min crit dimension for the current birth time,...
					minindex = 0;
				}

				prunelists.push_back(make_pair(curbirth,curdim));
			}
		}
	}

	// now loop over prunelists and prune them!
	typename vector<pair<BT,num> >::iterator pliter;
	for(pliter = prunelists.begin(); pliter != prunelists.end(); ++pliter)
	{
		this->pruneList(pliter->first,pliter->second);
	}
}


// removes a cell from the complex, returns an iterator
template <typename C, typename BT>
bool MComplex<C,BT>::removeCell(Cell<C,BT>* tokill, typename CLIST::iterator& pos, bool asCrit=false)
{
	if (tokill == NULL) return false;
	bool toret;
	// extract birth and dimension of cell to kill
	//BT birth = tokill->birth;
	//num dim = tokill->getDim();

	// extract old min dimension...
	//pair<num,bool> oldmd = getMinDim(birth);

	if (asCrit) // if removed as critical...
	{
		toret = this->Complex<C,BT>::removeCell(tokill, pos, false, true, true);
	}
	else // ordinary removal, not critical
	{
		toret = this->Complex<C,BT>::removeCell(tokill, pos, true, true, true);
	}
	//setMinIndex(tokill->birth);
	return toret; // no such cell was removed
}

// recomputes minindex after cell removal.
template <typename C, typename BT>
bool MComplex<C,BT>::setMinIndex(const BT birth)
{
	// extract new minimum dimension
	pair<num,bool> mdret = isred ? getMaxDim(birth) : getMinDim(birth);
	num newdim = mdret.first;


	//cout<<"\n ----- reset min index... BIRTH = "<<birth<<" and DIM = "<<newdim; cin.get();

	bool tracer = false;


	// set traces here:
	//if (birth == 1 && newdim == 2) tracer = true;

	if (mdret.second == false) // bad mindim?!!??!!
	{
		cout<<"no mindim for birth "<<birth<<", try reset!";
		cin.get();
		minindex = 0;
		return false;
	}

	if (tracer)
	{
		this->printSizeInfo(true);
		this->printSizeStruct();
	}

	// store current list of cells:
	CLIST* curl;
	// and iterate over it
	typename CLIST::iterator liter;
	// pointer to the cells encountered while iterating
	Cell<C,BT>* curcell = NULL;
	// extract sizes for this dimension
	pair<num,num> sz = getSizeInfo(birth,newdim);
	// and the list in question
	curl = (*this->clist[birth])[newdim];
	// is the final element critical?
	if (curl->at(sz.second)->isCrit)
	{
		cout<<"\n min dim "<<newdim<<" has final cell at ind "<<sz.second<<" critical!!?";
		cin.get();
		minindex = 0;
		return false;
	}
	// ok, final element is NOT critical, so...
	// if the last cell of this new dimension is more than halfway into its size,
	// we iterate backwards, otherwise forwards.
	bool iterback = (2*sz.second < (num)curl->size()) ? true : false;
	if (iterback)
	{
		if (tracer) cout << "\nback iteration through list..., starting " << sz.second;
		liter = curl->begin() + sz.second;
		for (; liter != curl->begin(); --liter)
		{
			curcell = *liter;
			// pass over nulls and removeds
			if (curcell == NULL) continue;
			if (!(curcell->isIn)) continue;
			// found critical!!
			if (curcell->isCrit) break;
		}
		minindex = curcell->ind + 1; // best guess for next uncritical cell!
		return true;
	}
	else // iterate forwards
	{
		if (tracer) cout <<"\nfront iteration...";
		if (minindex > sz.second && tracer)
		{
			cout<<"\n birth: "<<birth<<" and dimension "<<newdim<<" has index "<<minindex<<" and finalind "<<sz.second;
			minindex = 0;
			cin.get();
		}
		for (liter = curl->begin() + minindex; liter <= curl->begin() + sz.second; ++liter)
		{
			curcell = *liter;
			if (curcell == NULL) {minindex++; continue;} // move over nulls
			if (!(curcell->isIn)) {minindex++; continue;} // and removeds
			if (curcell->isCrit) // and already critical cells
			{
				//cout<<"\n pre-existing crit cells of birth "<<birth<<" in new min dim "<<newdim<<" at index "<<curcell->ind;
				minindex++;
				//cin.get();
				//minindex = 0;
				//return false;
				continue;
			}
			// if you reach here, this is the first non critical, non-removed, non-null
			// cell in this list. congratulations, set minindex and get the hell out.
			minindex = curcell->ind;
			return true;
		}
		// reached here??? somehow there were no candidates for min index cells
		// even though mindim suggested that there are
		//cout<<"\n min dim "<<newdim<<" lied at birth "<<birth;
	}

	// we really should never reach here!
	return false;
}

// checks if there are any uncritical cells remaining in this frame
template <typename C, typename BT>
bool MComplex<C,BT>::hasUncritCells(const BT& birth) const
{
    // if there are uncrit cells either in the top or bottom available dim,
    // then we have uncrit cells.
	return isred ? hasUncritCells(birth, this->getBotDim(birth).first) :
                   hasUncritCells(birth, this->getTopDim(birth).first);
}

// checks if there are uncritical cells in this frame at the given dim
template <typename C, typename BT>
bool MComplex<C,BT>::hasUncritCells(const BT& birth, const num& dim) const
{
	if (this->lsize(birth,dim) == 0) return false; // no cells of this birth and dim!

	// if you reached here then the frame exists and cells of the
	// required dim exist. so now,... a simple comparison of sizeinfo
	// and critinfo should do the trick!


	// we are now checking the "largest index cell still present" for this birth
	// and dimension to see if it is critical...


	typename COMPLEX::const_iterator findb = this->clist.find(birth);
	if (findb == this->clist.end()) return false; // no cells at all for this birth, so no uncrit cells
	FRAME* curf = findb->second;
	typename FRAME::const_iterator findd = curf->find(dim);
	if (findd == curf->end()) return false; // no cells of this dim, so no uncrit cells

	CLIST* curlist = findd->second; // got the list for this dimension!
	if (curlist->size() == 0) return false; // no cells, so no uncrit cells...

	// extract index of last cell in this list!
	num highestIndex = this->getSizeInfo(birth,dim).second;
	//cout<<" using highest index "<<highestIndex; cin.get();

	if (curlist->at(highestIndex)->isCrit) return false; // all cells critical!
	// reached here? the last cell is uncritical!
	return true;
}


#endif /* MORSE_HPP_ */

