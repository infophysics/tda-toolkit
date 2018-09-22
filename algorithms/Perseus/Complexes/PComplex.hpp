/*
 * PComplex.hpp
 *
 *  Created on: Jan 10, 2011
 *      Author: Vidit
 */

#ifndef PCOMPLEX_HPP_
#define PCOMPLEX_HPP_

# include "PComplex.h"

template <typename C, typename BT>
num PComplex<C,BT>::makeFromComplex(const Complex<C,BT>& topers, bool truncate = false)
{
	Cell<C,BT>* curcell; // cell from complex topers
	PCell<C,BT>* percell; // persistent cell, to be made

	PCMAP corr; // correspondence between cell and persistent cell


	// loop over the complex: this requires iterators
	typename COMPLEX::const_iterator cit;
	FRAME* curf;
	typename FRAME::const_iterator fit;
	CLIST* curl;
	typename CLIST::const_iterator lit;

    num topdim = topers.getTopDim().first;
    if (topdim == 0) truncate = false;

	// loop over all frames
	for (cit = topers.clist.begin(); cit != topers.clist.end(); ++cit)
	{
		curf = cit->second; // extract current frame
		// loop over current frame
		for (fit = curf->begin(); fit != curf->end(); ++fit)
		{
			curl = fit->second; //extract current list
			// finally, loop over list
			for (lit = curl->begin(); lit != curl->end(); ++lit)
			{
				curcell = *lit; // extract current cell
				if (curcell == NULL) continue; // check nullity
				if (!(curcell->isIn))
				{
				    cout<<" \n wtf, cell missing??? : "<<*curcell;
				    cin.get();
				    continue; // and removedness
				}
				// check truncation! if cell has max dim...
				if (truncate && curcell->getDim() == topdim)
				{
				    // and empty boundary...
				    if (curcell->getBD().size()==0)
				    continue; // then ignore it.
				}

				// okay, use cell to make percell
				percell = new PCell<C,BT>;
				percell->makeFromCell(curcell);
				corr[curcell] = percell;
				klist.push_back(percell);
				percell->kindex = klist.size()-1;
				//cout<<"\n made "<<*percell<<" from "<<*curcell;
			}
		}
	}
	PCVEC(klist).swap(klist);
	makeBDChains(corr);
	return klist.size();
}

// destroy persistent cells
template <typename C, typename BT>
void PComplex<C,BT>::destroyCells()
{
	typename PCVEC::iterator pciter;
	for (pciter = klist.begin(); pciter != klist.end(); ++pciter)
	{
		delete *pciter;
		*pciter = NULL;
	}
	// wipe out list...
	klist.clear();
	PCVEC(klist).swap(klist);
}


// use the correspondence map to build boundaries of all the percells
template <typename C, typename BT>
bool PComplex<C,BT>::makeBDChains(PCMAP& corr, bool makez2 = true)
{
	// first, declare the current cell
	Cell<C,BT>* curcell;
	// and its boundary candidate
	Cell<C,BT>* curcellBD;
	// similarly the persistent cell
	PCell<C,BT>* percell;
	//and boundary
	PCell<C,BT>* percellBD;
	C bdcoeff; // incidence coefficient

	// boundary iterator!
	typename CSTRUCT::const_iterator bditer;

	// to loop over correspondence map:
	typename PCMAP::const_iterator corrit;
	for (corrit = corr.begin(); corrit != corr.end(); ++corrit)
	{
		// current cell:
		curcell = corrit->first;
		// current pcell:
		percell = corrit->second;
		// iterate over current cell's boundary!
		for (bditer = curcell->getBD().begin(); bditer != curcell->getBD().end(); ++bditer)
		{
			curcellBD = bditer->first;
			bdcoeff = bditer->second;
			if (makez2) bdcoeff %= 2;
			percellBD = corr[curcellBD];
			percell->addBdryLink(percellBD,bdcoeff);
		}
		//cout<<"\n made bdry "<<percell->getBD()<<" from "<<curcell->getBD();
	}
	//cout<<" done all!!! "; cin.get();
	return true;
}

template <typename C, typename BT>
void PComplex<C,BT>::initPersData(const Complex<C,BT>& other)
{
	// first we initialize our interval and betti number structures:
	ints.clear();
	betti.clear();	// betti numbers, for each birth-stage

	// initialize persistence intervals
	for (num dim = 0; dim <= other.getTopDim().first; dim++)
	{
		//cout<<"ints vec populated to: "<<this->getTopDim().first;
		ints[dim] = new INTVEC;
	}

	// initialize betti number junk
	typename COMPLEX::const_iterator i;
	for(i=other.clist.begin(); i != other.clist.end() ; ++i)
	{
		// for this frame, create betti numbers initialized to 0
		//cout<<"\nmade betti ptr for "<<i->first;
		betti[i->first] = new vector<BT>;
		betti[i->first]->insert(betti[i->first]->begin(),other.getTopDim().first + 1, 0);   // fill with zeros for now
	}
}


template <typename C, typename BT>
void PComplex<C,BT>::destroyPersData()
{
	// destroy persistence intervals
	typename INT_STR::iterator iit;
	for (iit = ints.begin(); iit != ints.end(); ++iit)
	{
		iit->second->clear();
		delete iit->second;
	}

	// destroy betti number junk
	typename BETTI_STR::iterator i;
	for(i=betti.begin(); i != betti.end() ; ++i)
	{
		i->second->clear();
		delete i->second;
	}
	ints.clear();
	betti.clear();
}

template <typename C, typename BT>
void PComplex<C,BT>::showInts(ostream& out = cout)
{
	//out<<"\n\n:::::INTERVALS::::::";
	typename INT_STR::const_iterator it;
	INTVEC* curints;
	typename INTVEC::const_iterator iit;
	for(it = ints.begin(); it != ints.end(); ++it)
	{
		out<<"\nDIM: ["<<it->first<<"]\n";
		curints = it->second;
		for (iit = curints->begin(); iit!=curints->end(); ++iit)
		{
			out<<"["<<iit->first<<","<<iit->second<<"]\n";
		}
	}
}

template <typename C, typename BT>
void PComplex<C,BT>::showBetti(ostream& out = cout)
{
	//out<<"\n\n:::::BETTI NUMBERS BY FRAME:::::\n\n"; //cin.get();
	typename BETTI_STR::const_iterator bit;
	typename vector<BT>::const_iterator vit;
	vector<BT>* curvec;

	for (bit = betti.begin(); bit != betti.end(); ++bit)
	{
		curvec = bit->second;
		// output birth time of frame:
		out<<"\n "<<bit->first<<" ";
		for(vit = curvec->begin(); vit != curvec->end(); ++vit)
		{
			out<<*vit<<" ";
		}
	}
}

template <typename C, typename BT>
map<num,vector<pair<BT,BT> > > PComplex<C,BT>::getInts() const
{
    // return interval structure
    map<num,vector<pair<BT,BT> > > toret;
    typename INT_STR::const_iterator it;
	INTVEC* curints;
	typename INTVEC::const_iterator iit;

    vector<pair<BT, BT> > newints; // intervals for current dimension...
    for (it = ints.begin(); it != ints.end(); ++it)
    {
        newints.clear();
        curints = it->second;
        toret.insert(make_pair(it->first,newints));
        for (iit = curints->begin(); iit != curints->end(); ++iit)
        {
            toret[it->first].push_back(*iit);
        }
    }
    return toret;
}

template <typename C, typename BT>
bool PComplex<C,BT>::makeOutputFiles(const string& filename)
{
	std::ostringstream dimst; // string representation for dimension
	// first we make interval files for each dimension
	typename INT_STR::const_iterator it;
	INTVEC* curints;
	typename INTVEC::const_iterator iit;
	string fname; // file name for writing
	for(it = ints.begin(); it != ints.end(); ++it)
	{
		// create file name by appending dimension
		dimst.str("");
		dimst<<it->first;
		fname = filename + "_" + dimst.str() + ".txt";
		ofstream outf(fname.c_str());

		if (!outf.good())
		{
			cout<<"\nUnable to create output file "<<fname;
			return false;
		}
		curints = it->second;
		for (iit = curints->begin(); iit!=curints->end(); ++iit)
		{
			//outf<<"["<<iit->first<<","<<iit->second<<"]\n";
			outf<<iit->first<<" "<<iit->second<<"\n";
		}
		outf.close();
	}

	// now make betti number file!
	fname = filename + "_betti.txt";
	ofstream bfile(fname.c_str());
	// if the file doesn't get made, exit gracefully...
	if (!bfile.good())
	{
		cout<<"\nUnable to create output file "<<fname;
		return false;
	}
	// otherwise, write betti numbers to file
	showBetti(bfile);
	bfile.close();
	return true;
}

#endif /* PCOMPLEX_HPP_ */
