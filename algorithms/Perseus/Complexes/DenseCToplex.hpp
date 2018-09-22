/*
 * DenseCToplex.hpp
 *
 *  Created on: Dec 23, 2010
 *      Author: Vidit
 */

#ifndef DENSECTOPLEX_HPP_
#define DENSECTOPLEX_HPP_

#include "DenseCToplex.h"

// makes face links and inherits birth times for lower cubes once top cubes are populated!
template <typename C, typename BT>
bool DenseCToplex<C,BT>::makeFaces(bool incbirth = false)
{
    if (faced) return false; // don't make faces again, please.

	// iterate backwards through cube structure, adding boundary links to lower addins
	CVEC* curvec; // current cube vector corresponding to this addin
	//num topdim = cube.size();
	vector<num> curaddin, curanch, curext;
	num cursz; // size of current cube vector
	num i; num csz = cube.size();
	// backward iteration
	bool istop;

	typename CGRID::reverse_iterator rit;
	for (rit = cube.rbegin(); rit != cube.rend(); ++rit)
	{
	    istop = (rit == cube.rbegin()) ? true : false;
		i = csz - (rit - cube.rbegin()) - 1;
		// get the addin vector for this position
		getAddin(i,curaddin);
		if (DCFACETALK) {cout<<"\n***** Facing: "; c_print<vector<num> >(curaddin);}
		// extract extent vector for this addin
		getExtForAddin(curaddin,curext);
        //cout<<"\n     gots addin"; cin.get();
		curvec = cube.at(i); // extract current vector of cubes
		if (curvec == NULL) {continue;} // and ignore if NULL
		// iterate over current vector of cubes!
		cursz = curvec->size();
		for (num j = 0; j < cursz; ++j)
		{
		    //cout<<"\ntesting position ["<<i<<","<<j<<"]";
			if (cube.at(i)->at(j) == NULL) {continue;} // exclude NULL cubes
			if (istop && incbirth) cube.at(i)->at(j)->birth += (INITBT + 1);

			// get anchor for this cube
			getAnchorFast(j, curext, curanch);
			//cout<<"\n*************************anch: "; c_print<vector<num> >(curanch);
			// make face connections for this cell!
			 makeFaceLinks(cube.at(i)->at(j), curaddin, curanch);
		}
	}
	faced = true;
	return faced;
}

template <typename C, typename BT>
num DenseCToplex<C,BT>::getAddinPosFast(const vector<num>& advec) const
{
	// first we add all C(n,k) for all k < size(anch)...
	num toret = 0;
	for (num i=0; i<(num)advec.size();++i)
	{
		toret += nChoosek(extent.size(),i);
	}
	return toret + lexicoPosRev(extent.size(), advec);
}

// writes cubical toplex to morse cell complex
template <typename C, typename BT>
void DenseCToplex<C,BT>::writeComplex(MComplex<C,BT>& mycomp, bool killself = true)
{
	CVEC* curvec;
	for (num i = 0; i < (num)cube.size(); ++i)
	{
		curvec = cube.at(i);
		for (num j = 0; j < (num) curvec->size(); ++j)
		{
			if (curvec->at(j) == NULL) continue; // ignore null cubes
			mycomp.insertCell(curvec->at(j));
		}
		if (killself) curvec->clear(); // remove cell vector from toplex to save memory
	}
	cube.clear();
}



template <typename C, typename BT>
void DenseCToplex<C,BT>::makeFaceLinks(Cell<C,BT>* curcb, const vector<num>& addin, const vector<num>& anch)
{
	vector<num> fcaddin, fcanch; // face addin and anchor
	num cubeDim = addin.size();
	Cell<C,BT>* face;
	C coeff; // face incidence
	num fcaddpos, fcancpos;

	if (DCFACETALK)
	{
		cout<<"\n\nFACING "<<*curcb<<" with addin ";
		c_print<vector<num> >(addin);
		cout<<" and anchor ";
		c_print<vector<num> >(anch);
		cin.get();
	}

	for (num cdim = 0; cdim < 2*cubeDim; cdim+=2)
	{
		coeff = ((cdim/2) % 2 == 0) ? 1:-1;

		// ADDIN VECTOR:
		fcaddin = addin;
		// and remove cdim/2 th entry
		fcaddin.erase(fcaddin.begin()+(cdim/2));
		fcaddpos = getAddinPosFast(fcaddin);

		// FIRST FACE
		fcanch = anch;
		fcancpos = getAnchPos(fcanch,fcaddin);
		if (DCFACETALK)
		{
			cout<<"\n>>>>> lower face, coeff "<<-coeff<<" has addin ";
			c_print<vector<num> >(fcaddin);
			cout<<" and anchor ";
			c_print<vector<num> >(fcanch);
			cout<<" at pos ["<<fcaddpos<<","<<fcancpos<<"]";
			cin.get();
		}
		face = cube.at(fcaddpos)->at(fcancpos);

		// create face if null and add boundary link
		if (face == NULL)
		{
			face = new Cell<C,BT>(cubeDim-1);
			cube.at(fcaddpos)->at(fcancpos) = face;
		}

		curcb->addBDLink(face,-coeff);

		// SECOND FACE
		fcanch = anch;
		fcanch.at(addin.at(cdim/2)) += ES;
		fcancpos = getAnchPos(fcanch,fcaddin);

		if (DCFACETALK)
		{
			cout<<"\n>>>>> upper face, coeff "<<coeff<<" has addin ";
			c_print<vector<num> >(fcaddin);
			cout<<" and anchor ";
			c_print<vector<num> >(fcanch);
			cout<<" at pos ["<<fcaddpos<<","<<fcancpos<<"]";
			cin.get();
		}
		face = cube.at(fcaddpos)->at(fcancpos);
		// create face if null and add boundary link
		if (face == NULL)
		{
			face = new Cell<C,BT>(cubeDim-1);
			cube.at(fcaddpos)->at(fcancpos) = face;
		}

		curcb->addBDLink(face,coeff);
	}
	//if (curcb->getBD().size() != 2*cubeDim) cout<<"\ncube: "<<*curcb<<" and bd "<<curcb->getBD();
}

template <typename C, typename BT>
void DenseCToplex<C,BT>::Destroy()
{
	num addnum = cube.size();
	//cout<<"Destroying with size: "<<addnum; cin.get();
	for (num i = 0; i < (num)addnum; ++i)
	{
		// deallocate vector of cubes corresponding to this addin vector!
		cube.at(i)->clear();
		CVEC(*(cube.at(i))).swap(*(cube.at(i)));
		delete cube.at(i);
	}
	cube.clear();
	CGRID(cube).swap(cube);
}

template <typename C, typename BT>
void DenseCToplex<C,BT>::Init(const vector<num>& ext)
{
	num topdim = ext.size();

	for (num i=0; i<topdim; i++)
	{
		extent.push_back(ext.at(i)+1);
	}


	// first make as many inner vectors as the number of addins... \sum_{j=0}^n C(n,j) = 2^n
	double addnum = pow((float)2,(float)topdim);

	if (DCTTALK)
	{
		cout<<"\n\n *** Init DCT with ext vec "; c_print<vector<num> >(ext);
		cout<<"\n *** Extent: "; c_print<vector<num> >(extent); cout<<'\n';
		cout<<"\n *** adnum: " << addnum;
	}

	vector<num> curaddin, curext;
	CVEC* curvec;

	for (num i = 0; i < (num)addnum; ++i)
	{
		// extract addin vector for this dimension
		if(!getAddin(i,curaddin)) continue;
		// determine extents vector for this addin vector
		getExtForAddin(curaddin,curext);

		// now allocate appropriate sized cell vector, populate with NULLS
		curvec = new CVEC(vProd<num>(curext),(Cell<C,BT>*)NULL);

		if (DCTTALK)
		{
			cout<<"\n#"<<i<<", addin: "; c_print<vector<num> >(curaddin);
			cout<<" and extent: "; c_print<vector<num> >(curext);
			cout<<"\n&&&& alloc vec of size "<<vProd<num>(curext)<<" for addin "; c_print<vector<num> >(curaddin); cin.get();
		}
		// and add curvec to "img" of addin vector
		cube.push_back(curvec);
	}
}

template <typename C, typename BT>
bool DenseCToplex<C,BT>::getAddin(num pos, vector<num>& add) const // recover addin vector from lexico position
{
	pos++;
	num dimpos = 0; num k = 0;
	num n = extent.size();

	if (DCTTALK) cout<<"\nExtracting addin for position "<<pos;

	// first determine the eventual size from position
	do
	{
		dimpos += nChoosek(n,k);
		if (DCTTALK) cout<<"\n>>>>>>>>>>>>> n "<<n<<" k "<<k<<" dimpos "<<dimpos;
		if (pos <= dimpos) break;
		k++;
	} while (k <= n);
	// check validity of k
	if (k > n) return false;
	// adjust dimpos to be less than pos again
	dimpos -= (nChoosek(n,k));
	num p = pos - dimpos;
	if (p==0) k--;
	// finally call lexicopos
	lexicoPos (n,k,p,add);

	if (DCTTALK)
	{
		cout<<"\n ... getting Addin, lexicopos call with n "<<n<<" k "<<k<<" p "<<p<<" returning ";
		c_print<vector<num> >(add);
	}

	return true;
}


// recovers anchor point from given position "pos" and addin vector "advec".
// stores in "ans".
template <typename C, typename BT>
bool DenseCToplex<C,BT>::getAnchor(num pos,  const vector<num>& advec, vector<num>& ans) const // recover anchor point from lexico position
{
	vector<num> dext;
	// first, determine the extents for this addin vector!
	getExtForAddin(advec,dext);

	// now we have the extents, and can decompose to recover the anchor point!
	decompose<num,num>(pos, dext, ans);
	return true;
}

// recovers anchor point from given position "pos" and extents vector "ext".
// stores in "ans".
template <typename C, typename BT>
bool DenseCToplex<C,BT>::getAnchorFast(num pos, const vector<num>& ext, vector<num>& ans) const // recover anchor point from lexico position
{
	// now we have the extents, and can decompose to recover the anchor point!
	decompose<num,num>(pos, ext, ans);
	return true;
}

// recovers anchor point from given position "pos" and addin vector index "addind".
// stores in "ans".
template <typename C, typename BT>
bool DenseCToplex<C,BT>::getAnchor(num pos,  const num addind, vector<num>& ans) const // recover anchor point from lexico position
{
	// extract addin vector from index!
	vector<num> advec;
	if (!getAddin(addind,advec)) return false;

	return getAnchor(pos,advec,ans);
}

template <typename C, typename BT>
void DenseCToplex<C,BT>::getExtForAddin(const vector<num>& advec, vector<num>& ext) const
{
	ext = extent; // inherit max possible extent...

	//go over addin vector, decrementing extent at each dimension found
	num dim = advec.size();
	for (num i=0; i<dim; i++)
	{
		// if this dimension shows up in the addin vector, decrement the extent!
		ext.at(advec.at(i))--;
	}
}

template <typename C, typename BT>
num DenseCToplex<C,BT>::getAnchPos(const vector<num>& anch, const num addpos) const
{
	// recover addin vector for this addin position
	vector<num> advec;
	if(!getAddin(addpos,advec)) return 0;
	// now call the function to recover anchor position from addin vector
	return getAnchPos(anch,advec);
}

template <typename C, typename BT>
num DenseCToplex<C,BT>::getAnchPos(const vector<num>& anch, const vector<num>& advec) const
{
	num toret = 0;
	// now use this addin to generate extents vector
	vector<num> dext;
	getExtForAddin(advec,dext);
	// finally, use this extent vector to recompose the answer
	toret = recompose<num>(anch,dext);
	return toret;
}

template <typename C, typename BT>
num DenseCToplex<C,BT>::getAnchPosFast(const vector<num>& anch, const vector<num>& ext) const
{
	return recompose<num>(anch,ext);
}

template <typename C, typename BT>
pair<num,bool> DenseCToplex<C,BT>::makeFromFile (ifstream& infile)
{
	if (CUBTOPTALK)
	{
		cout<<"\n file reading "; cin.get();
	}
	num topdim=0;
	// get top dimension of cubical set from the first line of input file
	if(!infile.eof()) infile >> topdim; else return make_pair(0,false);
	if(topdim<0) return make_pair(0,false); // negative dimension count??? i quit!

	vector<num> ext;
	num extval; // stores extent value
	num totnum = 1;

	// build dimension-bounding vector
	for(int i=0; i<topdim; i++)
	{
		// add every dimension to addins since cubes are maximal
		if(infile.good())
		{
			infile >> extval;
			if (extval <= 0) // ignore non-positive dimension!
			{
				topdim--;
				if (topdim == 0)
				{
					return make_pair(0,false);
				}
				continue;
			}
			totnum*=extval;
			ext.push_back(extval);
		}
		else
		{
			// if file finishes prematurely, then the declared topdim is wrong
			return make_pair(0,false);
		}
	}

	// now initialize with this extent vector!
	Init(ext);
	bool incb = false; // increment births while making faces in case INITBT is encountered...

	// store top cube information into the appropriate spot of cube structure!
	num count = 0, actct = 0;
	num topadd = (num) pow ((float)2, (float)topdim);
	CVEC* topvec = cube.at(topadd-1);
	BT birthtime;

	//cout<<"\ntotnum = "<<totnum;


	while(infile.good() && count < totnum)
	{

		count++;
		// get current birth-time from file
		infile >> birthtime;

		//if (birthtime == 510) cout<<" \n 510 encountered at "<<count;


		//if (infile.eof()) break;
		if (birthtime == BANBT) continue; // ignore "banned" birthtime, as defined in Global/Pflags.h


        // kelly test...
		//if (birthtime > 513) continue;

		if (birthtime == INITBT) incb = true; // if INITBT is found, set flag to increment all
		// create new top cube and set its birth time
		topvec->at(count-1) = new Cell<C,BT>(topdim);
		topvec->at(count-1)->birth = birthtime;
		actct++;
	}
	if(CUBTOPTALK)
	{
		cout<<" File read... ";
		cin.get();
	}
	// with top cubes loaded, make faces!
	makeFaces(incb);
	//cout<<*this;
	//cin.get();
	return make_pair(actct,true);
}

// warning: enabling savemem kills original toplex while copying over cells to complex
template <typename C, typename BT>
bool DenseCToplex<C,BT>::quickWriteComplex(MComplex<C,BT>& mycomp, bool savemem = true)
{

	if (cube.size()==0) return false; // nothing to write

	// obtain top dimension:
	num topdim = extent.size();
	if (DCTSPEEDTALK){cout<<"\n\n *** quicksert with topdim = "<<topdim; cin.get();}

	set<BT> births;
	// now extract set of all birthtimes:
	CVEC* topvec = cube.at(cube.size()-1);

	if (topvec == NULL) return false; // no top cubes!

    if (DCTSPEEDTALK){cout<<"\n  inserting births "<<topdim; cin.get();}

	// iterate over top dimensional cubes, make a set of all birth times!
	typename CVEC::const_iterator topit;
	for (topit = topvec->begin(); topit != topvec->end(); ++topit)
	{
		if(*topit != NULL) births.insert((*topit)->birth);
	}
	// and get total number of birth times encountered:
	// create birth vector for position lookup
	vector<BT> bvec;
	bvec.assign(births.begin(), births.end());
	num bsize = bvec.size();
	births.clear();


	// create 3d array of cubes, a vector for each by birth and dimension:
	vector<vector<CVEC*>*> insertArray;
	make3Dvector<Cell<C,BT>*>(bsize,topdim,insertArray); // allocate array space...
	if (DCTSPEEDTALK)
    {
        cout<<"\n 3d vec called with "<<bsize<<" and "<<topdim; cin.get();
    }
	// now loop over cube structure, and copy over cells to the insert Array at the appropriate spots
	typename CGRID::const_iterator grit; // grid iterator
	typename CVEC::const_iterator vit; // vector iterator
	CVEC* curvec; // current cube vector
	Cell<C,BT>* curcell;
	//cout<<"\nReordering..."; cin.get();
	for(grit = cube.begin(); grit != cube.end(); ++grit)
	{
		curvec = *grit;
		if (curvec == NULL) continue;
		for (vit = curvec->begin(); vit != curvec->end(); ++vit)
		{
			curcell = *vit;
			// add current cube (if non null) to insertion structure:
			if (curcell == NULL) continue;

			// okay, this is complicated: the first dim of insertArray is the birth time, and to
			// find the position corresponding to the current birth we must look up the sorted
			// vector of births, i.e. bvec. the second dim of insertArray is just the cell dimension

			//print bdry of cube in case of 2 coefficients!
//			if (curcell->getDim() == 2)
//			{
//				//cout<<"\n bd of insert: "<<*curcell<<" is "<<curcell->getBD();
//			}
			//cout<<"\n    insert of "<<*curcell<<" attempted at: "<<getVPos<BT>(bvec,curcell->birth)<<", "<<curcell->getDim(); cin.get();
			insertArray.at(getVPos<BT>(bvec,curcell->birth))->at(curcell->getDim())->push_back(curcell);
		}
		if (DCTSPEEDTALK) {cout<<" end of line, clearing "; cin.get();}
		if (savemem)
		{
			curvec->clear();
			CVEC(*curvec).swap(*curvec);
		}
	}
	//cout<<"Done!"; cin.get();
	if( DCTSPEEDTALK)
	{
		cout<<"\n Insert to cell complex.... "; cin.get();
	}
	// at this point, we have the cell info stored in insertArray.

	// lastly, insert using insertArray!
	for (num btime = 0; btime < bsize; ++btime)
	{
		if (DCTSPEEDTALK) cout<<"\n Inserting birth: "<<bvec.at(btime);
		for (num curdim = 0; curdim < topdim; ++curdim)
		{
			if (DCTSPEEDTALK)
			{
				cout<<"\n        dim: "<<curdim<<", size: "<<insertArray.at(btime)->at(curdim)->size(); //cin.get();
				cin.get();
			}
			mycomp.quicksert(*(insertArray.at(btime)->at(curdim)));
			insertArray.at(btime)->at(curdim)->clear();
		}
		//insertArray.at(btime)->clear();
	}
	// deallocate insertArray structure!
	vector<CVEC*>* cur2dvec;
	CVEC* cur1dvec;
	for (num i = 0; i < (num) insertArray.size(); ++i)
	{
		cur2dvec = insertArray.at(i);
		for (num j = 0; j < (num) cur2dvec->size(); ++j)
		{
			cur1dvec = cur2dvec->at(j);
			delete cur1dvec;
		}
		cur2dvec->clear();
		vector<CVEC*>(*cur2dvec).swap(*cur2dvec);
		delete cur2dvec;
	}
	insertArray.clear();
	vector<vector<CVEC*>*>(insertArray).swap(insertArray);
	return true; // done
}

// adds a single top cube via its coordinates and birth time. Call AFTER init()
template<typename C, typename BT>
void DenseCToplex<C,BT>::addTopCube(const vector<num>& coords, const BT& btime)
{
    BT birth = btime;
	// store top cube information into the appropriate spot of cube structure!
	num topdim = coords.size();
	num addinpos = (num) pow ((float)2, (float)topdim) - 1;

	CVEC* topvec = cube.at(addinpos);
	Cell<C,BT>* toins;

    if (birth != BANBT)
    {
		// create new top cube and set its birth time
		toins = new Cell<C,BT>(topdim);
        toins->birth = birth;
		//cout<<" ins pos: "<<getAnchPos(coords,addinpos)<<" vs sz "<<topvec->size();
		//cin.get();
		topvec->at(getAnchPos(coords,addinpos)) = toins;
	}
}

// call after Init() and AddCube()
template<typename C, typename BT>
void DenseCToplex<C,BT>::ComputePersistence(string fname = "output")
{
    // make faces from top cubes: since only miro uses this, and he can not
    // stop putting in zeros, increment all the birth times just in case
    makeFaces(true);

    // write to complex
    MComplex<C,BT>* cptr = new MComplex<C,BT>;
    this->writeComplex(*cptr);


    //cout<<" about to cored: "; cin.get();
    // perform morse coreductions
    AlternateAndUpdate(cptr, false, 0.2);
    //cubcomp.MorseWrapper_Cored(false);

    // compute persistence intervals of reduced complex
    PComplex<C,BT> pcomp;
    pcomp.COMPUTE_INTERVALS(*cptr);

    // generate output
    // pcomp.showBetti();
    pcomp.makeOutputFiles(fname);
    delete cptr;

}

#endif /* DENSECTOPLEX_HPP_ */
