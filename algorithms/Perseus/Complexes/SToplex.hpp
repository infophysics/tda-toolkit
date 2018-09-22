/*
 * Toplex.hpp
 *
 *  Functions to create Cell complexes out of simplicial toplex information!
 */

#ifndef TOPLEX_HPP_
#define TOPLEX_HPP_

# include "SToplex.h"

// inserts a single simplex into the complex, preferring bottom-up construction
// in terms of dimension. Useless for toplexes in general, but useful for Rips
// complexes etc.
template <typename C, typename PS, typename BT>
bool SToplex<C,PS,BT>::addCoface(Simplex<C,PS,BT>* toin)
{
	if (toin == NULL) return false;
	PT_SET verts = toin->getVerts(); // extract set of vertices of this simplex.

	typename PT_SET::iterator hat;
	// generate vertex sets for all boundary simplices and check for inclusion
	// so as to update boundary/coboundary links suitably.
	C coeff; int facect = 0; // incidence coefficient!
	Cell<C,BT>* curbd; // current boundary simplex

	//cout<<"  add attempt: "<<*toin; cin.get();

	// simplex already exists in complex, we exit with failed status
	if (isIn(verts) != NULL)
	{
		//cout<<"   Doubling Insert???? "<<*toin;
		//cin.get();
		return false;
	}
	//cout<<"start comp insert"; cin.get();

	// ADD TO TOPLEX STRUCTURE
	typename SIMP_MAP::iterator finddim = mymap.find(toin->getDim());
	SIMP_RG* rg;

	// possibly we don't have this dimension yet!
	if (finddim == mymap.end())
	{
		// so make the range,
		rg = new SIMP_RG;
		// and insert it into the map, probably at the highest level!
		mymap.insert(mymap.begin(),make_pair(toin->getDim(),rg));
	}
	else // found dim, so extract range
	{
		rg = finddim->second;
	}

	// and insert at appropriate hash!
	rg->insert(make_pair(toin->getHash(getHashMod()),toin));


	// MAKE BOUNDARY AND COBOUNDARY LINKS
	verts = toin->getVerts();
	//cout<<" start making cb links "; cin.get();
	for (hat = toin->begin(); hat != toin->end(); ++hat)
	{
		//cout<<" vertex io: "; cin.get();
		verts.erase(*hat); // hat out current vertex
		coeff = (facect++ % 2 == 0) ? 1 : -1; // set coefficient

		// now locate the boundary simplex corresponding to the remaining vertices!
		curbd = isIn(verts);
		if (curbd == NULL)
		{
            //cout<<"\n               fail bdry #"<<facect<<" for "<<*toin; cin.get();
			// cptr_print<PT_SET >(verts); cin.get();
			return false; // could not find a boundary simplex!!!
		}
		// otherwise, boundary simplex is found! we make boundary links...
		curbd->addCBLink(toin,coeff);
		// and restore hatted vertex for next time
		verts.insert(*hat);
	}
	return true;
}


// checks to see if simplex determined by tofind lies in given list of simplex pointers
// returns pointer of found simplex, or NULL
template <typename C, typename PS, typename BT>
Simplex<C,PS,BT>* SToplex<C,PS,BT>::isIn (const PT_SET& tofind)
{
	num curd = tofind.size() - 1; // dimension we are seeking
	typename SIMP_MAP::iterator finddim = mymap.find(curd);
	if (finddim == mymap.end()) return NULL; // dimension not found


	// otherwise, dimension found at iterator position finddim!

	Simplex<C,PS,BT>* toret = NULL; // initial return value
	num hash = getPtSetHash<PS>(tofind, getHashMod());

	SIMP_RG* myrange = finddim->second; // range of simplex
	typename SIMP_RG::iterator cursimp;
	pair<typename SIMP_RG::iterator, typename SIMP_RG::iterator> ancrange;
	// bounds of iterators for equal range
	ancrange = myrange->equal_range(hash);
	// search across this range...
	for (cursimp = ancrange.first; cursimp != ancrange.second; ++cursimp)
	{
		if (cursimp->second->getVerts() == tofind) // compare vertices!
		{
			toret = cursimp->second; // found! set return value appropriately
			break; // and get out of loop, we are done here!
		}
	}
	return toret;
}

// make faces of a given simplex
template <typename C, typename PS, typename BT>
bool SToplex<C,PS,BT>::makeFacesOf(Simplex<C,PS,BT>* simp, bool inheritBirths = true)
{
	if (simp == NULL) return false; // can't handle null simplex
	if (simp->getDim()==0) return true; // no faces for 0 dimensional/already faced simplices!

	PT_SET verts = simp->getVerts(); // boundary vertex set
	typename PT_SET::iterator hat; // points to vertex which is being ignored when making face
	C coeff; // coefficient of this face
	num facect = 0;


	typename SIMP_MAP::iterator finddim;
	typename SIMP_RG::iterator findhash;
	SIMP_RG* currg;

	Simplex<C,PS,BT>* sloc; // face simplex possibly located in complex

	// iterate over vertices of this simplex.
	for (hat = simp->begin(); hat != simp->end(); ++hat)
	{
		verts.erase(*hat); // hat out this vertex
		coeff = (facect % 2 == 0) ? 1 : -1; // and its incidence coefficient

		// now check if current bdry simplex is in its anchor list
		//cout<<"isin call: "; cin.get();
		sloc = isIn(verts); // check if face is in corresponding list
		//cout<<" back with "<<sloc; cin.get();


		if (sloc == NULL) // no such face exists
		{
			// so add this new simplex to the map structure
			sloc = new Simplex<C,PS,BT>(verts);
			//cout<<"\n    attempting insertion of: "<<*sloc; cin.get();
			// it is possible that there  tDim());
			finddim = mymap.find(sloc->getDim());
			if (finddim == mymap.end())
			{
				currg = new SIMP_RG;
				mymap.insert(mymap.end(),make_pair(sloc->getDim(),currg));
				//cout<<"inserted!"; cin.get();
			}
			else
			{
				currg = finddim->second;
			}
			//cout<<"range insert!"; cin.get();

			currg->insert(make_pair(sloc->getHash(getHashMod()),sloc));
			//cout<<"\n     added: "<<*bdsimp<<" to "<<*firstbdvert<<" of dim "<<curd-1<<" at "<<bdsimp;
		}
		simp->addBDLink(sloc,coeff,true,inheritBirths);
		facect++; // update number of faces handled
		//cout<<"\n^^^^^^^^^^^^^^^^^added face: "<<*sloc; cin.get();
		// restore boundary for next face, idiot:
		verts.insert(*hat);
	}

//	// now make faces of faces etc...
//	typename list<Simplex<C,PS,BT>*>::const_iterator face;
//	for (face = facelist.begin(); face != facelist.end(); ++face)
//	{
//		makeFacesOf(*face); // recursive call
//	}
	return true;
}

// adds faces of given simplex to the existing top-dim map vector of simplices
template <typename C, typename PS, typename BT>
bool SToplex<C,PS,BT>::addSimpFaces(bool inheritBirths = true)
{
	typename SIMP_MAP::const_iterator map_iter; // iterates over map
	//num curdim = 0; // current dim being iterated over

	typename SIMP_RG::const_iterator rg_iter; // iterates over range of current point
	SIMP_RG* currange;

	num count = 0;
	//num rgcount = 0;
	// iterate over the map:
	for (map_iter = mymap.begin(); map_iter != mymap.end(); ++map_iter)
	{
		count++;
		//cout<<"\n facing for dim: "<<map_iter->first; cin.get();
		//cout<<*this; cin.get();
		//curdim = map_iter->first; // extract current dim

		//if (curdim == 0) continue; // don't make faces of 0 dimensional stuff!


		currange = map_iter->second; // and the range
		//cout<<"\n     of size: "<<currange->size(); cin.get();
		//cout<<"\n\n Point: "<<*curpt; cin.get();
		//rgcount = 0;
		// iterate over the range!
		for (rg_iter = currange->begin(); rg_iter != currange->end(); ++rg_iter)
		{
			//cout<<"\n        face#: "<<rgcount++; cin.get();
			//cout<<".... handling: "<<rgcount; cin.get();
			//cout<<"\n\nprocessing simplex of dim: "<<(rg_iter->second)->getDim(); cin.get();
			// and make faces of this simplex.
			//cout<<" call going: "; cin.get();
			makeFacesOf(rg_iter->second, inheritBirths);
		}
	}
	return true;
}



// now store all the points encountered as keys in a map, where
// the data for any point is (pointers to) ALL the simplices whose first
// vertex is that point.
template <typename C, typename PS, typename BT>
bool SToplex<C,PS,BT>::buildTopSimpMap (const vector<Simplex<C,PS,BT>*>& topsimps, bool makeverlist = true)
{
    verlist.clear();
	//cout<<" top simp size "<<topsimps.size(); cin.get();
	if (topsimps.size()==0 || topsimps.front()==NULL) return false;  // no simplices to add

	num topd;

	Simplex<C,PS,BT>* cursimp;

	set<Point<PS>*, ptcomplex<PS> > vertices; // store all vertices
	typename PT_SET::const_iterator viter; //

	// iterate over top simplices and extract vertices:
	typename vector<Simplex<C,PS,BT>*>::const_iterator curtop;
	for (curtop = topsimps.begin(); curtop != topsimps.end(); ++curtop)
	{
		cursimp = *curtop;
		topd = cursimp->getDim(); // get top simplex's dimension
		// update top dim of toplex!
		if (topd > topdim) topdim = topd;

		// now extract vertices!
		for (viter = cursimp->begin(); viter != cursimp->end(); ++viter)
		{
			vertices.insert(*viter);
		}
		if (makeverlist) verlist.insert(vertices.begin(),vertices.end());
	}

    if (SIMTOPTALK)
    {
        cout<<"\n built vertex list";
        if (MAKEBPS) cin.get();
    }

	// set number of distinct vertices:
	psize = vertices.size();

	typename SIMP_MAP::iterator mit;
	SIMP_RG* currg;
	num hmod = getHashMod();

	// now we can compute the hashes using topdim and psize, and hence
	// store the top simplices. sweep over them again...
	for (curtop = topsimps.begin(); curtop != topsimps.end(); ++curtop)
	{
		cursimp = *curtop;
		mit = mymap.find(cursimp->getDim()); // seek this hash
		if(mit == mymap.end()) // not found!
		{
			currg = new SIMP_RG; // define range, and...
			mymap.insert(mymap.end(),make_pair(cursimp->getDim(),currg));
		}
		else currg = mit->second;
		currg->insert(currg->end(),make_pair(cursimp->getHash(hmod),cursimp));
	}

    if (SIMTOPTALK)
    {
        cout<<"\n top map size: "<<mymap.size();
        if (MAKEBPS) cin.get();
    }
	//cout<<" top map size: "<<mymap.size(); cin.get();
	return true;
}

template <typename C, typename PS, typename BT>
num SToplex<C,PS,BT>::getHashMod() const
{
	return 1+(1+topdim)*psize*psize;
}


template <typename C, typename PS, typename BT>
void SToplex<C,PS,BT>::writeComplex(MComplex<C,BT>& simpcomp) const
{
	// now we fill up the complex with the simplices in mymap!
	typename SIMP_MAP::const_iterator mypair; // iterates over map
	typename SIMP_RG::const_iterator cursimp;
	SIMP_RG* currg;

	for (mypair = mymap.begin(); mypair != mymap.end(); ++mypair)
	{
		currg = mypair->second;
		for (cursimp = currg->begin(); cursimp != currg->end(); ++cursimp)
		{
			// insert this cell
			simpcomp.insertCell((Cell<C,BT>*)cursimp->second);
			//cout<<"\ninserting: "<<**cursimp;
		}
	}
}


// routine to make a "Complex<C,BT> simpcomp" out of a vector of top-dimensional
// simplices... the dimension of the first simplex in the vector is assumed to be
// universal, and simplices in the vector with other dimensions are ignored!
template <typename C, typename PS,typename BT>
bool SToplex<C,PS,BT>::makeCellComplex(const vector<Simplex<C,PS,BT>*>& topsimps, MComplex<C,BT>& simpcomp, bool inheritBirths = true)
{
	int topd = populateMap(topsimps, inheritBirths);
	if (topd < 0) return false;
	simpcomp.Init();
	writeComplex(simpcomp); // make cell complex!
	return true;
}

template <typename C, typename PS,typename BT>
num SToplex<C,PS,BT>::populateMap(const vector<Simplex<C,PS,BT>*>& topsimps, bool inheritBirths = true)
{
	if (topsimps.size()==0 || topsimps.front()==NULL) return -1;  // no simplices to add
	int topd = topsimps.front()->getDim(); // and the top dimension
	if(buildTopSimpMap(topsimps))
	{
		addSimpFaces(inheritBirths); // add all faces to map structure
	}
	return topd;
}

// wipes out the map structure;
// if deep is enabled, all pointers to points and simplices are deleted
template <typename C, typename PS,typename BT>
void SToplex<C,PS,BT>::Destroy(bool deep=false)
{
    //cout<<"\n :: SToplex destroyer \n\n    vertices... "; cin.get();

    //PT_SET vertices; vertices.clear();
    //getVerts(vertices);

    //cout<<" gotverts "; cin.get();
    // delete vertices
    typename PT_SET::iterator vit;
    for (vit = verlist.begin(); vit != verlist.end(); ++vit)
    {
       delete *vit;
    }

    //cout<<"     done! \n\n stop struct... "; cin.get();

    // now for simplices...
	typename SIMP_RG::iterator cl_it; // iterate over multimap range
	// iterate over the map
	typename SIMP_MAP::iterator m_iter;
	SIMP_RG* currg;

	for (m_iter = mymap.begin(); m_iter != mymap.end(); ++m_iter)
	{
		currg = m_iter->second;

		if (deep)
		{
			for (cl_it = currg->begin(); cl_it != currg->end(); ++cl_it)
			{
				// kill simplex
				delete cl_it->second;
			}
			//delete curpt; // and also the point
		}
		// but certainly we should clear the multimap
		currg->clear();
		delete currg;
	}
	// and now clean up the map
	mymap.clear();
	//cout<<" \n TERMINATED!"; cin.get();
}

// extract vertex points!
template <typename C, typename PS, typename BT>
void SToplex<C,PS,BT>::getVerts(PT_SET& allverts) // extracts all vertices;
{
	allverts.clear();
	// iterate over map
	typename SIMP_MAP::const_iterator miter;
	typename SIMP_RG::const_iterator riter;
	SIMP_RG* currg;
	Simplex<C,PS,BT>* cursimp;
	//PT_SET::const_iterator pit; // iterate over current simps vertices...

    //cout<<"\n getting 0 dim... "; cin.get();

    miter = mymap.find((num) 0); //extract 0-dim simplices

	if (miter != mymap.end()) // if there are some...
	{
	    //cout<<" getting rng "; cin.get();
	    currg = miter->second;
	    for (riter = currg->begin(); riter != currg->end(); ++riter)
	    {
	        //cout<<"\n getting simplex... "; cin.get();
	        cursimp = riter->second;
	        //cout<<" got... insert?"; cin.get();
	        allverts.insert(cursimp->verts.begin(),cursimp->verts.end()); // first and only vertex: dimension is zero
            //cout<<" ins suc @"; cin.get();
	    }
	}
}

// creates simplicial toplex from file input
template <typename C, typename PS, typename BT>
pair<num,bool> SToplex<C,PS,BT>::makeFromFile(ifstream& infile, bool hasbirths = true,
                                              bool mfd = true, bool hasdimcap = false)
{
	num count = 0; // total number of top simplices read
	num topdim=0; // top dimension of simplices encountered
	num ptdim = 0; // dimension of the space from which points are being picked
	// get top dimension of simplices (i.e. #vertices-1)
	if(infile.good())
	{
		if(mfd) infile >> topdim; // get uniform top dimension from file in case of manifold triangulation
	}
	else return make_pair(0,false);

	if(infile.good()) infile >> ptdim; else return make_pair(0,false);

	// extract cap on dimension, if needed!
	num dimcap = 0;
	if (hasdimcap)
	{
		if (infile.good()) infile >> dimcap;
		else return make_pair(0,false);
	}

	map<Point<PS>*, num, ptcomplex<PS> > pts; // set of all points and corresponding indices
	typename map<Point<PS>*, num, ptcomplex<PS> >::iterator finder;

	// map from set of vertices to birth times...
	map<PT_SET, BT> vertbirth;
	pair<typename map<PT_SET, BT>::iterator, bool> insres; // result of inserting into map
	typename map<PT_SET, BT>::iterator curptsetpos; // position of current point-set in vertbirth


	if((topdim < 0) || ptdim==0) return make_pair(count,false); // negative dimension count??? i quit!

	Point<PS>* curpt; // current point being read from file
	PT_SET curptset; // vertex set for current simplex
	PS curcoord; // current coordinate of point
	BT curborn; // current birth time

	srand(time(NULL));

	// now we are prepared to receive simplex data...
	// here if ptdim = m and topdim = n, each simplex is a collection of (m+1)n coordinates,
	// the first n are for vertex 1, the next n for vertex 2,... etc.

	num curptind = 0;
    curborn = INITBT; // initialize current birth time
    bool breakout = false;

	while (infile.good())
	{
		if (!mfd) infile >> topdim; // get top dim of each simplex in case of non manifoldness
		if (topdim < 0) continue;
		curptset.clear();
		// loop over simplex dimension
		for (num j = 0; j <= topdim; j++)
		{
		    if (infile.eof() && j < topdim)
            {
                breakout = true;
                break;
			}
		    curpt = new Point<PS>;

			// loop over point dimension
			for(num i = 0; i < ptdim; i++)
			{
				infile >> curcoord; // get current coordinate
				if (infile.eof() && i < ptdim - 1)
				{
				    breakout = true;
				    break;
				}
				curpt->push(curcoord); // and add to current point
			}

			if (breakout)
			{
			    delete curpt;
			    curpt = NULL;
			    break;
			}

			finder = pts.find(curpt);
			// if not found, add suitably
			if (finder == pts.end())
			{
				// add to this simplex's vertex list
				curptset.insert(curpt);
				// and to set of all vertices with its index
				pts.insert(pts.end(),make_pair(curpt,curptind++));
				//cout<<"\n adding pt: "<<*curpt;
			}
			else // otherwise we already have this point located by finder, so...
			{
				delete curpt; // no need to store new copy
				curptset.insert(finder->first); // insert old copy into vertex set
			}
		}
		if (breakout)
		{
		    break;
		}

		if (hasbirths)
        {
            infile >> curborn;
        }

		// check to see if we already have this point-set:
		// that is, are we being fed a duplicate simplex???
		insres = vertbirth.insert(make_pair(curptset,curborn));
		curptsetpos = insres.first; // extract location in vertbirth of curptset

		if (insres.second == false) // point set already existed, non unique!
		{
			// if the new birth time is smaller, update!
			// WHEN THE LIST OF TOP SIMPLICES IN THE FILE IS NOT UNIQUE, ONLY
			// THE OLDEST (SMALLEST BIRTH TIME) WILL SURVIVE!
			if (curborn < curptsetpos->second) curptsetpos->second = curborn;
		}
	}


	// if we are restricting dimension,...
	if (hasdimcap)
	{
		capDimension(vertbirth, dimcap);
	}
	// now extract number of top simplices!
	count = (num)vertbirth.size(); // total number of unique simplices found
	Simplex<C,PS,BT>* cursimp;
	// now make top simplex vector from the capped vertbirth structure
	vector<Simplex<C,PS,BT>*> topsimps;
	// loop over vertbirths
	//cout<<"szinfo0"; cin.get();

	typename map<PT_SET,BT>::const_iterator vset;
	for (vset = vertbirth.begin(); vset != vertbirth.end(); ++vset)
	{
		cursimp = new Simplex<C,PS,BT>(vset->first); // make new simplex
		cursimp->birth = vset->second; // and feed it a birth time
		// now store it
		topsimps.push_back(cursimp);
	}

	if (SIMTOPTALK)
	{
	    cout<<"\n processed "<<topsimps.size()<<" file simplices... ";
	    if (MAKEBPS) cin.get();
	}

	// now build top simplices
	buildTopSimpMap(topsimps,false); // build map for top simplices...

	if (SIMTOPTALK)
	{
	    cout<<"\n top simp map built! ";
	    if (MAKEBPS) cin.get();
	}

	addSimpFaces(); // and add faces recursively.
	if (SIMTOPTALK)
	{
	    cout<<"\n faces made ";
	    if (MAKEBPS) cin.get();
	}
	//if (SIMTOPTALK) cout<<*this;
	return make_pair(count,true);
}

// given a map from unique point sets to births... representing simplices
// if any simplex in the vector has a dimension
// larger than "capdim" then this routine will replace that simplex with a
// collection of its capdim-dimensional faces. This ensures that the highest
// dimension simplex has dimension <= capdim
template <typename C, typename PS, typename BT>
num SToplex<C,PS,BT>::capDimension(map<PT_SET,BT>& vertbirth, num dimcap)
{
	num n, k = dimcap+1; // for n choose k
	bool trace = false;

	// to iterate over the top simplices!
	typename map<PT_SET,BT>::iterator vset;
	typename map<PT_SET,BT>::iterator finder; // result of insert attempt!
	vector<Point<PS>*> ptvec;
	PT_SET facepts;
	vector<num> faceverts;
	map<PT_SET,BT> toadd; // pairs to add to vertbirth
	vector<typename map<PT_SET,BT>::iterator > torem; // pairs to remove from vertbirth
	num count = 0; // how many simplices get capped?

    if (CDTALK)
    {
        cout<<"\n\n Capping Dimension to "<< dimcap <<" on "<<vertbirth.size()<<" Top Simplices!";
    }

	for (vset = vertbirth.begin(); vset != vertbirth.end(); ++vset)
	{
		// INSERT TRACES HERE
		//if (vset->second == 38) trace = true;

		// the current vertex set is vset->first. Check its size, the dimension of
		// the corresponding simplex has dimension (n - 1)
		n = (num)vset->first.size();
		// if the dimension is less than the cap, ...
		if (n <= k)
		{
		    if (CDTALK)
		    {
		        cout<<"\n     ignoring simplex ";
		        cptr_print<PT_SET>(vset->first);
		    }
			// ignore! this dimension is fine.
			continue;
		}
		count++;
		// IF WE GET HERE, WE NEED TO CONSTRUCT FACES AND INSERT!
		if(CDTALK || trace)
		{
			cout<<"\n\n "<<count<<". DECOMPOSING: ["<<n<<"]: "; cptr_print<PT_SET>(vset->first);
			cout<<" BORN AT "<<vset->second;
			cin.get();
		}
		// first, extract the current point set into a vector
		ptvec.assign(vset->first.begin(), vset->first.end());

		// loop over the number of size-(k)-subsets of [n]
		for (num p = 1; p <= nChoosek(n,k); ++p)
		{
			facepts.clear();
			lexicoPos(n,k,p,faceverts);
			// now faceverts contains indices into ptvec of the k-face we care about
			// make pt-set corresponding to these.
			for (num i = 0; i < (num)faceverts.size(); ++i)
			{
				facepts.insert(ptvec.at(faceverts.at(i)));
			}
			if (CDTALK || trace)
			{
			    cout<<"\n ----- RECOMPOSING: ["<<k<<"]: ";
			    cptr_print<PT_SET>(facepts);
			}
			// finally, try to insert facepts into vertbirth with the given birth
			// value from the top simplex which we are decomposing
			finder = toadd.find(facepts);
			if (finder == toadd.end()) // this face is not already present...
			{
				if (CDTALK || trace)
				{
					cout<<"\n ----- Keeping Birth Time of "<<vset->second<<", New add!";
				}
				// then insert it into toadd:
				toadd.insert(make_pair(facepts,vset->second));
			}
			else // if found, update birth time to be minimum
			{
				if (vset->second < finder->second)
				{
					if (CDTALK || trace)
					{
						cout<<"\n ----- Altering Birth Time to "<<vset->second<<" from "<<finder->second;

					}
					finder->second = vset->second;
				}
			}
		}
		torem.push_back(vset); // remove this vertex set, the dimension is too high.
	}
	// now we want to remove all the elements of vertbirth corresponding to
	// iterators in torem
	typename vector<typename map<PT_SET,BT>::iterator >::iterator killme;
	for (killme = torem.begin(); killme != torem.end(); ++killme)
	{
		vertbirth.erase(*killme); // erase over-dimensional simplex!
	}
	//finally, we loop over the elements in toadd and add to vertbirths
	typename map<PT_SET,BT>::iterator addme;
	pair <typename map<PT_SET,BT>::iterator, bool> insres; // result of insertion attempt
	for (addme = toadd.begin(); addme != toadd.end(); ++addme)
	{
		insres = vertbirth.insert(*addme);
		// in case this element already exists, compare birth times and set to minimum
		if (insres.second == false)
		{
		    //cout<<" \n weird ins dupe??";
			if (insres.first->second > addme->second)
			{
				insres.first->second = addme->second;
			}
		}
	}
	return (num)vertbirth.size();
}

// top level wrapper: goes from vector of simplices to generating
// persistence output directly. birth times are NOT inherited.
template <typename C, typename PS, typename BT>
bool SToplex<C,PS,BT>::topWrapper(const vector<Simplex<C,PS,BT>*>& topsimps, string fname = "output")
{
	MComplex<C,BT> mycomp;
	if (! makeCellComplex(topsimps, mycomp, false) ) return false;
	// loaded up all simplices into complex mycomp, so reduce!
	mycomp.MorseWrapper(true);
	// reductions complete, compute persistence intervals!
	PComplex<C,BT> pcomp;
	pcomp.COMPUTE_INTERVALS(mycomp);
	// persistence intervals computed, make output files
	pcomp.makeOutputFiles(fname);
	return true;
}

// barycentric subdivision, resulting stoplex stored in "sd"
template <typename C, typename PS, typename BT>
bool SToplex<C,PS,BT>::barySub(SToplex<C,PS,BT>& sd)
{
    // map each simplex sigma of this complex to a vertex in sd.
    Simplex<C,PS,BT>* sigma;

    // store barycenter of sigma...
    Point<PS>* bcpt; // the vertex,...
    Simplex<C,PS,BT>* bcsigma; // the actual 0-simplex.

    // iterate over simplicial complex...
    typename SIMP_MAP::iterator mit;
    typename SIMP_RG::iterator rit;
    SIMP_RG* currg;


    map<Simplex<C,PS,BT>*,Simplex<C,PS,BT>*> bcent;

    num simpcount = 0; // count simplices in this complex, ie vertices in sd.

    srand(time(NULL)); // initialize randomizer for hashes

    // BUILD MAP:: SIMPLEX TO BARYCENTER-VERTEX (0-Simplex in SD)
    for (mit = mymap.begin(); mit != mymap.end(); ++mit)
    {
        currg = mit->second;
        for (rit = currg->begin(); rit != currg->end(); ++rit)
        {
            sigma = rit->second; // extract simplex
            // generate point.
            bcpt = new Point<PS>; bcpt->push((PS)simpcount);

            // make zero simplex...
            bcsigma = new Simplex<C,PS,BT>(bcpt);
            // and build the map

            if (BARYTALK) cout<<"\n -- Map "<<*sigma<<" to "<<simpcount;

            bcent.insert(make_pair(sigma,bcsigma));
            simpcount++;
        }
    }

    // number of vertices in sd is known: will help generate hashes
    sd.psize = simpcount;
    vector<Simplex<C,PS,BT>*> sdsimplist;

    // now loop over simplices ASCENDING in dimension, to build subdivided simplices.
    typename SIMP_MAP::reverse_iterator revit;
    set<PT_SET> sdresults;

    // iteration over sdresults:
    typename set<PT_SET>::iterator setit;
    PT_SET cursdptset; // current subdiv point set
    Simplex<C,PS,BT>* sdsimp; // current simplex resulting from subdiv of sigma


    for (revit = mymap.rbegin(); revit != mymap.rend(); ++revit)
    {
        currg = revit->second;
        for (rit = currg->begin(); rit != currg->end(); ++rit)
        {
            sigma = rit->second; // current simplex
            //cout<<"\n    proc: "<<*sigma<<" inserting "<<*(bcent[sigma]);
            //sd.insertVertex(bcent[sigma]); // add corresponding barycenter to sd

            // go home for dimension zero:
            if (sigma->getDim() == 0) continue;

            sdresults.clear();
            // collect the results of subdividing sigma:
            subdivide(sigma, bcent, sdresults);

            //cout<<"\n    sdsize: "<<sdresults.size();
            // iterate over sdresults...
            for (setit = sdresults.begin(); setit != sdresults.end(); ++setit)
            {
                cursdptset = *setit; // extract current point set...

                // now make a simplex out of cursdptset!!
                sdsimp = new Simplex<C,PS,BT>(cursdptset);
                if (BARYTALK) cout<<" \n  contemplating coface: "<<*sdsimp;

                // and add to the sd complex as a coface: all the faces should be there already
                sdsimplist.push_back(sdsimp);
            }
        }
    }
    sd.buildTopSimpMap(sdsimplist);
    sd.addSimpFaces();
    return true;

}

// takes a vertex zerosimp and a simplex tosus.
// returns a SET of point-sets, representing all possible simplices arising from
// the suspension of tosus at zerosimp. In particular, all these simplices contain
// "zerosimp"...

template <typename C, typename PS, typename BT>
bool SToplex<C,PS,BT>::subdivide (Simplex<C,PS,BT>* sigma, map<Simplex<C,PS,BT>*,
                                  Simplex<C,PS,BT>*>& bc, set<PT_SET>& result)
{
    if (sigma == NULL) return false;
    // local variables... these are NASTY.
    PT_SET curptset;

    result.clear();


    if (sigma->getDim()==0) // base case of recursion: just return bc[sigma] in set...
    {
        // insert PT_SET of sigma into result:
        result.insert((bc[sigma])->getVerts());
        return true;
    }

    // For extracting the faces tau of sigma:
    Simplex<C,PS,BT>* tau;
    typename CSTRUCT::const_iterator cit;
    set<PT_SET> taures; // stores result of subdivide(tau)
    typename set<PT_SET>::iterator tit;
    PT_SET curset;

    for (cit = sigma->getBD().begin(); cit != sigma->getBD().end(); ++cit)
    {
        tau = (Simplex<C,PS,BT>*) cit->first; // extract tau...
        taures.clear(); // wipe the result structure...
        subdivide (tau, bc, taures); // recursive call...

        if (BARYTALK) cout<<"\n                 face "<<*tau<<" returns list of size "<<taures.size();

        // want to add bc[sigma] to EACH ELEMENT OF TAURES!
        for (tit = taures.begin(); tit != taures.end(); ++tit)
        {
            curset = *tit;
            curset.insert(bc[sigma]->getVerts().begin(),bc[sigma]->getVerts().end());
            result.insert(curset);
        }
    }
    return true;
}



template <typename C, typename PS, typename BT>
bool SToplex<C,PS,BT>::insertVertex (Simplex<C,PS,BT>*& tau)
{
    if (tau->getDim()!=0)
    {
        return false;
    }
    typename SIMP_MAP::iterator finddim = mymap.find(0);
    SIMP_RG* currg;

    // Insert point directly into the toplex structure
	if (finddim == this->mymap.end()) // found dimension 0?
	{
		// no, make range
		currg = new SIMP_RG;
		mymap.insert(this->mymap.end(),make_pair(0,currg));
	}
	else
	{
		// yes, get range
		currg = finddim->second;
	}

	currg->insert(make_pair(tau->getHash(this->getHashMod()),tau));
    return true;
}

// for ricci flow data
// creates simplicial toplex from file input
template <typename C, typename PS, typename BT>
pair<num,bool> SToplex<C,PS,BT>::makeRicci(ifstream& tetfile, ifstream& edgefile)
{
	num count = 0; // total number of top simplices read
	num topdim=0; // top dimension of simplices encountered
	num ptdim = 1; // dimension of the space from which points are being picked
	// get top dimension of simplices (i.e. #vertices-1)
	if(tetfile.good())
	{
		tetfile >> topdim; // get uniform top dimension from file in case of manifold triangulation
	}
	else return make_pair(0,false);

	map<Point<PS>*, num, ptcomplex<PS> > pts; // set of all points and corresponding indices
	typename map<Point<PS>*, num, ptcomplex<PS> >::iterator finder;

	// map from set of vertices to birth times...
	map<PT_SET, BT> vertbirth;
	pair<typename map<PT_SET, BT>::iterator, bool> insres; // result of inserting into map
	typename map<PT_SET, BT>::iterator curptsetpos; // position of current point-set in vertbirth


	if((topdim < 0) || ptdim==0) return make_pair(count,false); // negative dimension count??? i quit!

	Point<PS>* curpt; // current point being read from file
	PT_SET curptset; // vertex set for current simplex
	PS curcoord; // current coordinate of point

	srand(time(NULL));

	// now we are prepared to receive simplex data...
	// here if ptdim = m and topdim = n, each simplex is a collection of (m+1)n coordinates,
	// the first n are for vertex 1, the next n for vertex 2,... etc.

	num curptind = 0;
    bool breakout = false;

	while (tetfile.good())
	{
		curptset.clear();
		// loop over simplex dimension...
		for (num j = 0; j <= topdim; j++)
		{
		    if (tetfile.eof() && j < topdim)
            {
                breakout = true;
                break;
			}
		    curpt = new Point<PS>;

			// get the single coordinate, point dimension = 1 !
			for(num i = 0; i < ptdim; i++)
			{
				tetfile >> curcoord; // get current coordinate
				if (tetfile.eof() && i < ptdim - 1)
				{
				    breakout = true;
				    break;
				}
				curpt->push(curcoord); // and add to current point
			}

			if (breakout)
			{
			    delete curpt;
			    curpt = NULL;
			    break;
			}

			finder = pts.find(curpt);
			// if not found, add suitably
			if (finder == pts.end())
			{
				// add to this simplex's vertex list
				curptset.insert(curpt);
				// and to set of all vertices with its index
				pts.insert(pts.end(), make_pair(curpt,curptind++));
				//cout<<"\n adding pt: "<<*curpt;
			}
			else // otherwise we already have this point located by finder, so...
			{
				delete curpt; // no need to store new copy
				curptset.insert(finder->first); // insert old copy into vertex set
			}
		}
		if (breakout)
		{
		    break;
		}

		// uniquely populate point structure, all births are INITBT right now
		insres = vertbirth.insert(make_pair(curptset,INITBT));
		curptsetpos = insres.first; // extract location in vertbirth of curptset

	}

	// now extract number of top simplices!
	count = (num)vertbirth.size(); // total number of unique simplices found
	Simplex<C,PS,BT>* cursimp;
	// now make top simplex vector from the capped vertbirth structure
	vector<Simplex<C,PS,BT>*> topsimps;
	// loop over vertbirths
	//cout<<"szinfo0"; cin.get();

    // populate top simplices...
	typename map<PT_SET,BT>::const_iterator vset;
	for (vset = vertbirth.begin(); vset != vertbirth.end(); ++vset)
	{
		cursimp = new Simplex<C,PS,BT>(vset->first); // make new simplex
		cursimp->birth = vset->second; // and feed it INIT birth time
		// now store it
		topsimps.push_back(cursimp);
	}

	if (SIMTOPTALK)
	{
	    cout<<"\n processed "<<topsimps.size()<<" file simplices... ";
	    if (MAKEBPS) cin.get();
	}

	// now build top simplices
	buildTopSimpMap(topsimps,false); // build map for top simplices...

	addSimpFaces(); // and add faces recursively.

    // inherit all births from edge information!
    if(!readEdgeBirths(edgefile, pts))
    {
        cout<<"\n error: unable to read edge birth file!";
        cin.get();
        return make_pair(0,false);
    }
    inheritBirths(1);

	//if (SIMTOPTALK) cout<<*this;
	return make_pair(count,true);
}

// given an already populated simplicial complex of dim > 1,
// reads in a file with edge birth information and assigns
// that information to the edges. Then use InheritBirths(1)
// to propagate birth information to ALL simplices
template <typename C, typename PS, typename BT>
bool SToplex<C,PS,BT>::readEdgeBirths(ifstream& edgefile,
                        map<Point<PS>*, num, ptcomplex<PS> >& pts)
{
    // read edgefile to extract edge-birth pairs, each edge
    // being a PT_SET of size 2...
    Point<PS>* curpt = new Point<PS>;
    PT_SET curedge;
    PS coord;
    BT edgebirth;

    bool breakout = false;

    Simplex<C,PS,BT>* foundedge;

    // locate point in pts structure using this:
    typename map<Point<PS>*, num, ptcomplex<PS> >::iterator finder;

    if (edgefile.bad())
    {
        cout<<"\n unable to read edge file";
        return false;
    }
    //cout<<" made it ";

    // traverse edgefile:
    while(edgefile.good())
    {
        curedge.clear();
        // read in edge vertices:
        // this "2" changes to d+1 if we want arbitrary dimensions
        // rather than edges (d=1):
        for (num i = 0; i < 2; ++i)
        {
            if (edgefile.eof())
            {
                breakout = true;
                break;
            }

            curpt->clear();
            // read in first point
            edgefile >> coord;
            curpt->push(coord);

            // now we have a point... use finder:
            finder = pts.find(curpt);
            if (finder == pts.end())
            {
                cout<<"\n Edge location error: Vertex "<<*curpt<<" undefined!";
                cin.get();
                return false;
            }
            else
            {
                curedge.insert(finder->first);
            }
        }

        if (breakout) break;
        // now read in the birth value which we want to assign to the located edge
        edgefile >> edgebirth;

        // okay, curedge is what we want to search for in the
        // simplex structure
        foundedge = isIn(curedge);
        // WARN IF EDGE NOT FOUND!!
        if (foundedge == NULL)
        {
            cout<<"\n Warning: birth time supplied for edge which is not the face of any known\
                    tetrahedron... ignoring!";
            cin.get();
            continue;
        }
        // OTHERWISE UPDATE BIRTH ACCORDING TO FILE INFORMATION!
        else
        {
            if (edgebirth != INITBT)
            {
                foundedge->birth = edgebirth;
            }
        }
    }
    delete curpt;
    return true;
}

// uses the standard max/min approach to inherit
// birth information up or down the simplicial complex
// assuming the existence of birth info at dimension dim!
template <typename C, typename PS, typename BT>
bool SToplex<C,PS,BT>::inheritBirths(const num dim)
{
    //cout<<"\n\n inherit dim "<<dim;

    typename SIMP_MAP::const_iterator miter;
    typename SIMP_RG::const_iterator liter;
    Simplex<C,PS,BT>* cursimp;
    typename CSTRUCT::const_iterator chainit;

	typename SIMP_MAP::const_iterator finddim = mymap.find(dim);
	// if dimension not found, no way to access birthed simplices!
	if (finddim == mymap.end()) return false;

	// otherwise, having found the dimension, start miter above:

    // propagate DOWN dimension using min-of-coboundary
	for (miter = finddim; miter != mymap.end(); ++miter)
	{
	    //cout<<"\n\n --- BDIM-DN "<<miter->first; cin.get();
	    if (miter->first == dim) continue;
		// iterate over range map
		for (liter = miter->second->begin(); liter != miter->second->end(); ++liter)
		{
		    // extract simplex:
		    cursimp = liter->second;
		    // and loop over boundary, seeking maximal birth:
		    for (chainit = cursimp->getCB().begin(); chainit != cursimp->getCB().end(); ++chainit)
		    {
		        // if something in the coboundary is smaller or if you
		        // are uninitialized, inherit the birth!
		        if (chainit->first->birth < cursimp->birth || cursimp->birth == INITBT)
		        {
		            if (chainit->first->birth != INITBT)
		            {
		                cursimp->birth = chainit->first->birth;
		            }
		        }
		    }
		    if (cursimp->birth == INITBT)
		    {
		        cout<<"\n init error downwards "<<*cursimp;
		        cin.get();
		    }
        }
	}

	// now get reverse iterator...
	miter = mymap.find(dim);
	typename SIMP_MAP::const_reverse_iterator biter(miter);

	// Propagate DOWN using min of coboundary!
	for (; biter != mymap.rend(); ++biter)
	{
	    //cout<<"\n\n --- BDIM-UP "<<biter->first;
	    if (biter->first == dim) continue;
	    // iterate over range map
		for (liter = biter->second->begin(); liter != biter->second->end(); ++liter)
		{
		    // extract simplex:
		    cursimp = liter->second;
		    //cout<<" cur simplex: "<<*cursimp<<" bdry: "<<cursimp->getBD();
		    // and loop over boundary, seeking maximal birth:
		    for (chainit = cursimp->getBD().begin(); chainit != cursimp->getBD().end(); ++chainit)
		    {
		        // if something in the boundary is larger or if you are
		        // uninitialized, INHERIT THE BIRTH!
		        if (chainit->first->birth > cursimp->birth || cursimp->birth == INITBT)
		        {
		            if (chainit->first->birth != INITBT)
		            {
		                cursimp->birth = chainit->first->birth;
		            }
		        }
		    }
            if (cursimp->birth == INITBT)
		    {
		        cout<<"\n init error upwards "<<*cursimp;
		        cin.get();
		    }
        }
	}
    return true;
}



#endif /* TOPLEX_HPP_ */
