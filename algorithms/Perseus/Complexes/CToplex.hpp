/*
 * CToplex.hpp
 *
 *  Member functions for cubical toplexes
 */

#ifndef CTOPLEX_HPP_
#define CTOPLEX_HPP_

#include <fstream>
#include "CToplex.h"

// inserts parallel faces represented by add and <anc1, anc2>, faces of cur
template <typename C, typename PS, typename BT>
void CToplex<C,PS,BT>::insertFacePairOf(Cell<C,BT>* cur, vector<num>* add, Point<PS>* anc1,
		                                Point<PS>* anc2, const C coeff)
{
	// locates addin
	typename CUBE_MAP::iterator addfind;
	typename CUBE_RGMAP::iterator ancfind;

	// range map
	CUBE_RGMAP* rgmap;
	Cell<C,BT>* face1;
	Cell<C,BT>* face2; //parallel face pointers

	addfind = mymap.find(add);

	if (CUBFACETALK)
	{
		cout<<"\n    Search addin: < ";
		for(unsigned int i=0; i<add->size(); i++) cout<<add->at(i)<<" ";
		cout<<">";
	}

	if (addfind == mymap.end()) // no such addin exists
	{
		if (CUBFACETALK) cout<<"... not found!";
		rgmap = new CUBE_RGMAP;
		// so we should create the addin and its ranges to add to mymap

		// make new face with bdry/cobdry links
		face1 = new Cell<C,BT>(add->size()); // face with (add, anc1);
		//face1->addCBLink(cur,-1);
		cur->addBDLink(face1,coeff);
		// set up range map with anchor->face
		(*rgmap)[anc1] = face1;

		// now the parallel face:
		face2 = new Cell<C,BT>(add->size()); // face with (add, anc2);
		//face2->addCBLink(cur,1);
		cur->addBDLink(face2,-coeff);
		(*rgmap)[anc2] = face2;

		if (CUBFACETALK) cout<<"\n          Inserting: "<<*face1<<" and "<<*face2;
		// okay, insert this range map into mymap corresponding to addin
		mymap.insert(mymap.begin(), pair<vector<num>*,CUBE_RGMAP*>(add,rgmap));;
	}
	else // okay, this addin already exists at addfind!
	{
		delete add; // so we delete it...
		add = addfind->first; // and set the pointer to the mapped value
		//check if the first anchor exists as well
		ancfind = addfind->second->find(anc1);
		if (ancfind == addfind->second->end()) // not found!
		{
			// so make a new face
			face1 = new Cell<C,BT>(add->size()); // face with (add, anc1);
			// and add to range map
			(*(addfind->second))[anc1] = face1;

		}
		else // found at ancfind already
		{
			delete anc1;
			// face already exists
			face1 = ancfind->second;
		}

		// okay, now update bd/cobd links from face to cur etc.
		//face1->addCBLink(cur,-1);
		cur->addBDLink(face1,coeff);

		ancfind = addfind->second->find(anc2);
		if (ancfind == addfind->second->end()) // not found!
		{
			// so make a new face
			face2 = new Cell<C,BT>(add->size()); // face with (add, anc2);
			// and add to range map
			(*(addfind->second))[anc2] = face2;

		}
		else // found at ancfind already
		{
			delete anc2;
			// face already exists
			face2 = ancfind->second;
		}

		// okay, now update bd/cobd links from face to cur etc.
		//face2->addCBLink(cur,1);
		cur->addBDLink(face2,-coeff);
	}
}

// inserts faces for a single cube
template <typename C, typename PS,typename BT>
bool CToplex<C,PS,BT>::singleCubeFaces(typename CUBE_MAP::reverse_iterator additer, typename CUBE_RGMAP::iterator rgiter)
{
	vector<num>* curaddin = additer->first;
	Point<PS>* curanc = rgiter->first;

	Cell<C,BT>* curcube = rgiter->second;

	vector<num>* faceaddin; // addin of face
	Point<PS>* faceanch1; // anchor of face
	Point<PS>* faceanch2; // anchor of parallel face

	C coeff;

	// generate faces:
	for (num cdim = 0; cdim < 2*curcube->getDim()-1; cdim+=2)
	{
		coeff = ((cdim/2) % 2 == 0) ? -1:1;
		// CURRENT FACE
		faceaddin = new vector<num>(*(curaddin));
		// and remove cdim/2 th entry
		faceaddin->erase(faceaddin->begin()+(cdim/2));
		// and its anchor point
		faceanch1 = new Point<PS>(*(curanc));

		// make anchor point for parallel face
		faceanch2 = new Point<PS>(*(curanc));
		(*faceanch2)[curaddin->at(cdim/2)]+= ES;

		insertFacePairOf(curcube, faceaddin, faceanch1, faceanch2, coeff);
	}
	return true;
}



// takes in mymap populated by top dimensional cubes and makes all
// faces which are then inserted into the map
template <typename C, typename PS, typename BT>
void CToplex<C,PS,BT>::addCubeFaces()
{

	typename CUBE_RGMAP::iterator rgiter;
	typename CUBE_MAP::reverse_iterator additer;
	vector<num>* curaddin;

	CUBE_RGMAP* currg; // current range map

	// iterate over addin vectors backwards!
	for (additer = mymap.rbegin(); additer != mymap.rend(); ++additer)
	{
		curaddin = additer->first; // extract addin vector
		currg = additer->second; // and corresponding map

		if (curaddin->size()==0) break;

		if (CUBFACETALK)
		{
			cout<<"\n  Handling Addin: < ";
			for (unsigned int i=0; i<curaddin->size();i++) cout<<curaddin->at(i)<<" ";
			cout<<">";
			cin.get();
		}


		// now loop over the range of the map corresponding to curaddin
		for (rgiter = currg->begin(); rgiter != currg->end(); ++rgiter)
		{
			if (CUBFACETALK)
			{
				//cout<<"\n      making faces of: "<<*(rgiter->second);
				//cin.get();
			}
			// and make faces for each cube in the range of this addin
			if(!(rgiter->second->isFaced))
			{
				singleCubeFaces(additer,rgiter);
				rgiter->second->isFaced = true;
			}
		}
	}
}

template <typename C, typename PS, typename BT>
bool CToplex<C,PS,BT>::buildTopCubeMap(const num topdim, const vector<Point<PS>*>& anchs, const vector<BT>* births = NULL)
{
	// build top dim cubical map!
	if (topdim < 0) return false; // negative dimensions disallowed!

	bool dobirths = false;
	if (births != NULL)
	{
		//cout<<" bs "<<births->size()<<" vs as "<<anchs.size();
		//cin.get();
		if (births->size() >= anchs.size())
		{
			//cout<<"\n>>>>birthing!!! ";
			dobirths = true;
		}
	}

	// holds current (topdim) cube
	Cell<C,BT>* curcube;
	CUBE_RGMAP* rgmap;
	num count = 0;
	// iterate over provided points:
	typename vector<Point<PS>*>::const_iterator curpt;

	// if we are starting from scratch...
	if (mymap.size()==0)
	{
		//cout<<">>>> NEW MAP START!!"; cin.get();
		// populate first addin vector, [0,1,...,topdim-1] of size topdim
		vector<num>* topadd = new vector<num>;
		for (num dim=0; dim<topdim; dim++) topadd->push_back(dim);

		// make a map structure to feed the range of topadd:
		rgmap = new CUBE_RGMAP;
		// and finally, assign the top dimensional addin vector this range map we've built!
		mymap[topadd] = rgmap;
	}
	else // the top-dim addin etc exists, we just need to add to the map!
	{
		//cout<<">>>> OLD MAP REMAINS!!";
		// extract the (end-1)^th iterator!
		rgmap = (--mymap.end())->second;
		//cout<<" \n target addin = "; c_print<vector<num> >(*((--mymap.end())->first));
		//cin.get();
	}

	for (curpt = anchs.begin(); curpt != anchs.end(); ++curpt)
	{
		curcube = new Cell<C,BT>(topdim);
		if (dobirths) curcube->birth = births->at(count);
		if (CUBTOPTALK)
		{
			cout<<"\n added top cube: "<<*curcube;
			//cin.get();
		}
		// attempt to store anchor->cube link...
		if (!(rgmap->insert(make_pair(*curpt,curcube)).second))
		{
			// if it already exists, delete cube please...
			cout<<"\n >>>> already exists: "<<*(*curpt);
			delete curcube;
			continue;
		}
		count++;
	}
	return true;
}

template <typename C, typename PS, typename BT>
bool CToplex<C,PS,BT>::makeCellComplex(const num topd, const vector<Point<PS>*>& ancs, MComplex<C,BT>& mycomp)
{
	if (ancs.size() == 0) return false;
	buildTopCubeMap(topd, ancs);
	addCubeFaces();
	writeComplex(mycomp);
	return true;
}

template <typename C, typename PS, typename BT>
void CToplex<C,PS,BT>::writeComplex(MComplex<C,BT>& mycomp, bool clearself = false)
{
	//long count = 0;
	// iterate over addins
	typename CUBE_MAP::const_iterator myit;
	CUBE_RGMAP* currg;
	// iterate over range maps
	typename CUBE_RGMAP::const_iterator rgit;

	Cell<C,BT>* curcell;
	// loop over map, addin -> range
	for (myit = mymap.begin(); myit != mymap.end(); ++myit)
	{
		currg = myit->second;
		// loop over range, ie anchor -> cube
		for (rgit = currg->begin(); rgit != currg->end(); ++rgit)
		{
			curcell = rgit->second;
			mycomp.insertCell(curcell);

			if (curcell->getBD().size() != 2*curcell->getDim())
			{
			    cout<<"\n bdsz error: "<<*curcell;
			    cin.get();
			}
//			if (curcell->birth == 4 && curcell->dim == 0 && curcell->ind == 16)
//			{
//				cout<<"\ntraced: "<<*curcell<<" addin ";
//				c_print<vector<num> >(*(myit->first));
//				cout<<" anchor "<<*(rgit->first);
//				cin.get();
//			}
			// tracer... wtf is it?
		}
		// clean out range if done...
		if (clearself)
		{
			currg->clear();
			delete currg;
		}
	}
	if (clearself) mymap.clear();
	//cout<<"in complex..."; cin.get();
}

template <typename C, typename PS, typename BT>
void CToplex<C,PS,BT>::writeComplexExcept(MComplex<C,BT>& mycomp, const BT& except, bool clearself = false)
{
	//long count = 0;
	// iterate over addins
	typename CUBE_MAP::const_iterator myit;
	CUBE_RGMAP* currg;
	// iterate over range maps
	typename CUBE_RGMAP::const_iterator rgit;

	Cell<C,BT>* curcell;
	// loop over map, addin -> range
	for (myit = mymap.begin(); myit != mymap.end(); ++myit)
	{
		currg = myit->second;
		// loop over range, ie anchor -> cube
		for (rgit = currg->begin(); rgit != currg->end(); ++rgit)
		{
			curcell = rgit->second;
			if (curcell->birth == except) continue; // ignore banned birth time
			mycomp.insertCell(curcell);
//			if (curcell->birth == 4 && curcell->dim == 0 && curcell->ind == 16)
//			{
//				cout<<"\ntraced: "<<*curcell<<" addin ";
//				c_print<vector<num> >(*(myit->first));
//				cout<<" anchor "<<*(rgit->first);
//				cin.get();
//			}
			// tracer... wtf is it?
		}
		// clean out range if done...
		if (clearself)
		{
			currg->clear();
			delete currg;
		}
	}
	if (clearself) mymap.clear();
	//cout<<"in complex..."; cin.get();
}



// Destruction Helper, wipes out map structure!
template <typename C, typename PS,typename BT>
void CToplex<C,PS,BT>::Destroy()
{
	CUBE_RGMAP* rgmap;
	// iterate over range
	typename CUBE_RGMAP::iterator r_iter;
	// iterate over map
	typename CUBE_MAP::reverse_iterator m_iter;

	for (m_iter = mymap.rbegin(); m_iter != mymap.rend(); m_iter++)
	{
		rgmap = m_iter->second;
		// wipe out range map
		if (rgmap != NULL)
		{
			for (r_iter = rgmap->begin(); r_iter != rgmap->end(); ++r_iter)
			{
				delete r_iter->first; // delete anchor point
			}
		}
		// in any case, we clear out the map:
		rgmap->clear();

		delete rgmap;
		rgmap = NULL;
		//remove addin vector as well
		delete m_iter->first;
	}
	// finally, wipe the map clean
	mymap.clear();
}

// read in file struct, return top dimension and an anchor pt->birth time map
template <typename C, typename PS, typename BT>
pair<num,bool> CToplex<C,PS,BT>::readFileInfo(ifstream& infile, num& topdim, PB_MAP& anchors)
{
	num count = 0;
	// get top dimension of cubical set from the first line of input file
	if(!infile.eof()) infile >> topdim; else return make_pair(0,false);
	if(topdim<0) return make_pair(0,false); // negative dimension count??? i quit!
	BT birthtime;
	pair<typename PB_MAP::iterator, bool> insres;

	// POINT STUFF
	PS curcoord; // current coordinate being read from file
	vector<PS> ptcoords; // coordinates of anchor point
	Point<PS>* anchor; // anchor point
	bool breakout = false;

	while(infile.good())
	{
		ptcoords.clear(); // clear coordinates
		// get anchor point coordinates, updating extent as needed!
		for(num j = 0; j < topdim; ++j)
		{
			//cout<<"\nptrdr!";
			// get current coordinate
			infile >> curcoord;
			if (infile.eof() && j < topdim-1 )
			{
			    breakout = true;
			    break;
			}
			ptcoords.push_back(curcoord);
		}

		if (breakout)
		{
		    break;
		}

		// get current birth-time from file
		infile >> birthtime;
		if (birthtime == BANBT) continue; // ignore "banned" birthtime, as defined in Global/Pflags.h

        // kelly test
        //if (birthtime > 513) continue;

		anchor = new Point<PS>(ptcoords);
		// keep point only if unique point!
		insres = anchors.insert(make_pair(anchor,birthtime));
		if (!insres.second)
		{
			if(FLOWTALK) cout<<"\nignoring duplicated point..."<<*anchor;
			delete anchor;
		}
	}
	count = anchors.size();
	return make_pair(count,true);
}


// reads in a cubical toplex from file. the file structure is like:
// if dense,...
// <top dimension = n>
// <n coordinates determining lexico-minimal corner of cube> <birth of cube>
// <n coordinates....> <birth> etc.
// returns pair: <#top cubes read in, success?>
template <typename C, typename PS, typename BT>
pair<num,bool> CToplex<C,PS,BT>::makeFromFile(ifstream& infile)
{
	num topdim;
	PB_MAP anchors;
	// read point-birth info from the file!
	pair<num,bool> toret = readFileInfo(infile, topdim, anchors);
	//cout<<" read "; cin.get();
	if (toret.second == false) return toret;
	//cout<<"count: "<<count; cin.get();
	vector<Point<PS>*> ancvec;
	vector<BT>* births = new vector<BT>;

	typename PB_MAP::const_iterator pb;
	for(pb = anchors.begin(); pb != anchors.end(); ++pb)
	{
		ancvec.push_back(pb->first); // store point
		births->push_back(pb->second); // and its birth
	}
	// have everything, now build map!
	//cout<<" top map... "; cin.get();
	buildTopCubeMap(topdim,ancvec,births);
	// and add faces
	//cout<<" faces "; cin.get();
	addCubeFaces();
	//cout<<*this; cin.get();
	births->clear();
	delete births;
	return toret;
}

// removes all cubes except the ones with given birth time
template <typename C, typename PS, typename BT>
void CToplex<C,PS,BT>::removeAllBut(const BT& except)
{
	// iterate over cube structure...
	typename CUBE_MAP::iterator mit;
	CUBE_RGMAP* currg;
	typename CUBE_RGMAP::iterator rit;
	Cell<C,BT>* curcell; // current cell being considered

	// store iterator to rgmaps to be removed!
	vector<typename CUBE_RGMAP::iterator> tokill;
	// and an iterator over these
	typename vector<typename CUBE_RGMAP::iterator>::iterator killme;

	// loop over addin vectors
	for (mit = mymap.begin(); mit != mymap.end();++mit)
	{
		tokill.clear(); // clear out the anchor-cube pair iterator list
		currg = mit->second; // extract range map
		// loop over range map...
		for (rit = currg->begin(); rit != currg->end(); ++rit)
		{
			// get current cell
			curcell = rit->second;
			// if it has the exceptional birth time, leave it alone!
			if (curcell != NULL)
			{
				if (curcell->birth == except) continue;
			}
			// store this iterator for deletion!
			tokill.push_back(rit);
		}
		// now loop over tokill,...
		for (killme = tokill.begin(); killme != tokill.end(); ++killme)
		{
			// delete the anchor point! presumably, the actual cube is already deleted!
			delete (*killme)->first;
			// and remove the offensive elements from the current range!
			currg->erase(*killme);
		}
		// if we were being really nice, we would check addin vectors for
		// empty ranges and kill them too, but this should never be needed.
	}
}



// makes movie file... same input format as makeFromFile, but we use
// the constant MOVFRM from Global/PFlags.h to handle only so many
// frames at a time for memory conservation. Also, we want the last
// coordinate of the point to be birthtime + constant
template <typename C, typename PS, typename BT>
MComplex<C,BT>* CToplex<C,PS,BT>::readMovie(ifstream& infile)
{
	num topdim;
	PB_MAP anchors;
	// read point-birth info from file...
	pair<num,bool> toret = readFileInfo(infile, topdim, anchors);
	if (toret.second == false) return NULL;

	if (FLOWTALK || MOVTALK)
	{
		cout<<"\n Read "<<toret.first<<" top dimensional cubes from file!";
		if (MAKEBPS) cin.get();
	}

	// okay... now we have kosher top cubes, make a new map
	// from birth times to top cubes. inherit from pb-map!
	BM_MAP bmap;
	// also store uniquely sorted birth times
	set<BT> births;
	// populate new map from old, also set of births
	typename PB_MAP::const_iterator pbit;
	for (pbit = anchors.begin(); pbit != anchors.end(); ++pbit)
	{
		bmap.insert(make_pair(pbit->second,pbit->first));
		births.insert(pbit->second);
	}

	// now populate the cube structure using the alternate map...
	// in fact, we can clear out the old map first for memory reasons.

	anchors.clear();

	// declare complex
	MComplex<C,BT>* curcomp = new MComplex<C,BT>;
	MComplex<C,BT>* temp; // holds complex while copying...

	// now we feed the frames in chunks to the toplex...
	// then write and reduce, copy over, and then write and
	// reduce...
	num iterations = (num)((double)births.size() / (double)MOVFRM);
	num cur_it = 0; // current iteration level, goes from 0 to iterations...

	//if (MOVTALK) cout<<" iters: "<<iterations;

	// iterate over births
	typename set<BT>::const_iterator bit;

	// extract all points for this birth time
	typename BM_MAP::iterator low, high, trav;

	BT lowval = INITBT, highval = INITBT;
	// vectors to call buildtopcubes
	vector<Point<PS>*> ancvec;
	vector<BT>* bvec = new vector<BT>;

	// iterate over this chunk of frames...
	for (cur_it = 0; cur_it <= iterations; ++cur_it)
	{
		// get bounds on birth times for this given chunk
		lowval = highval+1; // 1, MOVFRM+1, 2*MOVFRM+1,...
		highval = (BT)(lowval + MOVFRM); // MOVFRM + 1, 2*MOVFRM + 1, ...

		if (FLOWTALK || MOVTALK)
		{
			cout<<"\n Iteration Number: "<<cur_it<<" of "<<iterations;
			cout<<"\n from "<<lowval<<" to "<<highval;
			if (MAKEBPS) cin.get();
		}

		low = bmap.lower_bound(lowval);
		high = bmap.upper_bound(highval);

		// now traverse between these values, feeding anchors and births
		bvec->clear(); vector<BT>(*bvec).swap(*bvec);
		ancvec.clear(); vector<Point<PS>*>(ancvec).swap(ancvec);

		for (trav = low; trav != high; ++trav)
		{
			bvec->push_back(trav->first); // populate birth
			ancvec.push_back(trav->second); // populate anchor
			//if (MOVTALK) cout<<"\n .... added: b="<<trav->first<<", anc="<<*(trav->second);
			//cin.get();
		}
		// memory savers:
		vector<BT>(*bvec).swap(*bvec);
		vector<Point<PS>*>(ancvec).swap(ancvec);

		// at this point, we have birth and anchors, so build top cube map!
		buildTopCubeMap(topdim,ancvec,bvec);
		//cout<<*this;
		//cin.get();
		// and make faces...
		addCubeFaces();
		// have the CToplex filled for the current chunk of frames, so...
		// 1. write it to complex, but leave out lowval - 1, which already exists (or is empty)
		writeComplexExcept(*(curcomp),lowval-1);

		if (FLOWTALK)
		{
			cout<<"\n------ Complex Size: "<<curcomp->size();
		}
		// remove all but the highval-birthtimed cells from CToplex for next round of insertion...
		removeAllBut(highval);
		// reduce the complex shamelessly, except leave the highval-1^th frame alone!
		curcomp->MorseCoreduceExcept(lowval-1,highval,false);
		// now copy over to temp...
		temp = new MComplex<C,BT>;
		curcomp->moveOver_Sort(*temp,false); // copy to temp... use coreduce
		// at this point temp holds all the info we want, so eliminate curcomp without
		// deleting cells
		curcomp->clear();
		delete curcomp;
		// now make curcomp point to temp, and we are set for the next iteration!
		curcomp = temp;
		if (FLOWTALK)
		{
			cout<<"\n------ Reduced Size: "<<curcomp->size();
		}
//		// first, delete point information bvec except at the last level
//		for (num i = 0; i < (num)bvec->size(); ++i)
//		{
//			if (bvec->at(i) != highval) delete ancvec.at(i);
//		}
		//curcomp->printSizeInfo();
		// also clear out births and anchors!
	}
	bvec->clear(); vector<BT>(*bvec).swap(*bvec);
	delete bvec; // delete memory allocated for birth vectors
	// finally, return reduced complex of movie!
	return curcomp;
}
#endif /* CTOPLEX_HPP_ */
