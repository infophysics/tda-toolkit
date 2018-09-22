/*
 * RIPS.hpp
 *
 *      Author: Vidit
 */


#ifndef RIPS_HPP_
#define RIPS_HPP_

#include <algorithm> // for set intersections!
#include <cmath> // for ceil and floor etc.
#include "RIPS.h" // for class definition


// point-birth pairs from file:
# define PB_PAIR map<Point<PS>*, pair<double,BT>, PT_ORDER >


// inserts vertex into rips complex and returns the index i of that vertex.
// to access this vertex, please use vert[i].
template <typename C, typename PS, typename BT>
num RIPS<C,PS,BT>::addVertex(const BT& btime)
{
    vector<PS> coord;
    coord.push_back((PS)vert.size());
	Point<PS>* mypt = new Point<PS>(coord); // create point
	//cout<<"\n adding "<<mypt;
	vert.push_back(mypt); // insert to vert structure
	birth.push_back(btime); // allocate birth time
	//cout<<"\n size of vert is "<<vert.size();
	return (num)(vert.size()-1);
}

// allocates neighbor upper triangular matrix. call after all the
// vertices have inserted using addVertex()
template <typename C, typename PS, typename BT>
void RIPS<C,PS,BT>::allocateNbrMatrix()
{
	this->psize = vert.size();
	//this->topdim = 1;
	//vector<BT>* drow;
	for (num i = 1; i < (num)vert.size(); i++)
	{
		isnbr.push_back(NULL);
	}
	if (RIPSTALK || RIPSNBRTALK) cout<<"\n     nbr matrix nulled with "<<vert.size()-1<<" rows ";

}

// inserts an edge of birth time "btime" between vertices of given
// index pt1Index and pt2index into "vert"
template <typename C, typename PS, typename BT>
bool RIPS<C,PS,BT>::addEdge(num pt1, num pt2, const BT& btime)
{
    BT born = btime;
	if (pt1 == pt2) return false; // no diagonal entries in adjacency matrix!
	// first order such that pt2index is LESS than pt1index
	if (pt2 > pt1)
	{
		// swap!
		num temp = pt1;
		pt1 = pt2;
		pt2 = temp;
	}
	// now check bounds!
	 if (pt2 < 0 || pt1 >= (num)vert.size()) return false;

	// if you get here then you can just change isnbr entry as needed!
	if (isnbr.at(pt1-1) == NULL)
	{
	    //cout<<"\n null row encountered for vertex "<<pt1;
	    isnbr.at(pt1-1) = new RROW;
	}

	BT vertmax = max(birth.at(pt1), birth.at(pt2));
	if (born < vertmax) born = vertmax;

	// add pt2, btime to notate the edge (pt1, pt2)
	isnbr.at(pt1-1)->push_back(make_pair(pt2,born));
	//cout<<"\n added edge from "<<pt1<<" to "<<pt2<<" born "<< born;
	return true;
}

// returns persistence intervals. call AFTER allocatenbrmatrix.
template <typename C, typename PS, typename BT>
void RIPS<C,PS,BT>::ComputePersistence(num upto, map<num, vector<pair<BT,BT> > >& pints)
{

	incRIPS(upto); // create rips complex and feed to stoplex...

	MComplex<C,BT> ccomp;
	this->writeComplex(ccomp);


	ccomp.MorseWrapper_Cored(false);
	PComplex<C,BT> pcomp;
	pcomp.COMPUTE_INTERVALS(ccomp);

	pints = pcomp.getInts();
}


// once neighbor matrix has been filled by calling addEdge(),
// call this function! We restrict the dimension of the complex
// using ptdim.
template <typename C, typename PS, typename BT>
void RIPS<C,PS,BT>::ComputePersistence(num upto, string fname = "output", bool printcomp=false)
{
	//cout<<*this;

	incRIPS(upto); // create rips complex and feed to stoplex...

	MComplex<C,BT>* ccomp = new MComplex<C,BT>;
	this->writeComplex(*ccomp);

    if(printcomp)
    {
        ofstream cfile("complex.txt");
        cfile << *ccomp;
        cfile.close();
    }

    //ccomp.printSizeInfo();
	// kill point-allocated memory!
    AlternateAndUpdate(ccomp);

	//cout<<ccomp;
	// reduce the complex
	//ccomp.MorseWrapper_Cored(false,0.2);
	// make persistence complex!
	PComplex<C,BT> pcomp;
	pcomp.COMPUTE_INTERVALS(*ccomp); // compute persistence intervals
	//pcomp.showBetti();
	pcomp.makeOutputFiles(fname); // write output!
	delete ccomp;
}

// once neighbor matrix has been filled by calling addEdge(),
// call this function! We restrict the dimension of the complex
// using ptdim.
template <typename C, typename PS, typename BT>
void RIPS<C,PS,BT>::ComputePersistence_GrowBalls(num upto, string fname = "output")
{
	//cout<<*this;

	incRIPS_GrowBalls(upto); // create rips complex and feed to stoplex...
	// create morse complex from rips complex
	//cout<<*this;

	MComplex<C,BT> ccomp;
	this->writeComplex(ccomp);

	// kill point-allocated memory!

	//cout<<ccomp;
	// reduce the complex
	ccomp.MorseWrapper_Cored(false,0.0);
	// make persistence complex!
	PComplex<C,BT> pcomp;
	pcomp.COMPUTE_INTERVALS(ccomp); // compute persistence intervals
	//pcomp.showBetti();
	pcomp.makeOutputFiles(fname); // write output!
}

// generates random rips complex by edge shuffling the edge births of
// the complete graph on n vertices. stops short when a claw is
// detected!
template <typename C, typename PS, typename BT>
void RIPS<C,PS,BT>::makeRandom_EdgeShuffle(num n, bool makeEdgefile=false)
{
        // add n vertices: labeled 0 through n-1, born at time 1.
        for (num i=0; i<n; ++i)
        {
            addVertex(1);
        }
        // allocate space for distances
        allocateNbrMatrix();

        num ecount = n*(n-1)/2; // get total edge count...
        // initialize random number generator
        srand(time(NULL));
        // insert edge labels:
        vector<num> edges;
        for (num i=0; i<ecount; ++i) edges.push_back(i);
        vector<num> coords;
        ofstream edgefile("edgeloc.txt");

        std::random_shuffle(edges.begin(),edges.end()); // permute edges...

        vector<num> degree(n,0); // n vertices get degree 0
        num eins = 0; // actual number of edges inserted
        // now add edges to rips complex:
         for (num pos = 0; pos < ecount; ++pos)
        {
            eins++;
            lexicoPos(n,2,edges.at(pos)+1,coords);
            addEdge(coords.at(0),coords.at(1),pos+2);
            // update and check degrees to see if we have a claw
            if (++(degree.at(coords.at(0))) == n-1) break;
            if (++(degree.at(coords.at(1))) == n-1) break;
            if (makeEdgefile)
            {
                edgefile <<pos+2 <<". ("<< coords.at(0)<<", "<<coords.at(1)<<")\n";
            }
        }
 /*       cout<<"\n Density: "<<((double)(eins))/((double)(ecount))<<"\n Continue?";
        char ans;
        cin >> ans;
*/
        edgefile.close();
        this->topdim = n;
}


// makes RIPS Complex out of a file containing the following point cloud data:
// <point dimension, n>
// <initial radius multiplier, \epsilon_0> <step size \delta\epsilon> <# steps>
// x_1 x_2 ... x_n r_1 where x = (x_1,...,x_n) is the point, r its initial radius
// y_1, y_2 ... etc.
// if comred is enabled then all points have a common initial radius at the beginning
// of the file
template <typename C, typename PS, typename BT>
pair<num,bool> RIPS<C,PS,BT>::makeFromFile_GrowBalls(ifstream& infile,
                            bool comred = false, bool capped=false, bool witness=false)
{
	map<Point<PS>*, double, ptcomplex<PS> > prmap; // point-radius pairs
	numsteps = -1; // test file input

    // get point dimension
	num ptdim = 0;
	if (infile.good())
	{
		infile >> ptdim;
	}
	else
	{
		return make_pair(0,false); // error, couldn't open file for reading
	}

	if (ptdim <= 0)
	{
		return make_pair(0,false); // bad dimension!
	}

	double radius = 0; num dimcap = ptdim;

	if (infile.good())
	{
		infile >> epsilon >> stepsize >> numsteps;
		if (comred) infile >> radius;
		if (capped) infile >> dimcap;
	}

	// verify positivity of radius multiplier and steps etc:
	if (epsilon < 0 || stepsize <= 0 || numsteps < 0)
	{
		//cout<<"epsposfail: "<<epsilon<<" "<<stepsize<<" "<<numsteps<<"\n";
		return make_pair(0,false);
	}

	vector<PS> ptcoords; // coordinates of current point
	PS coord;  // current coordinate and radius
	Point<PS>* curpt; // current point
	num flines = 0;

	srand(time(NULL));
    bool breakout = false;
	// read in point and radius!

	while (infile.good())
	{
		flines++;
		ptcoords.clear();
		for (num i=0; i<ptdim; i++)
		{
			infile>>coord;
			//cout<<"\n coordinate read: "<< coord;
			if (infile.eof() && i != ptdim - 1)
			{
			    breakout = true; // premature exit!
			    break;
			}
			ptcoords.push_back(coord);
		}

        if (breakout) break;

		// input radius for this point!
		if (!comred)
		{
			infile >> radius;
			// ensure positivity of radius
			if (radius < 0)	continue;
		}

		// create a point with these coordinates
		curpt = new Point<PS>(ptcoords);

		// try to insert current point and radius
		if (! (prmap.insert(make_pair(curpt, radius))).second )
		{
			// point already exists, so delete
			delete curpt;
		}
	}

	// now iterate over map of point-radius pairs and populate vert, rads, birth...
	typename map<Point<PS>*, double, ptcomplex<PS> >::const_iterator pr;
	for(pr = prmap.begin(); pr != prmap.end(); ++pr)
	{
		vert.push_back(pr->first);
		rads.push_back(pr->second);
	}

    if (FLOWTALK)
	{
	    cout<<"\n #points read: "<<vert.size();
	    if (MAKEBPS) cin.get();
	}

    // witness with delta = stepsz/3
    if (witness) this->witness(stepsize*4);

    if (FLOWTALK && witness)
    {
        cout<<", and #witnesses: "<<vert.size();
        if (MAKEBPS) cin.get();
    }

    this->psize = vert.size();
	this->topdim = ptdim;

	makeNbrMatrix_GrowBalls();

	if (FLOWTALK)
	{
	    cout<<"\nnbr matrix made: ";
	    //showNbrMatrix();
	    if (MAKEBPS) cin.get();
	}
	//cout<<"\n nbr matrix filled! making rips complex"; cin.get();
	incRIPS_GrowBalls(dimcap, 1+INITBT);
	//cout<<*this;
	return make_pair(vert.size(), true);
}

// modifies the vertex set so that no two points are less than delta apart!
template <typename C, typename PS, typename BT>
num RIPS<C,PS,BT>::witness (double delta)
{
    if (FLOWTALK)
    {
        cout<<"\n    witnessing with radius: "<<delta;
        if (MAKEBPS) cin.get();
    }
    // make sure no pruning is attempted for bad delta values
    if (delta <= 0) return 0;

    VERTS witness; // create witness vertex list... initially empty
    Point<PS>* curvert; // current vertex
    Point<PS>* curwit; // current witness

    vector<bool> tokeep; // the n-th flag decides whether the n-th vertex is a witness.

    // iterate over the vertices...
    bool quit;
    typename VERTS::const_iterator vit, wit;
    for (vit = vert.begin(); vit != vert.end(); ++vit)
    {
        if (RIPSWITTALK && (vit - vert.begin()) % 1000 == 0)
        {
            cout<<"\n    >>>> test vert num: "<<vit - vert.begin();
            cout<<" against witlist: "<<witness.size();
            if (MAKEBPS) cin.get();
        }
        quit = false;
        // extract current vertex
        curvert = *vit;
        // iterate over witness vertices
        for (wit = witness.begin(); wit != witness.end(); ++wit)
        {
            curwit = *wit;
            //cout<<"\n dist: "<<dist(*curvert,*curwit); cin.get();
            if (dist(*curvert,*curwit) < delta) // if the current vertex is witnessed,...
            {
                // flag for deletion!
                tokeep.push_back(false);
                // and move on to the next vertex.
                quit = true;
                break;
            }
        }
        if (quit) continue;
        // if we reached here, the current vertex was not witnessed.
        witness.push_back(curvert);
        tokeep.push_back(true);
    }

    vector<double> witrads; // radii of witness points
    // now, we should delete all those vertices which are not witnesses.
    for (num i=0; i < (num)tokeep.size(); ++i)
    {
        if (! tokeep.at(i)) delete vert.at(i);
        else witrads.push_back(rads.at(i));
    }
    // finally, replace vert by witness...
    vert.assign(witness.begin(),witness.end());
    rads.assign(witrads.begin(),witrads.end());

    return (num)vert.size();
}

template <typename C, typename PS, typename BT>
pair<num,bool> RIPS<C,PS,BT>::makeFromDistMatrixFile(ifstream& infile, bool iscormat=true, bool isuptr = false)
{
    num numpts;
    if (!infile.eof())
	{
		infile >> numpts;
	}
	else
	{
		//cout<<"initopen\n";
		return make_pair(0,false); // error, couldn't open file for reading
	}

    num maxdim; // how many dimensional simplices to build?
	// read starting distance epsilon, step size and number of steps
	if (!infile.eof())
	{
	    infile >> epsilon >> stepsize >> numsteps >> maxdim;
        if (epsilon < 0 || stepsize <= 0 || numsteps < 0 || maxdim < 0) return make_pair(0,false);
	}

	if (DISTMATTALK)
    {
        cout<<"\n eps: "<<epsilon<<" stepsize: "<<stepsize<<" #steps: "<< numsteps<<" dimcap: "<<maxdim;
        if (MAKEBPS) cin.get();
    }

	// a row of birth times
	//vector<BT>* drow;
	double curdist;
    num jend = numpts; // ending value for inner loop, change in case of upper triangular matrix!

    isnbr.clear();

	// loop over file
	for (num i = 0; i < numpts; i++)
	{
	    // create a vertex born at time 0:
	    addVertex(INITBT);
        // and null out its neighbor row
        isnbr.push_back(NULL);

        num cursteps;

	    if (isuptr) jend = numpts - i;
	    for (num j = 0; j < jend; j++)
	    {
	        // read distance from file
	        infile >> curdist;
	        //cout << "\n read: "<<curdist<<" at "<<i<<", "<<j;

	        if (!isuptr && j>=i) continue; // if not upper triangular, this is useless

            //cout<<" KEPT!!";
	        // if these are correlations, we probably want the distance to be 1 - (file val)!
	        if (iscormat) curdist = 1 - curdist;

	        // okay, now we have a distance. how does this relate to birth time?
            // well, if dist is less than initial epsilon... consider this edge born at time 0
            if (curdist < epsilon)
            {
                if (isnbr.at(i) == NULL)
                {
                    isnbr.at(i) = new RROW;
                }
                if (DISTMATTALK) cout<<" dist: "<<curdist<<"at i="<<i<<" j="<<j<<" and eps: "<<epsilon;
                //cin.get();
                addEdge(i,j,INITBT);
                //isnbr.at(i)->push_back(make_pair(j,INITBT));
            }
            else // if dist is larger than epsilon,
            {
                //compute number of steps needed to make it equal epsilon
                cursteps = (num)ceil((curdist - epsilon)/stepsize);
                if (DISTMATTALK) cout<<"\n   cursteps: "<<cursteps;

                if (cursteps <= numsteps)
                {
                    if (isnbr.at(i) == NULL)
                    {
                        isnbr.at(i) = new RROW;
                    }
                    //isnbr.at(i)->push_back(make_pair(j,cursteps + INITBT));
                    addEdge(i,j,cursteps + INITBT);
                    if (DISTMATTALK) cout<<"\n success nbr: "<<j<<" with birth "<<cursteps + INITBT;
            /*      if (cursteps == 0)
                    {
                        cout<<"\n cstep zero with eps"<<epsilon<<" curd "<<curdist<<" at (i,j) = "<<i<<" "<<j;
                    }
            */
                }
                //else drow->push_back(BANBT);
            }

	    }
	    //if (i > 0) isnbr.push_back(drow);
	}
	this->psize = numpts;
	// arbitrary topdim!
	this->topdim = maxdim;
	//cout<<"\n #file points: "<<flines; cin.get();
    if (DISTMATTALK) showNbrMatrix();
        //cout<<"\n nbr matrix filled! making rips complex"; cin.get();
	incRIPS(maxdim);
	if (DISTMATTALK) cout<<*this;
	return make_pair(numpts,true);
}

template <typename C, typename PS, typename BT>
void RIPS<C,PS,BT>::makeNbrMatrix_GrowBalls()
{
	RROW* drow; // row of point-birth pairs for current point

	double dratio;
	BT cursteps;

	//num minind1, minind2;
	//double mindist = 0;


	for (num i = 1; i < (num)vert.size(); i++)
	{
		drow = NULL;
		for(num j=0; j < i; j++)
		{
			cursteps = 1+INITBT;
			if (RIPSNBRTALK)
			{
				cout<<"\n dist: "<<*(vert[i])<<" to "<<*vert[j]<<" is "
				<<dist(*(vert[i]),*(vert[j]))<<" vs. eps-rad-sum: "
				<<(rads[i]+rads[j])*epsilon;
				//cin.get();
			}
/*			// initial seed values
			if (i == 1 && j == 0)
			{
				minind1 = j;
				minind2 = i;
				mindist = dist(*(vert[i]),*(vert[j]));
			}
			else // finding the minimum
			{
				if (dist(*(vert[i]),*(vert[j])) < mindist)
				{
					//cout<<"\n new min dist: "<<dist(*(vert[i]),*(vert[j]))<<" vs "<<mindist;
					minind1 = j;
					minind2 = i;
					mindist = dist(*(vert[i]),*(vert[j]));
				}
			}
*/
			// compute #iterations needed to make i^th and j^th vertices neighbors...

			dratio = (0.5*(1/stepsize))*((dist(*(vert.at(i)),*(vert.at(j))) / epsilon)
                      - (rads.at(i)+rads.at(j))) ;
            cursteps = (BT) (ceil(dratio));
            if (cursteps < INITBT) cursteps = INITBT;

			if (RIPSNBRTALK)
			{
				cout<<"\n dratio: "<<dratio<<" cursteps: "<<cursteps;
				cin.get();
			}
			// is edge birth time larger than the max allowed?
			if (cursteps <= numsteps)
			{
				if (RIPSNBRTALK) {cout<<"... YES, born at ["<<cursteps + INITBT<<"]\n";}
				if (drow == NULL)
				{
				    drow = new RROW;
				}
				drow->push_back(make_pair(j,cursteps + INITBT));
			}
/*			else
			{
				if (RIPSNBRTALK) {cout<<"... NO!\n";}
				drow->push_back(BANBT);
			}
*/		}

		isnbr.push_back(drow);
	}

	//cout<<"\n\n MIN RIPS DIST = "<<mindist<<" BETWEEN "<<*(vert.at(minind1))<<" at "<<minind1
	//				<<" and "<<*(vert.at(minind2))<<" at "<<minind2; cin.get();

		// display matrix structure
	if (RIPSNBRTALK) showNbrMatrix();
}

template <typename C, typename PS, typename BT>
void RIPS<C,PS,BT>::showNbrMatrix(ostream& out = cout) const
{
    num ptind;
    BT bir;
	out<<"\n\n++++++ RIPS Distance matrix ++++++++++++";
	for (num i = 0; i < (num)isnbr.size(); i++)
	{
		out<<"\n";
		//cout<<"\n i"; cin.get();
		if (isnbr.at(i) == NULL)
		{
		    cout<<"\n no nbrs for vertex "<<i+1;
		    continue;
		}
		cout<<"\n from vertex "<<i+1<<" to: ";
		for(num j = 0; j < (num)isnbr.at(i)->size(); j++)
		{
		    ptind = isnbr.at(i)->at(j).first;
		    bir = isnbr.at(i)->at(j).second;
			out<<"\n ~~~~~ vertex "<<ptind<<", birth "<<bir;
		}
	}
}

template <typename C, typename PS, typename BT>
void RIPS<C,PS,BT>::incRIPS_GrowBalls(const num topdim, const BT& start = INITBT)
{
	INSRIP dimsimps; // simplices sorted by dimension (>= 1)

	vector<num> nblist; // list of lower neighbors
	vector<BT> ebirths; // and their edge birth times

	Simplex<C,PS,BT>* zerosimp; // 0-simplex corresponding to this point!
	for(num cur = 0; cur < (num)vert.size(); ++cur)
	{

		if (BIGRIPSTALK || RIPSTALK )
		{
			cout<<"\nPOINT: "<<*(vert.at(cur))<<" at "<<cur;
			if (MAKEBPS) cin.get();
		}

		getLowerNbrs_GrowBalls(cur, nblist, ebirths, true); // extract neighbor list and edge births
		zerosimp = new Simplex<C,PS,BT>(vert.at(cur)); // make 0-simplex out of current vertex...
		zerosimp->birth = start;

		AddRIPSCofaces_GrowBalls(zerosimp, nblist, ebirths, topdim, dimsimps); // now add all cofaces
	}
	// inserting into toplex structure from rips info
	//cout<<" final toplex insertion... "; cin.get();
	insertToToplex(dimsimps);
}

template <typename C, typename PS, typename BT>
bool RIPS<C,PS,BT>::AddRIPSCofaces_GrowBalls(Simplex<C,PS,BT>* tau, const vector<num>& nblist,
		                                    const vector<BT>& ebirths, num topdim, INSRIP& dimsimps)
{
	if (tau == NULL) return false;
	if (RIPSTALK)
	{
		cout << "\n                     coface adder";
		cout<<" with simp: "<<*tau;
		cout<<"\n nblist ";
		c_print<vector<num> >(nblist);
		//showPoints(nblist);
		cin.get();
	}
	typename SIMP_MAP::const_iterator finddim;
	// add tau to complex... including all boundary links
	if (tau->getDim() == 0) // if this simplex is just a point...
	{
		typename SIMP_MAP::iterator finddim = this->mymap.find(0);
		SIMP_RG* currg;

		// Insert point directly into the toplex structure
		if (finddim == this->mymap.end()) // found dimension 0?
		{
			// no, make range
			currg = new SIMP_RG;
			this->mymap.insert(this->mymap.end(),make_pair(0,currg));
		}
		else
		{
			// yes, get range
			currg = finddim->second;
		}
		// make new cell corresponding to tau's vertex list!
		currg->insert(make_pair(tau->getHash(this->getHashMod()),tau));
	}
	else // higher dimension, add tau to the dimsimps structure!
	{
		storeCoface(tau, dimsimps);
	}
	// if tau was top-dimensional or has no neighbors, exit!!
	if (tau->getDim() >= topdim || nblist.empty()) return true;
	// Simplex<C,PS,BT>* sigma; // co-face of tau

	// otherwise, iterate over nblist, get vertex v
	typename vector<num>::const_iterator niter;
	typename vector<num>::iterator intend;

	Simplex<C,PS,BT>* sigma;
	PT_SET siglist; // store vertices of tau in siglist
	vector<num> intersects, vnbrs;
	vector<BT> dummy; // to call getlowernbrs

	Point<PS>* v; // current vertex under consideration
	num nbcount = 0; // neighbor count

	bool getbirth = (tau->getDim() == 0) ? true : false;

	for (niter = nblist.begin(); niter != nblist.end(); ++niter)
	{
		intersects.clear();
		//cout<<" ind: "<<*niter; cin.get();
		siglist = tau->getVerts();
		v = vert.at(*niter); //extract vertex v from nblist
		if ((siglist.insert(v)).second) // augment vertices of tau with this
		{
			sigma = new Simplex<C,PS,BT>(siglist);
			getLowerNbrs_GrowBalls(*niter,vnbrs,dummy,false);
			if (getbirth) sigma->birth = ebirths.at(nbcount++);
			// store intersection of vnbrs and nblist in intersects
			set_intersection(nblist.begin(),nblist.end(),vnbrs.begin(),vnbrs.end(),
					         inserter(intersects,intersects.begin()));

			if (INTRIPSTALK)
			{
				cout<<"\n***********intersect***************\n";
				c_print<vector<num> >(nblist); cout<<" and "; c_print<vector<num> >(vnbrs); cout<<" gives "; c_print<vector<num> >(intersects);
				cout<<"\n*************************************";
			}
			// and recursively call addcofaces!
			AddRIPSCofaces_GrowBalls(sigma, intersects, dummy, topdim, dimsimps);
		}
		else return false;
	}
	return true;
}



template <typename C, typename PS, typename BT>
void RIPS<C,PS,BT>::getLowerNbrs_GrowBalls(num index, vector<num>& nbinds, vector<BT>& ebirths, bool dobirths = true) const
{
	nbinds.clear(); // make the list of lower neighbors empty!
	ebirths.clear();
	if (index >= 1 && index < (num)vert.size())
	{
		//num lowind = 0;
		//cout<<"\n\nLOWER nbr at "<<index; cin.get(); cout<<"of "<<*(vert[index]);
		//cout<<"\n"<<isnbr.at(index)->size(); vprint<bool>(*(isnbr.at(index)));
        RROW* currow = isnbr.at(index-1);

		//typename RIPMAT::const_iterator lowerpt; // iterate over preceding points
		typename RROW::const_iterator rit;

		//cout<<"\nnbr cons: "<<lowind+1; cin.get(); cout<<" i.e. "<<*(vert[lowind]);
		if(currow != NULL)
		{
		    for (rit = currow->begin(); rit != currow->end(); ++rit)
		    {
                   nbinds.push_back(rit->first); // neighboring point
                   if (dobirths) ebirths.push_back(rit->second); // birth time of corresponding edge to lower nbr
		    }
        }
	}
}


template <typename C, typename PS, typename BT>
void RIPS<C,PS,BT>::makeNbrMatrix()
{
	RROW* drow;
	for (num i = 1; i < (num)vert.size(); i++)
	{
		//cout<<"\ni: "<<i;
		drow = NULL;
		for(num j=0; j < i; j++)
		{
			if (RIPSNBRTALK)
			{
				cout<<"\n   dist: "<<*(vert[i])<<" to "<<*vert[j]<<" is "
					<<dist(*(vert[i]),*(vert[j]))<<" vs. rad-sum: "
					<<(rads[i]+rads[j]);
				//cin.get();
			}

			if (dist(*(vert.at(i)),*(vert.at(j))) < (rads.at(i)+rads.at(j)))
			{
				if (RIPSNBRTALK) {cout<<"... YES!\n";}
				// set birth to default value... this will change when we insert
				// to the simplicial toplex if a boundary vertex happens to
				// have a higher birth time!
				if (drow == NULL)
				{
				    drow = new RROW;
				}
				drow->push_back(make_pair(j,INITBT));
			}
			else
			{
				if (RIPSNBRTALK) {cout<<"... NO!\n";}
				//drow->push_back(BANBT);
			}
		}
		isnbr.push_back(drow);
	}
	// display matrix structure
	if (RIPSNBRTALK)
	{
		showNbrMatrix();
	}
}

// creates rips complex from file with the following format:
// <dimension of points> n
// x_1 x_2 ... x_n r_1 b_1 where x = (x_1,...,x_n) is the point, r its radius, and b its birth-time
// y_1 y_2 etc

template <typename C, typename PS, typename BT>
pair<num,bool> RIPS<C,PS,BT>::makeFromFile(ifstream& infile)
{
	// get point dimension
	num ptdim = 0;
	if (!infile.eof()) infile >> ptdim;
	else return make_pair(0,false); // error, couldn't open file for reading

	if (ptdim <= 0) return make_pair(0,false); // bad dimension!

	PB_PAIR ptbirth; // point - birth map

	// okay... now prepare to read point info, insert birthtime, and create radius vectors!
	vector<PS> ptcoords;
	PS coord; double radius;
	BT curborn;
	Point<PS>* curpt;
	while(infile.good())
	{
		ptcoords.clear();
		for (num i=0; i<ptdim; i++)
		{
			infile>>coord;
			ptcoords.push_back(coord);
		}
		//cout<<"\n... "<<*curpt;
		// input radius for this point!
		infile >> radius;
		// ensure positivity of radius
		if (radius < 0)	continue;
		infile >> curborn;
		if (curborn == BANBT) continue;
		// create a point with these coordinates
		if (infile.eof()) break;
		curpt = new Point<PS>(ptcoords);

		// try to insert...
		if (! (ptbirth.insert(make_pair(curpt,make_pair(radius,curborn)))).second )
		{
			// if point already exists, delete current version
			delete curpt;
		}

	}

	// iterate over map and feed into vert, birth, rads etc.
	typename PB_PAIR::const_iterator pbit;
	for (pbit = ptbirth.begin(); pbit != ptbirth.end(); ++pbit)
	{
		vert.push_back(pbit->first);
		rads.push_back(pbit->second.first);
		birth.push_back(pbit->second.second);
    }
	// eradicate ptbirth:
	ptbirth.clear();
	this->psize = vert.size();
	this->topdim = ptdim;

	// populate RIPMAT structure...
	makeNbrMatrix();
	//cout<<" incrips call"; cin.get();
	// build RIPS complex incrementally.
	incRIPS(ptdim);

	//cout<<"ending make from file!"; cin.get();

	return make_pair(vert.size(),true);
}

// returns all neighbors of a given point (determined by
// its index in "vert") that precede it in the ordering!
template <typename C, typename PS, typename BT>
void RIPS<C,PS,BT>::getLowerNbrs(num index, vector<num>& nbinds) const
{
	//BT maxborn = BANBT;
	nbinds.clear(); // make the list of lower neighbors empty!
	if (index >= 1 && index < (num)vert.size())
	{
		RROW* currow = isnbr.at(index-1);
        typename RROW::const_iterator rit;

		if(currow != NULL)
		{
		    for (rit = currow->begin(); rit != currow->end(); ++rit)
		    {
                   nbinds.push_back(rit->first); // neighboring point
                   //if (dobirths) ebirths.push_back(rit->second); // birth time of corresponding edge to lower nbr
		    }
		}
	}
	// try to sort:
	sort(nbinds.begin(),nbinds.end());
}




// builds RIPS complex incrementally, stopping at dimension D
template <typename C, typename PS, typename BT>
void RIPS<C,PS,BT>::incRIPS(const num upto)
{
	INSRIP dimsimps; // simplices sorted by dimension (> 1)

	vector<num> nblist; // list of lower neighbors
	Simplex<C,PS,BT>* zerosimp; // 0-simplex corresponding to this point!
	for(num cur = 0; cur < (num)vert.size(); ++cur)
	{

		if (BIGRIPSTALK || RIPSTALK )
		{
			cout<<" showing pt: "; //cin.get();
			cout<<"\nPOINT: "<<(vert.at(cur))<<" at "<<cur;
			//cin.get();
			if (MAKEBPS) cin.get();
		}

		getLowerNbrs(cur, nblist); // extract neighbor list
		zerosimp = new Simplex<C,PS,BT>(vert.at(cur)); // make 0-simplex out of current vertex...
		zerosimp->birth = birth.at(cur);
		AddRIPSCofaces(zerosimp, nblist, upto, dimsimps, cur); // now add all cofaces
	}
	//cout<<" done recursive crap!"; cin.get();
	insertToToplex(dimsimps);
}

// DEFUNCT: THE ISNBR STRUCTURE IS NOW UNFRIENDLY FOR SEARCHING.
// uses the neighbor matrix to extract birth time for edge between two vertices!
template <typename C, typename PS, typename BT>
BT RIPS<C,PS,BT>::getEdgeBirth(num pt1, num pt2) const
{
	if (pt1 == pt2) return BANBT; // no diagonal entries in adjacency matrix!
	// first order such that pt2index is LESS than pt1index
	if (pt2 > pt1)
	{
		// swap!
		num temp = pt1;
		pt1 = pt2;
		pt2 = temp;
	}
	// now check bounds!
	if (pt2 < 0 || pt1 >= (num)vert.size()) return BANBT;
    if (isnbr.at(pt1-1) == NULL) return BANBT;

	// loop over the nbrs of pt1...
	for (num j=0; j < (num)isnbr.at(pt1-1)->size(); j++)
	{
	    if (isnbr.at(pt1-1)->at(j).first == pt2)
	    {
	        return isnbr.at(pt1-1)->at(j).second;
	    }
	}

	// if you get here then you can just return isnbr entry
	return (BT) BANBT;
}


// builds RIPS complex incrementally, stopping at dimension D
// assigns the same birth time epsilon to all points, and uses
// the output of getlowernbrs to assign birth times to
// edges.

template <typename C, typename PS, typename BT>
void RIPS<C,PS,BT>::insertToToplex (INSRIP& dimsimps)
{
	// insert into simplicial toplex from dimsimps and clear memory!
	typename INSRIP::iterator dim;
	SLIST* curlist;
	typename SLIST::iterator simp;
	num count = 0;

	// iterate over dimensions
	for (dim = dimsimps.begin(); dim != dimsimps.end(); ++dim)
	{
		//cout<<" insert dimension: "<<dim->first; cin.get();
		curlist = dim->second; // list of simplices of this dimension
		// iterate over this list:
		for (simp = curlist->begin(); simp != curlist->end(); ++simp)
		{
			if (DIMSIMPSTALK)
			{
				cout<<"\n::::::: INSERT ATTEMPT: "<<**simp;
				if (MAKEBPS) cin.get();
			}
			// add to complex!
			//cout<<"adding as coface..."; cin.get();
			this->addCoface(*simp);
			//cout<<"added!"; cin.get();
			if (RIPSTALK)
			{
				cout<<"\n inserted at "<<(*simp)->birth;
			}
			count++;
		}
		// now we have inserted everything from this list, so...
		curlist->clear();
		delete curlist;
	}
	// finally, clear out dimsimps:
	dimsimps.clear();
}


template <typename C, typename PS, typename BT>
void RIPS<C,PS,BT>::storeCoface(const Simplex<C,PS,BT>* tau, INSRIP& dimsimps) const
{
	if (RIPSTALK) cout<<"\n --------- storing: "<<*tau;

	num insdim = tau->getDim();
	// do we already have simplices of this dimension?
	typename INSRIP::iterator finddim = dimsimps.find(insdim);
	SLIST* inslist;
	// if not, make the corresponding list:
	if (finddim == dimsimps.end())
	{
		inslist = new SLIST;
		// and insert into dimsimps
		dimsimps.insert(dimsimps.end(),make_pair(insdim, inslist));
	}
	// if so, find appropriate list:
	else
	{
		inslist = finddim->second;
	}
	// at this point, we want to insert tau in inslist:
	inslist->push_back((Simplex<C,PS,BT>*)tau);

	// show the list to which this simplex has just been added:
	// cptr_print<SLIST>(*inslist);
}


template <typename C, typename PS, typename BT>
bool RIPS<C,PS,BT>::AddRIPSCofaces(Simplex<C,PS,BT>* tau, const vector<num>& nblist,
		                           num upto, INSRIP& dimsimps, num callerind)
{
	if (tau == NULL) return false;
	if (RIPSTALK)
	{
		cout << "\n                     coface adder, dimcap "<<upto;
		cin.get();
		cout<<" with simp: "<<*tau;
		cout<<"\n and nblist ";
		c_print<vector<num> >(nblist);
		cin.get();
		//showPoints(nblist);
	}

	// add tau to complex... including all boundary links
	if (tau->getDim() == 0) // if this simplex is just a point...
	{
		this->insertVertex(tau);
	}
	else // higher dimension, add tau to the dimsimps structure!
	{
		storeCoface(tau, dimsimps);
	}

	// if tau was top-dimensional or has no neighbors, exit!!
	if (tau->getDim() >= upto || nblist.empty()) return true;
	// Simplex<C,PS,BT>* sigma; // co-face of tau

	// otherwise, iterate over nblist, get vertex v
	typename vector<num>::const_iterator niter;
	typename vector<num>::iterator intend;

	Simplex<C,PS,BT>* sigma;
	PT_SET siglist; // store vertices of tau in siglist
	vector<num> intersects, vnbrs;
	Point<PS>* v; // current vertex under consideration

	for (niter = nblist.begin(); niter != nblist.end(); ++niter)
	{
		intersects.clear();
		//cout<<" ind: "<<*niter; cin.get();
		siglist = tau->getVerts();
		v = vert.at(*niter); //extract vertex v from nblist
		if ((siglist.insert(v)).second) // augment vertices of tau with this
		{

			sigma = new Simplex<C,PS,BT>(siglist);

			if (tau->getDim()==0) // if we are augmenting a point, we must have an edge...
			{
				// so get its birth time from nbr matrix!
				sigma->birth = getEdgeBirth(callerind,*niter);
			}
			getLowerNbrs(*niter,vnbrs);


			// store intersection of vnbrs and nblist in intersects
			set_intersection(nblist.begin(),nblist.end(),vnbrs.begin(),vnbrs.end(),
					         inserter(intersects,intersects.begin()));

			if (INTRIPSTALK)
			{
				cout<<"\n***********intersect***************\n";
				c_print<vector<num> >(nblist); cout<<" and ";
				c_print<vector<num> >(vnbrs); cout<<" gives ";
				c_print<vector<num> >(intersects);
				cout<<"\n*************************************";
			}
			// and recursively call addcofaces!
			AddRIPSCofaces(sigma, intersects, upto, dimsimps, callerind);
		}
		else return false;
	}
	return true;
}


// makes intersection true/false lower - triangular matrix!


// returns all neighbors of a given point (determined by
// its index in "vert") that precede it in the ordering!


// shows all points in rips complex:
template <typename C, typename PS, typename BT>
void RIPS<C,PS,BT>::showPoints(const vector<num>& toshow) const
{
	typename vector<num>::const_iterator it;
	cout<<"{";
	for (it = toshow.begin(); it != toshow.end(); ++it)
	{
		cout<<*(vert.at(*it));
	}
	cout<<"}\n";
}

// clear point allocated memory
template <typename C, typename PS, typename BT>
void RIPS<C,PS,BT>::clearPoints()
{
   //iterate over vertices!
   typename VERTS::iterator vit;
   for(vit = vert.begin(); vit != vert.end(); ++vit)
   {
       delete *vit;
   }
   vert.clear();
}



template <typename C, typename PS, typename BT>
void RIPS<C,PS,BT>::clear()
{
    //clearPoints(); // NOW HANDLED BY STOPLEX::DESTROY()
	// just destroy RIPMAT structure!
	typename RIPMAT::iterator rit;
	RROW* curvec;
	for (rit = isnbr.begin(); rit != isnbr.end(); ++rit)
	{
		curvec = *rit;
		if (curvec == NULL) continue;
		curvec->clear();
		(*curvec).swap(*curvec);
		delete curvec;
	}
}

template <typename C, typename PS, typename BT>
pair<num,bool> RIPS<C,PS,BT>::makeFromTimeSeriesFile(ifstream& infile)
{
    num wsize; // window size for delay reconstruction
    if (!infile.eof())
	{
		infile >> wsize;
	}
	else
	{
		//cout<<"initopen\n";
		return make_pair(0,false); // error, couldn't open file for reading
	}

    epsilon = 1;

    num maxdim; // how many dimensional simplices to build?
    double radius; // initial radius around point cloud
	// read in the step size and max number of steps
	if (!infile.eof())
	{
	    infile >> stepsize >> numsteps >> radius >> maxdim;
        if (epsilon < 0 || stepsize <= 0 || numsteps < 0 || maxdim < 0 || radius < 0)
        {
            return make_pair(0,false);
        }
	}
	else return make_pair(0,false);


    // store time series entries in a vector.
    PS cur;
    vector<PS> tser;
    while(!infile.eof())
    {
        infile >> cur;
        tser.push_back(cur);
    }
    // close file
    infile.close();

    num numpts = tser.size() - wsize;
    // window larger than time series!?
    if (numpts < 0)
    {
        return make_pair(0,false);
    }

    // make the points?!
    // current point from time series window
    vector<PS> ptcoords;
    Point<PS>* curpt;

    set<Point<PS>*, ptcomplex<PS> > prmap;

    //cout<<"\n making pts frm ts"; cin.get();

	typename vector<PS>::const_iterator titer;
	for(titer = tser.begin(); titer != tser.end() - wsize + 1; ++titer)
	{
		ptcoords.assign(titer,titer+wsize);
		// create a point with these coordinates
		curpt = new Point<PS>(ptcoords);
        //cout<<"\n... "<<*curpt; cin.get();

		// try to insert current point and radius
		if (! (prmap.insert(curpt).second ))
		{
			// point already exists, so delete
			delete curpt;
		}
	}

    //cout<<"\n making uniq "; cin.get();

	// now iterate over map of points and populate vert, rads, birth...
	typename set<Point<PS>*, ptcomplex<PS> >::const_iterator pr;
	for(pr = prmap.begin(); pr != prmap.end(); ++pr)
	{
		vert.push_back(*pr);
		//cout<<"\n rad add: "<<pr->second;
		rads.push_back(radius);
	}


	this->psize = vert.size();
	// arbitrary topdim!
	this->topdim = maxdim;
	//cout<<"\n #file points: "<<flines; cin.get();

	makeNbrMatrix_GrowBalls();
    //cout<<"\n nbr matrix filled! making rips complex"; cin.get();
	incRIPS_GrowBalls(maxdim, 1+INITBT);
	//cout<<*this;
	return make_pair(vert.size(),true);
}

#endif /* RIPS_HPP_ */
