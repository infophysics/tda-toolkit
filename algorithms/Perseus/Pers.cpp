 /*
 * Pers.cpp
 *
 * Contains main() function for Persistent Homology
 * Using Discrete Morse Theory
 */
#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>

#include "Cells/All.h"
#include "Complexes/All.h"
#include "Algos/All.h"

// class from which the birth times arise, usually = "int"
typedef int BC;
// class defining the ring over which we perform chain algebra, usually = "int"
typedef int CC;

using namespace std;

// silly string functions for dealing with input arguments, change case
inline void lowerCase(std::string& s)
{std::transform(s.begin(), s.end(), s.begin(), ::tolower);}

inline void upperCase(std::string& s)
{std::transform(s.begin(), s.end(), s.begin(), ::toupper);}

// wrapper declaration for converting input file data to cell complexes;
// actual functions follow main():

// sparse cubical complex from top cube info
bool buildCToplexFromFile (ifstream&, MComplex<CC,BC>&, CToplex<CC,int,BC>&);
// dense cubical complex from top cube info
bool buildDenseCToplexFromFile (ifstream&, MComplex<CC,BC>&, DenseCToplex<CC,BC>&);
// uniform triangulation from top simplex info
bool buildSToplexFromFile (ifstream&, MComplex<CC,BC>&, SToplex<CC,double,BC>&, bool, bool);
// non-uniform triangulation from top simplex info
bool buildNMFSToplexFromFile (ifstream&, MComplex<CC,BC>&, SToplex<CC,double,BC>&, bool, bool);
// rips complex from points with non-uniform birth times
bool buildRIPSFromFile (ifstream&, MComplex<CC,BC>&, RIPS<CC,double,BC>&);
// typical rips complex, growing balls around points
bool buildBRIPSFromFile (ifstream&, MComplex<CC,BC>&, RIPS<CC,double,BC>&,bool,bool,bool);
// rips complex from point-point distance matrix
bool buildBRIPSFromDistMatrixFile(ifstream&, MComplex<CC,BC>&, RIPS<CC,double,BC>&,bool);
// rips complex from delay reconstruction of time series
bool buildBRIPSFromTimeSeriesFile(ifstream&, MComplex<CC,BC>&, RIPS<CC,double,BC>&);


// main, here are the currently supported complexes:
// cubtop, scubtop, rips, brips, brips_comrad, simtop, nmfsimtop, movie
int main(int argc, char* argv[])
{
	bool savegens = false; // store generator chains
	bool truncate = false; // ignore boundary-less cells of top dimension when
	                       // computing persistence: for rips complexes only.
	// error message for wrong number of program arguments:
	if (argc < 3 || argc > 5)
	{
		cout<<"\nIncorrect number of arguments!";
		cout<<"\nUsage: <Program Name> <Input Type> <Input Filename> (optional)<Output Filename> (optional)<Engine>\n";
		cout<<"Where: \n  1. The acceptable <Input Type> is SimTop, CubTop, Rips, etc... \n";
		cout<<"  2. If <Output Filename> is unspecified, \"output\" will be used with suffixes.\n";
		cout<<"  3. The legal <engine> values are R, C and A for reduction, coreduction and alternation.";
		return -1;

	}

	// obtain string representations of input options:
	string input_type = argv[1];
	string infile = argv[2];
	string outfile = (argc >= 4)? argv[3] : "output";
	string eng = (argc >= 5) ? argv[4] : "a";

	//cout<<"\neng is "<<eng; cin.get();

	// verify input file integrity
	ifstream inf;
	inf.open(infile.c_str(), ifstream::in);

	if (!inf.good())
	{
		cout<<"\nUnable to Open File "<<infile<<" for reading!";
		return -2;
	}

	// declare generic cell complex
	MComplex<CC,BC>* ccomp = new MComplex<CC,BC>;

    // now we make a cell complex that depends on the choice of input
	lowerCase(input_type);
    // and choose our engine: coreduction, reduction, or the amazing alternator!
	lowerCase(eng);


    // Use this when Miro complains about yet another bug:
    if (input_type == "debug")
    {
        // build kahle's random rips complex on n vertices, where
        num n = atoi(outfile.c_str());
        RIPS<num,num,num> randcomp;

        randcomp.makeRandom_EdgeShuffle(n, true); // bool: make edge file

       //if (ans == 'y')
       randcomp.ComputePersistence(n, eng, true);
       return 0;

    }

	// use for dense cubical grid of top cubes
	else if (input_type == "cubtop")
	{
		DenseCToplex<CC,BC> mytop;
		if (! buildDenseCToplexFromFile(inf, *ccomp, mytop) )
		{
			cout<<"\nFile "<<infile<<" exists but contains bogus data";
			return -3;
		}

	}
	//use for sparse cubical toplex
	else if (input_type == "scubtop")
	{
		CToplex<CC, int, BC> mytop;
		if (! buildCToplexFromFile(inf, *ccomp, mytop) )
		{
			cout<<"\nFile "<<infile<<" exists but contains bogus data";
			return -3;
		}
		//cout<<mytop;
		//ccomp->checkComplex();
	}
	// :::::::::::::RIPS COMPLEX POINT DATA::::::::::::::::::::::
	// point-radius pairs for rips complex
	else if (input_type == "rips")
	{
		RIPS<CC, double, BC> rcomp;
		if (! buildRIPSFromFile (inf, *ccomp, rcomp) )
		{
			cout<<"\nFile "<<infile<<" exists but contains bogus data";
			return -3;
		}
		// debug?
		truncate = true;
	}
	// rips complex, each point has different initial radius
	else if (input_type == "brips")
	{
		RIPS<CC, double, BC> rcomp;
		if (! buildBRIPSFromFile (inf, *ccomp, rcomp, false, false, false) )
		{
			cout<<"\nFile "<<infile<<" exists but contains bogus data";
			return -3;
		}
		truncate = true;
	}
    else if (input_type == "brips_capped" )
	{
		RIPS<CC, double, BC> rcomp;
		if (! buildBRIPSFromFile (inf, *ccomp, rcomp, false, true, false) )
		{
			cout<<"\nFile "<<infile<<" exists but contains bogus data";
			return -3;
		}
		truncate = true;
	}
	// rips complex, there is a common initial radius!
	else if (input_type == "brips2" || input_type == "brips_comrad")
	{
		RIPS<CC, double, BC> rcomp;
		if (! buildBRIPSFromFile (inf, *ccomp, rcomp, true, true, false) )
		{
			cout<<"\nFile "<<infile<<" exists but contains bogus data";
			return -3;
		}
		truncate = true;
	}
    // rips complex with dimension cap
	else if (input_type == "brips_comrad_capped" || input_type == "brips_capped_comrad")
	{
		RIPS<CC, double, BC> rcomp;
		if (! buildBRIPSFromFile (inf, *ccomp, rcomp, true, true, false) )
		{
			cout<<"\nFile "<<infile<<" exists but contains bogus data";
			return -3;
		}
		truncate = true;
	}
	// witnessed rips complex
	else if (input_type == "brips_witness")
	{
		RIPS<CC, double, BC> rcomp;
		if (! buildBRIPSFromFile (inf, *ccomp, rcomp, true, true, true) )
		{
			cout<<"\nFile "<<infile<<" exists but contains bogus data";
			return -3;
		}
		truncate = true;
	}
	// rips complex from distance matrix
	else if (input_type == "distmat")
	{
		RIPS<CC, double, BC> rcomp;
		if (! buildBRIPSFromDistMatrixFile (inf, *ccomp, rcomp, false) )
		{
			cout<<"\nFile "<<infile<<" exists but contains bogus data";
			return -3;
		}
		truncate = true;
	}
	// rips complex from correlation matrix
	else if (input_type == "corrmat")
	{
		RIPS<CC, double, BC> rcomp;
		if (! buildBRIPSFromDistMatrixFile (inf, *ccomp, rcomp, true) )
		{
			cout<<"\nFile "<<infile<<" exists but contains bogus data";
			return -3;
		}
		truncate = true;
	}

	// time series delay reconstruction
	else if(input_type == "timeser")
    {
        RIPS<CC, double, BC> rcomp;
        if (!buildBRIPSFromTimeSeriesFile (inf, *ccomp, rcomp) )
        {
            cout<<"\nFile "<<infile<<" exists but contains bogus data";
			return -3;
        }
        truncate = true;
    }


	// :::::::::::MANIFOLD TOP SIMPLEX DATA:::::::::::::::::::::
	// list of top simplices, all of the same dimension
	else if (input_type == "simtop")
	{
		SToplex<CC, double, BC> stop;
		if (! buildSToplexFromFile (inf, *ccomp, stop, true, false) )
		{
			cout<<"\nFile "<<infile<<" exists but contains bogus data";
			return -3;
		}
	}

    // manifold stoplex with no birth time information:
	else if (input_type == "simtop_static")
	{
		SToplex<CC, double, BC> stop;
		if (! buildSToplexFromFile (inf, *ccomp, stop, false, false) )
		{
			cout<<"\nFile "<<infile<<" exists but contains bogus data";
			return -3;
		}
	}
	// ::::::::::::NON MANIFOLD TOP SIMPLEX DATA::::::::::::::::
	// list of top simplices of different dimensions
	else if (input_type == "nmfsimtop")
	{
		SToplex<CC, double, BC> stop;
		if (! buildNMFSToplexFromFile (inf, *ccomp, stop, false, true) )
		{
				cout<<"\nFile "<<infile<<" exists but contains bogus data";
				return -3;
		}
		//truncate = true;
	}

    else if (input_type == "nmfsimtop_capped")
	{
		SToplex<CC, double, BC> stop;
		if (! buildNMFSToplexFromFile (inf, *ccomp, stop, true, true) )
		{
				cout<<"\nFile "<<infile<<" exists but contains bogus data";
				return -3;
		}
		//truncate = true;
	}

	else if (input_type == "nmfsimtop_static_capped")
	{
		SToplex<CC, double, BC> stop;
		if (! buildNMFSToplexFromFile (inf, *ccomp, stop, true, false) )
		{
				cout<<"\nFile "<<infile<<" exists but contains bogus data";
				return -3;
		}
		//truncate = true;
	}


    else if (input_type == "nmfsimtop_static")
	{
		SToplex<CC, double, BC> stop;
		if (! buildNMFSToplexFromFile (inf, *ccomp, stop, false, false) )
		{
				cout<<"\nFile "<<infile<<" exists but contains bogus data";
				return -3;
		}
		//truncate = true;
	}

    // to run: ./perseus ricci "simtop file" "edge info"
	else if (input_type == "ricci")
	{
	    MComplex<CC,double>* dcomp =  new MComplex<CC,double>;

	    SToplex<CC, num, double> stop;
	    ifstream edgef;
	    cout<<" trying to open: "<<argv[3];
	    edgef.open(argv[3], ifstream::in);

	    stop.makeRicci(inf,edgef);
	    stop.writeComplex(*dcomp);

        //dcomp->printSizeInfo();
        //cin.get();
        //dcomp->checkComplex();
        //dcomp->printList(0,1);

	    AlternateAndUpdate(dcomp,savegens,0.2);
	    PComplex<CC,double> pcomp;
        if (FLOWTALK) cout<<"\nComputing Persistence Intervals!";
        if (MAKEBPS) {cout<<"... ";    cin.get(); }
        pcomp.COMPUTE_INTERVALS(*dcomp, savegens, truncate);

        //ccomp.showInts();
        pcomp.makeOutputFiles(outfile);
        //pcomp.showBetti();

        delete dcomp;
        return 0;
	}

    else if (input_type == "bary")
    {
        num totcount = 1;
        vector<vector<num> > sizes;
        vector<num> svec; // current size vec...

        SToplex<CC, double, BC>* stop;
        MComplex<CC,BC>* mycomp;

        //stop.barySub(bcstop);
        // output file?
        outfile = infile.replace(infile.find("."),4,"_bmat.txt");
        ofstream outf(outfile.c_str());
        //ostream& outf = cout;

        // subdivide simplicial complex

        num maxtowerht = 0; // largest number of iterations needed.
        num tht = 0;

        for (num count = 0; count < totcount; count++)
        {
            svec.clear();
            //cout<<"\n count "<<count; cin.get();

            stop = new SToplex<CC,double,BC>;

            pair<num,bool> iores = stop->makeFromFile(inf,false,false,false);
            if (iores.second == false)
            {
                outf<<"\n File Error! \n";
                return -1;
            }

            mycomp = new MComplex<CC,BC>;

            stop->writeComplex(*mycomp);
            cout<<"\n\n******* Unreduced size of original: "<< mycomp->size(); //cin.get();

            mycomp->hyperq = true;

            if (eng == "r")
            {
                ReduceAndUpdate(mycomp,savegens,0.0,true);
                //CoreduceAndUpdate(mycomp,savegens,0.0,true);
                // when done, reduce again with traces
            }
            else if (eng == "c")
            {
                CoreduceAndUpdate(mycomp,savegens,0.0,true);
                //ReduceAndUpdate(mycomp,savegens,0.0,true);
            }
            else if (eng == "a2")
            {
                AlternateAndUpdate(mycomp,savegens,0.0,true,true);
            }
            else
            {
               tht = AlternateAndUpdate(mycomp,savegens,0.0,false,true);
            }

            // get max iteration count:
            if ((num)svec.size() > maxtowerht) maxtowerht = svec.size();
            // add to the sizes structure:
            sizes.push_back(svec);
/*
            mycomp->printSizeInfo(true);

            vector<vector<CC> >bsmat;

            cout<<"\n   Density: "<<mycomp->getBoundaryDensity(5);
            cout<<"\n   UnitFrac: "<<mycomp->getUnitFraction(5);
            //mycomp->printList(0,)
            mycomp->boundaryMatrix(3,bsmat,1);
            outf << bsmat.size() <<" ";
            outf << bsmat.at(0).size();
            for (num i=0; i<(num)bsmat.size(); i++)
            {
                outf<<"\n";
                c_print(bsmat.at(i), outf);
            }
*/
            //cout<<*mycomp;
            delete mycomp;
            delete stop;
            inf.clear();
            inf.seekg(0);
        }


        //outf << maxtowerht << endl;
        // generate output file using sizes structure:
        for (num i=0; i<(num)sizes.size();++i)
        {
            //outf<<sizes.at(i).size(); // write iteration count
            for (num j=0; j<(num)sizes.at(i).size();++j)
            {
                outf<<" "<<sizes.at(i).at(j);
            }
            // okay, but if there was a shortfall compared to maxtowerht, want
            // to pad with zeros:
            for (num k=(num)sizes.at(i).size(); k<maxtowerht;++k)
            {
                outf<<" "<<0;
            }
            outf<<"\n";
        }

        outf.close();


        return 0;
    }

	else if (input_type == "cubmov")
	{
		delete ccomp;
		CToplex<CC,double,BC> mtop;
		MComplex<CC,BC>* mycomp = mtop.readMovie(inf);
		if (mycomp == NULL)
		{
			cout<<"\nFile "<<infile<<" exists but contains bogus data";
			return -3;
		}
		mycomp->MorseWrapper_Cored(savegens);
		PComplex<CC,BC> mypc;
		if (FLOWTALK)
		{
		    cout<<"\nComputing Persistence Intervals!";
            if (MAKEBPS) {cout<<"... ";    cin.get(); }
		}
		mypc.COMPUTE_INTERVALS(*mycomp, savegens);
		//ccomp.showInts();
		mypc.makeOutputFiles(outfile);
		mypc.showBetti();
		if (FLOWTALK) cout<<"\n\nDone!!! Please consult ["<<outfile<<"*.txt] for results.\n\n";
		delete mycomp;
		return 0;
	}

	else
	{
		cout<<"\nThe Data format \""<<input_type<<"\" is currently unsupported.";
		return 0;
	}

	if(FLOWTALK)
    {
        cout<<" Complex stored with "<<ccomp->size()<<" cells!";
        if (MAKEBPS) cin.get();
    }

    // optimize!!
    ccomp->hyperq = true;
    // choose engine
    if (eng == "r")
    {
       ReduceAndUpdate(ccomp,savegens,0.2);
    }
    else if (eng == "c")
    {
       CoreduceAndUpdate(ccomp,savegens,0.2);
    }
    else if (eng == "s") // skip (co)reductions altogether, just compute persistence
    {

    }
    else
    {
       AlternateAndUpdate(ccomp,savegens,0.2, false, true);
    }

 // bdry matrix for snf
/*
    MComplex<CC,BC>* rcomp = new MComplex<CC,BC>;
    ccomp->moveOver(*rcomp);
    delete ccomp;
    ccomp = rcomp;

    vector<vector<CC> >bmat;
    ccomp->boundaryMatrix(6,bmat);
    for (int i=0; i<bmat.size(); i++)
    {
        cout<<"\n";
        c_print(bmat.at(i));
    }
*/
    // dimensional size info:
    //ccomp->printSizeInfo(true);
    //cin.get();

	PComplex<CC,BC> pcomp;
	if (FLOWTALK) cout<<"\nComputing Persistence Intervals!";
	if (MAKEBPS) {cout<<"... ";    cin.get(); }
	pcomp.COMPUTE_INTERVALS(*ccomp, savegens, truncate);

	//ccomp.showInts();
	pcomp.makeOutputFiles(outfile);
    //pcomp.showBetti();

	delete ccomp;
	//delete ccomp;

	if (FLOWTALK) cout<<"\n\nDone!!! Please consult ["<<outfile<<"*.txt] for results.\n\n";
	return 0;
}

// subroutine to construct cell complex from dense cubical toplex information. The text
// file infile has the following format:
// <top dimension n>
// <dim 1 extent, e_1>
// <dim 2 extent, e_2>
// :
// :
// <dim n extent, e_n>
// <birth time of cube at 0,0,...,0>
// <birth time of cube at 1,0,...,0>
// :
// <birth time of cube at (e_1 - 1),0,...,0>
// <birth time of cube at 0,1,...,0>
// :
// <birth time of cube at (e_1 - 1),1,...,0>
// :
// :
// <birth time of cube at (e_1 - 1),(e_2 - 1),...,(e_n - 1)>
bool buildDenseCToplexFromFile (ifstream& inf, MComplex<CC,BC>& ccomp, DenseCToplex<CC,BC>& mytop)
{
	pair<num,bool> iores = mytop.makeFromFile(inf);
	inf.close();

	if (iores.second == false)
	{
		return false;
	}

	if (FLOWTALK) cout<<"\nRead "<<iores.first<<" birth times from Input File";
	if (MAKEBPS) cin.get();
	//cout<<mytop;
	//mytop.quickWriteComplex(ccomp);
	mytop.writeComplex(ccomp);
	mytop.Destroy();
	if (FLOWTALK) cout<<"\nWritten to Cell Complex!";
	if (MAKEBPS) cin.get();

	return true;
}

bool buildCToplexFromFile (ifstream& inf, MComplex<CC,BC>& ccomp, CToplex<CC,int,BC>& mytop)
{
	pair<num,bool> iores = mytop.makeFromFile(inf);
	inf.close();

	if (iores.second == false)
	{
		return false;
	}
	//cout<<mytop;

	if (FLOWTALK) cout<<"\nRead "<<iores.first<<" birth times from Input File";
	if (MAKEBPS) cin.get();
	mytop.writeComplex(ccomp);
	//mytop.Destroy();
	//cout<<ccomp;
	if (FLOWTALK) cout<<"\nWritten to Cell Complex!";
	if (MAKEBPS) cin.get();

	return true;
}

// subroutine to construct cell complex from RIPS complex point information. The TEXT file inf
// should have the following format:
// <dimension of embedding space, n>
// <n coordinates of point 1> <radius of neighbor cutoff ball> <birth time>
// <n coordinates of point 2> <radius of neighbor cutoff ball> <birth time>
// :
// :
// <n coordinates of point k> <radius of neighbor cutoff ball> <birth time>

bool buildRIPSFromFile (ifstream& inf, MComplex<CC,BC>& ccomp, RIPS<CC, double, BC>& rcomp)
{
	pair<num,bool> iores = rcomp.makeFromFile(inf);
	inf.close();

	if (iores.second == false)
	{
		return false;
	}

	if (FLOWTALK) cout<<"\nRead "<<iores.first<<" point/radius pairs and birth times!";
	if (MAKEBPS) cin.get();


	if (FLOWTALK) cout<<"\nWriting Cell Complex From RIPS Complex";
	if (MAKEBPS) cin.get();

	rcomp.writeComplex(ccomp);

	//cout<<rcomp; cin.get();

	if (FLOWTALK) cout<<"\nDone!";
	if (MAKEBPS) cin.get();

	return true;
}

// subroutine to construct cell complex from RIPS complex point information. The TEXT file inf
// should have the following format:
// <dimension of embedding space, n>
// <radius multiplier epsilon> <step size delta> <number of steps k>
// <n coordinates of point 1> <radius of neighbor cutoff ball>
// <n coordinates of point 2> <radius of neighbor cutoff ball>
// :
// :
// <n coordinates of point k> <radius of neighbor cutoff ball>

bool buildBRIPSFromFile (ifstream& inf, MComplex<CC,BC>& ccomp, RIPS<CC, double, BC>& rcomp,
                          bool comred = false, bool capped=false, bool witness=false)
{
	pair<num,bool> iores = rcomp.makeFromFile_GrowBalls(inf,comred,capped,witness);
	inf.close();

	if (iores.second == false)
	{
		return false;
	}

	if (FLOWTALK) cout<<"\nRead "<<iores.first<<" point/radius pairs and birth times!";
	if (MAKEBPS) cin.get();

	if (FLOWTALK) cout<<"\nWriting Cell Complex From RIPS Complex";
	if (MAKEBPS) cin.get();

	rcomp.writeComplex(ccomp);

	//cout<<rcomp; cin.get();

	if (FLOWTALK) cout<<"\nDone!";
	if (MAKEBPS) cin.get();

	return true;
}

// subroutine to build rips complex from symmetric distance matrix file. the format is:
// <number of points = n>
// <base distance> <increment step size> <number of increments>
// <dist to p1> <dist to p2> <dist to p3> .... <dist to pn>
// :
// :
bool buildBRIPSFromDistMatrixFile (ifstream& inf, MComplex<CC,BC>& ccomp, RIPS<CC, double, BC>&
                                   rcomp, bool iscorrmat = true)
{
    pair<num,bool> iores = rcomp.makeFromDistMatrixFile(inf,iscorrmat,false);
	inf.close();

	if (iores.second == false)
	{
		return false;
	}

	if (FLOWTALK) cout<<"\nRead "<<iores.first<<" point/distance pairs!";
	if (MAKEBPS) cin.get();

    //rcomp.showNbrMatrix();


	//if (FLOWTALK) cout<<"\nWriting Cell Complex From RIPS Complex";
	//if (MAKEBPS) cin.get();

	rcomp.writeComplex(ccomp);

	//cout<<rcomp; cin.get();

	if (FLOWTALK) cout<<"\nDone!";
	if (MAKEBPS) cin.get();

	return true;
}

bool buildBRIPSFromTimeSeriesFile (ifstream& inf, MComplex<CC,BC>& ccomp, RIPS<CC, double, BC>& rcomp)
{
    pair<num,bool> iores = rcomp.makeFromTimeSeriesFile(inf);
	inf.close();

	if (iores.second == false)
	{
		return false;
	}

	if (FLOWTALK) cout<<"\nRead "<<iores.first<<" time series points!";
	if (MAKEBPS) cin.get();

    //rcomp.showNbrMatrix();


	//if (FLOWTALK) cout<<"\nWriting Cell Complex From RIPS Complex";
	//if (MAKEBPS) cin.get();

	rcomp.writeComplex(ccomp);

	//cout<<rcomp; cin.get();

	if (FLOWTALK) cout<<"\nDone!";
	if (MAKEBPS) cin.get();

	return true;
}


// subroutine to construct cell complex from simplicial toplex information. the TEXT file inf
// should have the following format:
// <simplex dimension m>
// <ambient space dimension n>
// <x_11,...,x_1n> <x_21,...,x_2n> ... <x_(m+1)1,...x_(m+1)n> <birth time> { i.e. vertex coords and birth time }
// : [same as above for next simplex]
// : [ etc ]
bool buildSToplexFromFile (ifstream& inf, MComplex<CC,BC>& ccomp,
                           SToplex<CC,double,BC>& stop, bool births = true,
                           bool hascap=false)
{
	pair<num,bool> iores = stop.makeFromFile(inf, births, true, hascap);
	inf.close();

	if (iores.second == false)
	{
		return false;
	}
	if (FLOWTALK) cout<<"\nRead "<<iores.first<<" top simplices from Input File";
	if (MAKEBPS) cin.get();
	if (FLOWTALK) cout<<"\nWriting Cell Complex From Simplicial Complex";
	if (MAKEBPS) cin.get();

	stop.writeComplex(ccomp);

	if (FLOWTALK) cout<<"\nDone!";
	if (MAKEBPS) cin.get();
	return true;
}

// build cell complex for non-manifold simplicial toplex from file. the input TEXT file inf must have
// the following format:
// <dimension of embedding space n>
// <dim of simplex, m> <x_11,...,x_1n> <x_21,...,x_2n> ... <x_(m+1)1,...x_(m+1)n> <birth time>
// : [ similarly for other simplices, each with its own dim m ]

bool buildNMFSToplexFromFile (ifstream& inf, MComplex<CC,BC>& ccomp,
                              SToplex<CC,double,BC>& stop, bool tocap = false,
                              bool hasbirths = true)
{
	pair<num,bool> iores = stop.makeFromFile(inf,hasbirths,false,tocap);
	inf.close();

    //cout<<stop;

	if (iores.second == false)
	{
		return false;
	}
	if (FLOWTALK) cout<<"\nRead "<<iores.first<<" top simplices from Input File ";
	if (MAKEBPS) cin.get();

	if (FLOWTALK) cout<<"\nWriting Cell Complex From Non-manifold simplicial Complex";
	if (MAKEBPS) cin.get();

	//cout<<stop;
	stop.writeComplex(ccomp);

	//ccomp.checkComplex();
	if (FLOWTALK) cout<<"\nDone!";
	if (MAKEBPS) cin.get();
	return true;
}

