/*
 * Simplex.h
 *
 * Simplex Class definitions
 */

#ifndef SIMPLEX_H_
#define SIMPLEX_H_

#include <set>
#include "Point.hpp"
#include "Cell.hpp"

// set of points... see Point.h for definition of ptcomp<> etc.
#define PT_SET set<Point<PS>*, PT_ORDER >

template <typename C, typename PS, typename BT = int>
class Simplex : public Cell<C,BT>
{
public:
	PT_SET verts; // vertices of simplex stored as set of (pointers to) points

	// silly constructor
	Simplex():Cell<C,BT>()
	{
		verts.clear();
	}

	// constructor as set of points!!
	Simplex(const PT_SET& pts):Cell<C,BT>(pts.size()-1, 0)
	{
		verts = pts;
	}

	// copy constructor
	Simplex(const Simplex<C,PS,BT>& other):Cell<C,BT>((Cell<C,BT> &)other)
	{
		verts = other.verts;
	}

	// construct 0-dim simplex, just a single point!
	Simplex(const Point<PS>* onlyVert)
	{
		verts.insert((Point<PS>*)onlyVert);
	}

/*
    void print (ostream& out) const // print helper
	{
    	//this->Cell<C,BT>::print(out);
		typename PT_SET::const_iterator i;
		out<<"[ ";

		for (i = begin(); i != end(); ++i)
		{
			out<<**i<<" ";
		}
		out<<"]";
		if (this->isCrit) out<<" (*)";
	}
*/

    friend ostream& operator << (ostream& out, const Simplex<C,PS,BT>& toprint)
    {
    	toprint.print(out);
    	return out;
    }

    Point<PS>* getAnchor() const { return *(begin()); } // returns anchor point
    const PT_SET& getVerts() const {return verts;}

    bool addVert(Point<PS>*); // add a vertexS
	void setVerts(const PT_SET&); // set vertices to this point set
	typename PT_SET::const_iterator begin() const; // begin iter over verts
	typename PT_SET::const_iterator end() const; // end iter over verts
	num getDim () const; // get simplex dimension

	// simplex operators
	bool operator == (const Simplex<C,PS,BT>&) const; // check equality
	bool operator < (const Simplex<C,PS,BT>&) const; // check vert-lexico order
	Simplex<C,PS,BT>& operator = (const Simplex<C,PS,BT>&); // assignment
	num getHash(num) const;
	void ptHashPrint(ostream&) const;
	// destructor, will not erase actual points!
	~Simplex()
	{
		//cout<<" del simp ";
		verts.clear();
	}
};

// returns the hash of a set of points, summing and modding by modby
template <typename PS>
num getPtSetHash(const PT_SET& tohash, num modby)
{
	num toret = 0;
	typename PT_SET::const_iterator viter;
	for (viter = tohash.begin(); viter != tohash.end(); ++viter)
	{
		toret += (*viter)->getHash();
		toret %= modby;
	}
	return toret;
}



#endif /* SIMPLEX_H_ */
