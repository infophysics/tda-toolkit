/*
 * Cube.h
 *
 *  Created on: Jul 15, 2010
 *      Author: administrtor
 */

#ifndef CUBE_H_
#define CUBE_H_

// this is the size of the cube's edge
#define ES 1

#include <map>
#include <vector>

#include "Point.hpp"
#include "Cell.hpp"

// C is the class of ring coefficients for considering chains over
// cubes, PS is the point space. It is important that for any p in PS,
// p + ES is also in PS. Canonical example: PS = Integers, ES = 1 etc.
// Obviously, requiring the edge-size ES to lie in PS works.

template <typename C, typename PS, typename BT = int>
class Cube : public Cell<C,BT>
{
public:
	/* DATA */

	Point<PS>* anchor;	// lowest lexicographical vertex
	vector<num>* addin;	// dimensions to increment by ES

	// the cube being stored here is given by an anchor point
	// (p,q,r) and an addin vector with a set of dimensions.
	// if addin = [0,1] then the cube spans (p,q,r) and
	// (p+ES,q+ES,r) etc.

	/* FUNCTIONS */

	// constructors
	Cube():Cell<C,BT>()
	{
		anchor = NULL;
		addin = NULL;
	}
	// constructor with given addin vector and anchor point
	Cube(vector<num>* add, Point<PS>* anch)
	:Cell<C,BT>((add==NULL || add->size()==0)? 0:add->size()-1,0)
	{
		anchor = anch;
		addin = add;
	}
	// copy constructor
	Cube(const Cube<C,PS,BT>& other):Cell<C,BT>((Cell<C,BT>&) other)
	{
		anchor = other.anchor;
		addin = other.addin;
	}
	// destructor
	~Cube()
	{
		// nothing to do here, allocated memory
		// will be handled elsewhere
	}

	void print (ostream& out) const // print helper
	{
		typename vector<num>::iterator add_iter;
		if (anchor != NULL)
		{
			out<<"[ "<<*anchor<<" ; { ";
			if (addin != NULL)
			{
				for (add_iter = addin->begin(); add_iter != addin->end(); ++add_iter)
				{
					out<<*(add_iter)<<" ";
				}
				out<<"} ]";
				if (this->isCrit) out<<"(*)";
			}
		}
	}
	void setAnchor(const Point<PS>*); // set vertices to this point set
	void setAddin(const vector<num>*); // set addin vector
	num getDim() const; // get cube dimension

	// cube operators
	bool operator == (const Cube<C,PS,BT>& other) const; // check equality
	bool operator < (const Cube<C,PS,BT>& other) const; // check vert-lexico order


};

#endif /* CUBE_H_ */
