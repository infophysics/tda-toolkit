/*
 * Point.h
 *
 *  Created on: Jul 15, 2010
 *      Author: administrtor
 */

#ifndef POINT_H_
#define POINT_H_

#include <vector>
#include <cmath>
// for random hashing
#include <cstdlib>
#include <ctime>

#include "../Global/Global.h"

using namespace std;

// choose how to order points when (for instance) building simplices...
// set to ptcomp<PS> to order by memory address (fast!) or to ptcomplex<PS>
// to sort by lexico ordering on coordinates (slow, but easier to read and
// debug and ensure uniqueness!)
# define PT_ORDER ptcomp<PS>

template <typename PS>
class Point
{
public:
	vector<PS> coords; // vector of coordinates
	num hash; // hash value
		// silly constructor
	Point()
	{
		coords.clear();
		hash = (num)rand();
	}

	Point (const vector<PS>& co)
	{
		coords = co;
		hash = (num)rand();
	}

	Point (const Point<PS>& other)
	{
		coords.assign(other.coords.begin(), other.coords.end());
		hash = other.hash;
	}

	friend ostream& operator << (ostream& out, const Point<PS>& toprint)
	{
		typename vector<PS>::const_iterator coord;
		out<<"( ";
		for (coord = toprint.begin(); coord != toprint.end(); ++coord)
		{
			out<<(*coord)<<" ";
		}
		out<<")";
		return out;
	}

	void setCoords(const vector<PS>&); // sets coordinates to given vector
	int getDim() const; // returns dimension
	void clear(); // wipes out coordinates
	void push(const PS&); // adds dimension, (1,0) pushed 4 -> (1,0,4) etc.
	typename vector<PS>::const_iterator begin() const; // begin coordinate iter
	typename vector<PS>::const_iterator end() const;  // end coordinate iter
	void erase(const int); // delete a coordinate
	// point operations
	bool operator == (const Point<PS>& other) const; // equality check
	bool operator != (const Point<PS>& other) const; // inequality check
	bool operator < (const Point<PS>& other) const; // lexico ordering
	Point<PS>& operator = (const Point<PS>& other); // assignment
	PS& operator [] (const int); // to do pt[3] = 7, set coordinates individually
	double dist(const Point<PS>& other, double) const; // returns distance from other point
	num getHash() const;

	// destroys
	~Point()
	{
		// clears vertices!
		clear();
	}

};

// point comparison for map-making
// uses pointer address for quick comparison!!
template <typename PS>
struct ptcomp
{
	bool operator ()(Point<PS>* p, Point<PS>* q)
	{
		if (p == q) return false;
		if (p == NULL && q != NULL) return true;
		if (q == NULL) return false;
		return (p < q);
	}
};

// lexico ordering, uses actual points rather than
// addresses in memory for comparison. slower, but
// more useful and definitely prettier than ptcomp
template <typename PS>
struct ptcomplex
{
	bool operator ()(Point<PS>* p, Point<PS>* q)
	{
		if (p == q) return false;
		if (p == NULL && q != NULL) return true;
		if (q == NULL) return false;
		//cout<<"\n lex compare: "; cin.get();
		//cout<<p<<" and "<<q;
		bool toret = (*p < *q);
		//cout<<" ans "<<toret;
		return toret;
	}
};
#endif /* POINT_H_ */
