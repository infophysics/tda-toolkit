/*
 * Point.hpp
 *
 *  Created on: Jul 15, 2010
 *      Author: administrtor
 */

#ifndef POINT_HPP_
#define POINT_HPP_

# include "Point.h"

template <typename PS>
void Point<PS>::setCoords(const vector<PS>& co)
{
	coords = co;
}

template <typename PS>
Point<PS>& Point<PS>::operator = (const Point<PS>& other) // to do pt[3] = 7, sets coordinates!
{
	coords.assign(other.begin(), other.end());
	return *this;
}


template <typename PS>
PS& Point<PS>::operator [] (const int cdim) // to do pt[3] = 7, sets coordinates!
{
		return coords.at(cdim);
}

template <typename PS>
void Point<PS>::erase(const int pos) // erases coordinate at pos^th position!
{
	coords.erase(coords.begin()+pos);
}

// returns hash value, currently just address
template <typename PS>
num Point<PS>::getHash() const
{
	return hash;
}

template <typename PS>
int Point<PS>::getDim() const
{
	return coords.size();
}

template <typename PS>
void Point<PS>::clear()
{
	coords.clear();
	// release memory!
	vector<PS>(coords).swap(coords);
}

template <typename PS>
void Point<PS>::push(const PS& coor)
{
	coords.push_back(coor);
}

// iterate over coords!
template <typename PS>
typename vector<PS>::const_iterator Point<PS>::begin() const
{
	return coords.begin();
}

template <typename PS>
typename vector<PS>::const_iterator Point<PS>::end() const
{
	return coords.end();
}

// checks for equality
template <typename PS>
bool Point<PS>::operator == (const Point<PS>& other) const
{
	if (getDim() != other.getDim()) return false;
		for (int i=0; i<getDim(); i++)
	{
		if (coords[i] != other.coords[i]) return false;
	}
	return true;
}
	// checks for not-equality
template <typename PS>
bool Point<PS>::operator != (const Point<PS>& other) const
{
	return (! (*this == other));
}

// less than operator, lexicographical ordering
template <typename PS>
bool Point<PS>::operator < (const Point<PS>& other) const
{
	//cout<<"\n compare: "<<*this<<" at "<<this<<" to "<<other<<" at "<<&other;
	if (getDim() < other.getDim()) return true;
	if (getDim() > other.getDim()) return false;

	int td = getDim();
	for (int i=0; i<td; i++)
	{
		if (coords[i] == other.coords[i]) continue;
		if (coords[i] < other.coords[i])
		{
			return true;
		}
		else return false;
	}
	return false; // if equal
}

// returns distance to other point, will
// be negative in case of dimension mismatch
// computes l^p norm with default p = 2.
template <typename PS>
double Point<PS>::dist (const Point<PS>& other, double normtype = 2) const
{
	double toret = -1; // assumed distance is invalid, of course.
	num n = coords.size();

	if (normtype >= 1 && n == (num)other.coords.size())
	{
		toret = 0;
		for (int i=0; i<n; i++)
		{
			// sum up differences of coordinates raised to normtype
			toret += pow(abs((double)(coords[i] - other.coords[i])),(double)normtype);
		}
		// and now take the appropriate root
		toret = pow((float)toret,(float)(1/normtype));
	}
	return toret;
}

// generic distance function
template <typename PS>
double dist (const Point<PS>& P, const Point<PS>& Q, float norm=2)
{
	return P.dist(Q,norm);
}


#endif /* POINT_HPP_ */
