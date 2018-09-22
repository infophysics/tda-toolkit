/*
 * Simplex.hpp
 *
 *  Created on: Jul 15, 2010
 *      Author: administrtor
 */

#ifndef SIMPLEX_HPP_
#define SIMPLEX_HPP_

# include "Simplex.h"

// set vertices to given point set
template <typename C, typename PS, typename BT>
void Simplex<C,PS,BT>::setVerts(const PT_SET& pts)
{
	verts = pts;
}

// print out vertices in hash form
template <typename C, typename PS, typename BT>
void Simplex<C,PS,BT>::ptHashPrint(ostream& out = cout) const
{
	typename PT_SET::const_iterator vit;
	out<<"[ ";
	for (vit = verts.begin(); vit != verts.end(); ++vit)
	{
		out << (*vit)->getHash() <<" ";
	}
	out<<"]";
}

template <typename C, typename PS, typename BT>
bool Simplex<C,PS,BT>::addVert(Point<PS>* toadd)
{
	if (verts.find(toadd) == verts.end()) // if we don't have this vertex,...
	{
		// add it
		verts.insert(toadd);
		return true;
	}
	// vertex was already there, we did nothing
	return false;
}

// iterators to loop over vertices: begin
template <typename C, typename PS, typename BT>
typename PT_SET::const_iterator Simplex<C,PS,BT>::begin() const
{
	return verts.begin();
}

// iterators to loop over vertices: end
template <typename C, typename PS, typename BT>
typename PT_SET::const_iterator Simplex<C,PS,BT>::end() const
{
	return verts.end();
}

// dimension... returns -1 in case there are no vertices!
template <typename C, typename PS, typename BT>
num Simplex<C,PS,BT>::getDim () const
{
	return verts.size()-1;
}

// comparison: uses the sorting of sets using < on points!
template <typename C, typename PS, typename BT>
Simplex<C,PS,BT>& Simplex<C,PS,BT>::operator = (const Simplex<C,PS,BT>& other)
{
	verts = other.verts;
	return *this;
}

// comparison: uses the sorting of sets using < on points!
template <typename C, typename PS, typename BT>
bool Simplex<C,PS,BT>::operator == (const Simplex<C,PS,BT>& other) const
{
	if (getDim() != other.getDim()) return false;
	typename PT_SET::const_iterator i = begin(), j = other.begin();
	do
	{
		if (!(PT_ORDER(*i,*j)) && PT_ORDER(*j,*i)) return false;
		i++; j++;
	}while (i != end());
	return true;
}

// returns hash value, summing over point hashes... modded by argument!
template <typename C, typename PS, typename BT>
num Simplex<C,PS,BT>::getHash(num modby = 1) const
{
	return getPtSetHash<PS>(verts,modby);
}

// comparison: by dimension, and then lexicographical in point addresses
template <typename C, typename PS, typename BT>
bool Simplex<C,PS,BT>::operator < (const Simplex<C,PS,BT>& other) const
{
	if (getDim() < other.getDim()) return true;
	if (getDim() > other.getDim()) return false;

	// iterate over vertices
	typename PT_SET::const_iterator i = begin(), j = other.begin();
	do
	{
		if (PT_ORDER(*i, *j)) return true; // less than satisfied
		else if (PT_ORDER(*j, *i)) return false; // greater than satisfied
		i++; j++; // current points equal, go on!
	}while (i != end());
	return false;
}

#endif /* SIMPLEX_HPP_ */
