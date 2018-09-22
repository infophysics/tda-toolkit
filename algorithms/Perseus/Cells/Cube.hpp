/*
 * Cube.hpp
 *
 * Contains functions for handling the Cube class
 */

#ifndef CUBE_HPP_
#define CUBE_HPP_

# include "Cube.h"

// set vertices to this point
template <typename C, typename PS, typename BT>
void Cube<C,PS,BT>::setAnchor(const Point<PS>* anch)
{
	anchor = anch;
}

// set addin vector to this addin
template <typename C, typename PS, typename BT>
void Cube<C,PS,BT>::setAddin(const vector<num>* add) // set addin vector
{
	addin = add;
}


template <typename C, typename PS, typename BT>
num Cube<C,PS,BT>::getDim() const // get cube dimension, -1 if addin is NULL or empty
{
	return (addin == NULL) ? -1 : addin->size();
}

// check for equality
template <typename C, typename PS, typename BT>
bool Cube<C,PS,BT>::operator == (const Cube<C,PS, BT>& other) const // check equality
{
	return (addin == other.addin && anchor == other.anchor);
}

// check for ordering (lexico of anchor points, then dimension)
template <typename C, typename PS, typename BT>
bool Cube<C,PS,BT>::operator < (const Cube<C,PS,BT>& other) const // check vert-lexico order
{
	if (anchor == NULL || other.anchor == NULL) return false;
	if (*anchor < *(other.anchor)) return true;
	return false;
}

#endif /* CUBE_HPP_ */
