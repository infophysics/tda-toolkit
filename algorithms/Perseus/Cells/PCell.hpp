/*
 * PCell.hpp
 *
 *  Created on: Jan 10, 2011
 *      Author: Vidit
 */

#ifndef PCELL_HPP_
#define PCELL_HPP_

# include "PCell.h"

template <typename C, typename BT>
void PCell<C,BT>::makeFromCell(const Cell<C,BT>* other)
{
	dim = other->getDim();
	birth = other->birth;
	marked = false;
	TLIST = new PCCHAIN;

	kgen = other->kgen;
	//this->k_index = other->k_index;
	// we can change chain coefficients to field ones here maybe?
}

template <typename C, typename BT>
void PCell<C,BT>::print(ostream& out = cout) const
{
	out<<"<"<<birth<<";"<<dim<<","<<kindex<<">";
	if (isMarked()) out<<"(*)";
}

template <typename C, typename BT>
num PCell<C,BT>::getDim() const
{
	return dim;
}

template <typename C, typename BT>
bool PCell<C,BT>::addBdryLink(PCell<C,BT>*& toadd, const C& coeff)
{
	return bdry.addLink(toadd, coeff);
}

template <typename C, typename BT>
const PCCHAIN& PCell<C,BT>::getBD() const
{
	return bdry;
}

// DUMMY: DOES NOT ACTUALLY RETURN A COBOUNDARY. AVOID!
template <typename C, typename BT>
const PCCHAIN& PCell<C,BT>::getCB() const
{
	return bdry;
}

#endif /* PCELL_HPP_ */
