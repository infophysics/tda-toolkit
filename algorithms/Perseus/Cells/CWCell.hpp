/*
 * CWCell.hpp
 *
 * A single cell in a CW complex, along with attaching maps, which are
 * words in lower dimensional cells
 */

#ifndef CWCELL_HPP_
#define CWCELL_HPP_
/*
#include "CWCell.h"

//****************************************************************
//                            CW Cell Methods
// ****************************************************************

// output helper
template <typename C, typename BT>
void CWCell<C,BT>::print(ostream& out)
{
	out<<"[ "<<this->dim<<" , "<<this->ind<<" ] ";
	out<<" att: "<<attmap;
}

//****************************************************************
//                              WORD METHODS
// ***************************************************************

// return size of word, ignoring coefficients. <A,1><B,2><A,-1> has size 3
template <typename C>
int Word<C>::getSize() const
{
	return word.size();
}

// sum of  coefficients, <A,1><B,2><A,-1> returns 1 + 2 + -1 = 2
template <typename C>
C Word<C>::getCSum() const
{
	C toret = 0;
	typename WORD::const_iterator link;
	for (link = begin(); link != end(); ++link)
	{
		toret += link->second;
	}
	return toret;
}

// add a letter (with coefficient to the back of a word
template <typename C>
void Word<C>::addLetterBack(const Cell<C>* letter, const C& coeff)
{
	if (getSize() == 0)
	{
		word.push_back(pair<Cell<C>*,C>((Cell<C>*)letter,coeff));
	}
	else if (letter ==  back().first)
	{
		back().second += coeff;
		if (back().second == 0) // then we should remove this element
		{
			word.pop_back();
		}
	}
	else // just add it
	{
		word.push_back(pair<Cell<C>*,C>((Cell<C>*)letter,coeff));
	}
}

// add a letter with a coefficient to the front of a word
template <typename C>
void Word<C>::addLetterFront(const Cell<C>* letter, const C& coeff)
{
	if (getSize() == 0)
	{
		word.push_front(pair<Cell<C>*,C>(letter,coeff));
	}
	else if (letter ==  front().first)
	{
		front().second += coeff;
		if (front().second == 0) // then we should remove this element
		{
			word.pop_front();
		}
	}
	else // just add it
	{
		word.push_front(pair<Cell<C>*,C>(letter,coeff));
	}
}

// iterator to the first letter
template <typename C>
typename WORD::const_iterator Word<C>::begin() const
{
	return word.begin();
}

// iterator past the last letter
template <typename C>
typename list<pair<Cell<C>*,C> >::const_iterator Word<C>::end() const
{
	return word.end();
}

// iterator to the last element
template <typename C>
typename WORD::const_reverse_iterator Word<C>::rbegin() const
{
	return word.rbegin();
}

// iterator before the first element
template <typename C>
typename list<pair<Cell<C>*,C> >::const_reverse_iterator Word<C>::rend() const
{
	return word.rend();
}

// reference to the first letter
template <typename C, typename BT>
pair<Cell<C,BT>*,C>& Word<C>::front() const
{
	return (pair<Cell<C,BT>*,C>&)word.front();
}

// reference to the last letter
template <typename C, typename BT>
pair<Cell<C,BT>*,C>& Word<C>::back() const
{
	return (pair<Cell<C,BT>*,C>&)word.back();
}


// post-composes and assigns
template <typename C>
Word<C>& Word<C>::operator *= (const Word<C>& other)
{
	typename WORD::const_iterator othlink;
	for (othlink = other.begin(); othlink!=other.end();++othlink)
	{
		addLetterBack(othlink->first,othlink->second);
	}
	return *this;
}

// returns result of post composition without changing either operand
template <typename C>
Word<C> Word<C>::operator * (const Word<C>& other) const
{
	Word<C> toret = *this;
	toret *= other;
	return toret;
}

// abelianizes into chain, stored in result!
template <typename C>
void Word<C>::abel(Chain<C>& result)
{
	result.clear();

	// iterate over word:
	typename WORD::const_iterator letter;
	for (letter = begin(); letter != end(); ++letter)
	{
		result.addLink(letter->first,letter->second);
	}
}

*/
#endif /* CWCELL_HPP_ */
