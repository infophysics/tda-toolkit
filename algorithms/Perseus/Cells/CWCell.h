/*
 * CW.h
 *
 * CW cell class definitions
 */

#ifndef CWCELL_H_
#define CWCELL_H_
/*
#define WORD list<pair<D,C> >

#include <list>
#include <utility>
#include "Cell.h"

#define CWORD Word<C,CWCell<C,BT>*>

template <typename C, typename D> class Word;

// C is an Abelian group,
template <typename C, typename BT>
class CWCell : public Cell<C,BT>
{
public:

	// DATA
	CWORD attmap; // attaching map to lower cells

	// FUNTIONS
	void print(ostream& out); // output helper


};


// this class is a word in Cell pointers collected with coefficients
// from integral domain C. We intend to use C = integers to represent
// the degree of each attaching map
template <typename C, typename D>
class Word
{
public:

	// DATA

	// word data stored as a list of <X,C> pairs, with each X getting
	// paired with its coefficient in C. Note we do not Abelianize, so
	// for example the list <A,3><B,-2><A,-3> will not simplify to
	// <B,-2>, but <A,3><A,-2> is in fact <A,1>.
	WORD word;

	// FUNCTIONS

	// default constructor
	Word()
	{
		word.clear();
	}
	// copy constructor
	Word(const Word<C,D>& other)
	{
		word = other.word;
	}

	int getSize() const;
	C getCSum() const;

	// iterators for word structure
	typename WORD::const_iterator begin() const;
	typename WORD::const_iterator end() const;
	typename WORD::const_reverse_iterator rbegin() const;
	typename WORD::const_reverse_iterator rend() const;

	// access first and last links
	pair<C,D>& front() const;
	pair<C,D>& back() const;


	// print function
	friend ostream& operator << (ostream& out, const Word<C,D>& toprint)
	{
		typename WORD::const_iterator link;
		for (link = toprint.begin(); link != toprint.end(); ++link)
		{
			out<<"<"<<*(link->first)<<","<<link->second<<">";
		}
		return out;
	}

	void addLetterFront(const D&, const C&);
	void addLetterBack(const D&, const C&);

	// composition operators, NOT ABELIAN
	Word<C> operator * (const D&) const; // compose
	Word<C>& operator *= (const D&); // compose and assign
	bool operator == (const D&) const; // equality check


	// abelianize into chain!
	void abel(CCHAIN&);

	//destructor
	~Word()
	{
		word.clear();
	}
};

*/
#endif /* CWCELL_H_ */
