/*
 * Combinatorics.h
 *
 *  Created on: Dec 31, 2010
 *      Author: Vidit
 */

#ifndef COMBINATORICS_H_
#define COMBINATORICS_H_

# include "Size.h"
# include "Debuggers.h"
# include <algorithm>

# define dtype Cell<C,BT>*


// order nums increasingly
struct incnum
{
	bool operator() (num n1, num n2)
	{
		return (n1 < n2);
	}
};


// order nums decreasingly
struct decnum
{
	bool operator() (num n1, num n2)
	{
		return (n1 > n2);
	}
};



// create a 3d vector containing "datatypes" with allocated sizes for 1st and 2nd dims
template <typename datatype>
void make3Dvector(num dim1_size, num dim2_size, vector<vector<vector<datatype>*>*>& ans)
{
	for (num i = 0; i < dim1_size; ++i)
	{
		ans.push_back(new vector<vector<datatype>*>);
		for(num j=0; j < dim2_size; ++j)
		{
			ans.at(i)->push_back(new vector<datatype>);
		}
	}
	// now, ans.at(i)->at(j) is a vector of datatypes for every 0 <= i < dim1_size and 0 <= j < dim2_size

}

// given a SORTED vector, returns the position of the argument
template <typename datatype>
num getVPos(const vector<datatype>& svec, const datatype toget)
{
	typename vector<datatype>::const_iterator locate = lower_bound(svec.begin(),svec.end(),toget);
	return (num)(locate - svec.begin());
}


// returns C(n,k)
num nChoosek(const num n, num k)
{
	if (k > n || n <= 0) return 0;
	double ans = 1;
	//double numer = 1, denom = 1; // numerator and denominator

	// use symmetry!
	if (k > n/2) k = n-k;

	for (num i = 1; i <= k; i++)
	{
		ans *= (n-i+1);
		ans /= i;
		//cout<<"\nn-i+1 "<<n-i+1<<" div i "<<i<<" runans "<<runans;
	}
	return (num)(ans);
}


// Let [n] = {0,...,n-1}, and let [n]_k be subsets of [n] of size k.
// there are n-choose-k distinct subsets. pick any 1 <= p < n-choose-k.
// arrange [n]_k lexicographically. What is the p-th element? This
// solution comes from Buckles' algorithm.
bool lexicoPos (const num n, const num k, const num p, vector<num>& ans)
{

	if (p > nChoosek(n,k) || p < 0) return false;

	ans.clear();
	num combs, curpos = 0;

	if (n <= 0 || k > n) return false;

	if (k==0) return true;
	if (k==1)
	{
		//cout<<"one!"; cin.get();
		ans.push_back(p-1);
		return true;
	}


	for (num i=0; i < k-1; i++)
	{
		//cout<<"\n i: "<<i; cin.get();
		ans.push_back((i==0) ? -1 : ans.at(i-1));
		do
		{
			ans.at(i)++;
			combs = nChoosek(n-ans.at(i)-1,k-(i+1));
			//cout<<"\n   n "<<n-ans.at(i)<<" ch k "<<k-(i+1);
			curpos += combs;
			//cout<<"\n     combs: "<<combs<<" curpos "<<curpos; cin.get();
		} while(curpos < p);
		curpos -= combs;
	}
	ans.push_back(ans.at(k-2) + p - curpos);
	// fix start at 1 instead of 0:
	//for(num i=0; i<(num)ans.size(); i++) ans.at(i)--;
	return true;
}

// given a vector <a,b,c,...>, returns the product abc...
template <typename T>
T vProd (const vector<T> tomult)
{
	T ans = 1;
	for (num i = 0; i < (num)tomult.size(); ++i)
	{
		ans *= tomult.at(i);
	}
	return ans;
}


// reverses lexicopos. Given [n] = {0,...,n-1} and a subset
// [n]_k of size k, what is the rank of [n]_k among all k-sized
// subsets of [n] arranged lexicographically?
num lexicoPosRev(num n, const vector<num> subset)
{
	num position = 0; // default answer, signals error if size of subset is > 0!
	num k = subset.size();
	if (k > n || k == 0) return position; // can't have too many elements
	else if (k == 1) return subset.at(0);

	num j = 0;
	// otherwise, parse the subset,... increment position appropriately

	num enn, kay, first;

	first = subset.at(0);
	while (first > 0)
	{
		position += nChoosek(n-first,k-1);
		first--;
	}


	for (num i=0; i < k-2; ++i)
	{
		j = subset.at(i+1)-subset.at(i)-1;
		while (j>0)
		{
			enn = n - (subset.at(i) + j + 1);
			kay = k - i - 2;
			//cout<<"\n       diff i "<<i<<" j "<<j;
			//cout<<"\n               p = "<<position<<", adding ["<<enn<<"] C ["<<kay<<"] = "<<nChoosek(enn,kay);
			position += nChoosek(enn,kay);
			j--;
		}
	}
	position += subset.at(k-1) - subset.at(k-2) - 1;

	return position;
}

// extents is the dimension boundary vector,
// n is a single number. Let D be the size
// of extents. If we now lexicographically
// fill up a tensor of dimension extents(1)x
// extents(2) x...x extents(D) with numbers
// 0,1,...,n, what are the coordinates of n?
// Ask decompose: class T should be integral
template<class T, class S> void decompose(T n, const vector<T>& extents, vector<S>& ans)
{
	ans.clear();

	for (int i=0; i<(int)extents.size(); i++)
	{
		ans.push_back(S(n%extents.at(i)));
		n /= extents.at(i);
	}
}

// feed this (10,5,2) and it will spit out (1,10,50), basically
// ans_i = \prod [input_{j<i}]. this is useful when recomposing
// vector coordinates into a single index.

template <typename T>
void makeProd (const vector<T>& exts, vector<T>& ans)
{
	num d = exts.size();
	ans.assign(d,(T)1); // fill with ones
	for(num i = 0; i < d; i++)
	{
		for (num j = i+1; j < d; j++)
		{
			ans.at(j) *= exts.at(i);
		}
	}
}




// inverts decompose provided the same
// extents vector is used to make products.
// T must be an integer or long
template <typename  T>
T recompose(const vector<T>& coords, const vector<T>& exts)
{
	// get product:
	vector<T> prods;
	makeProd(exts,prods);

	// use product to recompose
	num dim = prods.size(); num ans = 0;
	// error!
	//cout<<"\n\n pdim: "<<dim<<" cdim "<<coords.size();
	//cout<<" prod: "; vprint<num>(prods);

	if (dim != (num)coords.size()) return ans;
	for (num i = 0; i < dim; i++)
	{
		ans += coords.at(i) * prods.at(i);
	}
	return ans;
}


#endif /* COMBINATORICS_H_ */
