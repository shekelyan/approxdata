/*
 * ApproxData Library
 * Copyright (c) 2018 Michael Shekelyan <michael.shekelyan@gmail.com>
 */

/*
Permission is hereby granted, free of charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions 
of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING 
BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef SHEKELYAN_UTILCOLLECTION_H
#define SHEKELYAN_UTILCOLLECTION_H

#include <approxdata/utils/utils.h>




class ReservoirSampling{

	private:
		long k;
		
		std::unique_ptr<std::mt19937> g;
		std::seed_seq seed;

	public:
		const long sampleSize;
		
		inline ReservoirSampling(string seedStr, long size) : seed(seedStr.begin(), seedStr.end() ), sampleSize(size){
		
			g.reset(new std::mt19937(seed));
			this->k = 0;
		}
		
		inline void reset(){
	
			k = 0;
		}
		
		inline long pos(){
		
			if (k < sampleSize){
			
				const long j = k;
				++k;
				return j;
			}
			
			std::uniform_int_distribution<long> d(0, k);
			const long j = d( (*g) );
			
			++k;
			return j < sampleSize ? j : -1;
		} 

};

namespace UtilCollection{



	inline shared_ptr< vector<long> > getRandomIndices(long srcSize, long dstSize, string seed){
	
		assert (dstSize < srcSize);
	
		shared_ptr< vector<long> > v(new vector<long>(dstSize, -1) );
	
		ReservoirSampling rs(seed, dstSize);
		
		rs.reset();
		
		for (long k = 0; k < srcSize; k++){
		
			const long j = rs.pos();
			
			if ( (j >= 0) && (j < dstSize) )
				v->at(j) = k;
		}
		
		std::sort(v->begin(), v->end() );
		
		assert (v->at(0) != -1);
		
		return v;
	}
	
	
	/*
	template<typename E>
	inline long binarySearchGreaterOrEqual(vector<E>* vec, E val){
	
		if (vec->size() == 0)
			return 0L;
	
		auto it = std::lower_bound(vec->begin(), vec->end(), val);
		return it == vec->end() ? vec->size() : it-vec->begin();
	}
	
	template<typename E>
	inline long binarySearchGreater(vector<E>* vec, E val){
	
		if (vec->size() == 0)
			return 0L;
	
		auto it = std::upper_bound(vec->begin(), vec->end(), val);
		return it == vec->end() ? vec->size() : it-vec->begin();
	}
	
	template<typename E>
	inline long binarySearchLessOrEqual(vector<E>* vec, E val){
	
		return binarySearchGreater(vec, val)-1L;
	}
	
	template<typename E>
	inline long binarySearchLess(vector<E>* vec, E val){
	
		return binarySearchGreaterOrEqual(vec, val)-1L;
	}*/
	
};

struct CompObj { 
  
	double key;
    long ind;

	
	bool operator<(const CompObj& other) const{
		return key < other.key;
    }
};

class Sorter{
public:
	vector<CompObj> vec;
	
	bool sorted = false;
	
	inline void clear(){
	
		vec.clear();
		sorted = true;
	}
	
	inline void reserve(long size){
	
		vec.reserve(size);
	}
	
	inline void add(long e, double key){
		
		CompObj obj;
		obj.ind = e;
		obj.key = key;
		vec.push_back(obj);
		
		sorted = false;
	}
	
	inline void sort(){
	
		if (sorted)
			return;
	
		std::sort(vec.begin(), vec.end() );
		sorted = true;
	}
	
	inline long getSortedIndex(long ind){
	
		sort();
		
		return vec[ind].ind;
	}
	
	template<typename E>
	inline shared_ptr<vector<E>> getSorted(vector<E>& v){
	
		sort();
		
		shared_ptr<vector<E>> ret(new vector<E>() );
		
		ret->reserve(v.size() );
		
		for (auto it = vec.begin(); it != vec.end(); it++)
			ret->push_back( v[it->ind] );
		
		return ret;
	}
	
	template<typename E>
	inline void sort(vector<E>& v){
		
		sort();
		
		shared_ptr<vector<E>> sorted = sorted(v);
		v = (*sorted);
	}
};

template<typename E>
class SkipLUT{

private:

	long nonEmpty;

public:
	vector<long> l;
	vector<long> u;
	
	vector<long> rev;
			
	
	inline SkipLUT(vector<E>& v) : l(v.size()), u(v.size() ){
		
		long ind = -1;
		
		for (long j = 0; j < v.size(); j++){
			
			if (v[j] == 0){
			
				l[j] = ind;
				u[j] = ind+1;
				
			} else {
			
				ind++;
				l[j] = ind;
				u[j] = ind;
			}
		}
		
		nonEmpty = ind+1;
		
		rev.reserve(nonEmpty);
		
		for (long j = 0 ; j < v.size(); j++){
			
			if (!isEmpty(j))
				rev.push_back(j);
		}
	}
	
	inline bool isEmpty(long k){
	
		return l[k] != u[k];
	}
	
	inline long getStart(long x){
	
		return u[x];
	}
	
	inline long getSmallIndex(long x){

		if (u[x] != l[x])
			return -1;
	
		return l[x];
	}
	
	inline long getLargeIndex(long x){
	
		return rev[x];
	}
	
	inline long getEnd(long x){
	
		return l[x];
	}
	
	inline long nonEmptyNum() const{
	
		return nonEmpty;	
	}
	

	

};


#endif