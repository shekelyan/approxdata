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

#ifndef SHEKELYAN_BOXINDEX_H
#define SHEKELYAN_BOXINDEX_H

#include <approxdata/data/data.h>

class BoxIndex{

protected:

	vector< vector<double> > discLut;
	vector< vector< vector<int> >> intersectedLut;
	
	inline void prepare(){
	
		for (int i = 0; i < DIMS; i++){
		
			vector<double> v;
	
			for (auto it = boxes.begin(); it != boxes.end(); it++){
	
				v.push_back(it->min[i]);
				v.push_back(it->max[i]);
			}
		
			sort( v.begin(), v.end() );
			v.erase( unique( v.begin(), v.end() ), v.end() );
			
			discLut.push_back(v);
			
			vector<vector<int>> lutDim;
			
			for (long j = 0; j <= (v.size() << 1); j++){
			
				vector<int> b;
				
				long k = 0;
				for (auto it = boxes.begin(); it != boxes.end(); it++){
			
					const long minIndex = getIndex( it->min[i], i );
					const long maxIndex = getIndex( it->max[i], i );
					
					if (minIndex <= j && maxIndex >= j)
						b.push_back(k);
						
					k++;
				}
				
				lutDim.push_back(b);
			}
			
			intersectedLut.push_back(lutDim);	
		}
	}

	inline long largestIndex(int dim) const{
	
		return discLut[dim].size() << 1;
	}
	
	inline long getIndex(double val, int dim) const{
	
		const vector<double>& v = discLut[dim];
	
		if (v.size() == 0)
			return -1;
	
		if (val < v[0])
			return 0;
		
		if (val == v[0])
			return 1;
		
		if (val == v[v.size()-1])
			return (v.size() << 1)-1;
		
		if (val > v[v.size()-1])
			return v.size() << 1;
		
		auto low = std::lower_bound (v.begin(), v.end(), val);  // first >= 
		
		const long lb = low-v.begin();
		
		if ( (*low) == val)
			return 2*lb+1;
		else
			return 2*lb;
	}

public:
	
	const int DIMS;
	
	vector<Box>& boxes;
	
	const bool fast;
	
	inline BoxIndex(int dims, vector<Box>& boxesPtr, int mode) : boxes(boxesPtr), fast(mode != 0), DIMS(dims){
		
		if (fast)
			prepare();
			
		for (auto it = boxes.begin(); it != boxes.end(); it++)
			it->count = 0;
	}
	
	
	inline void print(){
		
	}
	
	inline void count(Point& p){
	
		if (!fast){
			
			for (auto it = boxes.begin(); it != boxes.end(); it++){
		
				if ( it->contains(p) )
					it->count++;
			}
			
			return;
		}
		
		vector<int>* b = NULL;
		
		for (int i = 0; i < DIMS; i++){
		
			vector<int>* c = &(intersectedLut[i][getIndex(p[i], i)]);
			
			if (b == NULL || c->size() < b->size() )
				b = c;
		}
		
		if (b && b->size() > 0)
		for (auto it = b->begin(); it != b->end(); it++){
			
			const long i = (*it);
			
			if ( boxes[i].contains(p) )
				boxes[i].count++;
		}
	}
	

};







#endif