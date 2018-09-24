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

#ifndef SHEKELYAN_UTILGRID_H
#define SHEKELYAN_UTILGRID_H

#include <approxdata/utils/utils.h>


class LongCoords{

	private:
	
		const int DIMS;
		
		vector<long> dimMasks;
		vector<int> dimOffsets;
		
		long oneCoord;
		

	public:
		
		inline LongCoords(int dims) : DIMS(dims), dimMasks(dims), dimOffsets(dims){
		
			int pos = 0;
			
			const int w = 63/DIMS;
			
			for (int i = 0; i < DIMS; i++){
			
				dimMasks[i] = UtilBits::mask(pos, w);
				dimOffsets[i] = pos;
 				
 				pos += w;
			}
			
			oneCoord = 0;
			
			for (int i = 0; i < DIMS; i++)
				oneCoord = setDimCoord(oneCoord, i, 1);
		}
		
		inline int getDims(){
			
			return DIMS;
		}
		
		inline long decAllCoords(long x){
		
			for (int i = 0; i < DIMS; i++)
				if ( (x & dimMasks[i]) == 0)
					return -1;
			
			return x-oneCoord;
		}
		
		inline long incAllCoords(long x){
			
			for (int i = 0; i < DIMS; i++)
				if ( (x & dimMasks[i]) == dimMasks[i])
					return -1;
			
			return x+oneCoord;
		}
		
		inline bool contains(long imin, long imax, long x){
		
			if (imin == -1)
				return false;
				
			if (imax == -1)
				return false;
		
			for (int i = 0; i < DIMS; i++){
			
				const long m = dimMasks[i];
				
				const long xm = x & m;
				
				if ( xm < (imin & m) )
					return false;
				
				if ( xm > (imax & m) )
					return false;
			}
			
			return true;
		}
		
		inline long getDimCoord(long x, int i){
		
			return (x & dimMasks[i]) >> dimOffsets[i];
		}
		
		inline long setDimCoord(long x, int i, long v){
		
			return (v << dimOffsets[i]) |UtilBits::setZeroInMask(x, dimMasks[i]);
		}
		
		inline void getCoords(long x, vector<long>& ind){
		
			for (int i = 0; i < DIMS; i++)
				ind[i] = getDimCoord(x, i);
		}
		
		inline long getIndex(long x, long y){
		
			return (x << dimOffsets[0]) | (y << dimOffsets[1] ); 
		}
		
		inline long getIndex(long x, long y, long z){
		
			return (x << dimOffsets[0]) | (y << dimOffsets[1]) | (z << dimOffsets[2]); 
		}
		
		inline long getIndex(vector<long>& ind){
		
			long ret = 0;
		
			for (int i = 0; i < DIMS; i++){
			
				ret |= ind[i] << dimOffsets[i];
			}	
			
			return ret;
		}


};

template<bool rowMajor=true, bool bitshift=false>
class MultiDimIndex{

private:

	vector<long> lens;	
	vector<int> shifts;	
	vector<long> muls;
	
	long totalLength = 1;
	
	inline void prepare(){
		
		if (bitshift){
			for (int i = 0; i < lens.size(); i++)
				lens[i] = UtilMath::powerOfTwoUB(lens[i]);
		}
		
		long mul = 1;
		
		if (rowMajor){
			
			for (int i = lens.size()-1; i >= 0; i--){
		
				muls[i] = mul;
				mul *= lens[i];
			}
			
		} else {
			
			for (int i = 0; i < lens.size(); i++){
		
				muls[i] = mul;
				mul *= lens[i];
			}
		}
		
		totalLength = mul;
		
		if (bitshift){
			for (int i = 0; i < lens.size(); i++)
				shifts[i] = UtilMath::getPowerOfTwo(muls[i]);
		}
		
		//cout << UtilString::listToString(muls) << endl;
	}

public:

	inline MultiDimIndex(vector<long>& lengths): lens(lengths.size()), muls(lengths.size() ), shifts(lengths.size()){
	
		lens = lengths;
		
		prepare();
	}
	
	inline int getDims(){
	
		return lens.size();
	}	
	
	inline long getLength(int i){
	
		return lens[i];
	}
	
	inline long getLength(){
	
		return totalLength;
	}
	
	inline MultiDimIndex(int dims): lens(dims), muls(dims), shifts(dims){
	
		for (int i = 0; i < dims; i++){
		
			lens[i] = 1L << (62/dims);
		}
		
		prepare();
	}
	
	inline MultiDimIndex(long x): lens(1), muls(1), shifts(1){
	
		lens[0] = x;
		
		prepare();
	}
	
	inline MultiDimIndex(long x, long y): lens(2), muls(2), shifts(2){
	
		lens[0] = x;
		lens[1] = y;
		
		prepare();
	}
	
	inline MultiDimIndex(long x, long y, long z): lens(3), muls(3), shifts(3){
	
		lens[0] = x;
		lens[1] = y;
		lens[2] = z;
		
		prepare();
	}
	
	
	inline bool contains(long max, long x){
		
		if (rowMajor){
		
			if (bitshift){
			
				for (int i = 0; i < lens.size(); i++){
		
					const long maxi = max >> shifts[i];
					const long xi = x >> shifts[i];
					
					if (xi > maxi)
						return false;
					
					max -= maxi << shifts[i];
					x -= xi << shifts[i];
				}
				
			} else {
			
				for (int i = 0; i < lens.size(); i++){
		
					const long maxi = max/muls[i];
					const long xi = x/muls[i];
					
					if (xi > maxi)
						return false;
					
					max -= maxi * muls[i];
					x -= xi * muls[i];
				}
			}
		
		} else {
		
			if (bitshift){
			
				for (int i = lens.size()-1; i >= 0; i--){
		
					const long maxi = max >> shifts[i];
					const long xi = x >> shifts[i];
					
					if (xi > maxi)
						return false;
					
					max -= maxi << shifts[i];
					x -= xi << shifts[i];
				}
				
			} else {
				
				for (int i = lens.size()-1; i >= 0; i--){
		
					const long maxi = max/muls[i];
					const long xi = x/muls[i];
					
					if (xi > maxi)
						return false;
					
					max -= maxi * muls[i];
					x -= xi * muls[i];
				}
			}	
		}
		
		return true;
	}
	
	inline bool contains(long min, long max, long x){
		
		if (rowMajor){
		
			if (bitshift){
			
				for (int i = 0; i < lens.size(); i++){
		
					const long mini = min >> shifts[i];
					const long maxi = max >> shifts[i];
					const long xi = x >> shifts[i];
					
					if (xi > maxi || xi < mini)
						return false;
					
					min -= mini << shifts[i];
					max -= maxi << shifts[i];
					x -= xi << shifts[i];
				}
				
			} else {
			
				for (int i = 0; i < lens.size(); i++){
		
					const long mini = min/muls[i];
					const long maxi = max/muls[i];
					const long xi = x/muls[i];
					
					if (xi > maxi || xi < mini)
						return false;
					
					min -= mini * muls[i];
					max -= maxi * muls[i];
					x -= xi * muls[i];
				}
			}
		
		} else {
		
			if (bitshift){
			
				for (int i = lens.size()-1; i >= 0; i--){
		
					const long mini = min >> shifts[i];
					const long maxi = max >> shifts[i];
					const long xi = x >> shifts[i];
					
					if (xi > maxi || xi < mini)
						return false;
					
					min -= mini << shifts[i];
					max -= maxi << shifts[i];
					x -= xi << shifts[i];
				}
				
			} else {
				
				for (int i = lens.size()-1; i >= 0; i--){
		
					const long mini = min/muls[i];
					const long maxi = max/muls[i];
					const long xi = x/muls[i];
					
					if (xi > maxi || xi < mini)
						return false;
					
					min -= mini * muls[i];
					max -= maxi * muls[i];
					x -= xi * muls[i];
				}
			}	
		}
		
		return true;
	}
	
	inline long getIndex(long x){
	
		return x;
	}
	
	inline long getIndex(long x, long y){
	
		return rowMajor ? (bitshift ? ( y+(x << shifts[0]) ) : (y+x*muls[0]) ) : 
			(bitshift ? ( x+(y << shifts[1]) ) : (x+y*muls[1]) );
	}
	
	inline long getIndex(long x, long y, long z){
		
		return rowMajor ? (bitshift ? ( z+(y << shifts[1])+(x << shifts[0]) ) : (z+y*muls[1]+x*muls[0]) ):
			(bitshift ? ( x+(y << shifts[1])+(z << shifts[2]) ) : (x+y*muls[1]+z*muls[2]) );
	}
	
	inline long getIndex(LongCoords& lc, long coords){
	
		long ret = 0;
		
		if (bitshift){
		
			for (int i = 0; i < lens.size(); i++)
				ret += lc.getDimCoord(coords, i) << shifts[i];
		
		} else {
		
			for (int i = 0; i < lens.size(); i++)
				ret += lc.getDimCoord(coords, i)*muls[i];
		}
		
		return ret;
		
	}
	
	inline long getIndex(vector<long>& coords){
		
		long ret = 0;
		
		if (bitshift){
		
			for (int i = 0; i < lens.size(); i++)
				ret += coords[i] << shifts[i];
		
		} else {
		
			for (int i = 0; i < lens.size(); i++)
				ret += coords[i]*muls[i];
		}
		
		return ret;
	}

	inline long getCoords(long index, vector<long>& coords){
	
		coords = lens;
		
		if (rowMajor){
		
			if (bitshift){
			
				for (int i = 0; i < lens.size(); i++){
		
					coords[i] = index >> shifts[i];
					index -= coords[i] << shifts[i];
				}
			} else {
			
				for (int i = 0; i < lens.size(); i++){
		
					coords[i] = index/muls[i];
					index -= coords[i]*muls[i];
				}
			}
		
		} else {
		
			if (bitshift){
			
				for (int i = lens.size()-1; i >= 0; i--){
		
					coords[i] = index >> shifts[i];
					index -= coords[i] << shifts[i];
				}
				
			} else {
			
				for (int i = lens.size()-1; i >= 0; i--){
		
					coords[i] = index/muls[i];
					index -= coords[i]*muls[i];
				}
			}
		}
	}
};

template<typename E, bool FAST=false>
class MultiArray{
public:
	MultiDimIndex<true, FAST> m;
	
	vector<E> vec;
	long max;
	
	long wx;
	long wy;
	long wz;
	
	vector<int> lengths;
	
	inline MultiArray(long wx) : m(wx), wx(wx), wy(1), wz(1){
	
		vec.reserve(m.getLength());
		max = m.getIndex(wx-1);
		
		lengths.push_back(wx);
	}
	
	inline MultiArray(long wx, long wy) : m(wx,wy), wx(wx), wy(wy), wz(1){
		
		vec.reserve(m.getLength());
		max = m.getIndex(wx-1, wy-1);
		
		lengths.push_back(wx);
		lengths.push_back(wy);
	}
	
	inline MultiArray(long wx, long wy, long wz) : m(wx,wy,wz), wx(wx), wy(wy), wz(wz){
		
		vec.reserve(m.getLength());
		max = m.getIndex(wx-1, wy-1, wz-1);
		
		lengths.push_back(wx);
		lengths.push_back(wy);
		lengths.push_back(wz);
	}
	
	
	inline long getLengthX(){
	
		return wx;
	}
	
	inline long getLengthY(){
	
		return wy;
	}
	
	inline long getLengthZ(){
	
		return wz;
	}
	
	
	inline E& atIndex(long ind){
	
		return vec.at(ind);
	}
	
	inline void clear(){
	
		vec.clear();
	}
	
	
	
	inline void fill(E e){
		
		std::fill(vec.begin(), vec.end(), e);
		
		while (vec.size() < m.getLength() )
			vec.push_back(e);
	}
	
	inline long getLength(int i){
	
		return m.getLength(i);
	}
	
	inline E& at(long x){
	
		return vec.at( m.getIndex(x) );
	}
	
	inline E& at(long x, long y){
		
		return vec.at( m.getIndex(x,y) );
	}
	
	inline E& at(long x, long y, long z){
	
		return vec.at( m.getIndex(x,y,z) );
	}
	
	inline void push_back(E e){
		
		if (FAST){
		
			while (!m.contains(max, vec.size() ))
				vec.push_back(e);
		}
		
		vec.push_back(e);
	}
	
};


class MultipleChoiceKnapsack{

private:
	MultiArray<double> profits;
	MultiArray<double> weights;
	
	const long knapsacks;
	const long items;

public:

	inline MultipleChoiceKnapsack(long knapsacks, long items) : knapsacks(knapsacks), items(items), profits(knapsacks, items),  weights(knapsacks, items){
	
		profits.fill( std::numeric_limits<double>::lowest() );
		weights.fill( std::numeric_limits<double>::max() );
	}
	
	
	
	inline void add(int knapsack, int item, double profit, double weight){
	
		profits.at(knapsack, item) = profit;
		weights.at(knapsack, item) = weight;
	}
	
	inline vector<int> maximizeProfit(double weightConstraint){
	
		if (knapsacks <= 5 && items <= 100)
			return maximizeProfitBruteForce(weightConstraint);
	
		assert (false);
	}
	
	inline double getProfit(vector<int> v){
	
		double ret = 0;
	
		for (int i = 0; i < v.size(); i++)
			ret += profits.at(i, v[i]);
			
		return ret;
	}
	
	inline double getWeight(vector<int> v){
	
		double ret = 0;
	
		for (int i = 0; i < v.size(); i++)
			ret += weights.at(i, v[i]);
			
		return ret;
	}
	
	inline vector<int> maximizeProfitBruteForce(double weightConstraint){
	
		vector<int> ret(knapsacks, -1);
		
		double maxProfit = std::numeric_limits<double>::lowest();
		
		for (int j1 = 0; j1 < items; j1++){
	
			const double p1 = profits.at(0, j1);
			const double w1 = weights.at(0, j1);
		
			if (w1 > weightConstraint)
				continue;

			if (knapsacks == 1){
					
				if (p1 > maxProfit){
			
					ret[0] = j1;
					maxProfit = p1;
				}
				continue;
			}

			for (int j2 = 0; j2 < items; j2++){
		
				const double p2 = p1+profits.at(1, j2);
				const double w2 = w1+weights.at(1, j2);
			
				if (w2 > weightConstraint)
					continue;
				
				if (knapsacks == 2){
					
					if (p2 > maxProfit){
				
						ret[0] = j1;
						ret[1] = j2;
						maxProfit = p2;
					}
					continue;
				}
			
				for (int j3 = 0; j3 < items; j3++){
	
					const double p3 = p2+profits.at(2, j3);
					const double w3 = w2+weights.at(2, j3);
				
					if (w3 > weightConstraint)
						continue;
					
					if (knapsacks == 3){
					
						if (p3 > maxProfit){
					
							ret[0] = j1;
							ret[1] = j2;
							ret[2] = j3;
							maxProfit = p3;
						}
						continue;
					}
				
					for (int j4 = 0; j4 < items; j4++){
				
						const double p4 = p3+profits.at(3, j4);
						const double w4 = w3+weights.at(3, j4);
					
						if (w4 > weightConstraint)
							continue;
					
						if (knapsacks == 4){
							
							if (p4 > maxProfit){
					
								ret[0] = j1;
								ret[1] = j2;
								ret[2] = j3;
								ret[3] = j4;
								maxProfit = p4;
							}
							continue;
						}
						
						for (int j5 = 0; j5 < items; j5++){
				
							const double p5 = p4+profits.at(4, j5);
							const double w5 = w4+weights.at(4, j5);
					
							if (w5 > weightConstraint)
								continue;
							
							if (knapsacks == 5){
							
								if (p5 > maxProfit){
					
									ret[0] = j1;
									ret[1] = j2;
									ret[2] = j3;
									ret[3] = j4;
									ret[4] = j5;
									maxProfit = p5;
								}
								continue;
							}
						}//5
					}//4
				}//3	
			}//2
		}//1
		
		return ret;
	}
};


#endif