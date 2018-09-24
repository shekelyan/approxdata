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
#ifndef SHEKELYAN_QUANTILESUMMARY_H
#define SHEKELYAN_QUANTILESUMMARY_H

#include <approxdata/external/gk.h>
#include <approxdata/datasummaries/datasummaries.h>

template<typename E, typename F=std::less<E>>
class QuantileSummary{

	private:
	
		unique_ptr< GK<E,F> > gk;
		vector<E> borders;
		
		long count = 0;
		
		E min;
		E max;
		
		F comp;
		
		vector<double> ranks;
		
		const int dim;
		const double eps; // maximal quantile selectivity
		
		const double rankError;
		
	public:
		
		inline QuantileSummary(double eps, int dim, const E& domainMin, const E& domainMax, const F& comp) : rankError(eps/3), dim(dim), comp(comp), eps(eps){
		
			if (rankError < 1)
				gk.reset( new GK<E,F>(rankError, domainMin, domainMax, comp) );
			
			min = domainMax;
			max = domainMin;
			count = 0;
		}
		
		inline const vector<E> getBorders() const{
		
			return borders;
		}
		
		inline long getCount() const{
		
			return count;
		}
		
		inline double getEps() const{
		
			return eps;
		}
		
		inline long getQuantile(const E& p) const{
		
			assert (borders.size() > 1);
			
			const long ret = UtilHist::getBucket(borders, p, comp);
			
			assert (ret >= 0);
			assert (ret < getQuantileNum());
			
			return ret;
		}
		
		inline const E& getQuantileMin(long q) const{
		
			return borders[q];
		}
		
		inline const E& getQuantileMax(long q) const{
		
			return borders[q+1];
		}
		
		inline long getQuantileNum() const{
		
			return borders.size()-1;
		}
		
		
		inline double getRankUB(long q, QueryMode mode) const{
			
			return UtilMath::makeBetween<double>(0, getCount(), ranks[q+1]+rankError);
		}
		
		
		inline double getQuantileCount(long q, QueryMode mode) const{
			
			const double lb = ranks[q];
			
			if (mode == QueryMode::LB)
				return UtilMath::makeBetween<double>(0, getCount(), lb-rankError*getCount() );
			
			const double ub = ranks[q+1];
			
			if (mode == QueryMode::UB)
				return UtilMath::makeBetween<double>(0, getCount(), ub+rankError*getCount() );
			
			return lb*0.5+ub*0.5;
		}
		
		inline double getRank(const E& p, QueryMode mode) const{
			
			const long q = getQuantile(p);
			
			//cout << "getRank p " << p << " q " << q << endl;
			
			if (true)
				return getQuantileCount(q, mode);
			
			const double a = borders[q][dim];
			const double b = borders[q+1][dim];
			
			const double countA = ranks[q];
			const double countB = ranks[q+1];
			
			if (mode == QueryMode::LB)
				return UtilMath::makeBetween<double>(0, getCount(), countA-rankError );
			
			if (mode == QueryMode::UB)
				return UtilMath::makeBetween<double>(0, getCount(), countB+rankError );
			
			const double v = p[dim];
			
			if (v < a)
				return countA;
				
			if (v >= b)
				return countB;
				
			const double frac = (v-a)/(b-a); 
			
			// countA+(countB-countA)*frac <=> countA*(1-frac)+countB*frac
				
			return countA*(1-frac)+countB*frac;
			
		}
		
		inline double intervalCount(const E& qmin, const E& qmax, QueryMode mode) const{
			
			//cout << "intervalCount qmin " << qmin << " " << qmax << " qmax" << endl;
			
			if ( mode == QueryMode::UB)
				return UtilMath::makeBetween<double>(0, getCount(), getRank(qmax, QueryMode::UB)-getRank(qmin, QueryMode::LB) );
		
			if ( mode == QueryMode::LB)
				return UtilMath::makeBetween<double>(0, getCount(), getRank(qmax, QueryMode::LB)-getRank(qmin, QueryMode::UB) );
			
			return UtilMath::makeBetween<double>(0, getCount(), getRank(qmax, QueryMode::EST)-getRank(qmin, QueryMode::EST) );
			
			/*
			const double rank1 = getRank(qmin, mode);
			const double rank2 = getRank(qmax, mode);
			
			//cout << "rank1 " << rank1 << " " << rank2 << " rank2" << endl;
			
			return UtilMath::makeBetween<double>(0, getCount(), rank2-rank1);*/
		}
		
		inline void add(const E& p){
		
			//cout << min << " min " << max << " max " << p << " p" << endl;
		
			if (comp(p,min))
				min = p;
			
			if (comp(max,p))
				max = p;
			
			if (gk)
				gk->feed(p);
			
			count++;
		}
		
		inline void finalize(){
			
			if (gk)
				gk->finalize();
			
			borders.push_back(min);
			ranks.push_back(0);
			
			if (gk){
			
				const double step = rankError*0.5;
			
				for (double rank = step; rank < (1-step); rank += step){
			
					const Point p = this->gk->query_for_value(rank);
				
					if (comp(p, min))
						continue;
					
					assert(p.size() == min.size());
			
					borders.push_back(p);
				
					ranks.push_back(rank*getCount() );	
				}
			}
			
			borders.push_back(max);
			ranks.push_back(getCount() );
			
			if (gk)
				gk.reset();
		}
		
		inline void print() const{
		
			for (auto it = borders.begin(); it != borders.end(); it++)
				cout << (*it) << endl;
		}
};

#endif