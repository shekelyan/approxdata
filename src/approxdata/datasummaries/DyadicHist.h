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

#ifndef SHEKELYAN_DYADICHIST_H
#define SHEKELYAN_DYADICHIST_H

#include <approxdata/datasummaries/datasummaries.h>

class DyadicHistCalc{

private:
	class DyadicStorageFunction : public RealFunction{
	public:
		const int DIMS ;
	
		inline DyadicStorageFunction(int dims): DIMS(dims){
		
		}
	
		inline double operator()(double eps1) const override{

			assert (DIMS != -1);

			const long e = ceil(log2(1/eps1));
			return UtilMath::choose<long>( e+DIMS-1, e)*ceil(1/eps1);
		}
	};
	
	class DyadicPrecisionFunction : public RealFunction{
	public:
		const int DIMS ;
	
		inline DyadicPrecisionFunction(int dims): DIMS(dims){
		
		}
		
		inline double operator()(double eps1) const override{

			assert (DIMS != -1);

			long l = 2*ceil( log2(1.0/eps1) )-2;
			
			long ret = 0;
			
			for (int i = 0; i <= DIMS-2; i++)
				ret += 4*round(pow(l*1.0, i*1.0));
			
			ret += 2*pow(l, DIMS-1);
			
			return ret*eps1;
		}
	};

public:

	const int DIMS;

	inline DyadicHistCalc(int dims) : DIMS(dims){
	
	
	}

	inline double eps1ToKb(double eps1) const{
	
		DyadicStorageFunction f(DIMS);
		
		const double kbPerVal = 8.0/1024.0;
		
		return f(eps1)*kbPerVal;
	}
	
	
	inline double kbToEps1(double kb) const{
	
		DyadicStorageFunction f(DIMS);
		
		const double kbPerVal = 8.0/1024.0;
		
		return f.inv(0, 1, kb/kbPerVal, false);
	}
	
	inline double eps1ToEps(double eps1) const{
	
		DyadicPrecisionFunction f(DIMS);
		
		return f(eps1);
	}
	
	inline double epsToEps1(double eps) const{
	
		DyadicPrecisionFunction f(DIMS);
		
		return f.inv(0, eps, eps, true);
	}
	
	inline double epsToKb(double eps) const{
	
		return eps1ToKb(epsToEps1(eps) );
	}
	
	inline double kbToEps(double kb) const{
	
		return eps1ToEps(kbToEps1(kb) );
	}
};


typedef QuantileSummary<Point, PointComparator> GKSummary;





class DyadicNode{

typedef unique_ptr<DyadicNode> Node;	
	
private:

	unique_ptr<GKSummary> quantiles;
	
	const int dim;
	const int DIMS;
	const double eps1;
	
	vector<Node> children;
	
	int maxLevel = -1;
	
	inline Node& getChild(int level, int c){
	
		assert (children.size() > 0);
	
		//assert ( (c << (maxLevel-level)) < quantiles->getQuantileNum() );
	
		return children.at((1L << level)-1+c );
	}
	
	inline const Node& getChild(int level, int c) const{
	
		assert (children.size() > 0);
		
		//assert ( (c << (maxLevel-level)) < quantiles->getQuantileNum() );
		
		return children.at((1L << level)-1+c );
	}

public:
	
	inline DyadicNode(int dims, int mainDim, double eps1) : DIMS(dims), dim(mainDim), eps1(eps1){
		
		assert( dim < dims);
	}
	
	
	inline bool isEmpty() const{
		
		return quantiles ? false : true;
	}
	
	inline bool hasNoChildren() const{
	
		return children.size() == 0;
	}
	
	inline void print() const{
	
		if (isEmpty())
			return;
	
		cout << UtilString::repeat("\t", dim) << "eps1 " << eps1 << " count " << quantiles->getCount() << " quantiles " << quantiles->getQuantileNum() << endl;
		
		for (int level = 0; level <= maxLevel; level++){
		
			for (int j = 0; j < (1L << level); j++){
			
				const Node& n = getChild(level, j);
				
				cout << UtilString::repeat("\t", dim+1) << "child(" << level << "," << j << ") ";
				
				if (n){
				
					n->print();
				
				} else {
				
					//cout << UtilString::repeat("\t", dim+1) << "empty" << endl;
				}
				
				cout << endl;
			}
		}
	}
	
	inline double intervalCount(const Point& qmin, const Point& qmax, QueryMode mode=QueryMode::EST) const{
	
		if (isEmpty())
			return 0.0;
	
		return quantiles->intervalCount(qmin, qmax, mode);
	}
	
	inline long getSize() const{
		
		if (isEmpty())
			return 0;
		
		long ret = quantiles->getQuantileNum()+1;
		
		for (auto it = children.begin(); it != children.end(); it++){
		
			const Node& n = (*it);
			
			if (n)
				ret += n->getSize();
			
			ret++;
		}
	
		return ret;
	}
	
	inline void createChildren(){
		
		if (isEmpty())
			return;
			
		const long count = quantiles->getCount();
	
		assert (count > 0);
		
		assert (dim < (DIMS-1));
	
		const long qnum = quantiles->getQuantileNum();
		assert (qnum > 0);
		assert (qnum <= (1L << maxLevel) );
		assert (qnum >= (1L << (maxLevel-1) ) );
		
		for (int level = 0; level <= maxLevel; level++){
		
			//cout << "level " << level << endl;
		
			for (int j = 0; j < (1L << level); j++){
			
				const int shift = maxLevel-level;
			
				const long qmin = j << shift;
				const long qmax = ((j+1) << shift)-1;
				
				//cout << "j " << j << " qmin " << qmin << " qmax " << qmax << endl;
				
				Node n;
				
				const bool isReachable = qmin < qnum;
				
				if (isReachable){
				
					const double countUB = quantiles->getQuantileCount(qmax, QueryMode::UB)-quantiles->getQuantileCount(qmin, QueryMode::LB);
					
					// eps1*count == eps2*countUB <=> eps2 = eps1*count/countUB 
					
					const double eps2 = eps1*count/UtilMath::makeBetween<double>(0, quantiles->getCount(), countUB);
					
					n.reset(new DyadicNode(DIMS, dim+1, eps2) );
				}
				
				assert (children.size() == (j-1+(1L << level) ) );
				
				children.push_back( std::move(n) );
				
				const Node& m = getChild(level, j);
				
				assert (!isReachable || (m));
			}
		}
	}
	
	
	inline void test(long total){
	
		if (isEmpty() ){
		
			assert( total == 0);
			
			return;
		}
		
		assert (total == quantiles->getCount() );
		
		if (hasNoChildren() )
			return;
	
		for (int level = maxLevel; level >= 0; level--){
			
			long sum = 0;
			
			
			for (int q = 0; q < (1L << level); q++){
			
				Node& n = getChild(level, q);
				
				if (n){
					
					if (n->quantiles){
					
						const long count = n->quantiles->getCount();
					
						n->test(count);
					
						sum += count;
					}
				}
			}
			
			if (sum != total)
				cout << "level " << level << " sum " << sum << " != " << total << " total" << endl;
		}
		
	}	
	
	inline void add(int d, const Point& p){
	
		if (d == dim){
		
			if (isEmpty() ){
		
				const Point domainMin(DIMS, 0 );
				const Point domainMax(DIMS, 1+1e-20);
				const PointComparator domainComp(dim, true);
				
				quantiles.reset( new GKSummary(eps1, dim, domainMin, domainMax, domainComp) );
			}
			
			quantiles->add(p);
			
		} else if (d > dim){
		
			assert (!isEmpty());
			
			long q = quantiles->getQuantile(p);
			
			for (int level = maxLevel; level >= 0; level--){
				
				Node& n = getChild(level, q);
				
				assert (n);
				
				if (n)
					n->add(d, p);
				
				q >>= 1;
			}
		}
	}
	
	inline void finalize(int d){
	
		if (isEmpty() )
			return;
		
		if (dim == d){
		
			quantiles->finalize();
			
			const long qnum = quantiles->getQuantileNum();
		
			maxLevel = UtilMath::getPowerOfTwoUB(qnum);
		
			
			if (d < (DIMS-1))
				createChildren();
			
			return;
		}
		
		if (d > dim){
			
			for (auto it = children.begin(); it != children.end(); it++){
			
				Node& n = (*it);
				
				if (n)
					n->finalize(d);
			}
		}
	}
	
private:	
	inline double childBoxCount(int level, long q, const Point& qmin, const Point& qmax, QueryMode mode) const{
	
		if (isEmpty() )
			return 0.0;
		
		if (hasNoChildren() ){
		
			Point a = qmin;
			Point b = qmax;
		
			const int shift = maxLevel-level;
		
			const long x = q << shift;
			const long y = ((q+1) << shift)-1;
		
			a[dim] = UtilMath::maxVal<double>(a[dim], quantiles->getQuantileMin(x)[dim]);
			b[dim] = UtilMath::minVal<double>(b[dim], quantiles->getQuantileMax(y)[dim]);
			
			return intervalCount(a, b, mode);
		}
		
		const Node& n = getChild(level, q);
		
		if (n){
			
			return n->boxCount(qmin, qmax, mode);
			
		} else {
		
			Point a = qmin;
			Point b = qmax;
		
			const int shift = maxLevel-level;
		
			const long x = q << shift;
			const long y = ((q+1) << shift)-1;
		
			a[dim] = UtilMath::maxVal<double>(a[dim], quantiles->getQuantileMin(x)[dim]);
			b[dim] = UtilMath::minVal<double>(b[dim], quantiles->getQuantileMax(y)[dim]);
		
		
			return intervalCount(a, b, mode);
		}
	}

public:	
	inline double boxCount(const Point& qmin, const Point& qmax, QueryMode mode) const{
		
		if (isEmpty() )
			return 0.0;
		
		//cout << qmin << " " << qmax << endl;
		
		assert (quantiles->getCount() > 0);
		
		if (dim == (DIMS-1))
			return intervalCount(qmin, qmax, mode);
		
		const long q1 = quantiles->getQuantile(qmin);
		const long q2 = quantiles->getQuantile(qmax);
		
		//cout << "q1 " << q1 << " " << q2 << " q2 " << endl;
		
		if (q1 == q2)
			return childBoxCount(maxLevel, q1, qmin, qmax, mode);
		
		double ret = 0;
		
		ret += childBoxCount(maxLevel, q1, qmin, qmax, mode);
		ret += childBoxCount(maxLevel, q2, qmin, qmax, mode);
		
		/*
		
		vector<bool> check(1L << maxLevel);
		
		assert(!check[q1]);
		check[q1] = true;
		
		assert(!check[q2]);
		check[q2] = true;
		
		*/
		
		//cout << "[" << q1 << ":" << q1 << "] += " << boxCount(maxLevel, q1, qmin, qmax, mode) << endl;
		//cout << "[" << q2 << ":" << q2 << "] += " << boxCount(maxLevel, q2, qmin, qmax, mode) << endl;
		
		long s = numeric_limits<long>::max();
		long t = numeric_limits<long>::lowest();
		
		for (int level = 0; level <= maxLevel; level++){
			
			const int shift = maxLevel-level;
			
			const long a = q1 >> shift;
			const long b = q2 >> shift;
			
			const long x = a+1;
			const long y = b-1;
			
			if (x > y)
				continue;
			
			const long xmin = x << shift;
			const long xmax = ((x+1) << shift)-1;
			
			if (xmax < s){
			
				/*
					for (long j = xmin; j <= xmax; j++){
						assert(!check[j]);
						check[j] = true;
					}
				*/
				
				ret += childBoxCount(level, x, qmin, qmax, mode);
				
				//cout << "[" << xmin << ":" << xmax << "] += " << boxCount(level, x, qmin, qmax, mode) << endl;
				
				if (xmin < s)
					s = xmin;
					
				if (xmax > t)
					t = xmax;
			}
			
			if (x == y)
				continue;
			
			const long ymin = y << shift;
			const long ymax = ((y+1) << shift)-1;
			
			if (ymin > t){
				
				ret += childBoxCount(level, y, qmin, qmax, mode);
				
				/*
					for (long j = ymin; j <= ymax; j++){
				
						assert(!check[j]);
						check[j] = true;
					}
				*/
				
				//cout << "[" << ymin << ":" << ymax << "] += " << boxCount(level, x, qmin, qmax, mode) << endl;
				
				if (ymin < s)
					s = ymin;
					
				if (ymax > t)
					t = ymax;
			}
		}
		
		/*
		for (long j = 0; j < q1; j++)
			assert (!check[j]);
		
		for (long j = q1; j <= q2; j++)
			assert (check[j]);
		
		for (long j = q2+1; j < (1L << maxLevel); j++){
		
			if (check[j])
				cout << "j " << j << " q2 " << q2 << endl;
		
			assert (!check[j]);
		}*/
		
		//cout << "ret " << ret << endl;
		
		assert(quantiles);
		assert (ret <= quantiles->getCount());
		
		return ret;
	}

};


class DyadicHist : public DataSummary{

public:
	
	unique_ptr<DyadicNode> root;
	Params params;
	
	double eps1;
	double eps;
	
	const int DIMS;
	
	DyadicHistCalc calc;
	
	inline long getSize() const{
	
		if (root)
			return root->getSize();
		else
			return 0;
	}
	
	const long datasize;
	
	inline DyadicHist( shared_ptr<DataSet> data, const string paramStr) : datasize(data->getSize() ), DIMS(data->getDims() ), calc(DIMS), params(paramStr){
	
		eps = params.get("eps");
		const double kb = params.get("kb");
		
		eps1 = 0;
		
		if (eps > 0){
			eps1 = calc.epsToEps1(eps);
			eps = calc.eps1ToEps(eps1);
		}
		
		if (kb > 0){
			eps1 = calc.kbToEps1(kb);
			eps = calc.eps1ToEps(eps1);
		}
		
		root.reset( new DyadicNode(data->getDims(), 0, eps1) );
		
		for (int i = 0; i < DIMS; i++){
		
			for (auto it = data->begin(); it != data->end(); it++)
				root->add(i, (*it) );
			
			root->finalize(i);
		}
	}
	
	inline double boxCount(const Point& qmin, const Point& qmax, QueryMode mode=QueryMode::EST) const override{
	
		return root->boxCount(qmin, qmax, mode);
	}
	
	inline void print() const{
	
		return root->print();
	}
	
	inline void test() const{
	
		root->test(datasize);
	}
	
	inline const string toString() const{
	
		stringstream s;
		
		const double kb = getSize()*(8.0/1024.0);
	
		s << eps << " eps";
		s << " ";
		s << eps1 << " eps1";
		s << " ";
		s << kb  << " kb";
		s << " ";
		s << calc.eps1ToKb(eps1) << " expected kb";
		s << " ";
		s << calc.kbToEps(kb) << " expected eps";
		s << " ";
		
		return s.str();
	}
	
	friend inline std::ostream& operator<<(std::ostream &os, const DyadicHist& m) { 
    	
    	return os << m.toString();
	}
};

#endif