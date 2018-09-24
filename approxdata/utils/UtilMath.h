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

#ifndef SHEKELYAN_UTILMATH_H
#define SHEKELYAN_UTILMATH_H

#include <approxdata/utils/utils.h>



class Stats{

public:
	vector<double> v;

	bool sorted = false;

	const string name;
	const string unit;
	
	inline Stats(const string s="", const string u="") : name(s), unit(u){
	
		
	}

	inline void add(double p){
	
		v.push_back(p);
		sorted = false;
	}
	
	inline void sort(){
	
		if (sorted)
			return;
		
		std::sort(v.begin(), v.end() );
		sorted = true;
	}
	
	inline void reset(){
	
		v.clear();
		sorted = true;
	}
	
	inline double sum() const{
	
		double ret = 0;
	
		for (auto it = v.begin(); it != v.end(); it++)
			ret += (*it);
			
		return ret;
	}
	
	inline double max() const{
	
		assert (v.size() > 0);
	
		if (sorted)
			return v[v.size()-1];
		
		double ret = v[0];
		
		for (auto it = v.begin(); it != v.end(); it++)
			if ( (*it) > ret)
				ret = (*it);
		
		return ret;
	}
	
	inline double min() const{
	
		assert (v.size() > 0);
		
		if (sorted)
			return v[0];
	
		double ret = v[0];
		
		for (auto it = v.begin(); it != v.end(); it++)
			if ( (*it) < ret)
				ret = (*it);
		
		return ret;
	}
	
	inline double avg() const{
		
		return sum()/v.size();
	}
	
	inline long num() const{
		
		return v.size();
	}
	
	inline double median(double q = 0.5) const{
	
		assert (v.size() > 0);
		assert (sorted);
		
		return v[ round((v.size()-1)*q) ];
	}
	
	inline const string toString() const{
	
		stringstream s;
		
		s << name;
		s << " \t";
		s << "num " << num();
		s << " \t";
		s << "min " << min() << unit;
		s << " \t";
		s << "avg " << avg() << unit;
		s << " \t";
		s << "max " << max() << unit;
		
		return s.str();
	}
	
	friend inline std::ostream& operator<<(std::ostream &os, const Stats& m) { 
    	
    	return os << m.toString();
	}
	
	
	
	
};



class Interval{


public:
	double a;
	double b;
	
	bool aIncluded = true;
	bool bIncluded = true;
	
	inline Interval intersection(Interval& v){
	
		Interval ret;
		
		ret.a = a > v.a ? a : v.a;
		ret.b = b < v.b ? b : v.b;
		
		return ret;
	}
};


namespace UtilHist{

	// input: real "v" between 0 and 1
	// output: integer "discretize(v, b)" between 0 and (b-1)
	
	
	
	// find i in [0, v.size()-2] s.t. p >= v[i] && p < v[i+1]
	
	template<typename E, typename F = std::less<E>>
	inline long getBucket(vector<E> v, E p, F comp=F() ){
		
		auto it = std::upper_bound(v.begin(), v.end(), p, comp);
		
		if (it == v.begin())
			return 0;
		
		if (it == v.end() )
			return v.size()-2;
		
		return (it-v.begin() )-1;
	}
	
	inline long discretizePowerOfTwo(double v, int lod){
		
		if (v <= 0)
			return 0;
		
		if (v >= 1)
			return ((1L << lod)-1);
		
		//return ( (long) (v*(1L << lod)) );
		
		return ( (long) (v*(1L << 62)) ) >> (62-lod);
	}
	
	
	inline long discretize(double v, long buckets){
	
		assert (buckets > 0);
		
		if (v <= 0)
			return 0;
		
		if (v >= 1)
			return buckets-1;
		
		return (long) (v*buckets);
	}

	inline double bucketMin(double bucket, double buckets){
		
		return (0.0+bucket)/buckets;
	}
	
	inline double bucketMax(double bucket, double buckets){
		
		return (1.0+bucket)/buckets;
	}

};

class RealFunction{

public:
	inline virtual double operator()(const double x) const = 0;
	inline virtual ~RealFunction(){};
	
	
	inline double inv(double domainMin, double domainMax, double targetY, const bool inc=true, int steps=1000){
		
		const RealFunction& f = (*this);
		
		const double a = domainMin*0.9+domainMax*0.1;
		const double b = domainMin*0.1+domainMax*0.9;
		
		//const bool monotonicallyIncreasing = f(a) < f(b);
		
		double xmin = domainMin;
		double xmax = domainMax;
		
		//if (inc)
		//	cout << "f(" << a << ") = " << f(a) << " < " << f(b) << " = f(" << b << ")" << endl;
		//else
		//	cout << "f(" << a << ") = " << f(a) << " >= " << f(b) << " = f(" << b << ")" << endl;
		
		for (int i = 0; i < 1000; i++){
		
			double x = xmin*0.5+xmax*0.5;
			
			const double y = f(x);
			
			//cout << "f(" << x << ") = " << y << endl;
			
			if (y == targetY)
				return x;
				
			if (inc){
			
				if (y > targetY){
				
					xmax = x;
					
				} else {
				
					xmin = x;
				}
				
			} else {
			
				if (y < targetY){
				
					xmax = x;
					
				} else {
				
					xmin = x;
				}
			}
		}
		
		return inc ? xmin : xmax;
	}
};




enum class SpatialRelation : int{ NO_INTERSECTION = 0, CONTAINED=2, OVERLAP=3, CONTAINS=4, EQUAL=5};

namespace UtilMath{



	template <class T>
	inline T choose(T n, T k, T limit=0){
	
		if (k > n)
			return 0;
		
		if (k*2 > n)
			k = n-k;
		
		if (k == 0)
			return 1;
	
		T ret = n;
		
		for (int i = 2; i <= k; ++i){
		
			ret *= n-i+1;
			ret /= i;
			
			if (limit != 0 && ret > limit)
				return ret;
		}
		
		return ret;
	}

	inline bool isPowerOfTwo(long x){
	
		return x && !(x & (x - 1));
	}
	
	inline int getPowerOfTwo(long x){
    	
    	int ret = -1;
    	
		for ( ; x != 0; x >>= 1)
			ret++;
		
		return ret;
	}
	
	inline SpatialRelation intervalIntersection(double amin, double amax, double bmin, double bmax){
	
		if (amin == bmin && amax == bmax)
			return SpatialRelation::EQUAL;
	
		if (amin <= bmin && amax >= bmax)
			return SpatialRelation::CONTAINS;
	
		if (bmin <= amin && bmax >= amax)
			return SpatialRelation::CONTAINED;
			
		if (amin <= bmin && amin >= bmax)
			return SpatialRelation::OVERLAP;
			
		if (amax <= bmin && amax >= bmax)
			return SpatialRelation::OVERLAP;
			
		return SpatialRelation::NO_INTERSECTION;		
	}
	

	inline double integrate(RealFunction& f, double a, double b, int steps=20){
	
		if (b == a)
			return 0;
	
		const double step = (b-a)/(steps-1);
	
		const double x1 = a+step;
		const double x2 = b-step;
		
		double sum = 0.5*f(a);
		
		for (double x = x1; x <= x2; x += step)
			sum += f(x);
		
		sum += 0.5*f(b);
		
		return step*sum;
	}
	
	template<typename E>
	inline E minVal(E a, E b){
	
		return a < b ? a : b;
	}
	
	template<typename E>
	inline E maxVal(E a, E b){
	
		return a > b ? a : b;
	}
	
	template<typename E>
	inline E absVal(E a){
	
		return a < 0 ? -a : a;
	}
	
	
	template<typename E>
	inline E makeBetween(E a, E b, E c){
		
		if (c < a)
			return a;
		
		if (c > b)
			return b;
		
		return c;
		//return maxVal(a, minVal(b, c));
	}
	
	inline double makeMultipleOf(double v, long m){
	
		if (m == 0)
			return v;
	
		return (((long) v)/m)*m;
	}
	
	inline double truncate(double min, double max, double x){
	
		if (x < min)
			return min;
		
		if (x > max)
			return max;
		
		return x;
	}
	
	inline long powerOfTwoLB(long m){
	
		long ret = (1L << 62);	
		
		if (m <= 0)
			return 0;
		
		while (ret != 0){
		
			if (ret <= m)
				return ret;
				
			ret >>= 1;
		}
		
		return 0;
	}
	
	inline long powerOfTwoUB(long m){
		
		if (m <= 0)
			return 0;
		
		for (long ret = (1L << 62); ret != 0; ret >>= 1){
		
			if (ret == m)
				return ret;
			
			if (ret & m)
				return ret << 1;
		}
		
		return 0;
	}
	
	inline int getPowerOfTwoUB(long x){
    	
    	
		return getPowerOfTwo(powerOfTwoUB(x));
	}
	
	template<typename E>
	inline double sum(const vector<E>& vec){
	
		double ret = 0;
		for (auto it = vec.begin(); it != vec.end(); it++)
			ret += (*it);
		
		return ret;
	}
	
	
	
	
	inline void optimalRounding2(vector<double>& v){

		const double vsum = sum<double>(v);
		double sumOfFloors = 0;
		
		for (long j = 0; j < v.size(); j++)
			sumOfFloors += floor(v[j]);
		
		const long roundDown = v.size()-(round(vsum)-sumOfFloors);
		
		Sorter s;
		
		for (long j = 0; j < v.size(); j++)
			s.add(j, v[j]-floor(v[j]) );
		
		for (long k = 0; k < roundDown; k++){
			
			const long j = s.getSortedIndex(k);
			v[j] = floor(v[j]);
		}
		
		for (long k = roundDown; k < v.size(); k++){
	
			const long j = s.getSortedIndex(k);
			v[j] = ceil(v[j]);
		}
		
		const double a = round(sum<double>(v));
		const double b = round(vsum);
		
		if (a != b)
			cout << "optimalRounding2 " << a << " != " << b << endl;
		
		//assert ( sum(v) == vsum);
	}
	
	
	
	inline void optimalRounding(vector<double>& v){
	
		const double a = round(sum<double>(v));
	
		vector<double> v2 = v;
		
		for (long j = 0; j < v2.size(); j++)
			if (v2[j] == ((long) v2[j]) )
				v2[j] = 0;
		
		SkipLUT<double> skip(v2);
		
		vector<double> v3(skip.nonEmptyNum() );
		
		for (long j = 0; j < v.size(); j++){
		
			if (!skip.isEmpty(j))
				v3[skip.getSmallIndex(j)] = v[j];
		}
		
		optimalRounding2(v3);
		
		for (long k = 0; k < v3.size(); k++){
		
			v[skip.getLargeIndex(k) ] = v3[k];
		}
		
		const double b = round(sum<double>(v));
		
		if (a != b)
			cout << "optimalRounding " << a << " != " << b << endl;
		
	}
};



enum class QueryMode : int{ LB=-1, EST=0, UB=1};


class OneDimEstimator{
public:

	virtual ~OneDimEstimator(){};
	inline virtual double intervalCount(double qmin, double qmax, QueryMode mode=QueryMode::EST) const = 0;

protected:
	OneDimEstimator(){}
    OneDimEstimator(const OneDimEstimator&){}
    OneDimEstimator& operator=(const OneDimEstimator&){ return *this; }
};


class Histogram1d{

public:

	virtual ~Histogram1d(){};
	virtual long getBucketNum() const = 0;
	virtual double getBucketMin(long bucket) const = 0;
	
	virtual long getBucket(double d) const = 0;
	virtual double getBucketMax(long bucket) const = 0;
	virtual long& getBucketCount(long bucket) = 0;
	virtual void add(double v) = 0;

protected:
	Histogram1d(){}
    Histogram1d(const Histogram1d&){}
    Histogram1d& operator=(const Histogram1d&){ return *this; }

};

class EquiWidth1d : public Histogram1d, public OneDimEstimator{

private:

	template <typename T>
	inline const string listToString(const vector<T>& v) const{
	
		if (v.size() == 0)
			return "NULL";
	
		if (v.size() > 10){
		
			stringstream s;
	
			for (int i = 0; i < 5; i++){
				
				s << v.at(i);
				s << ", ";
			}
			
			s << "...";
			
			for (int i = 4; i >= 0; i--){
		
				s << ", ";
				s << v.at(v.size()-1-i);
			}
			
			return s.str();
		}
	
		stringstream s;
	
		for (int i = 0; i < v.size(); i++){
		
			if (i > 0)
				s << ", ";
				
			s << v.at(i);
		}
		
		return s.str();
	}

	inline long getSum(long a, long b) const{
		
		if (a > b)
			return 0;
	
		long ret = 0;
		
		if (!cumulative){
		
			for (long j = a ; j <= b; j++)
				ret += counts[j];
				
			return ret;		
		}
		
		if (a == 0)
			return counts[b];
			
		return counts[b]-counts[a-1];
	}
	

public:
	vector<long> counts;
	
	bool cumulative = false;
	
	double min;
	double max;
	
	double normMul;
	
	double bucketDiv;
	
	bool normalize;
	
	int lod = -1;
	
	inline EquiWidth1d(long buckets, double min=0, double max=1) : min(min), max(max), bucketDiv( buckets*1.0/(max-min) ), normMul(1.0/(max-min)), counts(buckets, 0){
	
		if (min != 0 || max != 1)
			normalize = true;	
			
		if (UtilMath::isPowerOfTwo(buckets)){
		
			lod = UtilMath::getPowerOfTwo(buckets);
			
		} else {
			
			lod = -1;
		}
	}
	
	inline EquiWidth1d(const EquiWidth1d& e){
		
		counts.reserve(e.counts.size() );
		
		for (auto it = e.counts.begin(); it != e.counts.end(); it++)
			counts.push_back( (*it) );
		
		cumulative = e.cumulative;
		min = e.min;
		max  = e.max;
		normMul = e.normMul;
		bucketDiv = e.bucketDiv;
		normalize = e.normalize;
		lod = e.lod;
	}
	
	
	~EquiWidth1d() override{}
	
	inline void print(){
	
		for (long j = 0; j < getBucketNum(); j++)
			if (getBucketCount(j) > 0)
				cout << "[" << getBucketMin(j) << "," << getBucketMax(j) << "]:{" << getBucketCount(j) << "}" << endl;
		
		cout << endl;
	}
	
	inline double intervalCount(double qmin, double qmax, QueryMode mode=QueryMode::EST) const override{
		
		const long imin = getBucket(qmin);
		const long imax = getBucket(qmax);
		
		if (imin == imax){
		
			if (mode == QueryMode::LB)
				return 0;
		
			const double c1 = getSum(imin, imin);
			
			if (mode == QueryMode::UB)
				return c1;
		
			return c1*(qmax-qmin)*bucketDiv;
		}
		
		const long c2 = getSum(imin+1, imax-1);
		
		if (mode == QueryMode::LB)
			return c2;
		
		const long c1 = getSum(imin, imin);
		const long c3 = getSum(imax, imax);
		
		if (mode == QueryMode::UB)
			return c1+c2+c3;
		
		const double r1 = UtilMath::absVal<double>( getBucketMax(imin)-qmin )*bucketDiv;
		const double r3 = UtilMath::absVal<double>( qmax-getBucketMin(imax) )*bucketDiv;		
		
		return (r1*c1)+c2+(r3*c3);
	}
	
	
	inline bool makeCumulative(bool b = true){
	
		const bool ret = cumulative;
	
		if (cumulative == b)
			return ret;
		
		if (b){
		
			for (long j = 1; j < counts.size(); j++)
				counts[j] += counts[j-1];
				
		} else {
		
			for (long j = counts.size()-1; j >= 1; j--)
				counts[j] -= counts[j-1];
		}
		
		cumulative = b;
		
		return ret;
	}
	
	
	inline long getBucket(double d) const override{
	
		if (normalize)
			return UtilHist::discretize( (d-min)*normMul, counts.size() );
		
		if (lod > 0)
			UtilHist::discretizePowerOfTwo(d, lod);
		
		return UtilHist::discretize(d, counts.size() );
		
	}
	
	inline void add(double v) override{
	
		assert (!cumulative);
	
		counts[getBucket(v)]++;
	}
	
	inline double getBucketMin(long bucket) const override{
	
		if (normalize)
			return min+( UtilHist::bucketMin(bucket, counts.size()) )*(max-min);
		
		return UtilHist::bucketMin(bucket, counts.size());
	}
	
	inline double getBucketMax(long bucket) const override{
	
		if (normalize)
			return min+( UtilHist::bucketMax(bucket, counts.size()) )*(max-min);
		
		return UtilHist::bucketMax(bucket, counts.size());
	}
	
	inline long getBucketNum() const override{
	
		return counts.size();
	}
	
	inline long getTotalBucketCount() const{
	
		if (cumulative)
			return counts[counts.size()-1];
		
		return UtilMath::sum<long>(counts);
	}
	
	inline long& getBucketCount(long bucket) override{
	
		if (cumulative)
			assert (false);
		
		return counts[bucket];
	}
	
	inline long& operator[](size_t bucket){
	
		return getBucketCount(bucket);
	}
	
	inline long& at(size_t bucket){
	
		return getBucketCount(bucket);
	}
	
	inline const string toString() const{
	
		stringstream s;
		
		s << "equiwidth " << getBucketNum() << " buckets between " << min << " and " << max;
		
		if (cumulative)
			s << " cumulative";
		
		s << " count " << getTotalBucketCount();
		
		s << " counts " << listToString<long>(counts);
		
		return s.str();
	}
	

};



#endif