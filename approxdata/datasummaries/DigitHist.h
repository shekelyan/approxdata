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

#ifndef SHEKELYAN_DIGITHIST_H
#define SHEKELYAN_DIGITHIST_H

#include <approxdata/utils/utils.h>
#include <approxdata/data/data.h>
#include <approxdata/datasummaries/datasummaries.h>

/*

Class for the digit histograms composing DigitHist summaries.

*/







class DigitHistogram : public DataSummary{

typedef unique_ptr<DigitHistogram> DH;
typedef unique_ptr<SparseGridHist> Hist;
typedef unique_ptr<EquiWidth1d> Marg;

protected:

	long mul = 1; // (Digit)-multiplier of all counts
	const int DIMS; // Number of dimensions
	bool AVI = false; // Intra-bucket distribution according to attribute value independence assumption?
	vector<OneDimEstimator*> oneDimEsts; // Intra-bucket distribution estimators (pointers to 1D histograms)
	
	inline DigitHistogram(int dims, bool avi=false) : DIMS(dims), AVI(avi){
		
		mul = 1;
	}
	
	inline DigitHistogram(int dims, int s, Hist& h, bool avi) : DigitHistogram(dims, avi){
		
		gridhist = std::move(h);
		mul = 1L << s;
	}
	
public:
	
	
	vector<Marg> marginals; // 1D histograms for intra-bucket distribution
	Hist gridhist; // multidimensional grid histogram
	
	// Constructor only used internally (by other class in the library)
	inline DigitHistogram(DataSet& data, Params& params) : DigitHistogram(data.getDims() ){
		
		summarize(data, params);
	}
	
	// Constructor receiving a dataset and string with parameters
	
	inline DigitHistogram(DataSet& data, string paramStr="-kb 100") : DigitHistogram(data.getDims() ){
	
		Params params(paramStr);
		summarize(data, params);
	}
	
	// All counts of the grid histogram are multiples of {mul} (they are stored in gridhist divided by mul)
	inline long getMultiplier() const{
	
		return mul;
	}
	
	~DigitHistogram() override {};
	
	inline void checkConsistency(){
		
		assert (AVI);
		assert (marginals.size() > 0);
		
		for (int i = 0; i < DIMS; i++){
		
			marginals[i]->makeCumulative(false);
			
			const int margLod = gridhist->z.dimMaskBits[i];
			const int srcLod = UtilMath::getPowerOfTwo(marginals[i]->counts.size() );
			
			const long buckets = gridhist->getBucketNum();
			const long margBuckets = 1L << margLod;
			const long srcBuckets = 1L << srcLod;
			
			vector<long>& src = marginals[i]->counts;
			
			if (src.size() != srcBuckets)
				cout << "src.size" << src.size() << " != " << srcBuckets << " srcBuckets" << endl;
			
			
			assert (src.size() == srcBuckets);
			
			
			for (long jsrc = 0; jsrc < srcBuckets; jsrc++)
				assert(src[jsrc] >= 0);
			
			if (margLod >= srcLod){
			
				const int lodDiff = margLod-srcLod;
				
				for (long j = 0; j < buckets; j++){
			
					const long jsrc = gridhist->z.getDimCoord(gridhist->elems[j], i) >> lodDiff;
					const long c = gridhist->counts[j]*mul;
					
					assert (jsrc >= 0);
					assert (jsrc < src.size() );
					
					assert (src[jsrc] >= 0);
					
					assert (c >= 0);
					assert (c <= src[jsrc]);
					
					src[jsrc] -= c;
				}
				
			} else {
				
				vector<long> tmp(margBuckets, 0);
				
				for (long j = 0; j < buckets; j++){
			
					const long x = gridhist->z.getDimCoord(gridhist->elems[j], i);
					const long y = gridhist->counts[j] * mul;
					
					assert (x >= 0);
					assert (x < tmp.size() );
					
					assert (y > 0);
					
					tmp[x] += y;
				}
				
				const int lodDiff = srcLod-margLod;
				
				vector<double> v( 1L << lodDiff );
			
				for (long jmarg = 0; jmarg < margBuckets; jmarg++){
					
					long b = tmp[jmarg];
					
					if (b == 0)
						continue;
					
					double vsum = 0;
					
					const long jsrc = jmarg << lodDiff;
					
					for (long jv = 0; jv < v.size(); jv++){
						v[jv] = src[jsrc+jv];
						vsum += v[jv];
					}
					
					if (vsum == 0)
						continue;
					
					//if (b > vsum)
						//b = vsum;
					
					const double sumDiv = b*1.0/vsum;
					
					assert(sumDiv > 0);
					assert(sumDiv <= 1);
					
					for (long jv = 0; jv < v.size(); jv++)
						v[jv] *= sumDiv;
					
					UtilMath::optimalRounding(v);
					
					for (long jv = 0; jv < v.size(); jv++){
					
						const long c = v[jv];
						
						if (c == 0)
							continue;
						
						if (c > src[jsrc+jv]){
							cout << "c " << c << " > " << src[jsrc+jv] << endl;
							cout << endl;
						}
						
						assert(c <= src[jsrc+jv]);
						
						src[jsrc+jv] -= c;
					}
				}
			}
			
			for (long jsrc = 0; jsrc < srcBuckets; jsrc++)
				assert(src[jsrc] >= 0);
			
			assert (UtilMath::sum<long>(src) == 0 );
		}
		
		cout << "check passed!" << endl;
	}
	
	inline void take(vector<Marg>& initMarginals){
		
		assert (AVI);
		assert (initMarginals.size() > 0);
		
		for (int i = 0; i < DIMS; i++){
		
			initMarginals[i]->makeCumulative(false);
			
			const int margLod = gridhist->z.dimMaskBits[i];
			
			assert (UtilMath::isPowerOfTwo( initMarginals[i]->counts.size() ) );
			
			const int srcLod = UtilMath::getPowerOfTwo(initMarginals[i]->counts.size() );
			const int dstLod = srcLod;
			
			const long buckets = gridhist->getBucketNum();
			const long margBuckets = 1L << margLod;
			const long srcBuckets = 1L << srcLod;
			const long dstBuckets = 1L << dstLod;
			
			Marg m(new EquiWidth1d(dstBuckets) );
			
			vector<long>& src = initMarginals[i]->counts;
			vector<long>& dst = m->counts;
			
			if (src.size() != srcBuckets)
				cout << "src.size" << src.size() << " != " << srcBuckets << " srcBuckets" << endl;
			
			
			assert (src.size() == srcBuckets);
			assert (dst.size() == dstBuckets);
			
			for (long jsrc = 0; jsrc < srcBuckets; jsrc++)
				assert(src[jsrc] >= 0);
			
			if (margLod >= srcLod){
			
				const int lodDiff = margLod-srcLod;
				
				for (long j = 0; j < buckets; j++){
			
					const long jsrc = gridhist->z.getDimCoord(gridhist->elems[j], i) >> lodDiff;
					const long c = gridhist->counts[j]*mul;
					
					assert (jsrc >= 0);
					assert (jsrc < dst.size() );
					assert (jsrc < src.size() );
					
					assert (src[jsrc] >= 0);
					
					assert (c >= 0);
					assert (c <= src[jsrc]);
					
					dst[jsrc] += c;
					src[jsrc] -= c;
				}
				
			} else {
				
				vector<long> tmp(margBuckets, 0);
				
				for (long j = 0; j < buckets; j++){
			
					const long x = gridhist->z.getDimCoord(gridhist->elems[j], i);
					const long y = gridhist->counts[j] * mul;
					
					assert (x >= 0);
					assert (x < tmp.size() );
					
					assert (y > 0);
					
					tmp[x] += y;
				}
				
				const int lodDiff = srcLod-margLod;
				
				vector<double> v( 1L << lodDiff );
			
				for (long jmarg = 0; jmarg < margBuckets; jmarg++){
					
					long b = tmp[jmarg];
					
					if (b == 0)
						continue;
					
					double vsum = 0;
					
					const long jsrc = jmarg << lodDiff;
					
					for (long jv = 0; jv < v.size(); jv++){
						v[jv] = src[jsrc+jv];
						vsum += v[jv];
					}
					
					if (vsum == 0)
						continue;
					
					if (b > vsum)
						b = vsum;
					
					const double sumDiv = b*1.0/vsum;
					
					assert(sumDiv > 0);
					assert(sumDiv <= 1);
					
					for (long jv = 0; jv < v.size(); jv++)
						v[jv] *= sumDiv;
					
					UtilMath::optimalRounding(v);
					
					for (long jv = 0; jv < v.size(); jv++){
					
						const long c = v[jv];
						
						if (c == 0)
							continue;
						
						if (c > src[jsrc+jv]){
							cout << "c " << c << " > " << src[jsrc+jv] << endl;
							cout << endl;
						}
						
						assert(c <= src[jsrc+jv]);
						
						src[jsrc+jv] -= c;
						dst[jsrc+jv] += c;
					}
				}
			}
			
			for (long jsrc = 0; jsrc < srcBuckets; jsrc++)
				assert(src[jsrc] >= 0);
			
			m->makeCumulative(true);
			
			marginals.push_back(std::move(m) );
			oneDimEsts.push_back( marginals[i].get() );	
		}
		
		assert (marginals.size() == DIMS);
	}
	
	inline void add(Point& p){
	
		gridhist->addPoint(p);
			
		if (AVI){
			for (int i = 0; i < DIMS; i++)
				marginals[i]->add(p[i] );
		}
	}
	
	bool operator<(const DigitHistogram& other) const{
	
		if (mul == other.mul)
			return gridhist->total < other.gridhist->total;
		
		return mul < other.mul;
    }
	
	inline void reduceLod(){
	
		gridhist->reduceLod();
	}
	
	inline long getLargestCount() const{
	
		long ret = 0;
		
		for (long j = 0; j < gridhist->counts.size(); j++)
			if (gridhist->counts[j] > ret)
				ret = gridhist->counts[j];
		
		return ret;
	}
	
	
	inline int getLod() const{
	
		return gridhist->z.getLod();
	}
	
	inline double getUError(int steps=30) const{
	
		return gridhist->getUError()*mul;
	}
	
	inline const string toString() const{
	
		stringstream s;
		
		s << gridhist->toString();
		
		s << " " << mul;
		
		return s.str();
	}
	
	
	
	
	inline DH createDigitHistogram(int s, int w, long maxBuckets, int lod=62){
		
		Hist retgrid( new SparseGridHist(gridhist->DIMS, maxBuckets, UtilMath::minVal<int>(lod, gridhist->getLod())) );
		
		for (long j = 0; j < gridhist->getBucketNum(); j++){
			
			const long c = (gridhist->counts[j] >> s) & UtilBits::mask(w);
			
			if (c == 0)
				continue;
			
			const int lodDiff = (gridhist->getLod()-retgrid->getLod() );
			
			assert (lodDiff >= 0);
			
			const long v = gridhist->elems[j] >> lodDiff;
			retgrid->addSortedZ(v, c);
		}
		
		retgrid->finalize();
		
		DH ret(new DigitHistogram(DIMS, s, retgrid, AVI) );
		
		return ret;
	}
	
	inline long count(){
	
		return UtilMath::sum<long>(gridhist->counts);
	}
	
	inline void summarize(DataSet& data, Params& params){
		
		const int initLod = params.get("lod") ? params.get("lod") : 62;
		
		const double maxSparse = params.get("kb")/( 8.0/(1024) );
		const double maxMarg = params.get("margkb")/( 8.0/(1024) );
		
		AVI = !params.get("nomarg") && maxMarg > 0;
		mul = 1;
		
		//cout << "create summary DigitHistogram " << params << endl;
		
		long marginalBuckets = UtilMath::powerOfTwoUB( maxMarg/DIMS );
		long gridBuckets =  maxSparse;
		
		/*
		if (params.get("check")){
		
			initLod = params.get("check")*DIMS;
			
			AVI = true;
			
			marginalBuckets = 1L << (10+ (int) params.get("check"));
			gridBuckets = 10000000;
		}*/
		
		assert (gridBuckets > 0);
		
		if (AVI)
			assert (marginalBuckets > 0);
		
		gridhist.reset(new SparseGridHist(DIMS, gridBuckets, initLod) );
		
		if (AVI){
			
			for (int i = 0; i < DIMS; i++){
			
				Marg m(new EquiWidth1d(marginalBuckets));
			
				marginals.push_back( std::move(m) );
				oneDimEsts.push_back( marginals[i].get() );
			}
		}
		
		for (auto it = data.begin(); it != data.end(); it++)
			add( (*it) );
		
		gridhist->finalize();
		
		cout << gridhist->toString() << endl;
		
		if (AVI)
		for (int i = 0; i < DIMS; i++){
		
			marginals[i]->makeCumulative(true);
			//cout << "i " << i << " " << marginals[i]->toString() << endl;
		}
	}
	
	inline unique_ptr<RangeQueries> getBuckets() const{
	
		unique_ptr<RangeQueries> ret(new RangeQueries() );
		
		for (auto it = gridhist->elems.begin(); it != gridhist->elems.end(); it++){
		
			Box b = gridhist->getBox( (*it) );
			ret->add( b );
		}
		
		return ret;
	}
	
	inline unique_ptr<RangeQueries> getBuckets(const Box& box) const{
	
		unique_ptr<RangeQueries> ret(new RangeQueries() );
		
		const long zmin = gridhist->getZ(box.min);
		const long zmax = gridhist->getZ(box.max);
		
		long k = 0;
		for (auto it = gridhist->elems.begin(); it != gridhist->elems.end(); it++, k++){
		
			if ( !gridhist->z.contains(zmin, zmax, (*it)) )
				continue;
			
			Box b = gridhist->getBox( (*it) );
			
			b.count = gridhist->counts[k];
			ret->add(b);
		}
		
		return ret;
	}
	
	inline double marginalBoxCountUpperBound(int i, double qmin, double qmax) const {
	
		if (AVI)
			return marginals[i]->intervalCount(qmin, qmax, QueryMode::UB);
		else
			return mul*gridhist->getTotalNum();
	}
	
	inline double boxCount(const Point& qmin, const Point& qmax, QueryMode mode = QueryMode::EST) const override{
		
		if (AVI){
		
			assert (oneDimEsts.size() == DIMS);
			return mul*gridhist->boxCount(qmin, qmax, mode, oneDimEsts);
			
		} else{
			
			return mul*gridhist->boxCount(qmin, qmax, mode);
		}
	}
	
	#if defined(SHEKELYAN_SELIMAGE)  
	
	inline void visualize(DataImage& img) const{
	
		array<float,3> black = {0,0,0};
		
		array<float,3> green = {0,100,0};
		array<float,3> violet = {255,0,255};
		array<float,3> red = {255,0,0};
		array<float,3> blue = {0,0,255};
		array<float,3> gray = {100,100,100};
		
		double sum = 1.0;
		
		shared_ptr<RangeQueries> buckets = getBuckets(img.getDataBox());
			
		for (auto it = buckets->ranges.begin(); it != buckets->ranges.end(); it++)
			sum += it->count;
		
		Box b;
		
		for (auto it = buckets->ranges.begin(); it != buckets->ranges.end(); it++){
			
			b = (*it);
			b.count /= sum;
				
			img.fillRect(b, black);	
		}
	}
	
	#endif
	
};





class DigitHistBinary{

public:
	const byte cont = 0x80; // == 128 == 0b10000000
	const byte sel = 0x7F;  // == 127 == 0b01111111
	
	const long FLAG_AVI = 1 << 0;
	const long FLAG_TWO_LEVEL = 1 << 1;
	
	byte* current = NULL;
	
	byte* mbytes = NULL;
	byte* mlod = NULL;
	
	long flags = 0;
	int DIMS = -1;
	int K = -1;
	
	vector<byte*> marginals;
	
	template<bool countOnly=false>
	inline long writeByte(const byte& b){
	
		if (countOnly)
			return 1;
			
		(*current) = b;
		current++;
		
		return 1;
	}
	
	inline long getZ(const ZCurveCoords& z, const Point& p) const{
	
		long ret = 0;
		
		for (int i = 0; i < DIMS; i++)
			ret = z.setDimCoord(ret, UtilHist::discretizePowerOfTwo(p[i], z.getLod(i) ), i);
		
		return ret;
	}
	
	
	template<bool countOnly=false>
	inline long writeInteger(long v){
	
	
		if (v == 0){
			
			if (countOnly)
				return 1;
			else
				return writeByte<countOnly>(0);
		}
	
		long ret = 0;
		
		while (v > 0){
		
			if (countOnly)
				ret += 1;
			else
				ret += writeByte<countOnly>( (v > sel) ? ( (v & sel)|cont ) : (v & sel) );
			v >>= 7;
		}
		
		return ret;
	}
	
	
	
	template<bool countOnly=false>
	inline long writeDouble(double d){
	
		if (countOnly)
			return 8;
	
		long ret = 0;
		vector<byte> arr(8);
		
		double* ds = reinterpret_cast<double*>(arr.data() );
		
		ds[0] = d;
		
		for (auto it = arr.begin(); it != arr.end(); it++)
			ret += writeByte<countOnly>( (*it) );
		
		return ret;
	}
	
	
	
	inline int getBytesPerValue(const EquiWidth1d& h){
	
		long max = 0;
			
		for (auto it = h.counts.begin(); it != h.counts.end(); it++){
		
			if ( (*it) > max)
				max = (*it);
		}
		
		int bytesPerValue = 0;
		
		while (max != 0){
			max >>= 8;
			bytesPerValue++;
		}
	
		return bytesPerValue;
	}
	
	
	
	template<bool countOnly=false>
	inline long writeMarginal(const EquiWidth1d& h){
	
		long ret = 0;
	
		const int bytesPerValue = getBytesPerValue(h);
			
		ret += writeInteger<countOnly>(h.counts.size()*bytesPerValue); // 8
			
		for (auto it = h.counts.begin(); it != h.counts.end(); it++){
			
			long val = (*it);
			
			for (int k = 0; k < bytesPerValue; k++){
				
				ret += writeByte<countOnly>(val & 255); // 9
				val >>= 8;
			}
		}
		
		return ret;
	}
	
	template<bool countOnly=false>
	inline long writeRootCell(const long& e2, const int& lodRed, long& k, const SparseGridHist& h1){
	
		long ret = 0;
	
		long last1 = (e2 << lodRed)-1;
			
		while (k < h1.elems.size() ){
	
			const long e1 = h1.elems[k];
		
			if ( (e1 >> lodRed) > e2)
				break;
			
			const long c1 = h1.counts[k];
		
			ret += writeInteger<countOnly>(e1-last1);
			ret += writeInteger<countOnly>(c1);
			
			k++;
			
			last1 = e1;
		}
		
		return ret;
	}
	
	template<bool countOnly=false>
	inline long writeGridHistogram(const SparseGridHist& h2, const SparseGridHist& h1, const int& lodRed){
		
		long ret = 0;
		
		long k1 = 0;
		long last2 = -1;
		
		for (long j = 0; j < h2.elems.size(); j++){
		
			const long e2 = h2.elems[j];
			const long c2 = h2.counts[j];
		
			ret += writeInteger<countOnly>( e2-last2 );
			ret += writeInteger<countOnly>( c2 );
			
			{
				long k1_ = k1;
				const long b2 = writeRootCell<true>( e2, lodRed, k1_, h1);
				ret += writeInteger<countOnly>(b2);
			}
			
			ret += writeRootCell<countOnly>(e2, lodRed, k1, h1);
			
			last2 = e2;
		}
		
		return ret;
	}
	
	template<bool countOnly=false>
	inline long writeDigitHistogram(const SparseGridHist& h1, const long mul){
		
		long ret = 0;
		
		//ret += writeInteger<countOnly>(123);
		
		ret += writeInteger<countOnly>(mul); // 101
		
		SparseGridHist h2 = h1;
		
		const ZCurveCoords& z1 = h1.z;
		
		const int lodRed = z1.getLod()/2;
		
		h2.reduceLod(lodRed);
		const ZCurveCoords& z2 = h2.z;
		
		const int lod1 = z1.getLod();
		const int lod2 = z2.getLod();
		
		ret += writeInteger<countOnly>(lod1); // 102
		ret += writeInteger<countOnly>(lod2); // 103
		
		//cout << "write lod1 " << lod1 << " lod2 " << lod2 << " lodRed " << lodRed << endl;
		
		{
			const long bytes = writeGridHistogram<true>( h2, h1, lodRed);
			ret += writeInteger<countOnly>(bytes); // 104
		}
		
		ret += writeGridHistogram<countOnly>( h2, h1, lodRed);
		
		return ret;
	}
	
	
	template<bool countOnly=false>
	inline long write2(const Box& boundingBox, const vector<shared_ptr<DigitHistogram> >& digitHistograms){
	
		DIMS = boundingBox.min.size();
		K = digitHistograms.size();
	
		long ret = 0;
		
		ret += writeInteger<countOnly>(DIMS);  // 1 DIMS
		
		for (int i = 0; i < DIMS; i++)
			ret += writeDouble<countOnly>(boundingBox.min[i]); // 2 MIN
		
		for (int i = 0; i < DIMS; i++)
			ret += writeDouble<countOnly>(boundingBox.max[i]); // 3 MAX
		
		ret += writeInteger<countOnly>(FLAG_AVI | FLAG_TWO_LEVEL); // 4 FLAG
		
		ret += writeInteger<countOnly>(K);  // 5 K
		
		// WRITE MARGINALS
		
		ret += writeInteger<countOnly>(123);
		
		for (int k = 0; k < K; k++){
		
			assert (digitHistograms.size() == K);
		
			for (int i = 0; i < DIMS; i++){
			
				assert (digitHistograms[k]->marginals.size() == DIMS);
				assert (digitHistograms[k]->marginals[i]);
				const EquiWidth1d& h = (*(digitHistograms[k]->marginals[i]) );
				
				ret += writeByte<countOnly>( h.lod); // 6
			}
		}
		
		for (int k = 0; k < K; k++){
			for (int i = 0; i < DIMS; i++){
			
				const EquiWidth1d& h = (*(digitHistograms[k]->marginals[i]) );
				ret += writeByte<countOnly>( getBytesPerValue(h)); // 7
			}
		}
		
		for (int k = 0; k < K; k++){
			for (int i = 0; i < DIMS; i++){
			
				const EquiWidth1d& h = (*(digitHistograms[k]->marginals[i]) );
				ret += writeMarginal<countOnly>(h); //8
			}
		}
		
		ret += writeInteger<countOnly>(123);
		
		// WRITE GRID HISTOGRAMS
		
		for (int k = 0; k < K; k++)
			ret += writeDigitHistogram<countOnly>( (*(digitHistograms[k]->gridhist)), digitHistograms[k]->getMultiplier() );
		
		return ret;
		
	}
	
	template<bool countOnly=false>
	inline long write(byte* data, const Box& boundingBox, const vector<shared_ptr<DigitHistogram> >& digitHistograms){
	
		if (!countOnly)
			current = data;
		
		long ret = 0;
	
		const long bytes = write2<true>(boundingBox, digitHistograms); 
		
		ret += writeInteger<countOnly>(bytes); // 0 BYTES
		ret += write2<countOnly>(boundingBox, digitHistograms);
		
		return ret;
		
	}
	
	
	inline int readByte(){
	
		const int b = (*current);
		current++;
			
		return b;
	}
	
	inline long readInteger(){
	
		long ret = 0;
	
		for (int j = 0; j < 64; j+=7){
		
			const int b = readByte();
			
			ret |= (b & sel) << j;
			
			if (! (b & cont) )
				return ret;
		}
		
		assert (false);
	}
	
	
	inline double readDouble(){
		
		vector<byte> arr(8);
		
		for (int i = 0; i < 8; i++)
			arr[i] = readByte();
		
		double* ds = reinterpret_cast<double*>(arr.data() );
		
		return ds[0];
	}
	
	
	inline long intervalSum(const int& k, const int& i, const long& a, const long& b) const{
		
		if (a > b)
			return 0;
		
		const int bytesPerVal = mbytes[k*DIMS+i];
		const byte* m = marginals[k*DIMS+i];
		
		const int pos1 = (a-1)*bytesPerVal;
		const int pos2 = b*bytesPerVal;
		
		long ret = 0;
		
		for (int x = 0; x < bytesPerVal; x++){
			
			ret += m[pos2+x] << (8*x);
			
			if (a > 0){
				ret -= m[pos1+x] << (8*x);
			}
		}
		
		return ret;
	}
	
	inline double intervalCount(int k, int i, double qmin, double qmax, QueryMode mode){
	
		const int lod = mlod[k*DIMS+i];
		
		const long imin = UtilHist::discretizePowerOfTwo(qmin, lod);
		const long imax = UtilHist::discretizePowerOfTwo(qmax, lod);
		
		const double bucketDiv = (1L << lod);
		
		if (imin == imax){
		
			if (mode == QueryMode::LB)
				return 0;
			
			const double c1 = intervalSum(k, i, imin, imin);
			
			if (mode == QueryMode::UB)
				return c1;
			
			return c1*(qmax-qmin)*bucketDiv;
		}
		
		const double c2 = intervalSum(k, i, imin+1, imax-1);
		
		
		if (mode == QueryMode::LB)
			return c2;
		
		const double c1 = intervalSum(k, i, imin, imin);
		const long c3 = intervalSum(k, i, imax, imax);
		
		if (mode == QueryMode::UB)
			return c1+c2+c3;
		
		const double r1 = UtilMath::absVal<double>( UtilHist::bucketMax(imin, 1L << (mlod[k]) )-qmin )*bucketDiv;
		const double r3 = UtilMath::absVal<double>( qmax-UtilHist::bucketMin(imax, 1L << (mlod[k]) ) )*bucketDiv;		
		
		return (r1*c1)+c2+(r3*c3);
	}
	
	inline double digitHistogramBoxCount(int k, const Point& qmin, const Point& qmax, QueryMode mode){
		
		//assert (readInteger() == 123);
		
		const long mul = readInteger();
		
		const int lod1 = readInteger();
		const int lod2 = readInteger();
		
		const int lodRed = lod1-lod2;
		
		//cout << "lod1 " << lod1 << " lod2 " << lod2 << " lodRed " << lodRed << endl;
		
		//const long vals = readInteger();
		
		const ZCurveCoords z1(DIMS, lod1);
		const ZCurveCoords z2(DIMS, lod2);
		
		const long zmin1 = getZ(z1, qmin);
		const long zmax1 = getZ(z1, qmax);
		
		const long zmin2 = getZ(z2, qmin);
		const long zmax2 = getZ(z2, qmax);
		
		assert (zmin1 != -1);
		assert (zmax1 != -1);
		
		assert (zmin2 != -1);
		assert (zmax2 != -1);
		
		const long zmin1c = z1.incAllCoords(zmin1);
		const long zmax1c = z1.decAllCoords(zmax1);
		
		const long zmin2c = z2.incAllCoords(zmin2);
		const long zmax2c = z2.decAllCoords(zmax2);
		
		Point tempmin(DIMS);
		Point tempmax(DIMS);
		
		for (int i = 0; i < DIMS; i++){
				
			const double bmin = UtilHist::bucketMin( z1.getDimCoord(zmin1, i), z1.getCellsPerDim(i) );
			const double bmax = UtilHist::bucketMax( z1.getDimCoord(zmin1, i), z1.getCellsPerDim(i) );
			
			const double imin = UtilMath::maxVal<double>( qmin[i], bmin);
			const double imax = UtilMath::minVal<double>( qmax[i], bmax);
			
			if (imin > imax){
			
				tempmin[i] = 0;
				
			} else {
			
				const double a = intervalCount(k, i, imin, imax, mode);
				const double b = intervalCount(k, i, bmin, bmax, mode);
				
				tempmin[i] =  b > 0 ? UtilMath::makeBetween<double>(0.0,1.0, (a/b)) : 0;
			}
		}
		
		for (int i = 0; i < DIMS; i++){
				
			const double bmin = UtilHist::bucketMin( z1.getDimCoord(zmax1, i), z1.getCellsPerDim(i) );
			const double bmax = UtilHist::bucketMax( z1.getDimCoord(zmax1, i), z1.getCellsPerDim(i) );
					
			const double imin = UtilMath::maxVal<double>( qmin[i], bmin);
			const double imax = UtilMath::minVal<double>( qmax[i], bmax);
			
			if (imin > imax){
			
				tempmax[i] = 0;
				
			} else {
			
				const double a = intervalCount(k, i, imin, imax, mode);
				const double b = intervalCount(k, i, bmin, bmax, mode);
				
				tempmax[i] = b > 0 ? UtilMath::makeBetween<double>(0.0,1.0, (a/b)) : 0;
			}
		}
		
		const long bytes3 = readInteger(); // 104
		
		byte* next3 = current+bytes3;
		
		long last2 = -1;
		
		double ret = 0;
		
		while (current < next3){

			
			const long e2 = last2+readInteger();
			last2 = e2;
			const long c2 = readInteger();
			
			const long bytes2 = readInteger();
			
			byte* next2 = current+bytes2;
			
			
			// query box contains root cell "e2"?
			if (z2.contains(zmin2c, zmax2c, e2)){
				
				ret += c2; // add count "c2" of root cell "e2"
				current = next2; // skip to next root cell
				continue;
			}
			
			// query box does not intersect root cell "e2"?
			if (!z2.contains(zmin2, zmax2, e2)){
				
				current = next2; // skip to next root cell
				continue; 
			}
			
			long last1 = (e2 << lodRed)-1;
			
			while (current < next2){
			
				const long e1 = last1+readInteger();
				last1 = e1;
				const long c1 = readInteger();
				
				
				// query box contains cell "e1"?
				if (z1.contains(zmin1c, zmax1c, e1)){
					ret += c1;
					continue;
				}
			
				if (mode == QueryMode::LB)
					continue;
				
				// query box does not intersect cell "e1"?
				if (!z1.contains(zmin1, zmax1, e1))
					continue; 
			
				if (mode == QueryMode::UB){
					ret += c1;
					continue;	
				}
			
				if (mode == QueryMode::EST){
			
					double frac = 1.0;
			
					for (int i = 0; i < DIMS; i++){
				
						if (z1.sameDimCoord(zmin1, e1, i) )
							frac *= tempmin[i];
						else if (z1.sameDimCoord(zmax1, e1, i))
							frac *= tempmax[i];
					}
				
					assert (frac >= 0);
					assert (frac <= 1);
				
					if (frac > 0)
						ret += frac*c1;
				}
			}
			
			assert (current == next2);
		}
		
		assert (current == next3);
		
		return mul*ret;
	}
	
	inline long getBytes( const SparseGridHist& h1, const long mul ){
	
		return writeDigitHistogram<true>(h1, mul);
	}
	
	inline double boxCount(byte* data, const Point& qmin, const Point& qmax, QueryMode mode = QueryMode::EST){
		
		current = data;
		
		const long bytes = readInteger(); // 0 BYTES
		
		byte* end = current+bytes;
		
		DIMS = readInteger(); // 1 DIMS
		
		current += DIMS*8*2; // 2 MIN, 3 MAX
		
		flags = readInteger(); // 4 FLAGS
		
		K = readInteger(); // 5 K
		
		{
		
			assert (readInteger() == 123);
		
			mlod = current; // 6
		
			current += DIMS*K; 
		
			mbytes = current; // 7
		
			current += DIMS*K;
		
			
		
			marginals.clear();
		
			marginals.reserve(K*DIMS);
			
			for (int k = 0; k < K; k++){
			
				for (int i = 0; i < DIMS; i++){
				
					const long bytes = readInteger(); // 8
				
					marginals.push_back(current);
					current += bytes;
				}
			}
			
			assert (readInteger() == 123);
		}
		
		
		
		double ret = 0;
		
		for (int k = 0; k < K; k++)
			ret += digitHistogramBoxCount(k, qmin, qmax, mode);
		
		assert (current == end);
		
		for (int i = 0; i < DIMS; i++){
		
			double ub = 0;
			
			for (int k = 0; k < K; k++)
				ub += intervalCount(k, i, qmin[i], qmax[i], QueryMode::UB);
			
			if (ub < ret)
				ret = ub;
		}
		
		marginals.clear();
		
		return ret;
	}

};



class DigitHistReader: public DataSummary{

	private:
		
		
		vector<byte> bytes;

	public:
	
		inline DigitHistReader(const string file){
		
			ifstream in(file, std::ios::binary | std::ios::in);
			
			in.seekg(0, std::ios::end);
			
			size_t len = in.tellg();
			
			bytes.reserve(len);
			
			for (long j = 0; j < len; j++)
				bytes.push_back(0);
			
			char* temp = reinterpret_cast<char*>(bytes.data() );
			
			
			in.seekg(0, std::ios::beg); 
			in.read(temp, len);
			
			in.close();
		}
		
		
		
		inline double boxCount(const Point& qmin, const Point& qmax, QueryMode mode=QueryMode::EST) const override {
		
			DigitHistBinary b;
			
			vector<byte> copy = bytes;
			
			return b.boxCount(copy.data(), qmin, qmax, mode);
		}
};


// DigitHist data summary as described in http://www.vldb.org/pvldb/vol10/p1514-shekelyan.pdf (except two-level representation)

// brief overview of terminology:
// lod = level of detail (binary logarithm of number of cells/buckets)

class DigitHist : public DataSummary{

	typedef shared_ptr<DigitHistogram> DH;
	typedef shared_ptr<DigitHistogram> Hist;
	
	json js;

public:
	
	~DigitHist() override {};
	
	vector<DH> digitHistograms;
	
	int DIMS;
	
	Params params;
	
	const string dataFile;
	
	Box dataBox;
	
	inline DigitHist(DataSet& data, string paramStr) : params(paramStr), DIMS(data.getDims() ), dataFile(data.getFile() ){
		
		summarize(data, params);
	}
	
	inline DigitHist(shared_ptr<DataSet> data, string paramStr) : params(paramStr), DIMS(data->getDims() ), dataFile(data->getFile() ){
		
		summarize( (*data), params);
	}
	
	inline DigitHist(unique_ptr<DataSet> data, string paramStr) : params(paramStr), DIMS(data->getDims() ), dataFile(data->getFile() ){
		
		summarize( (*data), params);
	}
	
	
	inline long getBytes(const DigitHistogram& dh) const{
	
		DigitHistBinary b;
		
		return b.writeDigitHistogram<true>( *(dh.gridhist), dh.getMultiplier());
	}
	
	inline void summarize(DataSet& data, Params& params){
		
		
		Point min(data.getDims(), 0);
		Point max(data.getDims(), 1);
		
		dataBox = Box(min, max, data.getSize() );//data.getDescription().getBoundingBox();
		
		js["data"]["file"] = data.getFile();
		js["data"]["dimensionality"] = data.getDims();
		//js["data"]["dimensions"] = data.getDimensionsString();
		js["data"]["size"] = data.getSize();
		
		js["summary"]["technique"] = "digithist";
		js["summary"]["parameters"] = params.toString();
		
		
		Params paramsInit;
		
		if (params.get("tkb") == 0)
			params.set("tkb", 100000);
		
		paramsInit.set("kb", params.get("tkb") );
		
		const bool AVI = !params.get("nomarg");
		
		if (!AVI)
			paramsInit.set("nomarg", 1);
		
		const long bytes = params.get("kb")*1024;
		
		const long margBytes = AVI ? bytes*0.25 : 0;
		
		paramsInit.set("margkb", margBytes/1024.0);
		
		cout << "create summary DigitHist " << params << endl;
		
		const int DIMS = data.getDims();
		
		
		// 5.1.1 Initial Histograms
		
		// 1. Create an initial d-dimensional histogram H and, for each
		// dimension i, a one-dimensional marginal histogram.
		
		DH init( new DigitHistogram(data, paramsInit) );
		// both H and M1, ... Md are stored in "init"
		
		cout << "init lod " << init->getLod() << endl;
		
		const int K = params.get("k") > 0 ? params.get("k") : 4;
		
		const long digitBytes = bytes-margBytes;
		
		const double minBytesPerElem = 2;
		const double maxBytesPerElem = 16;
		
		const long minElems = floor(digitBytes/maxBytesPerElem);
		const long maxElems = ceil(digitBytes/minBytesPerElem);
		
		cout << "maxElems " << maxElems << endl;
		
		// 5.1.2 Digit Histogram Compression
		
		// 2. Decompose and compress initial histogram into K digit histograms of maximum S bytes
		
		vector<int> lods = {init->getLod() };
		
		{
		
			if ( getBytes( (*init) ) <= digitBytes){
		
				digitHistograms.push_back( std::move(init) );
				return;
			}
		
			double maxProfit = std::numeric_limits<double>::lowest();
			DH bestInit = std::move( init->createDigitHistogram(0, 62, maxElems, init->getLod() ) );
			
			Timer t("digithist");
		
			t.start(10, 1, init->getLod()*init->getLod()*K );
		
			// Algorithm 1: DigitHistCompress, Line 2
		
			while (getBytes( (*init) ) > digitBytes ){
				
				const long ub = UtilMath::powerOfTwoUB( init->getLargestCount() );
				const int countLod = UtilMath::getPowerOfTwo(ub);
				const int step = ceil(countLod*1.0/K);
				
				// brute-force solver for multiple-choice knapsack problems
				MultipleChoiceKnapsack knap(K, init->getLod()+1);
			
				// Algorithm 1, Line 3: decompose H into (D[K1], . . . , D[0])
			
				for (int k = 0; k < K; k++){
					
					DH dh = std::move( init->createDigitHistogram(k*step, step, maxElems) );
					
					while (getBytes( (*dh) ) > digitBytes && dh->gridhist->getLod() > 1)
						dh->gridhist->reduceLod();
				
					while (dh->gridhist->getLod() > 1){
				
						t.tick();
				
						const int lod = dh->gridhist->getLod();
				
						const double profit = -dh->getUError();
						const double weight = getBytes( (*dh) );//dh->gridhist->getBytes();
				
						knap.add(k, lod, profit, weight );
						dh->gridhist->reduceLod();
				
					}
				}
			
				// Algorithm 1, Line 4 and 5
				
				vector<int> candLods = knap.maximizeProfit(digitBytes);
			
				const double profit = knap.getProfit(candLods);
				const double weight = knap.getWeight(candLods);
				
				assert (weight <= digitBytes);
			
				// Algorithm 1, Line 6
			
				if (profit >= maxProfit){
			
					// Algorithm 1, Line 7
					
					lods = candLods;
					maxProfit = profit;
					bestInit = std::move( init->createDigitHistogram(0, 62, maxElems, init->getLod() ) );	
				}
				
				// Algorithm 1, Line 8
				
				init->reduceLod();
			}
		
			t.end();
		
			
			//long total = 0;
		
			for (int k = 0; k < lods.size(); k++){
		
				const long ub = UtilMath::powerOfTwoUB( bestInit->getLargestCount() );
				const int countLod = UtilMath::getPowerOfTwo(ub);
				const int step = ceil(countLod*1.0/K);
			
				DH dh = std::move( bestInit->createDigitHistogram(k*step,step, maxElems, lods.at(k)) );
			
				//for (long j = 0; j < dh->gridhist->counts.size(); j++)
				//	total += dh->gridhist->counts[j] * dh->getMultiplier();
				
				digitHistograms.push_back( std::move(dh) );
			}
		}
		
		// 5.1.3 Marginal Histograms
		
		// 3. Distribute each initial marginal histogram over the marginal histograms
		
		if (AVI){
		
			Sorter s;
		
			for (int k = 0; k < lods.size(); k++)
				s.add(k, -lods.at(k) );
			
			for (int kk = 0; kk < lods.size(); kk++){
				digitHistograms[s.getSortedIndex(kk)]->take( init->marginals );
			}
		}
		
		
		js["digithist"]["k"] = digitHistograms.size();
		
		
		//cout << "total " << total << " data size " << data.getSize() << endl;
		
		/*
		assert (total <= data.getSize() );
		assert (total >= data.getSize() );
		
		if (true){
			
			const long siz = bestInit->gridhist->counts.size();
			
			map<long, long> c;
			
			for (int k = 0; k < K; k++){
			
				const long ub = UtilMath::powerOfTwoUB( bestInit->getLargestCount() );
				const int countLod = UtilMath::getPowerOfTwo(ub);
				const int step = ceil(countLod*1.0/K);
				
				DH dh = std::move( bestInit->createDigitHistogram(k*step, step, siz+1 ) );
				
				assert (dh->getLod() == bestInit->getLod() );
				
				for (long j = 0; j < dh->gridhist->counts.size(); j++){
			
					assert (dh->gridhist->counts[j] > 0);
 					
					c[dh->gridhist->elems[j]] += dh->gridhist->counts[j] * dh->getMultiplier();
				}
			}
			
			assert (c.size() <= bestInit->gridhist->counts.size() );
			assert (c.size() >= bestInit->gridhist->counts.size() );
			
			for (long j = 0; j < bestInit->gridhist->counts.size(); j++){
				
				if ( c[bestInit->gridhist->elems[j]] != bestInit->gridhist->counts[j])
					cout << c[bestInit->gridhist->elems[j]] << " != " << bestInit->gridhist->counts[j] << endl;
				
				assert (c[bestInit->gridhist->elems[j]] == bestInit->gridhist->counts[j]);
			}
		}
		
		
		
		if (AVI)
		for (int i = 0; i < DIMS; i++)
			cout << "i " << i << " strays " << UtilMath::sum(init->marginals[i]->counts) << endl;
		
		//std::sort(digitHistograms.begin(), digitHistograms.end() );
		
		cout << "finished" << endl;*/
	}
	
	// approximated count in box is equal to sum of component counts
	
	inline double boxCount(const Point& qmin, const Point& qmax, QueryMode mode=QueryMode::EST) const override {
		
		double ret = 0;
		
		for (int k = 0; k < digitHistograms.size(); k++)
			ret += digitHistograms[k]->boxCount(qmin, qmax, mode);
		
		for (int i = 0; i < DIMS; i++){
		
			double r = 0;
		
			for (int k = 0; k < digitHistograms.size(); k++)
				r += digitHistograms[k]->marginalBoxCountUpperBound(i, qmin[i], qmax[i]);
		
			if (r < ret)
				return r;
		}
		
		if (mode == QueryMode::EST){
		
			assert (ret <= boxCount(qmin, qmax, QueryMode::UB) );
			assert (ret >= boxCount(qmin, qmax, QueryMode::LB) );
		}
		
		return ret;
	}
	
	inline DigitHist(const string& jsonFile, const string& binaryFile){
		
		/*
		int K = 0;
		{
			ifstream in(jsonFile);
		
			//assert (in.open() );
		
			in >> js;
		
			DIMS = js["data"]["dimensionality"];
			
			//AVI = !js["digithist"]["nomarg"];
			
			K = js["digithist"]["k"];
			
			in.close();
		}
		
		ifstream in(binaryFile, std::ios::in | std::ios::binary);
		
		assert (in.open() );
		
		for (int k = 0; k < K; k++){
			
			unique_ptr<DigitHistogram> dh( new DigitHistogram(in) );		
			
			digitHistograms.push_back( std::move(dh) );
		}
		
		in.close();*/
	}
	
	inline void writeJSON(const string& jsonFile){
	
		vector<long> marginalBytes;
		vector<long> bytes;
		vector<long> gridResolutions;
		
		js["digithist"]["gridresolution"] = gridResolutions;
		js["digithist"]["bytes"] = bytes;
		//js["digithist"]["marginalBytes"] = ;
		js["digithist"]["k"] = digitHistograms.size();
		
		ofstream jsout(jsonFile);
		
		jsout << js;
		
		jsout.close();
	
	}
	
	inline void writeBinary(const string& binaryFile){
		
		ofstream out(binaryFile, std::ios::out | std::ios::binary);
		
		DigitHistBinary bin;
		
		const long byteNum = bin.write<true>(NULL, dataBox, digitHistograms);
		
		vector<byte> vec(byteNum, 0);
		
		bin.write(vec.data(), dataBox, digitHistograms);
		
		char* c = reinterpret_cast<char*>(vec.data() );
		
		out.write(c, byteNum);
		out.close();
	}

	#if defined(SHEKELYAN_SELIMAGE) 
	
	inline void visualize1(DataImage& img) const{
	
		array<float,3> black = {0,0,0};
		array<float,3> green = {0,100,0};
		array<float,3> violet = {255,0,255};
		array<float,3> red = {255,0,0};
		array<float,3> blue = {0,0,255};
		array<float,3> gray = {100,100,100};
		
		const int K = digitHistograms.size();
		
		double sum = 1.0;
		
		for (int k = 0; k < K; k++){
			
			unique_ptr<RangeQueries> buckets = std::move( digitHistograms[k]->getBuckets(img.getDataBox()) );
			
			for (auto it = buckets->ranges.begin(); it != buckets->ranges.end(); it++)
				sum += it->count;
		}
		
		for (int k = 0; k < K; k++){
			
			unique_ptr<RangeQueries> buckets = std::move( digitHistograms[k]->getBuckets(img.getDataBox()) );
			
			Box b;
			
			for (auto it = buckets->ranges.begin(); it != buckets->ranges.end(); it++){
			
				b = (*it);
				b.count /= sum;
				
				img.fillRect(b, black);	
			}
		}
	}
	
	inline void visualize(DataImage& img) const{
	
		array<float,3> black = {0,0,0};
		array<float,3> green = {0,100,0};
		array<float,3> violet = {255,0,255};
		array<float,3> red = {255,0,0};
		array<float,3> blue = {0,0,255};
		array<float,3> gray = {100,100,100};
		
		const int K = digitHistograms.size();
		
		Sorter s;
		
		for (int k = 0; k < K; k++)
			s.add(k, -digitHistograms[k]->getMultiplier() );
		
		for (int kk = 0; kk < K; kk++){
		
			for (int k = 0; k < K; k++){
			
				if (k != s.getSortedIndex(kk))
					continue;
			
				shared_ptr<RangeQueries> buckets = digitHistograms[k]->getBuckets(img.getDataBox());
				
				if (k == 3){
				
					for (auto it = buckets->ranges.begin(); it != buckets->ranges.end(); it++)
						img.drawRect( (*it), green, 1.0);
				}
				
				if (k == 2){
			
					for (auto it = buckets->ranges.begin(); it != buckets->ranges.end(); it++)
						img.drawRect( (*it), violet, 1.0);
				}
			
				if (k == 1){
			
					for (auto it = buckets->ranges.begin(); it != buckets->ranges.end(); it++)
						img.drawRect( (*it), blue, 0.2);
				}
			
				if (k == 0){
			
					for (auto it = buckets->ranges.begin(); it != buckets->ranges.end(); it++)
						img.drawRect( (*it), red, 1.0);
				}
			}
		}
	}
	
	#endif
	
};
#endif