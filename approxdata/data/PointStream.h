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

#ifndef SHEKELYAN_POINTSTREAM_H
#define SHEKELYAN_POINTSTREAM_H

#include <approxdata/data/data.h>

class PointStream{

public:

	virtual ~PointStream(){}
	inline virtual bool read(Point& p) = 0;
	
	inline virtual bool skip(long off) = 0;
	inline virtual void close() = 0;
};


class RandomGenerator{

public:

	std::seed_seq seed1;
	std::mt19937 generator;
	
	std::uniform_real_distribution<double> rnd;
	
	inline RandomGenerator(const string seed, double a = 0, double b = 1) : seed1(seed.begin(), seed.end()), generator(seed1), rnd(a,b){
	
	}
	
	inline double randVal(){
	
		return rnd(generator);
	}
};

class NormGenerator{

public:

	std::seed_seq seed1;
	std::mt19937 generator;
	
	std::normal_distribution<double> rnd;
	
	inline NormGenerator(const string seed, double mean = 0, double std = 1) : seed1(seed.begin(), seed.end()), generator(seed1), rnd(mean,std){
	
	}
	
	inline double randVal(){
	
		return rnd(generator);
	}
};

class ZipfPointStream : public PointStream{

public:
	
	const string seed1;
	const int DIMS;
	
	vector<long> clusterCards;
	vector<Point> clusterCenters;
	vector<double> clusterWidths;
	
	unique_ptr<NormGenerator> ng;
	
	long clusterIndex = -1;
	long pointIndex = -1;
	
	// seed1(seed.begin(), seed.end()),
	inline ZipfPointStream(const string seed, int dims, long clusters, long points) : seed1(seed.begin(), seed.end() ), DIMS(dims){
		
		addClusters(clusters, points);
	}
	
	inline ~ZipfPointStream() override{
	}
	
	inline void addCluster(const Point& p, long points, double std){
		
		clusterCenters.push_back(p);
		clusterCards.push_back(points);
		clusterWidths.push_back(std);
	}
	
	inline void addClusters(long clusters, long points){
	
		//cout << "addClusters(" << clusters << "," << points << ")" << endl;
	
		RandomGenerator rg(seed1, 0,1);
		
		Point p(DIMS);
		
		long cards = points;
		
		for (int k = 0; k < clusters; k++){
		
			const long card = (k >= (clusters-1)) ? cards : zipf(k+1, 1, points)*points;
			
			cards -= card;
			
			double w = zipfwidth(k+1, 1.0, clusters);
			
			for (int i = 0; i < DIMS; i++)
				p[i] = rg.randVal();
			 
			addCluster(p, card,  w);
		}
		
		clusterIndex = -1;
		pointIndex = -1;
		
	}
	
	inline double zipf(long k, double s, long N) const{
	
		double div = 0;
		
		for (long n = 1; n <= N; n++){
		
			div += (1/pow(n, s) );
		}
		
		return ( 1/pow(k,s) ) / div;
	} 
	
	inline double zipfwidth(long k, double s, long N) const{
	
		double div = 0;
		
		for (long n = 1; n <= N; n++){
		
			div += (1/pow(n, s) );
		}
		
		return 0.5/div;
	}
	
	inline bool skip(long off) override{
	
		Point p(DIMS);
	
		for (long i = 0; i < off; i++)
			read(p);
	}
	
	inline bool read(Point& p) override{
		
		if ( clusterIndex >= ((long) clusterCards.size()) )
			return false;
		
		if ( (pointIndex < 0) || (pointIndex >= clusterCards.at(clusterIndex)) ){
			
			clusterIndex++;
			
			if (clusterIndex >= ((long) clusterCards.size()) )
				return false;
			
			pointIndex = 0;
			
			//cout << "**" << clusterIndex << " " << pointIndex << endl;
			
			ng.reset( new NormGenerator(seed1+to_string(clusterIndex), 0, clusterWidths.at(clusterIndex)) );
		}
		
		assert (ng);
		
		const Point& cc = clusterCenters.at(clusterIndex);
		
		for (int i = 0; i < p.size(); i++)
			p[i] = cc[i]+ng->randVal();
		
		pointIndex++;
		
		//cout << "***" << clusterIndex << " " << pointIndex << endl;
		
		return true;
	}
	
	virtual inline void close() override{
		
		
	}
	
};

class VectorPointStream : public PointStream{

public:

	long pos = 0;
	
	shared_ptr<vector<Point>> vec;

	inline VectorPointStream(shared_ptr<vector<Point>> v){
		
		pos = 0;
		vec = v;
	}

	virtual ~VectorPointStream(){}
	
	inline bool skip(long off) override{
	
		if (off <= 0)
			return true;
	
		pos += off;
		
		if ((pos-1) >= vec->size() )
			return false;
		
		return true;
	}
	
	inline bool read(Point& p) override{
	
		if (pos >= vec->size() )
			return false;
		
		p = (*vec)[pos];
		pos++;
		
		return true;
	}
	
	inline void close() override{
	
		
	}
};


template<typename E>
class BufferedBinaryReader{

	protected:
	
		ifstream r;
		
		int vpos;
		
		vector<E> v;
		
		char* buf;
		
		long bufferSize;
	
	public:
		
		inline BufferedBinaryReader(const string f, int bufSize=8194) : r(f, std::ios::binary|std::ios::in), bufferSize(bufSize){
			
			const long valSize = bufSize/sizeof(E);
			v.reserve(valSize);
			
			for (int i = 0; i < valSize; i++)
				v.push_back(numeric_limits<E>::lowest() );
			
			buf = reinterpret_cast<char*>(v.data());
			
			vpos = v.size();
		}
		
		
		virtual inline void skip(long off){
			
			//assert (vpos <= v.size() );
			
			if (vpos == v.size() ){
			
				r.ignore(off);
				return;
			}
			
			const long bytesLeft = (v.size()-vpos)*sizeof(E);
			
			if (off <= bytesLeft){
			
			 	//assert ( (off/sizeof(E))*sizeof(E) == off);
				
				vpos += off/sizeof(E);
				return;
			}
			
			vpos += bytesLeft/sizeof(E);
			off -= bytesLeft;
			
			//assert (off >= 0);
			
			if (off > 0)
				r.ignore(off);
			
			vpos = v.size();
		}
		
		inline E read(){
		
			//assert (vpos <= v.size() );
			
			if (vpos == v.size() ){
				
				r.read( buf, bufferSize);
				vpos = 0;
			}
			
			const E ret = v[vpos];
			
			//assert (ret != numeric_limits<E>::lowest());
			
			vpos++;
			
			return ret;
		}
		
		inline void close(){
		
			r.close();
		}
};

class BinaryReader{

	protected:
	
		ifstream r;
		
		std::vector<char> vec;

	public:
		
		inline BinaryReader(string file, long buf=8192) : vec(buf){
		
			r.open(file, std::ios::binary|std::ios::in);
			
			if (!r.is_open())
				throw Exception("file not open! "+file);
			
			r.rdbuf()->pubsetbuf(&vec.front(), vec.size() );
		}
		
		virtual inline void skip(long off){
		
			r.rdbuf()->pubseekoff(off, std::ios_base::cur, std::ios_base::in);
		}
		
		virtual inline uint64_t readLongLong(){
		
			uint64_t ret;
			r.read( reinterpret_cast<char*>( &ret ), sizeof(ret) );
			return ret;
		}
		
		virtual inline uint64_t readUnsignedLongLong(){
		
			uint64_t ret;
			r.read( reinterpret_cast<char*>( &ret ), sizeof(ret) );
			return ret;
		}
		
		virtual inline long readLong(){
		
			long ret;
			r.read( reinterpret_cast<char*>( &ret ), sizeof(ret) );
			return ret;
		}
		
		virtual inline int readInt(){
		
			int ret;
			r.read( reinterpret_cast<char*>( &ret ), sizeof(ret) );
			return ret;
		}
		
		virtual inline string readString(){
		
			const int bytes = readInt();
			char* cs = new char[bytes];
			
			r.read(cs, bytes);
			string ret(cs);
			
			return ret;
		}
		
		virtual inline double readDouble(){
		
			double ret;
			r.read( reinterpret_cast<char*>( &ret ), sizeof(ret) );
			return ret;
		}
		
		virtual inline void close(){
			
			r.close();
		}
};

class BinaryPointStream : public PointStream{

	public:
	
		BufferedBinaryReader<double> r;
		//BinaryReader r;
		
		const int DIMS;
		
		inline BinaryPointStream(string file, int dims, long skipBytes=0, long buf=8192) : r(file, buf), DIMS(dims){
			
			if (skipBytes > 0)
				r.skip(skipBytes);
		}
		
		inline bool read(Point& p) override{
	
			for (int i = 0; i < DIMS; i++)
				p[i] = r.read(); // r.readDouble();
				
			return true;
			
		}
		
		inline bool skip(long off) override{
		
			if (off <= 0)
				return true;
			
			r.skip(DIMS*sizeof(double)*off);
			
			return true;
		}
		
		inline void close() override{
			
			r.close();
		}
};

class CSVPointStream : public PointStream{

	public:
		
		char sep = ',';
		
		bool skipFirstLine = true;
		
		ifstream r;

		int DIMS = 0;

		inline CSVPointStream(string file){
			
			r.open(file);
			assert (r.is_open());
			
			string line;
	
			if (getline (r,line))
				DIMS = std::count(line.begin(), line.end(), sep)+1;	
		}
		
		inline void close() override{
	
			r.close();
		}
		
		virtual inline bool read(Point& p) override{
	
			string line;
	
			if (!getline (r,line))
				return false;
			
			stringstream lineStream(line);
			
			string cell;
			
			for (int i = 0; i < DIMS; i++){
			
				if (!getline(lineStream, cell, sep )){
					
					assert (false);
					return false;
				}
				
				istringstream numStream(cell);
			
				double val;
			
				if (!(numStream >> val))
					return false;
				
				if (!(val == val))
					throw Exception("not a number!");
				
				p[i] = val;
			
			}
			
			return true;
		}
		
		virtual inline bool skip(long off){
		
			if (off <= 0)
				return true;
		
			string line;
			
			for (long j = 0; j < off; j++){
				
				if (!getline (r,line))
					return false;
			}
			
			return true;
		}
	
		
};


#endif