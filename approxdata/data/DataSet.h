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

#ifndef SHEKELYAN_DATASET_H
#define SHEKELYAN_DATASET_H

#include <approxdata/data/data.h>
#include <approxdata/external/json.h>

using json = nlohmann::json;


class CSVReader{
public:
	std::ifstream in;
	
	inline CSVReader(const string file): in(file){
	
	}
	
	bool isNumberChar(char c) const{

		if ((c >= '0') && (c <= '9'))
			return true;
	
		if (c == '-')
			return true;
		
		if (c == '.')
			return true;
		
		return false;
	}
	
	inline void writeDataFiles(const string out, const std::vector<int>& selectedDims){
	
		std::vector<double> p;
		
		std::vector<double> min;
		std::vector<double> max;
		
		ofstream dat;
		dat.open(out+".dat", std::ios::binary|std::ios::out);
	
		int dims = -1;
	
		long size = 0;
	
		if (in.is_open()) {
				
			std::string line;
		
			while (std::getline(in, line)){
		
				double d;
				unique_ptr<std::stringstream> s(new std::stringstream() );
				
				
				bool number = false;
			
				for (int i = 0; i < line.length(); i++){
				
					const char c = line[i];
			
					if (c == '"'){
					
						i++;
						while (line[i] != '"')
							i++;
					
						continue;
					}
			
					if ( isNumberChar(c) ){
				
						number = true;
						(*s) << c;
						continue;
				
					} else if (number){
				
						if (!((*s) >> d)){
							assert (false);
						}
				
						s.reset( new std::stringstream());
					
						p.push_back(d);
				
						number = false;
					}
	
				}
			
				if ((*s) >> d)
					p.push_back(d);
			
				s.reset( new std::stringstream());
			
				number = false;
			
				if (p.size() > 0){
					
					if (dims == -1){
						dims = p.size();
					
						for (int i = 0; i < selectedDims.size(); i++){
						
							min.push_back( std::numeric_limits<double>::max() );	
							max.push_back( std::numeric_limits<double>::lowest() );
						}	
					}
					
					for (int j = 0; j < selectedDims.size(); j++){
					
						const int i = selectedDims[j];
						
						const double v = p[i];
						
						if (v < min[j])
							min[j] = v;
							
						if (v > max[j])
							max[j] = v;
							
						dat << v;
					}
					size++;
					
					assert (p.size() == dims);
				}
			
				p.clear();
		
			}
		}
		
		dat.close();
		
		ofstream json(out+".json");
		
		json << "{" << endl;
		
		json << "\t\"name\": \"" << out << "\"," << endl;
		
		json << "\t\"data\": [" << endl;
		
		json << endl;
		json << "\t\t{" << endl;
		
		json << "\t\t\t\"file\": \"" << out <<".dat\"," << endl;
		
		json << "\t\t\t\"type\": \"multidimensional\"," << endl;
		json << "\t\t\t\"size\": " << size << "," << endl;
		json << "\t\t\t\"dimensionality\": " << selectedDims.size() << "," << endl;
		
		json << endl;
		
		json << "\t\t\t\"dimensions\":[" << endl;
		json << endl;
		
		for (int i = 0; i < selectedDims.size(); i++){
		
			json << "\t\t\t\t{\"domain\": \"real\", \"boundaries\":[" << min[i] << "," << max[i] << "]}";
			
			if (i != selectedDims.size()-1)
				json << ",";
				
			json << endl;
		}
		
		json << endl;
		json << "\t\t\t\]" << endl;;
		
		
		json << "\t\t}" << endl;
		
		json << "\t]" << endl;
		
		json << "}" << endl;
		
		
		/*
		ini << "[description]" << endl << endl;
		
		ini << "file = " << out << ".dat" << endl;
		
		ini << "type = real" << endl;
		
		ini << "dims = " << selectedDims.size() << endl;
		
		ini << "size = " << size << endl;
		
		ini << endl << "[boundaries]" << endl << endl;
		
		for (int i = 0; i < selectedDims.size(); i++){
		
			ini << "min" << i << " = " << min[i] << endl;
			ini << "max" << i << " = " << max[i] << endl;
		}
	
		ini.close();*/
		
		json.close();
		
		
		in.close();
	}
};


// meta information about dataset to help read it etc
class DataDescription{

	private:
		vector<double> domainAdd;
		vector<double> domainMul;
		vector<double> domainDiv;
		
		vector<string> name = {"noname"};
		vector<string> dataFile = {"nofile"};
		
		vector<long> domainDistinct;
		vector<string> types;
		
		long size = 0;
		
		int dims = 0;
		
		Box boundingBox;
		

	public:
		
		inline DataDescription(int dims, string type="real"): dims(dims){
			
			setName("noname");
		}
		
		inline long getSize() const{
		
			return size;
		}
		
		inline void setSize(long size){
		
			this->size = size;
		}
		
		inline int getDims() const{
		
			return dims;
		}
		
		inline void setName(const string name){
		
			this->name.clear();
			this->name.push_back(name);
		}
		
		inline void setDataFile(const string file){
		
			this->dataFile.clear();
			this->dataFile.push_back(file);
		}
		
		inline DataDescription(const string file){
			
			
			if (UtilString::endsWith(file, "json")){
				readJSON(file);
				return;
			}
			
			assert (false);
			/*
			setName("data");
			Point domainMin;
			Point domainMax;
			
			if (UtilString::endsWith(file, "ini")){
			
				IniReader r(file);
				dims = r.getInt("description", "dims");
			
				assert (dims > 0);
			
				size = r.getInt("description", "size");
				setDataFile( r.getString("description", "file") );
				
				
				for (int i = 0; i < dims; i++){
			
					domainMin.push_back(r.getDouble("boundaries", "min"+to_string(i)));
					domainMax.push_back(r.getDouble("boundaries", "max"+to_string(i)));
				
					domainDistinct.push_back(r.getInt("boundaries", "distinct"+to_string(i)));
				}
			}
			
			
			Box box;
			
			box.setMin(domainMin);
			box.setMax(domainMax);
			
			setBoundingBox(box);*/
		}
		
		inline void readJSON(const string file){
		
			ifstream in(file);
				
			json j;
	
			j << in;

			dims = j["data"][0]["dimensionality"];
			size = j["data"][0]["size"];
	
			setName(j["name"]);
	
			setDataFile( j["data"][0]["file"] );
	
			Point domainMin(dims);
			Point domainMax(dims);
	
			for (int i = 0; i < dims; i++){
	
				domainMin[i] = j["data"][0]["dimensions"][i]["boundaries"][0];
				domainMax[i] = j["data"][0]["dimensions"][i]["boundaries"][1];
		
				domainDistinct.push_back(numeric_limits<long>::max() );
			}
			
			setBoundingBox(Box(domainMin, domainMax, -1));
		}
		
		inline void writeJSON(const string file){
			
			assert (false); // not yet implemented
		}
		
		inline const string& getDataFile() const{
		
			return dataFile[0];
		}
		
		inline void setBoundingBox(const Box& b){
			
			boundingBox = b;
			
			domainAdd.clear();
			domainMul.clear();
		
			for (int i = 0; i < dims; i++){
			
				domainAdd.push_back(-boundingBox.min[i]);
				domainMul.push_back( 1.0/(boundingBox.max[i]-boundingBox.min[i]) );
			}	
		}
		
		inline const Box& getBoundingBox() const{
		
			return boundingBox;
		}
		
		inline const string& getName() const{
		
			return name[0];
		}
		
		inline double normalise(double d, int i) const{
			
			d += domainAdd[i];
			d *= domainMul[i];
			
			return d;
		}
		
		inline void normalise(Point& d) const{
		
			for (int i = 0; i < dims; i++){
			
				d[i] += domainAdd[i];
				d[i] *= domainMul[i];
			}
		}
		
		inline double denormalise(double d, int i) const{
			
			d /= domainMul[i];
			d -= domainAdd[i];
			
			return d;
		}
		
		inline void normalise(Point& d, const vector<int>& dimVec) const{
			
			int k = 0;
			for (auto it = dimVec.begin(); it != dimVec.end(); it++){
			
				const int i = (*it);
				d[k] += domainAdd[i];
				d[k] *= domainMul[i];
				k++;
			}
		}
		
		
		
		inline void print() const{
		
			cout << "dataFile " << getDataFile() << endl;
			cout << "dims " << getDims() << endl;	
			cout << "size " << getSize() << endl;	
			
		}
		
		inline void processTuple(Point& t){
		
			boundingBox.enclose(t);	
			size++;
		}
		
		inline string toString() const{
		
			std::stringstream ss;
	
			ss << "data \t file" << getDataFile() << " dims " << getDims() << " size " << getSize() << endl;
			
			ss << boundingBox << endl;
			
			return ss.str();
		}
	
		friend inline std::ostream& operator<<(std::ostream &os, const DataDescription& m) { 
			
			return os << m.toString();
		}
};


class IntSeq{



	vector<int> seq;

public:


	inline IntSeq(){
	
	
	}

	inline IntSeq(const IntSeq& s){
	
		seq = s.seq;
	}

	inline IntSeq(IntSeq& s){
	
		seq = s.seq;
	}
	
	inline const vector<int>& getSequence(){
	
		return seq;
	}

	inline IntSeq(const string s){
		
		stringstream commaSplit(s);
		string commaPart;
		
		while ( getline(commaSplit, commaPart, ',') ){
		
			stringstream rangeSplit(commaPart);
			string rangePart;
		
			int a = -1;
			int b = -1;
			
			if (getline(rangeSplit, rangePart, ':')){
			
				istringstream numStream(rangePart);
				numStream >> a;
				b = a;
				
			}
			
			if (getline(rangeSplit, rangePart, ':')){
			
				istringstream numStream(rangePart);
				numStream >> b;
			}
			
			assert (a != -1 && b != -1);
			
			if (a <= b){
			
				for (int i = a; i <= b; i++)
					seq.push_back(i);
				
			} else {
			
				for (long i = a; i >= b; i--)
					seq.push_back(i);
			}
		}
		
	}
	
	inline int get(int i) const{
	
		return seq[i];
	}
	
	inline bool operator==(const IntSeq& s) const{
	
		if (seq.size() != s.seq.size() )
			return false;
	
		vector<int> a = seq;
		vector<int> b = s.seq;
		
		sort(a.begin(), a.end());
		sort(b.begin(), b.end());
		
		for (int j = 0; j < a.size(); j++){
		
			if (a[j] != b[j])
				return false;
		}
		
		return true;
	}
	
	inline bool contains(int i) const{
	
		for (auto it = seq.begin(); it != seq.end(); it++)
			if ( (*it) == i)
				return true;
				
		return false;
	}
	
	
	inline int size() const{
	
		return seq.size();
	}
	
	inline string toString() const{
		
		std::stringstream ss;
	
		for (int i = 0; i < seq.size(); i++){
			if (i > 0)
				ss << ",";
			
			ss << to_string(seq[i]);
		}
		
		return ss.str();
	}
	
	friend inline std::ostream& operator<<(std::ostream &os, const IntSeq& m) { 
    	
    	return os << m.toString();
	}

};

class DimSelection{


	public:
	
	
		IntSeq dataDims;
		IntSeq selectedDims;
		
		vector<int> ind;
		
		int dims = 0;
		
		bool eq = false;
		
		// selected dimensions (first dimension is "1") given as comma-separated string, e.g., "1,4,5"
		
		
		inline void setDataDims(const string s){
		
			dataDims = IntSeq(s);
			
			selectDims(s);
			
			assert (dataDims == selectedDims);
		}
		
		inline int getDataDims() const{
		
			return dataDims.size();
		}
		
		inline int getDims() const{
		
			return selectedDims.size();
		}
		
		inline void selectDims(const string s){
		
			selectedDims = IntSeq(s);
			
			dims = selectedDims.size();
			
			ind.clear();
			
			int dim = 0;
			
			for (int j = 0; j < dataDims.size(); j++){
			
				if (selectedDims.contains( dataDims.get(j) ))
					ind.push_back(dim);
				
				dim++;
			}
			
			eq = (selectedDims == dataDims);
			
			assert (ind.size() == selectedDims.size() );
		}
		
		
		inline void print(){
		
			cout << "dims: ";
			UtilString::printList(ind);
		}
		
		inline string getFileString() const{
		
			std::stringstream ss;
		
			for (int i = 0; i < dims; i++){
				if (i > 0)
					ss << "_";
				
				ss << to_string(ind[i]);
			}
			
			return ss.str();
		}
		
		inline string getString() const{
		
			std::stringstream ss;
		
			for (int i = 0; i < dims; i++){
				if (i > 0)
					ss << ",";
				
				ss << to_string(ind[i]+1);
			}
			
			return ss.str();
		}
		
		inline void translate(Point& highdim, Point& lowdim) const{
			
			if (eq){
			
				lowdim = highdim;
				
			} else {
			
				for (int i = 0; i < dims; i++)
					lowdim[i] = highdim[ind[i]];
			}
			
			
			
		}
};





class DataSet{

	class iterator{
	
	private:
	
		std::unique_ptr<PointStream> rawDataStream;
		std::unique_ptr<Timer> t;
		
		const DataSet* ds;
		
		Point rawPoint;
		Point current;
		long pos = 0;
		const long datasize;
		
		inline void read(){
		
			assert(!sortingIterator);
			
			rawDataStream->read(rawPoint);
			ds->dimSel.translate(rawPoint, current);
			ds->getDescription().normalise(current, ds->dimSel.ind);
		}
	
		// sorting related
		
		const bool verbose;
		const bool sortingIterator;
		
		int currentBlock = -1;
		long blockPos = 0;
		
		vector<Point> block;
		vector<Point> borders;
		
		vector<int> lex;
	
		const int DIMS;
		
		bool blockSorted = false;
		
		const PointComparator comp;
	
	public:
		
		
		inline iterator() : sortingIterator(false), comp(0), DIMS(0), datasize(0), block(0), borders(0), lex(0), verbose(false){
		
			
		}
		
		inline iterator(const DataSet* ds, bool verb) : sortingIterator(false), comp(0), verbose(verb), block(0), borders(0), lex(0), ds(ds), pos(0), rawPoint(ds->getRawDims() ), DIMS(ds->getDims() ),current(ds->getDims() ), datasize(ds->getSize() ){
		
			if (verbose)
				t.reset(new Timer(ds->getDescription().getName() ));
		
			if (verbose)
				t->start(10, 1000, datasize);
			
			rawDataStream.reset(ds->getStream() );
			
			read();
			
			pos = 0;
		}
		
		inline iterator(const DataSet* ds, const PointComparator& comp, long blockSize) : verbose(ds->getVerbose() ), comp(comp), sortingIterator(true), ds(ds), DIMS(ds->getDims() ), datasize(ds->getSize() ){
		
			block.reserve(blockSize);
			
			if (verbose){
				t.reset(new Timer(ds->getDescription().getName()+"sorted" ));
				t->start(10, 100, datasize);
			}
			
			
			shared_ptr<vector<Point>> s = std::move( ds->getRandomSubsetPoints(100000, "abc") );
			
			vector<Point>& sample = (*s);
			
			const long sampleBlockSize = round(0.9*blockSize*1.0/datasize*sample.size() );
			
			sort(sample.begin(), sample.end(), comp);
			
			{
				Point first(DIMS);
				first.setAll(std::numeric_limits<double>::lowest() );
				borders.push_back(first);
			}
			
			
			for (long j = sampleBlockSize; j < (sample.size()-sampleBlockSize); j += sampleBlockSize){
			
				assert(sample.at(j).size() == DIMS);
			
				borders.push_back(sample.at(j));
				
			}
			
			{
				Point last(DIMS);
				last.setAll(std::numeric_limits<double>::max() );
				borders.push_back(last);
			}
			
			(*this) += 0;
		}
		
		inline Point& operator*(){
		
        	return sortingIterator ? block[blockPos] : current;
    	}
    	
    	inline void close(){
    	
    		if (pos < datasize)
    			if (verbose)
    				t->end();
    	}
    	
    	inline Point& getPoint(){
		
			return sortingIterator ? block[blockPos] : current;
    	}
    	
    	inline Point& getRawPoint(){
		
			assert (!sortingIterator);
			
        	return rawPoint;
    	}
    	
    	inline const long operator()() const{
    		
    		return 0;
    	}
    	
    	inline long getPos() const{
    	
    		return pos;
    	}
    	
		inline bool operator !=(const long b) const{
		
			if (pos < datasize){
			
				return true;
				
			} else {
		
				if (verbose)
					t->end();
				
				return false;
			}
			
		}
		
		inline void loadNextBlock(){
	
			currentBlock++;
			
			block.clear();
			
			if (currentBlock >= (borders.size()-1) ){
			
				blockSorted = true;
				return;
			}
			
			PointStream* ps = ds->getStream();
			
			const Point& min = borders.at(currentBlock);
			const Point& max = borders.at(currentBlock+1);
			
			const int i = comp.getMainDim();
			
			const double mini = min[i];
			const double maxi = max[i];
			
			{
				
				for (auto it = ds->beginUnsorted(false); it != ds->end(); it++){
				
					const Point& p = (*it);
					const double pi = p[i];
			
					if (pi > maxi)
						continue;
			
					if (pi < mini)
						continue;
					
					if (comp.lessOrEqual(min, p) && comp.less(p, max) )
						block.push_back(p);
				}
				
			}
			
			//cout << "mini " << mini << " maxi " << maxi << endl;
			//cout << "block.size() " << block.size() << endl;
			
			blockSorted = false;
		}
		
		inline iterator& operator +=(const long n){
	
			if (sortingIterator){
			
				blockPos += n;
				pos += n;
			
				if (pos < datasize){
			
					while (blockPos >= block.size()){
			
						blockPos -= block.size();
						loadNextBlock();
					}
				
					if (!blockSorted){
			
						//UtilPoint::sort(lex, block);
						
						sort(block.begin(), block.end(), comp);
						
						//assert (UtilPoint::isSorted(lex, block) );
						blockSorted = true;
					}
				}
				
				if (verbose)
					t->tick();
				
				return (*this);
			}
	
			if (n <= 0)
				return (*this);
	
			pos += n;
			
			if (pos >= datasize)
				return (*this);
			
			if (verbose)
				t->tick();
			
			if (n > 1)
				rawDataStream->skip(n-1);
			
			read();
			
			return (*this);
		}
		
		inline void print(){
			
			const Point p = getPoint();
		
			cout << p << endl;
		}
		
		inline iterator& operator++(){
			
			return ( (*this) += 1 );
		}
		
		inline iterator& operator++(int n){
			
			return ( (*this) += 1 );
		}
		
		
	};


protected:
	
	inline DataSet(const string s){// : t(s){
	
		
	}
	
	
	
	inline virtual PointStream* getStream() const = 0;

	bool sortingIterator = false;
		
	//vector<int> sortingLex = {};
	
	PointComparator comp = 0;
	
	long sortingIteratorBuffer = 0;
	
	
	bool verbose = false;

public:

	inline const string& getFile(){
		
		return getDescription().getDataFile();
	}
		

	virtual ~DataSet(){};
	

	DimSelection dimSel;
	
	
	
	iterator end;
	
	inline virtual void enableSorting(const PointComparator& comp1, long sb= 1000000){
	
		comp = comp1;
		sortingIterator = true;
		sortingIteratorBuffer = sb;
	}
	
	inline void disableSorting(){
	
		sortingIterator = false;
	}
	
	iterator begin() const{
	
		const DataSet* ds = this;
		
		if (sortingIterator)
			return iterator(ds, comp, sortingIteratorBuffer);
		
		return iterator(ds, getVerbose());
	}
	
	
	iterator beginUnsorted(bool verb=true) const{
	
		const DataSet* ds = this;
		
		return iterator(ds, getVerbose() && verb);
	}
	
	inline virtual const DataDescription& getDescription() const = 0;
	
	inline void setVerbose(bool verbose){
	
		this->verbose = verbose;
	}
	
	inline bool getVerbose() const{
		
		return verbose;
	}

	inline virtual long getSize() const = 0;
	
	inline void add(Point& p){
	
		throw "modification not supported!";
	}
	
	inline long boxCount(const Point& qmin, const Point& qmax) const{
	
		long ret = 0;
	
		for (auto it = beginUnsorted(); it != end(); it++){
		
			if (UtilPoint::contains(qmin, qmax, (*it) ))
				ret++;
		}
		
		return ret;
	}


	inline void selectDims(const string s){
	
		dimSel.selectDims(s);
		
	}
	
	inline void setDataDims(const string s){
	
		dimSel.setDataDims(s);
	}
	
	inline void selectAllDims(){
		
		selectDims("1:"+to_string(getDescription().getDims() ));
	}
	
	inline void print(){
	
		getDescription().print();
		dimSel.print();
		
		const bool oldVerbose = verbose;
		verbose = false;
		
		for (auto it = beginUnsorted(); it != end(); it++){
		
			it.print();
		
			if (it.getPos() == 4)
				it += getSize()-10;
		}
		
		verbose = oldVerbose;
	}
	
	inline int getDims() const{
	
		return dimSel.getDims();
	}
	
	inline int getRawDims() const{
	
		return dimSel.getDataDims();
	}
	
	inline virtual shared_ptr<DataSet> getSorted(const PointComparator& comp, long memory=1000000) const = 0;
	inline virtual shared_ptr<DataSet> getRandomSubset(long size, const string seed="abc123") const = 0;
	inline virtual shared_ptr<DataSet> getRandomPoints(long size, const string seed="abc123") const = 0;
	inline virtual shared_ptr<DataSet> getCopy(shared_ptr<vector<Point>> v) const = 0;
	
	inline shared_ptr<DataSet> getEmptyCopy(long size=1) const{
		
		shared_ptr<vector<Point>> v(new vector<Point>());
		
		if (size >= 0)
			v->reserve(size);
		
		return getCopy(v);
	}
	
	inline virtual shared_ptr< vector<Point> > getRandomSubsetPoints(long size, const string seed="123") const{
		
		
		shared_ptr<vector<Point>> vec(new vector<Point>() );
		
		assert (size <= getSize() );
		
		if (size >= getSize() ){
		
			for (auto it = beginUnsorted(false); it != end(); it++)
				vec->push_back( (*it) );
			
			std::random_shuffle(vec->begin(), vec->end() );
			
			return vec;
		}
		
		shared_ptr<vector<long>> inds = UtilCollection::getRandomIndices( getSize(), size, seed);
		
		vec->reserve(size);
		
		auto it1 = beginUnsorted(false);
		long last = 0;
		
		long k = 0;
		
		for (auto it = inds->begin(); it != inds->end(); it++){
		
			const long ind = (*it);
		
			it1 += ind-last;
			
			Point p = (*it1);
			
			vec->push_back(p);
			
			last = ind;
		}
		
		it1.close();
		
		std::random_shuffle(vec->begin(), vec->end() );
		
		
		return vec;
	}
	
	inline virtual shared_ptr< vector<Point> > getRandomUnifPoints(long size, string seed="123") const{
		
		shared_ptr<vector<Point>> vec(new vector<Point>() );
		
		vec->reserve(size);
		
		UniformPoints unif(getRawDims(), seed );
		
		for (long j = 0; j < size; j++){
		
			Point p;
			
			unif.randomPoint(getDescription().getBoundingBox(), p);
			
			vec->push_back(p);
		}
		
		return vec;
	}
	
	
	const virtual Point operator[](std::size_t idx) const{
		
		auto it = begin();
		
		it += idx;
		
		if (it != end() )
			return (*it);
		
		throw string("out of bounds!");
	}
	
	inline double getCount(double selPercent) const{
	
		return selPercent/100.0*getSize();
	}
	
};


////////////////////////////////////////////////////////////////////////////////////////////////

class VectorDataSet : public DataSet{

private:
		
	DataDescription desc;

protected:
	
	inline PointStream* getStream() const override{
	
		return new VectorPointStream(vec);
	}

public:
	
	shared_ptr<vector<Point>> vec;
	
	inline VectorDataSet(const Box& b, shared_ptr<vector<Point>> vec, int dims) : desc(dims), vec(vec), DataSet("vector"){
	
		setDataDims("1:"+to_string(desc.getDims() ) );
		
		desc.setSize(vec->size());
		desc.setBoundingBox(b);
	}
	
	inline void add(Point& p){
	
		vec->push_back(p);
	}
	
	inline const DataDescription& getDescription() const override{
	
		return desc;
	}
	
	inline long getSize() const override{
	
		return vec->size();
	}
	
	inline shared_ptr<DataSet> getCopy(shared_ptr<vector<Point>> v) const override{
	
		Point pmin(getDims(), 0);
		Point pmax(getDims(), 1);
	
		Box b(pmin, pmax);	
		
		shared_ptr<DataSet> ret( new VectorDataSet(b, v , getDims() ));
		
		ret->setDataDims(dimSel.selectedDims.toString() );
		
		ret->setVerbose(getVerbose() );
		return ret;
	}
	
	inline shared_ptr<DataSet> getSorted(const PointComparator& comp, long memory=10000000) const override{
	
		shared_ptr<vector<Point>> vecCopy(new vector<Point>() );
		
		vecCopy->reserve(vec->size() );
		
		for (auto it = vec->begin(); it != vec->end(); it++)
			vecCopy->push_back( (*it) );
		
		shared_ptr<DataSet> ret( new VectorDataSet(desc.getBoundingBox(), vecCopy, getDims() ));
		ret->setDataDims(dimSel.selectedDims.toString() );
		ret->setVerbose(getVerbose() );
		ret->enableSorting(comp);
		
		return ret;
	}
	
	inline void enableSorting(const PointComparator& comp, long memory=0) override{
	
		sort( vec->begin(), vec->end(), comp );
		sortingIterator = false;
		sortingIteratorBuffer = memory;
	}
	
	inline shared_ptr<DataSet> getRandomSubset(long size, const string seed="123") const override{
		
		return getCopy( getRandomSubsetPoints(size, seed) );
	}
	
	inline shared_ptr<DataSet> getRandomPoints(long size, const string seed="123") const override{
		
		return getCopy( getRandomUnifPoints(size, seed) );
	}	
};

////////////////////////////////////////////////////////////////////////////////////////////////

class FileDataSet : public DataSet{

private:
		
	DataDescription desc;
	const string descFile;
		
protected:
	
	inline PointStream* getStream() const override{
	
		if (UtilString::contains(desc.getDataFile(), ".csv") )
			return new CSVPointStream(desc.getDataFile());
			
		if (UtilString::contains(desc.getDataFile(), ".pdat") ){	
		
			return new BinaryPointStream(desc.getDataFile(), desc.getDims(), 16);
		}
			
		if (UtilString::contains(desc.getDataFile(), ".dat") ){
		
			return new BinaryPointStream(desc.getDataFile(), desc.getDims());
		}
		
		return NULL;
	}

public:

	inline FileDataSet(const string descFile) : descFile(descFile), desc(descFile), DataSet(descFile){
	
		setDataDims("1:"+to_string(desc.getDims() ) );
	
		
		
	}
	
	inline shared_ptr<DataSet> getSorted(const PointComparator& comp, long memory=10000000) const override{
	
		shared_ptr<DataSet> ret( new FileDataSet(descFile));
		
		ret->setDataDims(dimSel.dataDims.toString() );
		ret->selectDims(dimSel.selectedDims.toString() );
		ret->setVerbose(getVerbose() );
		ret->enableSorting(comp, memory);
		
		return ret;
	}
	
	
	inline const DataDescription& getDescription() const override{
	
		return desc;
	}
	
	inline long getSize() const override{
	
		return desc.getSize();
	}
	
	inline shared_ptr<DataSet> getCopy(shared_ptr<vector<Point>> v) const override{
	
		Point pmin(getDims(), 0);
		Point pmax(getDims(), 1);
	
		Box b(pmin, pmax);	
		shared_ptr<DataSet> ret( new VectorDataSet(b, v, getDims() ));
		
		ret->selectDims(dimSel.selectedDims.toString() );
		
		ret->setVerbose(getVerbose() );
		return ret;
	}
	
	inline shared_ptr<DataSet> getRandomSubset(long size, const string seed="123") const override{
		
		return getCopy( getRandomSubsetPoints(size, seed ) );
	}
	
	inline shared_ptr<DataSet> getRandomPoints(long size, const string seed="123") const override{
		
		return getCopy( getRandomUnifPoints(size, seed ) );
	}
};


////////////////////////////////////////////////////////////////////////////////////////////////

class ZipfDataSet : public DataSet{

private:
		
	DataDescription desc;
	const string descFile;
	
	const string seed;
	const int DIMS;
	const long points;
	const long clusters;
	
protected:
	
	inline PointStream* getStream() const override{
		
		return new ZipfPointStream(seed, DIMS, clusters, points);
	}

public:

	inline ZipfDataSet(int dims, long pointNum, long clusterNum, const string s="abc") : seed(s), DIMS(dims), clusters(clusterNum), points(pointNum), descFile(""), desc(dims), DataSet("zipf"){
	
		setDataDims("1:"+to_string(desc.getDims() ) );
		
		desc.setSize(pointNum);
		
		Box b1(DIMS);
		desc.setBoundingBox(b1);
		
		for (auto it = this->begin(); it != this->end(); it++)
			b1.enclose( (*it) );
		
		desc.setBoundingBox(b1);
	}
	
	
	inline const DataDescription& getDescription() const override{
	
		return desc;
	}
	
	inline long getSize() const override{
	
		return desc.getSize();
	}
	
	inline shared_ptr<DataSet> getSorted(const PointComparator& comp, long memory=10000000) const override{
	
		shared_ptr<DataSet> ret( new ZipfDataSet(DIMS, points, clusters, seed));
		
		ret->setDataDims(dimSel.dataDims.toString() );
		ret->selectDims(dimSel.selectedDims.toString() );
		ret->setVerbose(getVerbose() );
		ret->enableSorting(comp, memory);
		
		return ret;
	}
	
	
	inline shared_ptr<DataSet> getCopy(shared_ptr<vector<Point>> v) const override{
	
		Point pmin(getDims(), 0);
		Point pmax(getDims(), 1);
		
		Box b(pmin, pmax);	
		shared_ptr<DataSet> ret( new VectorDataSet(b, v, getDims() ));
		
		ret->selectDims(dimSel.selectedDims.toString() );
		
		ret->setVerbose(getVerbose() );
		return ret;
	}
	
	inline shared_ptr<DataSet> getRandomSubset(long size, const string seed="123") const override{
		
		return getCopy( getRandomSubsetPoints(size, seed ) );
	}
	
	inline shared_ptr<DataSet> getRandomPoints(long size, const string seed="123") const override{
		
		return getCopy( getRandomUnifPoints(size, seed ) );
	}
};

#endif