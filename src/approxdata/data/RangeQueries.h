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

#ifndef SHEKELYAN_RANGEQUERIES_H
#define SHEKELYAN_RANGEQUERIES_H

#include <approxdata/data/data.h>

enum QueryDistribution{

	BIASED, UNBIASED

};

class RangeQueries{

public:
	vector<Box> ranges;
	
	inline RangeQueries(){
	
	
	}
	
	inline RangeQueries(const string file){
		
		add(file);
	}
	
	inline void addLegacy(const string file){
	
		BinaryReader r(file);
	
		const long size = r.readLong();
		const int DIMS = r.readInt();
		
		assert (size > 0);
		assert (size < 1000000000);
		
		for (long j = 0; j < size; j++){
			
			Point qmin(DIMS);
			Point qmax(DIMS);
			
			for (int i = 0; i < DIMS; i++)
				qmin[i] = r.readDouble();
			
			for (int i = 0; i < DIMS; i++)
				qmax[i] = r.readDouble();
			
			Box b(qmin, qmax, r.readDouble() );
			
			ranges.push_back(b);
		}
		
		r.close();
	}
	
	inline RangeQueries(DataSet& data, long target, double selmin, double selmax, bool biased=true, const string seed0 = "abc"){
	
		add(data, biased, selmin, selmax, target, seed0);
	}
	
	inline void print(){
	
		cout << "number of ranges " << ranges.size() << endl;
		
		long k = 0;
		for (auto it = ranges.begin(); it != ranges.end(); it++){
		
			if (k < 5 || k > (ranges.size()-10)){
			
				it->print();
			}
			k++;
		}
	}
	
	inline void add(const nlohmann:json& j){
		
		const long size = j["queries"].size();
		
		for (long k = 0; k < size; k++){
			
			long count = j["queries"][k]["count"];
			
			Point min(j["queries"][k]["min"].size() );
			
			for (int i = 0; i < min.size(); i++)
				min[i] = j["queries"][k]["min"][i];
				
			Point max(j["queries"][k]["max"].size() );
			
			for (int i = 0; i < max.size(); i++)
				max[i] = j["queries"][k]["max"][i];
				
			Box b(min, max, count);
			
			ranges.push_back(b);
		}
	}
	
	inline void add(string file){
	
		ifstream in(file);
		
		nlohmann::json j;
		
		in >> j;
		
		in.close();
		
		add(j);
	}
	
	inline void write(string file){
	
		ofstream out(file);
		
		out << "{" << endl;
		out << "\"queries\":[" << endl << endl;
		
		for (long j = 0; j < ranges.size(); j++){
		
			if (j > 0)
				out << ",";
		
			out << "{" << endl;
			
			out << "\"count\":" << UtilString::doubleToString(ranges[j].count) << "," << endl;
			out << "\"min\": " << ranges[j].min << "," << endl;
			out << "\"max\": " << ranges[j].max << "" << endl;
			out << "}";
		}
		
		out << endl << endl << "]" << endl;
		out << "}" << endl;
		
		out.close();
	}
	
	inline void boxCount(DataSet& data){
		
		const double mode0ops = data.getSize()*size();
		const double mode1ops = (size()*size())+data.getSize()*(0.01*size());
		
		BoxIndex bi(data.getDims(), ranges, ((mode0ops < mode1ops) || (size() < 10)) ? 0 : 1);	
		
		for (auto it = data.begin(); it != data.end(); it++ )
			bi.count( (*it) );
		
	}
	
	inline void boxCount(shared_ptr<DataSet> data){
		
		boxCount(*data);
		
	}
	
	
	inline void add(shared_ptr<DataSet> data, bool getBiased, double sel, long target, const string seed0 = "abc"){
		
		const int DIMS = data->getDims();
		
		const bool oldVerbose = data->getVerbose();
		
		data->setVerbose(false);
		
		shared_ptr<DataSet> biased = data->getRandomSubset(target, seed0+" 1");
		shared_ptr<DataSet> unbiased = data->getRandomPoints(target, seed0+" 2");
		
		shared_ptr<DataSet> unbiased2 = data->getRandomPoints(target, seed0+" 3");
		shared_ptr<DataSet> unbiased3 = data->getRandomPoints(target, seed0+" 4");
		
		shared_ptr<DataSet> centers = getBiased ? biased : unbiased;
		
		Sorter sorter;
		
		sorter.reserve(getSize() );
		
		const long nn = UtilMath::truncate(2, data->getSize(), round( data->getCount(sel) ) );
		
		Timer t("createQueries tar("+to_string(target)+") nn("+to_string(nn)+")");
		
		t.start(5, 10, centers->getSize() );
		
		auto centersIt = centers->begin();
		auto unbias2It = unbiased2->begin();
		auto unbias3It = unbiased3->begin();
		
		for (; centersIt != centers->end(); centersIt++, unbias2It++, unbias3It++){
		
			Point p = (*centersIt);
			Point p1 = (*unbias2It);
			Point p2 = (*unbias3It);
			
			Box b;
			
			b.enclose(p1);
			b.enclose(p2);
			
			for (int i = 0; i < DIMS; i++){
			
				const double delta = p[i]-(b.max[i]*0.5+b.min[i]*0.5);
				
				b.min[i] = UtilMath::truncate(0, 1, b.min[i]+delta);
				b.max[i] = UtilMath::truncate(0, 1, b.max[i]+delta);
			}
			
			sorter.clear();
			
			{
			long pos = 0;
			for (auto it = data->begin(); it != data->end(); it++, pos++)
				sorter.add( pos, b.volumeDistance( (*it) ));
			}
			
			Box rbox;
			
			for (long j = 0; j < nn; j++){
			
				auto it = data->begin();
				
				it += sorter.getSortedIndex(j);
				
				Point pp = (*it);
				
				rbox.enclose( pp);
			}
			
			add(rbox);
			
			t.tick();
		}
		
		t.end();
		
		data->setVerbose(oldVerbose);
		
	}
	
	
	
	
	
	inline void add(DataSet& data, bool getBiased, double selmin, double selmax, long target, const string seed0 = "abc"){
	
		const int DIMS = data.getDims();
		
		shared_ptr< RangeQueries > ret2(new RangeQueries() );
		
		shared_ptr<DataSet> big = data.getRandomSubset(UtilMath::minVal<long>(data.getSize(), 1000000), seed0);
		
		big->setVerbose(false);
		
		for (long jj = 0; ( ret2->getSize() < target ) && ( jj < 100 ) ; jj++){
		
			string seed1 = seed0+" "+to_string(jj);
		
			shared_ptr< RangeQueries > ret(new RangeQueries() );
		
			shared_ptr<DataSet> smp = big->getRandomSubset(10000, seed1);
		
			for (long j = 0; ( ret->getSize() < target ) && ( j < 100 ) ; j++){
		
				string seed2 = seed1+" "+to_string(j);
		
				shared_ptr< RangeQueries > cand(new RangeQueries() );
			
				cand->add(smp, getBiased, selmin*0.5+selmax*0.5, target*2, seed2);
				
				cand->boxCount(big);
				
				shared_ptr< RangeQueries > toAdd = cand->getWithCounts( big->getCount(selmin), big->getCount(selmax));
				
				ret->add( toAdd->getRandomSubset(target, seed2) );
			}

			ret->boxCount(data);

			ret2->add( ret->getWithCounts( data.getCount(selmin), data.getCount(selmax)) );		
			//cout << "ret2->getSize() " << ret2->getSize() << endl;
		}
		
		add( ret2->getRandomSubset(target, seed0) );
		
	}
	
	inline void add(shared_ptr<DataSet> data, bool getBiased, double selmin, double selmax, long target, const string seed0 = "abc"){
	
		add( (*data), getBiased, selmin, selmax, target, seed0);
	}
	
	inline void print(long size){
	
		
		long k = 0;
		for (auto it = ranges.begin(); it != ranges.end(); it++){
		
			if (k < 5 || k > (ranges.size()-10)){
			
				it->print(size);
			}
			k++;
		}
	}
	
	
	const Box& operator[](std::size_t idx){
		
		return ranges[idx];
	}
	
	
	inline void add(Box& cb){
	
		ranges.push_back(cb);
	}
	
	/*
	inline void add(Box& b){
	
		Box cb(b.min, b.max, 0);
	
		ranges.push_back(cb);
	}*/
	
	inline void add(shared_ptr<vector<Box>> b){
	
		add( (*b) );
	}
	
	inline void setBox(long k, Box& b){
	
		ranges[k].setMin(b.min);
		ranges[k].setMax(b.max);
	}
	
	inline void add(vector<Box>& b){
	
		for (auto it = b.begin(); it != b.end(); it++){
		
			add( (*it) );
		}
	}
	
	
	inline void add(RangeQueries& rq){
	
		add(rq.ranges);
	}
	
	inline void add(shared_ptr<RangeQueries> rq){
	
		add( (*rq) );
	}
	
	
	inline long size(){
	
		return ranges.size();
	}
	
	inline long getSize(){
	
		return ranges.size();
	}
	
	inline void setCount(long ind, long count){
	
		ranges[ind].count = count;
	}
	
	inline shared_ptr<RangeQueries> getWithCounts(double min, double max){
	
		shared_ptr<RangeQueries> ret(new RangeQueries() );
		
		for (auto it = ranges.begin(); it != ranges.end(); it++)
			if ( it->count >= min && it->count <= max)
				ret->add( (*it) );
		
		return ret;
	}
	
	
	inline shared_ptr<RangeQueries> getRandomSubset(long s, const string seed="123"){
		
		if (s >= size() ){
		
			shared_ptr<RangeQueries> ret(new RangeQueries() );
			
			for (auto it = ranges.begin(); it != ranges.end(); it++)
				ret->add( (*it) );
			
			return ret;
		}
		
		vector<long> inds(s);
		
		ReservoirSampling rs(seed, s);
		
		for (long pos = 0; pos < size(); pos++){
		
			const long j = rs.pos();
			
			if ( (j >= 0) && (j < s))
				inds[j] = pos;
		}
		
		shared_ptr<RangeQueries> ret(new RangeQueries() );
		
		for (auto it = inds.begin(); it != inds.end(); it++)
			ret->add( ranges[(*it)] );
		
		return ret;
	}
};


class RangeQueryFile{
public:
	vector<Box> boxes;
	
	vector<double> lbs;
	vector<double> ests;
	vector<double> ubs;
	
	inline void read(const string file){
	
		ifstream in(file);
		
		json js;
		
		in >> js;
		
		in.close();
		
		const long qs = js["queries"].size();
		
		for (long j = 0; j < qs; j++){
		
			double count = 0;
			double est = -1;
			double lb = -1;
			double ub = -1;
			
			Box b;
			
			try{ 
				count = js["queries"][j]["count"];
			} catch (std::exception& e){
			
			}
			
			b.setCount(count);
			
			try{
				est = js["queries"][j]["estimate"];
			} catch (std::exception& e){
			
			}
			
			
			try{
				lb = js["queries"][j]["bounds"][0];
				ub = js["queries"][j]["bounds"][1];
			} catch (std::exception& e){
			
			}
			
			try{
			
				Point pmin;
				Point pmax;
				
				for (int i = 0; i < js["queries"][j]["min"].size(); i++)
					pmin.push_back(  js["queries"][j]["min"][i] );
				
				for (int i = 0; i < js["queries"][j]["max"].size(); i++)
					pmax.push_back(  js["queries"][j]["max"][i] );
				
				b.setMin(pmin);
				b.setMax(pmax);
			} catch (std::exception& e){
			
			}
			
			add(b, est, lb, ub);
		}
		
		in.close();
	}
	
	inline void write(const string file){
	
		ofstream out(file);
		
		out << "{ \"queries\":" << endl << "[" << endl;
		
		for (long j = 0; j < boxes.size(); j++){
		
			if (j != 0)
				out << "," << endl;
			
			out << "{";
			
			out << "\"min\": " << boxes[j].min;
			out << ", \"max\": " << boxes[j].max;
			out << ", \"count\": " << boxes[j].count;
			
			if ( (lbs.size() != 0) && (ubs.size() != 0))
				out << ", \"bounds\": [" << lbs[j] << "," << ubs[j] << "]";
			
			if ( (ests.size() != 0) )
				out << ", \"estimate\": " << ests[j];
			
			out << "}";
		}
		
		out << endl << "]" << endl;
		
		out << endl << "}" << endl;
	
		out.close();
	}
	
	inline void add(const Box& b, double lb=-1, double est=-1, double ub=-1){
	
		boxes.push_back(b);
			
		if (est != -1)
			ests.push_back(est);
			
		if (lb != -1)
			lbs.push_back(lb);
			
		if (ub != -1)
			ubs.push_back(ub);
	}


};




#endif