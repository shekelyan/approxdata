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

#ifndef SHEKELYAN_DATASUMMARY_H
#define SHEKELYAN_DATASUMMARY_H

#include <approxdata/datasummaries/datasummaries.h>



class Params{

public:

	map<string, double> values;
	
	const string s;
	
	string command;
	
	inline Params(Params& p) : s(p.s){
	
		values = p.values;
	}
	
	~Params(){}
	
	inline Params() : s(""){
		
		
	}
	
	
	inline Params(const string s_) : s(s_){
	
		const string s__ = UtilString::replaceAll(s_ , "--", "-");
		const string s = UtilString::replaceAll(s__ , "=", " ");
		
		std::stringstream ss(s);
    	std::string item;
    	std::string item2;
    	
    	bool first = true;
    	
    	while (std::getline(ss, item, '-')) {
        	
        	std::stringstream ss2(item);
        	
        	int k = 0;
        	
        	string cmd = "undefined";
        	double val = 1.0;
        	
        	while (std::getline(ss2, item2, ' ')) {
        	
        		if (first){
        			
        			command = item2;
        			break;
        		}
        	
        	
        		if (k == 0){
        		
        			cmd = item2;
        		}
        		
        		if (k == 1){
        		
        			string clean = item2;
        			
        			clean.erase(std::remove(clean.begin(), clean.end(), ' '), clean.end());
        			clean.erase(std::remove(clean.begin(), clean.end(), '\t'), clean.end());
					
        			try{
		
						std::string::size_type sz;
			 
						const double d = std::stod(clean, &sz);
			 
						 if (sz == clean.size() )
							val = d;
			
					} catch (...){
		
						
					}
        		}
        		
        		k++;
        	}
        	
        	if (first){
        	
        		first = false;
        		continue;
        	}
        	
        	if (UtilString::equals(cmd, "undefined"))
        		;
        	else
        		values[cmd] = val;
    	}
	}
	
	
	inline const string& getCommand() const{
	
		return command;
	}
	
	inline const string toString() const{
	
		stringstream s;
		
		s << getCommand() << " ";
		
		for (auto it = values.begin(); it != values.end(); it++){
		
			s << "--" << it->first << "=" << it->second << " ";
		}
		
		return s.str();
	}
	
	friend inline std::ostream& operator<<(std::ostream &os, const Params& m) { 
    	
    	return os << m.toString();
	}
	
	inline void print() const{
	
		cout << toString() << endl;
	}
	
	inline void set(const string s, double v){
	
		values[s] = v;
	}
	
	inline double get(const string s){
	
		if (values.count(s) == 0)
			return 0.0;
		
		return values[s];
	}
	
};

class DataSummary{
public:
	
	virtual ~DataSummary(){};
	
	inline virtual double boxCount(const Point& qmin, const Point& qmax, QueryMode mode=QueryMode::EST) const = 0;
	//inline virtual void summarize(DataSet& data, Params& params) = 0;
	
protected:
    DataSummary(){}
    DataSummary(const DataSummary&){}
    DataSummary& operator=(const DataSummary&){ return *this; }
};



class UnifOneDimEstimator : public OneDimEstimator{
public:
	inline double intervalCount(double qmin, double qmax, QueryMode mode=QueryMode::EST) const override{
	
		return qmax-qmin;
	}
	
	~UnifOneDimEstimator() override{}
};


inline void evaluate(DataSet& data, DataSummary& summary, RangeQueries& queries, int num=1000000000, double scale=1.0){
	
	Timer timer("runQueries");
	timer.start(0, queries.ranges.size()/10, queries.ranges.size() );
	
	Stats ests("est");
	
	Stats counts("tru");
	Stats relErr("err", "%");
	Stats relWidth("wid", "%");

	int c = 0;
	for (auto it = queries.ranges.begin(); it != queries.ranges.end(); it++){
	
		const double t = it->count;
	
		const double lb = summary.boxCount( it->min, it->max, QueryMode::LB )*scale;
		const double est = summary.boxCount( it->min, it->max, QueryMode::EST )*scale;
		const double ub = summary.boxCount( it->min, it->max, QueryMode::UB )*scale;
		
		if ( timer.tick() || (ub < t) || (lb > t))
			cout << "LB " << lb << " EST " << est << " UB " << ub << " T " << t << " qmin " << it->min << " qmax " << it->max << endl;
		
		//if ( (ub < t) || (lb > t) ){
		//	cout << data.boxCount(it->min, it->max) << endl;
		//}
		
		assert (lb <= t);
		assert (ub >= t);
		
		ests.add(est);
		counts.add(t);
	
		relErr.add( 100*UtilMath::absVal(est-t)/t );
		relWidth.add( 100*UtilMath::absVal(ub-lb)/t );
		
		c++;
		
		if (c >= num)
			break;
	}
	
	timer.end();
	cout << "/runQueries" << endl;
	cout << relErr << endl;
	cout << relWidth << endl;
	cout << ests << endl;
	cout << counts << endl;

}

inline void evaluate(shared_ptr<DataSet> ds, DataSummary& summary, RangeQueries& queries, int num=1000000000, double scale=1.0){

	evaluate( (*ds), summary, queries, num, scale);
}

#endif