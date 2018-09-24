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

#include <approxdata/datasummaries/datasummaries.h>

inline shared_ptr<DataSet> getData(const string s){

	shared_ptr<DataSet> ret;
	
	if (UtilString::startsWith(s, "synth")){
	
		Params p(s);
		
		assert ( UtilString::equals(p.get("type"), "zipf"));
		
		ret.reset( new ZipfDataSet(p.get("dims"), p.get("size"), p.get("clusters")) );
		
	} else {
		
		ret.reset( new FileDataSet(s) );
	}
	
	return ret;
}

int main (int argc, char **argv){
	
	stringstream s;
	
	s << "{";
	
	bool open = false;
	
	for (int i = 1; i < argc; i++){
		
		if (argv[i][0] == '-' && argv[i][1] == '-'){
		
			stringstream s1;
			stringstream s2;
		
			int j = 2;
		
			for (; (argv[i][j] != '=') && (argv[i][j] != '\0'); j++){
			
				s1 << argv[i][j];
			}
		
			if (argv[i][j] == '=')
				j++;
		
			for (; (argv[i][j] != '\0'); j++){
			
				s2 << argv[i][j];
			}
			
			if (i != 1)
				s << ",";
			
			s << "\"" << s1.str() << "\":";
		
			try {
				json j = json::parse(s2.str());
				s << s2.str();
				
			} catch(...){
				
				s << "\"" << s2.str() << "\"";
			}
			
			continue;
		}
		
		if (argv[i][0] == '-'){
			
			if (open)
				s << "[]";
			
			if (i != 1)
				s << ",";
				
			stringstream s1;
			
			for (int j = 1; (argv[i][j] != '\0'); j++){
			
				s1 << argv[i][j];
			}
			
			s << "\"" << s1.str() << "\":";
			open = true;
			
		} else {
			
			try {
				json j = json::parse(argv[i]);
				s << argv[i];
			} catch(...){
				
				s << "\"" << argv[i] << "\"";
			}
			
			open = false;
		}
	}
	
	s << "}";
	
	json j = json::parse(s.str() );
	
	cout << "command line interpreted as JSON: " << j << endl;
	
	if (!j["c"].is_null())
		j["convert"] = j["c"];
	
	if (!j["q"].is_null())
		j["querygen"] = j["q"];
	
	if (!j["i"].is_null())
		j["input"] = j["i"];
		
	if (!j["o"].is_null())
		j["output"] = j["o"];
	
	if (!j["o"].is_null())
		j["output"] = j["o"];
	
	if (!j["r"].is_null())
		j["runqueries"] = j["r"];
	
	if (!j["convert"].is_null()){
    
    	cout << "convert comma-separated-value (.csv) data to library data format (.json+.dat)" << endl;
    	
    	const string input = j["input"];
    	
    	CSVReader csv(input);
    	
    	const string dims = j["convert"];
    	
    	IntSeq seq(dims);
    	
    	csv.writeDataFiles(j["output"], seq.getSequence() );
    	
    	cout << "write meta data to " << j["output"] <<".json" << endl;
    	cout << "write binary data to " << j["output"] <<".dat" << endl;
    	
    	return 0;
    }
    
    if (!j["querygen"].is_null()){
    
    	const string input = j["input"];
    	const string output = j["output"];
    	const string paramStr = j["querygen"];
    	
    	
    	
    	
    	cout << "load data " << input << endl;
    	
    	shared_ptr<DataSet> data = getData(input);
    	
    	/*
    	try {
    	
    		json j = parse("{\"queries\":[
    		
    		
    	
    	
    	} catch(...) {
    	
    	
    	
    	}*/
    	
    	if (UtilString::endsWith(paramStr, ".json")){
    	
    		cout << "create queries " << paramStr << endl;
    	
    		RangeQueries q(paramStr);
    		
    		q.boxCount(data);
    		
    		RangeQueryFile f;
    	
			for (auto it = q.ranges.begin(); it != q.ranges.end(); it++)
				f.add( (*it) );
		
			cout << "write output to " << j["output"] << endl;
		
			f.write(output);
		
			return 0;
    	}
    	
    	
    	Params p(paramStr);
    	
    	cout << "create queries " << paramStr << endl;
    	
    	int type = -1;
    	
    	if ( UtilString::equals(p.getCommand(), "unbiased") )
    		type == 0;
    	
    	if ( UtilString::equals(p.getCommand(), "biased") )
    		type == 1;
    	
    	assert ( (type == 0) || (type == 1) );
    	
    	RangeQueries q( (*data), p.get("num"), p.get("selectivity")*0.5, p.get("selectivity")*1.5, type == 1 );
    	
    	RangeQueryFile f;
    	
    	for (auto it = q.ranges.begin(); it != q.ranges.end(); it++)
    		f.add( (*it) );
    	
    	cout << "write output to " << j["output"] << endl;
    	
    	f.write(output);
    	
    	return 0;
    }
    
    if (!j["summarize"].is_null()){
    	
    	const string input = j["input"];
    	
    	cout << "load data " << input << endl;
    	
    	shared_ptr<DataSet> data = getData(input);
    	
    	const string paramStr = j["summarize"];
    	
    	Params p(paramStr);
    	
    	if (UtilString::equals(p.getCommand(), "digithist")){
    	
    		cout << "create summary " << j["summarize"] << endl;
    		DigitHist dh((*data), paramStr);
    		
    		const string output = j["output"];
    		
    		cout << "write output to " << j["output"] << endl;
    		dh.writeBinary(output);
    		
    	} else {
    	
    		cout << "summary not supported! try --summary=digithist" << endl;
    	}
    	
    	return 0;
    }
    
    if (!j["runqueries"].is_null() ){
    
    	const string input = j["input"];
    	const string queries = j["runqueries"];
    
    	DigitHistReader r(input);
    	
    	RangeQueries rq(queries);
    	
    	RangeQueryFile f;
    	
    	cout << "run queries " << j["runqueries"] << endl;
    	
    	Stats selectivities("trueCounts");
	
		Stats errors("error");
		
		Stats addErrors("additiveError");
		Stats addBounds("additiveBounds");
	
		Stats mulErrors("multiplicativeError");
		Stats mulBounds("multiplicativeBounds");
	
		Stats relErrors("relativeError", "%");
		Stats relBounds("relativeBounds", "%");
    	
    	for (auto it = rq.ranges.begin(); it != rq.ranges.end(); it++){
    		
    		const double lowerBound = r.boxCount( it->min, it->max, QueryMode::LB);
    		const double estimate = r.boxCount( it->min, it->max, QueryMode::EST);
    		const double upperBound = r.boxCount( it->min, it->max, QueryMode::UB);
    		
    		const double trueValue = it->count;
		
			if (lowerBound > estimate)
				cout << "lowerBound " << lowerBound << " > " << estimate << " estimate" << endl;
		
			if (upperBound < estimate)
				cout << "upperBound " << upperBound << " < " << estimate << " estimate" << endl;
		
			if (lowerBound > trueValue)
				cout << "lowerBound " << lowerBound << " > " << trueValue << " trueValue" << endl;
		
			if (upperBound < trueValue)
				cout << "upperBound " << upperBound << " < " << trueValue << " trueValue" << endl;
    		
    		selectivities.add(trueValue);
		
			const double error = estimate-trueValue;
		
			errors.add(error);
		
			const double mulError = UtilMath::maxVal<double>(estimate/trueValue, trueValue/estimate);
			const double mulBound = upperBound/lowerBound;
		
			mulErrors.add(100.0*mulError );
			mulBounds.add(100.0*mulBound );
		
			const double addError = UtilMath::absVal<double>(estimate-trueValue);
			const double addBound = upperBound-lowerBound;
		
			addErrors.add(addError );
			addBounds.add(addBound );
		
			const double relError = addError/trueValue;
			const double relBound = addBound/trueValue;
		
			relErrors.add(100.0*relError );
			relBounds.add(100.0*relBound );
    		
    		f.add( (*it), lowerBound, estimate, upperBound);
    	}
    	
    	cout << selectivities << endl;
		cout << endl;
		cout << errors << endl;
		cout << endl;
		cout << mulErrors << endl;
		cout << mulBounds << endl;
		cout << endl;
		cout << addErrors << endl;
		cout << addBounds << endl;
		cout << endl;
		cout << relErrors << endl;
		cout << relBounds << endl;
    	
    	const string output = j["output"];
    	
    	f.write(output);
    	
    	cout << "write output to " << j["output"] << endl;
    	
    	return 0;
    }
	
}