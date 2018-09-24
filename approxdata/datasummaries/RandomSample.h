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

#ifndef SHEKELYAN_RANDOM_SAMPLE_H
#define SHEKELYAN_RANDOM_SAMPLE_H

#include <approxdata/datasummaries/datasummaries.h>

class RandomSample : public DataSummary{
	
	
private:

	shared_ptr<DataSet> sample;
	long dataSize;
	long sampleSize;
	
	double samplePointWeight;
	
public:
	
	const int DIMS;
	
	Params params;
	
	inline double getSamplePointWeight() const{
	
		return samplePointWeight;
	}
	
	inline RandomSample(DataSet& data, Params& params) : params(params), DIMS(data.getDims() ){
		
		summarize(data, params);
	}
	
	inline RandomSample(DataSet& data, string paramStr) : params(paramStr), DIMS(data.getDims() ){
		
		summarize(data, params);
	}
	
	inline RandomSample(shared_ptr<DataSet> data, string paramStr) : RandomSample( (*data), paramStr){
		
	}
	
	inline RandomSample(unique_ptr<DataSet> data, string paramStr) : RandomSample( (*data), paramStr){
		
	}
	
	~RandomSample() override {};
	
	inline void summarize(DataSet& data, Params& params){
		
		const int DIMS = data.getDims();
		
		dataSize = data.getSize();
		
		sampleSize = UtilMath::minVal<long>(dataSize, params.get("kb")/( DIMS*8.0/(1024) ) );
		
		cout << "sampleSize " << sampleSize << endl;
		
		sample = data.getRandomSubset(sampleSize, "abc");
		samplePointWeight = 1.0/sampleSize*data.getSize();
		sample->setVerbose(false);
	}	
	
	inline double boxCount(const Point& qmin, const Point& qmax, QueryMode mode = QueryMode::EST) const override{
		
		if (mode == QueryMode::LB)
			return sample->boxCount(qmin, qmax);
			
		if (mode == QueryMode::UB){
		
			const long samplePointsInside = sample->boxCount(qmin, qmax);
			const long samplePointsOutside = sampleSize-samplePointsInside;
		
			return dataSize-samplePointsOutside;
		}
		
		if (dataSize == 0)
			return 0;
		
		if (sampleSize == 0)
			return 0;
		
		return sample->boxCount(qmin, qmax)*samplePointWeight;
	}
};

#endif