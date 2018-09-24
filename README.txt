----------------------------------------------------------------------------
HOW MATURE IS THIS LIBRARY?
----------------------------------------------------------------------------

This is a preliminary version, whose sole purpose is to give interested researchers an idea what the library will be able to offer in the near future.

It contains a full implementation of DigitHist (http://www.vldb.org/pvldb/vol10/p1514-shekelyan.pdf), but there are still some minor issues that can impact its performance. I think it should be ok to get a rough idea, but I would not recommend it yet for serious experiments. For more information please contact me <michael.shekelyan@gmail.com>.

----------------------------------------------------------------------------
LICENSING?
----------------------------------------------------------------------------

This library is licensed under the MIT License <http://opensource.org/licenses/MIT>.

ApproxData Library
Copyright (c) 2018 Michael Shekelyan <michael.shekelyan@gmail.com>

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

----------------------------------------------------------------------------
EXTERNAL LIBRARIES?
----------------------------------------------------------------------------

Yes, it contains the external library "JSON for Modern C++ 3.20", which is also licensed under the MIT license.
Copyright (c) 2013-2018 Niels Lohmann <http://nlohmann.me>.
https://github.com/nlohmann/json

----------------------------------------------------------------------------
HOW CAN I COMPILE IT?
----------------------------------------------------------------------------

- 1. Go to directory where "src" is a sub-folder.
- 2. Compile using any modern (c++11) compiler the file "src/approxdata.cpp" and create the binary "bin/approxdata"

----------------------------------------------------------------------------
HOW CAN I RUN A SIMPLE DEMO?
----------------------------------------------------------------------------

- 1. Follow steps for "HOW CAN I COMPILE IT?"
- 2. Create queries by running the binary following parameters (remove line breaks):

--input="synth --type=zipf --dims=3 --size=1000000 --clusters=1000" 
--querygen="biased --num=1000 --selectivity=1" 
--output="queries/demo.json"

( Explanation:
- The input is a 3-dimensional synthetic dataset with 1000000 points.
- The action is the generation of 1000 range queries with around 1% selectivity, that are biased towards data points. 
- The generated queries are stored in a ".json"-file.
)

3. Create summary by running the binary with following parameters (remove line breaks):

--input="synth --type=zipf --dims=3 --size=1000000 --clusters=1000" 
--summarize="digithist --kb=1000"
--output="summaries/demo.summary"

( Explanation:
- The input is a 3-dimensional synthetic dataset with 1000000 points, the same one as used for the query generation.
- The action is to create a 1000kb digithist summary.
- The created summary is stored as a binary file.
)

4. Run queries by running the binary with the following parameters (remove line breaks):

--input="summaries/demo.summary"
--runqueries="queries/demo.json"
--output="results/demo.json"

----------------------------------------------------------------------------
HOW CAN I USE MY OWN QUERIES?
----------------------------------------------------------------------------

Proceed as in "HOW CAN I RUN A SIMPLE DEMO?", but use"--querygen=example.json" to a file "example.json" that specifies the queries (using the same JSON format as otherwise created queries!). The json-file only needs to specify the minimas and maximas of the queries, it will then compute the counts for the queries.

----------------------------------------------------------------------------
HOW CAN I USE MY OWN CSV DATASET?
----------------------------------------------------------------------------

- 1. Follow steps for "HOW CAN I COMPILE IT?"
- 2. Convert CSV file by running the binary following parameters (remove line breaks):

--input="data/example.csv"
--convert="0,1"
--output="data/example"

( Explanation:
- The input is a comma-separated-value-file.
- The action is to convert the CSV file into a meta and data file, selecting only the first two dimensions.
- The created summary is stored in the files "data/example.json" and "data/example.dat"
)

3. Proceed as in "HOW CAN I RUN A SIMPLE DEMO?"

