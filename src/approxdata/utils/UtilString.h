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

#ifndef SHEKELYAN_UTILSTRING_H
#define SHEKELYAN_UTILSTRING_H

#include <approxdata/utils/utils.h>

namespace UtilString{

	inline const string twoDigit(long l){
	
		if (l >= 10)
			return to_string(l);
		else
			return "0"+to_string(l);
	}

	inline const string dateTimeString(long l, bool printOrig = false){
	
		long y, m, d, h, min, sec;
		
		UtilTime::dateTimeExtract(l, y, m, d, h, min, sec);
		
		return "["+twoDigit(h)+":"+twoDigit(min)+"]["+twoDigit(sec)+"s]";
	
		//return to_string(y)+"-"+twoDigit(m)+"-"+twoDigit(d)+"_"+twoDigit(h)+":"+twoDigit(min)+"::"+twoDigit(sec)+(printOrig ? (" ["+to_string(l)+"]") : "");
	}
	
	
	inline const string bitString(long l){
	
		char* c = new char[64];
	
		for (int i = 0; i < 64; i++){
		
			c[63-i] = (l & (1L << i)) != 0 ? '1' : '0';
		}
		
		for (int i = 0; i < 64; i++){
		
			if (c[i] == '1')
				break;
				
			c[i] = '0';
		}
		
		string s(c);
		
		return s;
	
	}

	inline bool split(const string& str, const string& beginStr, string& a, string& b){
							
		std::size_t x = str.find_first_of(beginStr);
		
		if (x == std::string::npos)
			return false;
		
		a = str.substr(0,x);
		b = str.substr(x+1, str.size()-x-1);
		
		return true;
	}
	
	inline bool equals(string a, string b){
		
		return a.compare(b) == 0;
	}
	
	inline const string beforeSplit(const string& str, const string& split){
							
		std::size_t x = str.find_first_of(split);
		
		if (x == std::string::npos)
			return "";
		
		return str.substr(0,x);
	}
	
	//inline bool startsWith(const string& str, const string& split){
	//						
	//	return str.find_first_of(split) == 0;
	//}
	
	inline const string afterSplit(const string& str, const string& split){
							
		std::size_t x = str.find_first_of(split);
		
		if (x == std::string::npos)
			return str;
		
		return str.substr(x+1, str.size()-x-1);
	}
	
	inline const string afterLastSplit(const string& str, const string& split){
							
		std::size_t x = str.find_last_of(split);
		
		if (x == std::string::npos)
			return str;
		
		return str.substr(x+1, str.size()-x-1);
	}
	
	inline double strToDouble(const string& valstr){
	
		try{
		
			std::string::size_type sz;
			 
			const double d = std::stod(valstr, &sz);
			 
			 if (sz == valstr.size() )
				return d;
			
		} catch (...){
		
			
		}
		
		return std::numeric_limits<double>::infinity();
	}
	
	inline int strToInt(const string& valstr){
	
		try{
		
			std::string::size_type sz;
			 
			const int d = std::stoi(valstr, &sz);
			 
			 if (sz == valstr.size() )
				return d;
			
		} catch (...){
		
			
		}
		
		return std::numeric_limits<int>::min();
	}
	
	template <typename E>
	inline void getVal(const string& s, E& val){
	
		istringstream numStream(s);
			
		if (!(numStream >> val))
			assert(false);
	}
	
	
	template <typename E>
	inline bool splitGetVal(const string& str, const string& beginStr, E& a, E& b){
							
		std::size_t x = str.find_first_of(beginStr);
		
		if (x == std::string::npos)
			return false;
		
		string as = str.substr(0,x);
		string bs = str.substr(x+1, str.size()-x-1);
		
		getVal(as,a); 
		getVal(bs,b);
		
		return true;
	}
	
	
	inline bool split(const string& str,  const string& beginStr, const string endStr, string& a, string& b, string& c){
							
		std::size_t x = str.find_first_of(beginStr);
		
		if (x == std::string::npos)
			return false;
		
		std::size_t y = str.find_first_of(endStr);
		
		if (y == std::string::npos)
			return false;
		
		a = str.substr(0,x);
		b = str.substr(x+1, y-x-1);
		c = str.substr(y+1,str.size()-y-1);
		
		return true;
	}
	
	
	inline const string digitToStr(int i){
	
		if (i == 0)
			return "zero";
		
		if (i == 1)
			return "one";
	
		if (i == 2)
			return "two";
	
		if (i == 3)
			return "three";
			
		if (i == 4)
			return "four";
		
		if (i == 5)
			return "five";
	
		if (i == 6)
			return "six";
	
		if (i == 7)
			return "seven";
			
		if (i == 8)
			return "eight";
		
		if (i == 9)
			return "nine";
	
		if (i == 10)
			return "ten";
	
		if (i == 11)
			return "eleven";
			
		if (i == 12)
			return "twelve";
			
		if (i == 13)
			return "thirteen";
			
		if (i == 14)
			return "fourteen";
		
		if (i == 15)
			return "fifteen";
		
		if (i == 16)
			return "sixteen";
		
		if (i == 17)
			return "seventeen";
			
		if (i == 18)
			return "eighteen";
			
		if (i == 19)
			return "nineteen";
		
		if (i == 20)
			return "twenty";
		
		return "overtwenty";
	}
	
	
	inline bool endsWith(const string& value, const string& ending){
		
    	if (ending.size() > value.size())
    		return false;
    	
    	return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
    }
	
	
	inline bool startsWith(const string& value, const string& prefix){
	
		return value.substr(0, prefix.size()).compare(prefix) == 0;
    }
	
	
	template <typename E>
	inline string getDigitRepresentation(E n, int radix, int digitNum){
	
		stringstream s;
		
		for (int i = digitNum-1; i >= 0; --i){
		
			long coef = pow( radix, i);
		
			int digit = n /coef;
			
			s << digit;
			
			n -= digit*coef;
		}
		
		return s.str();	
	}
	
	template <typename T>
	inline const string listToString(const vector<T>& v){
	
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
	
	inline const string repeat(const string& s, int n){
	
		stringstream ss;
		
		for (int i = 0; i < n; i++){
		
			ss << s;
		}
		
		return ss.str();
	}
	
	template <typename T>
	inline const string horizontalLatexTable(const vector<T>& v){
	
		stringstream ss;
		
		ss << "\\begin{tabular}{|";
		
		for (int i = 0; i < v.size(); i++)
			ss << "c";
		
		ss << "|}\n";
		
		ss << "\\hline" << endl;
		
		for (int i = 0; i < v.size(); i++){
		
			if (i > 0)
				ss << " & ";
			
			ss << v.at(i);
		}
		ss << "\\\\" << endl;
		ss << "\\hline" << endl;
		
		ss << "\\end{tabular}" << endl;
	
		return ss.str();
	}
	
	/*
	inline vector<string>* newStringVector(string v1, string v2){
	
		vector<string>* ret = new vector<string>();
		
		ret->push_back(v1);
		ret->push_back(v2);
		
		return ret;
	}*/
	
	template <typename T>
	inline const string latexTable(const vector<vector<T>>& vs, const vector<string>& rownames, const string betweenRows=""){
	
		stringstream ss;
		
		ss << "\\begin{tabular}{|";
		
		ss << "l|";
		
		for (int i = 0; i < vs.at(0).size(); i++)
			ss << "c|";
		
		ss << "}\n";
		
		ss << "\\hline" << endl;
		
		for (int j = 0; j < vs.size(); j++){
		
			ss << rownames.at(j) << " & ";
			
			vector<T>& v = vs.at(j);
			
			for (int i = 0; i < v.size(); i++){
			
				if (i > 0)
					ss << " & ";
			
				ss << v.at(i);
			}
			ss << "\\\\" << endl;
			ss << "\\hline" << endl;
			ss << betweenRows << endl;
		}
		
		ss << "\\end{tabular}" << endl;
	
		return ss.str();
	}
	
	template <typename T>
	inline string listToStringJoin(const vector<T>& v, const string join=","){
	
		stringstream s;
	
		for (int i = 0; i < v.size(); i++){
		
			if (i > 0)
				s << join;
			
			s << v.at(i);
		}
		
		return s.str();
	}
	
	template <typename T>
	inline string arrayToStringJoin(int size, T* v, const string join=","){
	
		if (v == NULL)
			return "NULL";
		
		stringstream s;
	
		for (int i = 0; i < size; i++){
		
			if (i > 0)
				s << join;
			
			s << v[i];
		}
		
		return s.str();
	}
	
	inline void dateTimePrint(bool printOrig = true){
	
		cout << dateTimeString( UtilTime::dateTimeNow(), printOrig) << endl;
	}
	
	inline string dateTimeNowString(bool printOrig = true){
	
		return dateTimeString( UtilTime::dateTimeNow(), printOrig);
	}
	
	inline bool contains(const string& s1, const string& s2){
	
		return s1.find(s2) != string::npos;
	}
	
	/*
	inline bool containsAny(const string& s1, const vector<const string>& s2){
	
		for (auto it = s2.begin(); it != s2.end(); it++){
			
			if (contains(s1, (*it)))
				return true;
		}
		
		return false;
	}*/
	
	/*
	inline bool containsAll(const string s1, const vector<const string>& s2){
	
		for (auto it = s2.begin(); it != s2.end(); ++it){
		
			if (!contains(s1, (*it)))
				return false;
		}
		
		return true;
	}*/
	
	template <typename T>
	inline void printList(vector<T>& v){
	
		for (int i = 0; i < v.size(); i++){
		
			if (i > 0)
				cout << ", ";
				
			cout << v.at(i);
		}
		
		cout << endl;
	}
	
	
	/*
	inline void printList(Point& v){
	
		for (int i = 0; i < v.size(); i++){
		
			if (i > 0)
				cout << ", ";
				
			cout << v.at(i);
		}
		
		cout << endl;
	}*/
	
	
	template <typename T>
	inline void printList(shared_ptr<vector<T>> v){
	
		printList( (*v) );
	}
	
	
	
	inline void printException(std::exception& exception){
	
		cout << "exception " << exception.what() << endl;
		cerr << "exception " << exception.what() << endl;
	}
	template <typename T>
	inline void printArray(int len, T* v){
	
		for (int i = 0; i < len; i++){
		
			if (i > 0)
				cout << ", ";
			
			cout << v[i];
		}
		
		cout << endl;
	}
	
	
	inline void splitFirst(const string source, char split, string& first, string& rest){
	
		stringstream lineStream(source);
		
		string s;
			
		first = "";
		rest = "";
			
		if (!getline(lineStream, s, split))
			return;
		
		first = s;
		
		std::stringstream ss;	
		
		for (bool first = true; getline(lineStream, s, split); first = false){
		
			if(!first)
				ss << split;
			ss << s;
		}
		
		rest = ss.str();
	}
	
	
	inline const string replaceAll(const string& str, const string& from, const string& to) {
		
		string s = str;
		
		size_t start = 0;
		
		while((start = s.find(from, start)) != std::string::npos) {
			
			s.replace(start, from.length(), to);
			
			start += to.length();
		}
		
		return s;
	}

	inline const string numToStr(double d){
	
		if (d < 0)
			return "-"+numToStr(-d);
	
		if (d != d)
			return " ";
	
		if (d == ((double)((int)d)))
			return to_string((int) d);
		
		return replaceAll(to_string(d), ".", ",");
	}
	
	
	
	inline const string numToStr(double value, int decimals){
		
		std::ostringstream ss;
		ss << std::fixed << std::setprecision(decimals) << value;
		std::string s = ss.str();
		if(decimals > 0 && s[s.find_last_not_of('0')] == '.') {
			s.erase(s.size() - decimals + 1);
		}
		return s;
	}
	
	inline const string doubleToString(double value, int decimals=-1){
	
		
		if (value == floor(value))
			return to_string( ((long) value));
			
		std::ostringstream ss;
	
		if (decimals == -1)
			ss << std::fixed << std::setprecision(100) << value;
		else
			ss << numToStr(value, decimals);
		
		return ss.str();
	}
	
	
	
	
	inline void print(const string str, const vector<int>& coords){
	
		cout << str;
		
		for (int i = 0; i < coords.size(); i++)
			cout << " " << coords.at(i);
			
		cout << endl << std::flush;
	}
	
	
	
	inline void printErr(const string str, const vector<int>& coords){
	
		std::cerr << str;
	
		for (int i = 0; i < coords.size(); i++)
			std::cerr << " " << coords.at(i);
		
		std::cerr << endl << std::flush;
	}
	
	inline void print(const string str, const vector<double>& coords){
	
		cout << str;
	
		for (int i = 0; i < coords.size(); i++)
			cout << " " << coords.at(i);
		
		cout << endl << std::flush;
	}
	
	inline void print(std::string str, double d){
	
		cout << str << " " << d << endl << std::flush;
	}
	
	inline void printErr(std::string str, int dims, vector<double>& coords){
	
		std::cerr << str;
	
		for (int i = 0; i < coords.size(); i++)
			std::cerr << " " << coords.at(i);
			
		std::cerr << endl << std::flush;
	}
	
	
	inline const string msToString(double ms){
	
		if (ms <= 0)
			return "0 ms";
	
		stringstream s;
	
		bool first = true;
		
		for (int i = 0; i <= 4; i++){
		
			long mul = 0;
			
			if (i == 0)
				mul = 1000*60*60*24;
			if (i == 1)
				mul = 1000*60*60;
			if (i == 2)
				mul = 1000*60;
			if (i == 3)
				mul = 1000;
			if (i == 4)
				mul = 0;
			
			double val = UtilMath::makeMultipleOf(ms, mul);
			
			if (i == 3 && (val > 0 || (!first) ) ){
			
				if (!first)
					s << " ";
				
				s << UtilString::doubleToString(ms/mul, 3);
				
				s << " ";
				
				if (i == 0)
					s << "d";
				else if (i == 1)
					s << "hs";
				else if (i == 2)
					s << "mins";
				else if (i == 3)
					s << "secs";
				else if (i == 4)
					s << "ms";
				
				break;
			}
			
			if (val > 0){
			
				ms -= val;
			
				if (!first)
					s << " ";	
			
				first = false;
			
				s << UtilString::doubleToString(mul == 0 ? val : val/mul, 3);
				s << " ";
			
				if (i == 0)
					s << "d";
				else if (i == 1)
					s << "hs";
				else if (i == 2)
					s << "mins";
				else if (i == 3)
					s << "secs";
				else if (i == 4)
					s << "ms";
			}
		}
		
		return s.str();	
	}
	
	inline const string secsToString(long secs){
	
		return msToString(secs*1000);
	}
};


/*
class TextTable{

	public:

	vector<vector<const string>> rows;
	
	vector<int> columnMaxLengths;
	
	int colnum;
	const string sep = "  ";
	
	inline TextTable(int colnum) : columnMaxLengths(colnum, 0), colnum(colnum){
	
		
	}
	
	inline void addRow(vector<const string> cols){
	
		for (int i = 0; i < colnum; i++)
			columnMaxLengths.at(i) = UtilMath::maxVal<int>(columnMaxLengths.at(i), cols.at(i).size() );
			
		rows.push_back(cols);	
	}
	
	inline const string toString(){
	
		std::stringstream sbuf;
		
		for (int i = 0; i < rows.size(); i++){
		
			vector<const string>& row = rows.at(i);
		
			for (int j = 0; j < colnum; j++){
		
				const string s = row.at(j);
				
				sbuf << s;
				
				for (int k = s.size(); k < columnMaxLengths.at(j); k++)
					sbuf << ' ';
					
				sbuf << sep;
			}
			
			sbuf << '\n';
		}
		
		return sbuf.str();
	}
	
	
	
	
	

};*/



#endif