#ifndef useoption_h
#define useoption_h
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

class useoption {
public :
	
	useoption();
	useoption(vector<string> options_);
	~useoption();

	void            set_option(vector<string> lines);
	string          get_first(string s); // get first option , to be option name
	string          get_second(string s); // get second option , to be value
	// input string
	bool            call_bool(  string s);
	double          call_float( string s);
	int             call_int(   string s);
	string          call_string(string s);
	vector<string>  call_vector_string(string s);
	// input char
	bool            call_bool  ( char* s);
	double          call_float ( char* s);
	int             call_int   ( char* s);
	string          call_string( char* s);

	void            show_options();
	void            remove_options(string s);
	void            add_option(    string s);

private :
	vector<string>  options ;
};

useoption::useoption(){
	vector<string> options_ ;
	options = options_ ;
}

useoption::useoption(vector<string> options_){
	options = options_ ;
}
useoption::~useoption()
{
	options.clear();
}
void useoption::set_option(vector<string> lines){
	options = lines ;
}
string useoption::get_first(string s){
	size_t index = s.find_first_of("#");
	if (index != string::npos) s.erase(index,(s.length()-index));
	index = s.find_first_not_of(" ");
	if (index != string::npos) s.erase(0,index);
	index = s.find_first_of(" ");
	if (index != string::npos) return s.substr(0,index);
	else return s ;
}

string useoption::get_second(string s){
	size_t index = s.find_first_of("#");
	if (index != string::npos) s.erase(index,(s.length())-1);
	index = s.find_last_not_of(" ");
	if (index != string::npos) s.erase(index+1,(s.length())-1);
	index = s.find_first_of(" ");
	while (index != string::npos){
		s.erase(0,index+1);
		index = s.find_first_of(" ");
	}
	return s;
}
void useoption::show_options(){
	size_t maxL = 0 ;
	for (size_t i = 0; i < options.size(); i++){
		size_t len = get_second(options.at(i)).length() ;
		if (maxL < len ) maxL = len ;
	}
	
		cout << setw(20) << string("Option Name")  ;
		cout << setw(maxL+1) << string("Option value") << endl ;
		cout << "====================" ;
		for (size_t i = 0; i <= maxL+1; i++) cout << "=" ;
		cout << endl ;
	for (unsigned int i = 0; i < options.size(); i++){
		cout << setw(20) << get_first(options.at(i))  ;
		cout << setw(maxL+1) << get_second(options.at(i)) << endl ;
	}
}
bool useoption::call_bool(string s)
{
	for (vector<string>::iterator it = options.begin() ; it != options.end(); it++){	
		if ( s == get_first(*it) ) {
			s = get_second(*it) ;
			if (s == string("True") || s == string("true") || s == string("1") || s == string("yes") 
			||  s == string("Yes")  || s == string("YES")  || s == string("Ok") ) return true ; 
			break ;
		}
	}
	return false ;
}
int useoption::call_int(string s)
{
	for (vector<string>::iterator it = options.begin() ; it != options.end(); it++){	
		if ( s == get_first(*it) ) {
			return atoi(get_second(*it).data()) ;
			break ;
		}
	}
	return 0 ;
}
double useoption::call_float(string s)
{
	for (vector<string>::iterator it = options.begin() ; it != options.end(); it++){	
		if ( s == get_first(*it) ) {
			return atof(get_second(*it).data()) ;
			break ;
		}
	}
	return 0 ;
}

string useoption::call_string(string s)
{
	for (vector<string>::iterator it = options.begin() ; it != options.end(); it++){	
		if ( s == get_first(*it) ) {
			return get_second(*it) ;
			break ;
		}
	}
	return string("None") ;
}

vector<string> useoption::call_vector_string(string s)
{
	vector<string> value ;
	for (vector<string>::iterator it = options.begin() ; it != options.end(); it++){	
		if ( s == get_first(*it) ) {
			value.push_back(get_second(*it)) ;
		}
	}
	return value ;
}

void useoption::add_option(string s)
{
	options.push_back(s);
}

void useoption::remove_options(string s)
{
	while (this->call_string(s) != string("None")){
//	cout << this->call_string(s) << endl ;
		for (unsigned i = 0; i < options.size() ; i++){
			if( get_first(options[i]) == s ){
				options.erase(options.begin()+i);
				break ;
			}
		}
	}
}

bool useoption::call_bool(char* s){
	return this->call_bool(string(s));
}
double useoption::call_float(char* s){
	return this->call_float(string(s));
}
int useoption::call_int(char* s){
	return this->call_int(string(s));
}
string useoption::call_string(char* s){
	return this->call_string(string(s));
}
#endif

//#ifdef useoption_cxx
//#endif // #ifdef useoption_cxx
