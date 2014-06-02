#ifndef PARSE_INPUT
#define PARSE_INPUT

using namespace std;

#include <iostream>
using std::cout;
using std::endl;

#include <fstream>
using std::ifstream;

#include <cstring>

const int MAX_CHARS_PER_LINE = 512;
const int MAX_TOKENS_PER_LINE = 20;
const char* const DELIMITER = " :=";

const float EPS = 0.00001;

/***********************************/
/********** READ DATA SET **********/
/***********************************/

bool is_pow2(int atr_val) {
    int x = 1;

    for(int i = 0; i < 31; i++) {
        if(x == atr_val) return true;
        x = (x << 1);
    }
    return false;
}
#if 0
static bool is_pow2(unsigned int x) {
    while(x > 1) {
        if( x & 1 ) return false;
        x = x >> 1;
    }
    return true;
}
#endif

inline int parseIntAttribute(ifstream& fin, char* buf, const char* atr) {
    int ret_val;

    assert( (!fin.eof()) && "incomplete input.data file!" );
    fin.getline(buf, MAX_CHARS_PER_LINE);

    { // reading the atribute name
        char* token_atr = strtok(buf, DELIMITER);
        bool ok = (token_atr != NULL) && (strcmp(token_atr, atr) == 0);
        assert( ok && "illegal/missing attribute name in input.data file!");
        // the result of strtok is part of buf, 
        // i.e., no dealocation is necessary!
    }

    { // reading the attribute value
        char* token_val = strtok(NULL, DELIMITER); 
        assert( token_val && "missing attribute value in input.data file!");
        ret_val = atoi(token_val);
    }

    return ret_val;
}

void readDataSet(   unsigned int& outer, 
                    unsigned int& num_X,
                    unsigned int& num_Y,
                    unsigned int& num_T
) {
    char buf[MAX_CHARS_PER_LINE];
    int atr_val; bool atr_ok;

    // create a file-reading object
    ifstream fin;
    fin.open("input.data");
    assert(fin.good() && "Error opening input.data file!");

    { // seeting and checking legality of the outer loop count
        atr_val = parseIntAttribute(fin, buf, "OUTER");
        atr_ok  = atr_val > 0;
        assert(atr_ok && "outer loop count less than 0!");
        outer   = static_cast<unsigned int>(atr_val); 
    }

    { // seeting and checking legality of NUM_X
        atr_val = parseIntAttribute(fin, buf, "NUM_X");
        atr_ok  = (atr_val > 0) && (atr_val <= 512) && is_pow2(atr_val); 
        assert(atr_ok && "illegal NUM_X value!");
        num_X   = static_cast<unsigned int>(atr_val); 
    }

    { // seeting and checking legality of NUM_Y
        atr_val = parseIntAttribute(fin, buf, "NUM_Y");
        atr_ok  = (atr_val > 0) && (atr_val <= 512) && is_pow2(atr_val); 
        assert(atr_ok && "illegal NUM_Y value!");
        num_Y   = static_cast<unsigned int>(atr_val); 
    }

    { // seeting and checking legality of NUM_T
        atr_val = parseIntAttribute(fin, buf, "NUM_T");
        atr_ok  = atr_val > 0; 
        assert(atr_ok && "NUM_T value less or equal to zero!");
        num_T   = static_cast<unsigned int>(atr_val); 
    }

    fin.close();

    cout<<"Parameters: { OUTER = "<<outer<<" }, { NUM_X = "<<num_X<<" }, "
        <<"{ NUM_Y = "<<num_Y<<" }, NUM_T = { "<<num_T<<" } "<<endl;
}

bool validate( const REAL* res, const unsigned int N ) {
    bool  is_valid = true;
    float token_val;
    char buf[MAX_CHARS_PER_LINE];
    ifstream fin;
    fin.open("output.data");
    assert(fin.good() && "Error opening input.data file!");

    unsigned int i = 0;
    while( i < N && (!fin.eof()) ) {
        fin.getline(buf, MAX_CHARS_PER_LINE);
        char* token_atr = strtok(buf, DELIMITER);
        token_val = atof(token_atr);
        float err = fabs(token_val - res[i]);
        if ( err > EPS ) {
            is_valid = false;
            cout<<"INVALID: error["<<i<<"] = "<<err<<" > "<<EPS<<" ! "<<res[i]<<endl;
        }
        i++;
    }

    assert(i==N && "size of the legal output does not match size of result!");
    fin.close();
    return is_valid;
}

void writeResult( const REAL* res, const unsigned int N ) {
    //ofstream fout;
    //fout.open ("output.data");

    for( int k=0; k<N; ++k ) {
        cout<<res[k]<<endl;
    }

    //fout.close();
}

#endif // PARSE_INPUT

