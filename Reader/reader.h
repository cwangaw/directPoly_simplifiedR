#ifndef READER_INCLUDED
#define READER_INCLUDED

#ifndef NULL
# define NULL 0
#endif

#include <string>
#include <vector>
#include "dimen.h"

/* DEFINE STANDARD SIZES */
#define SIZE_LINE            121
#define SIZE_DIRECTORYNAME   121
#define SIZE_WORD             31

/* DEFINE ERROR AND RETURN CODES */
#define ERR_ON_LINE                     1
#define ERR_BEFORE_LINE                 2
#define ERR_MEMORY                      3
#define ERR_OPENING_FILE                4
#define ERR_NESTING_FILE                5
#define ERR_EOF                         6
#define ERR_UNKNOWN_CMD                 7
#define ERR_UNMATCHED_END_COMMENT_CMD   8
#define ERR_IN_CODE                     9
#define ERR_UNSUPPORTED_OPTION         10
#define ERR_IN_DATA_BLOCK              11
#define ERR_IN_GRID_ARRAY              12
#define ERR_KEY_NOT_READ               13
#define ERR_BEGIN_LIST_NOT_READ        14
#define ERR_END_LIST_NOT_READ          15
#define ERR_SECTION_READ               16
#define ERR_SUBSECTION_READ            17
#define ERR_SUBSUBSECTION_READ         18
#define ERR_SECTION_NOT_READ           19
#define ERR_SUBSECTION_NOT_READ        20
#define ERR_SUBSUBSECTION_NOT_READ     21
#define ERR_SECTION_KEY_NOT_READ       22
#define ERR_SUBSECTION_KEY_NOT_READ    23
#define ERR_SUBSUBSECTION_KEY_NOT_READ 24
#define READ_LITERAL_CMD               25
#define ERR_FILESYSTEM                 26

#define NO_UNITS                   100
#define ERR_UNITS_NEGATIVE         101
#define ERR_NONNESTED_PARENTHESES  102
#define ERR_ARITHMETIC_OP          103
#define ERR_UNITS_ADDITION_FOUND   104
#define ERR_UNITS_WRONG_DIMEN      105
#define ERR_UNITS_TOO_LONG         106
#define ERR_UNITS_IMPROPER_EXPR    107
#define ERR_UNKNOWN_UNIT           108
#define ERR_UNITS                  109

#define ERR_BAD_DATA               200
#define ERR_OUT_OF_RANGE           201

/* DECLARE USER FUNCTIONS */

int openReader(const char* inFileName0, const char* echoFileName0,
	       int length_fileName, int batch0, int echo0);
int openReader(const std::string& inFileName0, const std::string& echoFileName0,
	       int length_fileName, int batch0, int echo0);
void closeReader(int error);
int processReaderError(const int error);
int screenIn();

int readKey(char* key, int skip_echoWord);
int readKeys(char** keys, int nKeys, int& keyRead, int skip_echoWord);

int readSection(int skip);
int readSubSection(int skip);
int readSubSubSection(int skip);

int skipToSection(char* key);
int skipToSubSection(char* key);
int skipToSubSubSection(char* key);
int skipToSections(char** keys, int nKeys, int& keyRead);
int skipToSubSections(char** keys, int nKeys, int& keyRead);
int skipToSubSubSections(char** keys, int nKeys, int& keyRead);

int readScalar(int& var);
int readScalar(double& var);
int readScalar(double& var, int* dimen);
int readScalar(double& var, int* dimen, double& unitVal);
int readTemp(double& var);
int readScalar(float& var);
int readScalar(float& var, int* dimen);
int readScalar(float& var, int* dimen, double& unitVal);
int readTemp(float& var);

inline int readScalar(double& var, int l, int m, int t) {
  int dimen[NUMBER_DIMEN] = {l,m,t};
  return readScalar(var, dimen);
}
inline int readScalar(float& var, int l, int m, int t) {
  int dimen[NUMBER_DIMEN] = {l,m,t};
  return readScalar(var, dimen);
}

int readString(std::string& str);
int readString(char* str, int len);
int readLine(char* str, int len);

int readVector(std::vector<int>& array);
int readVector(std::vector<float>& array);
int readVector(std::vector<double>& array);
int readVector(std::vector<float>& array, int* dimen, double& unitVal);
int readVector(std::vector<double>& array, int* dimen, double& unitVal);

inline int readVector(std::vector<double>& array,
		     int l, int m, int t, double& unitVal) {
  int dimen[NUMBER_DIMEN] = {l,m,t};
  return readVector(array, dimen, unitVal);
}
inline int readVector(std::vector<float>& array,
		     int l, int m, int t, double& unitVal) {
  int dimen[NUMBER_DIMEN] = {l,m,t};
  return readVector(array, dimen, unitVal);
}

int readArray(int* array, int len);
int readArray(double* array, int len);
int readArray(float* array, int len);
int readArray(double* array, int len, int* dimen, double& unitVal);
int readArray(float* array, int len, int* dimen, double& unitVal);

inline int readArray(double* array, int len,
		     int l, int m, int t, double& unitVal) {
  int dimen[NUMBER_DIMEN] = {l,m,t};
  return readArray(array, len, dimen, unitVal);
}
inline int readArray(float* array, int len,
		     int l, int m, int t, double& unitVal) {
  int dimen[NUMBER_DIMEN] = {l,m,t};
  return readArray(array, len, dimen, unitVal);
}

/* DECLARE SPECIALIZED ROUTINES */

int readBlock(int* array, int len, int& nRead, int nToReadTotal);
int readBlock(double* array, int len, int& nRead, int nToReadTotal);
int readBlock(double* array, int len, int& nRead,
	      int nToReadTotal, double units);
int readBlock(float* array, int len, int& nRead, int nToReadTotal);
int readBlock(float* array, int len, int& nRead,
	      int nToReadTotal, double units);

int readBlock(std::vector<int>& array, int& nRead);
int readBlock(std::vector<double>& array, int& nRead);
int readBlock(std::vector<double>& array, int& nRead, double units);
int readBlock(std::vector<float>& array, int& nRead);
int readBlock(std::vector<float>& array, int& nRead, double units);

int readPWS(int stopEOL=0);
int skipWS(int stop_after_EOL=0, int stop_after_EOF=0);
int skipLine();
void skipWord(int stopNonPWS, int stopListSymbols, int echoWord);

char getc();
void ungetc(char c);

void reader_echo(char* str);
void reader_echo(char c);
void reader_echo(int i);
void reader_echo(double x);
void reader_echo(float x);

int readUnits(double& value, int* dimen);
int readTempUnits(double& temp, int echoConversion);

int beginList();
int endList();
int skipToEndList();

void beginSkip();
void endSkip();

#endif
