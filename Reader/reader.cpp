/******************************************************************************
Read Scalars (in serial)

This package reads scalar data.  It consists of the two files
  reader.c and reader.h .

Todd Arbogast
The University of Texas at Austin
August 1996
Modified July 1997, July 2005, June 2007, and June 2015

SYNTACTIC STRUCTURE OF GENERAL INPUT FILES

The general input file(s) contain:
  sectional commands
  keywords
  data items
  list delineators
  white space (WS), including proper white space (PWS)
  comments (explicit and implicit)
  commands
Within an input file, data are arranged in sections, subsections, and
  subsubsections and ordered item-by-item, possibly separated by
  list delineators, special commands, comments, and white-space.

Sections, Subsections, and Subsubsections are initiated by a
  SECTION, SUBSECTION, or SUBSUBSECTION command, respectively,
  and a specified string.  See below for command syntax.
Keywords are specified strings of characters.

Data items consist of:
  int
  float
  double
  word of characters (terminated by a SPACE, newline (EOL: NL or CR), or TAB)
  line of characters (including spaces and tabs, terminating before an EOL)
List delineators are:
  BEGIN_LIST '{' to indicate the beginning of a list;
  END_LIST '}' to indicate the end of a list.
    Lists are user specified numbers of items.  The deliminators check that
    the user has counted the number of list items correctly.

WS (white space) consists of SPACE ' ', TAB '\t',
  NL '\n' (end-of-line), CR '\r' (end-of-line),
  COMMA ',', SEMICOLON ';', and COLON ':'.  WS is ignored on input.
PWS (proper white space) consists of SPACE ' ', TAB '\t', NL '\n', CR '\r'.
  Command names and keywords are terminated by PWS.

Explicit comments appear
  After the COMMENT symbol '#' (the comment runs to an EOL symbol,
    and the '#' cannot be imbedded in a string).
  Between the BEGIN_COMMENT and END_COMMENT commands (see below).
Implicit comments can be allowed:
  Before the sectional commands.  This enables skipping of unnessessary
  sectional units.

Commands alter the way the data file(s) are read.  Commands take the form
  COMMAND symbol '$' followed immediately (no space)
    by a command name.  Note that the COMMAND symbol is ignored
    in any explicit comment (but not in implicit comments).
  Commands are
  $ is a null command
  $$ is a null command
  $INCLUDE, followed by PWS and a new file name.
    This causes subsequent reading to take place in the new file, until
    the EOF (end-of-file) is read.  Reading then returns to the original file.
    A finite levell of nesting is imposed.
  $BEGIN_COMMENT causes input to be ignored until the matching END_COMMENT
    command is read (such commands are nested) or the EOF is reached.
  $END_COMMENT terminates a BEGIN_COMMENT command.
  ${# is equivalent to BEGIN_COMMENT.
  $#} is equivalent to END_COMMENT.
  $LITERAL terminates reading WS and other symbols.  Used to allow special
    symbols in data strings (',', ';', '#', and '$').  (This command is
    not needed before INCLUDE file names.)
  $SECTION begins a section of the data file.  Reading may not pass a
    section command except by explicit action.    
  $S is equivalent to SECTION.
  $SUBSECTION begins a subsection of the data file.  Reading may not pass
    a section or subsection command except by explicit action.    
  $SS is equivalent to SUBSECTION.
  SUBSUBSECTION begins a subsubsection of the data file.  Reading
    may not pass any sectional command except by explicit action.    
  $SSS is equivalent to SUBSUBSECTION.
  $IGNORE_WORDS, followed by PWS and an integer.  Causes that number of
    words to be ignored from the input stream.  Note that the include
    command only is not ignored.  This command can be used to ignore
    preambles in included data files.
  $IGNORE_LINES, followed by PWS and an integer.  Causes that number of
    lines (i.e., up to an EOL) to be ignored from the input stream.  The line
    containing the integer is included in the count of lines ignored, and
    ignored from that point on.  Note that the include
    command only is not ignored.  This command can be used to ignore
    preambles in included data files.

REMARKS:
  (1) The commands readBlock<int>, readBlock<double[,units]>, and
      readBlock<float[,units]> read a block of data, and it may contain the
      '@' repetition character.  A set of "n" single data items, each with the
      value "value", may be indicated by "n@value".  PWS may precede "value",
      but only spaces may precede '@'.
  (2) Integers may not be followed by ".".
  (3) Arrays are read with readArray, using the syntax:
        1. constant <value> [<len>, for std::vector]
        2. linear <value1> <value2> [except int and std::vector]
        3. uniform <value1> <value2> (same as linear)
        4. { <value> ... } (no units allowed)
        5. { [<units>] <value> ... } (units apply to all numbers)
        6. ??? (return ERR_UNSUPPORTED_OPTION for possible further processing)
******************************************************************************/

#include <stdlib.h>
#include <cstring>
#include <iostream>
#include <fstream>
using namespace std;
#include "expr.h"
#include "reader.h"
#include "debug.h"

#define SIZE_INTBUF           4096
#define SIZE_DOUBLEBUF        4096
#define SIZE_FLOATBUF         4096
#define DEPTH_INCLUDES           5
#define MAX_LENGTH_UNITS_STR   121

/* DEFINE PROPER WHITESPACE */
#define NL        '\n'   /* end of line character */
#define CR        '\r'   /* end of line character */
#define SPACE      ' '   /* space character */
#define TAB        '\t'  /* tab character */

/* DEFINE IMPROPER WHITESPACE */
#define COMMENT    '#'   /* comment character */
#define COMMA      ','   /* comma separator */
#define SEMICOLON  ';'   /* semicolon separator */
#define COLON      ':'   /* colon separator */

/* DEFINE SPECIAL SYMBOLS */
#define REPETITION '@'
#define BEGIN_LIST '{'
#define END_LIST   '}'
#define BEGIN_UNIT '['
#define END_UNIT   ']'

/* DEFINE PROCESSOR COMMANDS */
#define COMMAND            '$' /* processor commands follow this symbol */
#define SIZE_COMMAND_NAME  21
static char beginCommentCmd[]       = "BEGIN_COMMENT";
static char endCommentCmd[]         = "END_COMMENT";
static char shortBeginCommentCmd[]  = "{#";
static char shortEndCommentCmd[]    = "#}";
static char includeCmd[]            = "INCLUDE";
static char literalCmd[]            = "LITERAL";
static char sectionCmd[]            = "SECTION";
static char shortSectionCmd[]       = "S";
static char subSectionCmd[]         = "SUBSECTION";
static char shortSubSectionCmd[]    = "SS";
static char subSubSectionCmd[]      = "SUBSUBSECTION";
static char shortSubSubSectionCmd[] = "SSS";
static char ignoreWordsCmd[]        = "IGNORE_WORDS";
static char ignoreLinesCmd[]        = "IGNORE_LINES";

static int inFileNo = -1;
static int SIZE_FILENAME;
static char* inFileName[DEPTH_INCLUDES];
static int inFileLineNo[DEPTH_INCLUDES];
static istream* inFile[DEPTH_INCLUDES];
static ifstream inFile_ifstream[DEPTH_INCLUDES];
#define instrm (*inFile[inFileNo])
static ostream* echoFile;
static ofstream echoFile_ofstream;

static int skipping = 0;
static int beginSkipFileLineNo;
static char* beginSkipFileName;

static int batch = 1;
static int echoing = 0;

static char* save_keyWord_1;
static char** save_keyWord;
static int save_nKeyWords;

/* STATIC FUNCTIONS */

static int beginFile();
static int finishFile();
static void getString(char* string, int len);
static int processCmd();
static int ignoreLW(int ignoreWords);
static void multByUnits(double* array, int n, double units);
static void multByUnits(float* array, int n, double units);

/******************************************************************************
Title: getc
Description: Get a character.  Essentially the usual getc C routine.
             This routine is dangerous, simce normal checking for
	     special characters is omitted.
*/
char getc() {
  char c;
  instrm.get(c);
  return c;
}

/******************************************************************************
Title: ungetc
Description: Unget a character.  Use only after getc, as in the usual
             getc and ungetc C routines.
*/
void ungetc(char c) {
  instrm.putback(c);
}

/******************************************************************************
Title: echo
Description: Echo a <type>.
*/
void reader_echo(const char* s) { if(echoing) *echoFile << s; }
void reader_echo(char c)        { if(echoing) *echoFile << c; }
void reader_echo(int i)         { if(echoing) *echoFile << i; }
void reader_echo(double x)      { if(echoing) *echoFile << x; }
void reader_echo(float x)       { if(echoing) *echoFile << x; }

/******************************************************************************
Title: beginSkip()
Description: Set the skip line number.
*/
void beginSkip() {
  skipping = 1;
  beginSkipFileLineNo = inFileLineNo[inFileNo];
  beginSkipFileName = inFileName[inFileNo];
}

/******************************************************************************
Title: endSkip()
Description: Set the skip line number.
*/
void endSkip() {
  skipping = 0;
}

/******************************************************************************
Title: readPWS
Description: Read PWS, up to non-punctuation white
             space or one less than the length of the WS, and add
	     the WS terminating character.
Input:  stopEOL  Switch to stop just after EOL is read.
Return: whether stopped at EOL.
*/
int readPWS(int stopEOL) {
  char c;

  do {
    c = getc();
    if(c == NL || c == CR) {
      inFileLineNo[inFileNo]++;
      if(stopEOL) return 1;
    }
  } while(c == SPACE || c == TAB || c == NL || c == CR);
  ungetc(c);
  return 0;
}

/******************************************************************************
Title: getString
Description: Read a string, up to proper WS
             (or one less than the length of the string), and add
	     the string terminating character.  Do not read initial
             proper WS.  Do not echo the string
Input:  len     The length of the string.
Output: string  The string read.
*/
static void getString(char* string, int len) {
  char c;
  int i=0;

  for(i=0; i<len-1; i++) {
    c = getc();
    if(c == NL || c == CR || instrm.eof() || c == SPACE || c == TAB ) {
      ungetc(c);
      break;
    }
    string[i] = c;
  }
  string[i] = '\0';
}

/******************************************************************************
Title: beginFile
Description: Begin reading a file.
Return: Error flag (see reader.h).
*/
static int beginFile() {
  if(inFileNo >= 0) {
    if(inFileNo >= DEPTH_INCLUDES-1) return ERR_NESTING_FILE;
    readPWS(0);
    getString(inFileName[inFileNo+1],SIZE_FILENAME);
  }
  if(!strcmp(inFileName[inFileNo+1],"")) {
    inFile[inFileNo+1] = &cin;
  } else {
    inFile_ifstream[inFileNo+1].open(inFileName[inFileNo+1],ios::in);
    inFile[inFileNo+1] = &inFile_ifstream[inFileNo+1];
    if(!inFile_ifstream[inFileNo+1] || !inFile[inFileNo+1]->good()) return ERR_OPENING_FILE;
  }
  inFileNo++;
  inFileLineNo[inFileNo] = 1;
  
  //ungetc(NL);

  for(int i=0; i<inFileNo; i++) cout << " ";
  if(screenIn()) {
    cout << "> Reading from screen" << endl;
  } else {
    cout << "> Reading from file " << inFileName[inFileNo] << endl;
  }
  return 0;
}

/******************************************************************************
Title: finishFile
Description: Finish reading a file.
Return: Error flag (see reader.h).
*/
static int finishFile() {
  if(inFile[inFileNo] != &cin) inFile_ifstream[inFileNo].close();
  inFileNo--;

  for(int i=0; i<inFileNo; i++) cout << " ";

  if(inFileNo < 0) {
    cout << "< Finished reading" << endl;
    if(batch) {
      return ERR_EOF;
    } else {
      int error;
      inFileNo = -1;
      strcpy(inFileName[0],"");
      error = beginFile();
      return error;
    }
  }

  if(screenIn()) {
    cout << "< Reading from screen" << endl;
  } else {
    cout << "< Reading from file " << inFileName[inFileNo] << endl;
  }
  return 0;
}

/******************************************************************************
Title: processCmd
Description: Process the commands.  Flag as an error reading the
             SECTIONAL commands.
Return: Error or return flag (see reader.h).
*/
static int processCmd() {
  char c;
  int error = 0;
  char cmd[SIZE_COMMAND_NAME];

#define CMD_IS(cmdName) !strcmp(cmd,cmdName)

  getString(cmd,SIZE_COMMAND_NAME);
  if(CMD_IS(includeCmd)) {
    error = beginFile(); if(error) return error;
  } else if( CMD_IS(beginCommentCmd) || CMD_IS(shortBeginCommentCmd) ) {
    int nBeginComments = 1;
    while(1) {
      do {
	c = getc();
	if(c == NL || c == CR) inFileLineNo[inFileNo]++;
      } while(c != COMMAND && !instrm.eof());
      if(instrm.eof()) {
	ungetc(c);
	return 0;
      }
      getString(cmd,SIZE_COMMAND_NAME);
      if( CMD_IS(endCommentCmd) || CMD_IS(shortEndCommentCmd) ) {
	nBeginComments--;
	if(nBeginComments == 0) return 0;
      } else if( CMD_IS(beginCommentCmd)
		|| CMD_IS(shortBeginCommentCmd) ) nBeginComments++;
    }
  } else if( CMD_IS(endCommentCmd) || CMD_IS(shortEndCommentCmd) ) {
    return ERR_UNMATCHED_END_COMMENT_CMD;
  } else if(CMD_IS(literalCmd)) {
    do {
      c = getc();
      if(c == NL || c == CR) inFileLineNo[inFileNo]++;
    } while(c == SPACE || c == TAB || c == NL || c == CR);
    ungetc(c);
    if(echoing) *echoFile << COMMAND << literalCmd << " ";
    return READ_LITERAL_CMD;
  } else if( CMD_IS(sectionCmd) || CMD_IS(shortSectionCmd) ) {
    if(echoing) *echoFile << COMMAND << shortSectionCmd << " ";
    return ERR_SECTION_READ;
  } else if( CMD_IS(subSectionCmd) || CMD_IS(shortSubSectionCmd) ) {
    if(echoing) *echoFile << COMMAND << shortSubSectionCmd << " ";
    return ERR_SUBSECTION_READ;
  } else if( CMD_IS(subSubSectionCmd) || CMD_IS(shortSubSubSectionCmd) ) {
    if(echoing) *echoFile << COMMAND << shortSubSubSectionCmd << " ";
    return ERR_SUBSUBSECTION_READ;
  } else if( CMD_IS(ignoreWordsCmd) ) {
    error = ignoreLW(1); if(error) return error;
  } else if( CMD_IS(ignoreLinesCmd) ) {
    error = ignoreLW(0); if(error) return error;
  } else {
    return ERR_UNKNOWN_CMD;
  }
  return 0;
}

/******************************************************************************
Title: skipWS
Description: Skip past white space and comments.
             Processes commands along the way.
	     Conditionally stop just after (the last of several) EOL is read.
Input:  stopEOL  Switch to stop just after EOL is read.
        stopEOF  Switch to stop just after EOF is read.
Return: Error flag (see reader.h).
*/
int skipWS(int stopEOL, int stopEOF) {
  char c;
  int error = 0;

  while(1) {
    c = getc();
    if(instrm.eof()) {
      error = finishFile(); if(error) return error;
      if(stopEOF) return 0;
    } else {
      switch(c) {
      case NL:
      case CR:
	if(stopEOL) {
	  do {
	    inFileLineNo[inFileNo]++;
	    c = getc();
	    if(instrm.eof()) {
	      error = finishFile(); if(error) return error;
	      c = getc();
	    }
	  } while(c == NL || c == CR);
	  ungetc(c);
	  return 0;
	} else {
	  inFileLineNo[inFileNo]++;
	}
      case SPACE:
      case TAB:
      case COMMA:
      case SEMICOLON:
      case COLON:
	break;
      case COMMENT:
	do { c = getc(); } while(c != NL && c != CR && !instrm.eof());
	ungetc(c);
	break;
      case COMMAND:
	error = processCmd();
	if(error) {
	  if(error == READ_LITERAL_CMD) return 0;
	  return error;
	}
	break;
      default:
	ungetc(c);
	return 0;
      }
    }
  }
}

/******************************************************************************
Title: skipLine
Description: Skip to the end of a line.
             Do not process commands along the way!
Return: Error flag (see reader.h).
*/
int skipLine() {
  char c;
  int error = 0;

  while(1) {
    c = getc();
    switch(c) {
    case NL:
    case CR:
      inFileLineNo[inFileNo]++;
      return 0;
    default:
      if(instrm.eof()) {
	error = finishFile(); if(error) return error;
	return 0;
      }
    }
  }
}

/******************************************************************************
Title: skipWord
Description: Skip past the next word (of non special characters).
Input:  stopNonPWS       Switch to stop just before non PWS or
                         COMMENT or COMMAND symbols.
        stopListSymbols  Switch to stop just before a BEGIN_LIST
                         or END_LIST symbol.
        echoWord         Switch to echo the word (if echoing is enabled).
Return: Error flag (see reader.h).
*/
void skipWord(int stopNonPWS, int stopListSymbols, int echoWord) {
  char c;
  int echoed = 0;

  echoWord = echoWord && echoing;

  while(1) {
    c = getc();
    switch(c) {
    case SPACE:
    case TAB:
    case NL:
    case CR:
    case COMMA:
    case SEMICOLON:
    case COLON:
    case COMMENT:
    case COMMAND:
      if(stopNonPWS) {
	ungetc(c);
	goto exit;
      }
      break;
    case BEGIN_LIST:
    case END_LIST:
      if(stopListSymbols) {
	ungetc(c);
	goto exit;
      }
      break;
    default:
      if(instrm.eof()) {
	ungetc(c);
	goto exit;
      }
    }
    if(echoWord) {
      *echoFile << c;
      echoed = 1;
    }
  }
exit:
  if(echoed) *echoFile << endl;
}

/******************************************************************************
Title: ignoreLW
Description: Ignore lines or words.  Process the INCLUDE command, but
             no other commands.
Input:  ignoreWords      Switch to ignore words vs. lines.
Return: Error flag (see reader.h).
*/
static int ignoreLW(int ignoreWords) {
  int n;
  char c;
  int error = 0;
  char cmd[SIZE_COMMAND_NAME];

  int* nRead;
  int nWordsRead = 0;
  int nLinesRead = 0;

  readPWS(0);

  instrm >> n; if(!instrm.good()) return ERR_ON_LINE;
  c = getc();
  if(c == '.') return ERR_ON_LINE;
  ungetc(c);

  nRead = (ignoreWords) ? &nWordsRead : &nLinesRead;

  while(*nRead < n) {
    c = getc();
    switch(c) {
    case NL:
    case CR:
      inFileLineNo[inFileNo]++;
      nLinesRead++;
      break;
    case SPACE:
    case TAB:
      break;
    case COMMAND:
      getString(cmd,SIZE_COMMAND_NAME);
      if(!strcmp(cmd,includeCmd)) {
	int nOld = inFileLineNo[inFileNo];
	error = beginFile(); if(error) return error;
	nLinesRead += inFileLineNo[inFileNo-1] - nOld - 1;
      } else {
	nWordsRead++;
      }
      break;
    default:
      if(instrm.eof()) {
	error = finishFile(); if(error) return error;
      }
      ungetc(c);
      skipWord(0,0,0);
      nWordsRead++;
    }
  }
  return 0;
}

/******************************************************************************
Title: readKeys
Description: Skip or read until one of several key words is read.
             Processes commands along the way.
             Returns the C-number of the key read in "keyRead".
             If reading, not skipping, the key must be read next from
             the input, except commands, WS, and explicit comments.
Input:  keyWord         Pointers to the key words to search for.
                        These key words should be statically defined,
			so that error conditions can be processed properly.
        nKeys           Number of key words.
        skip_echoWord   Flag:  1 = skip to the key (otherwise read to the key)
                              -1 = echo the word read (if echoing is enabled).
Output: keyRead  The C-index of the key read.
Return: Error flag (see reader.h).
*/
int readKeys(char** keyWord, int nKeys, int& keyRead, int skip_echoWord) {
  int k;
  char c,cc;
  int error = 0;
  int next;
  int* possibleKey;
  int stillPossible;

  int skip = (skip_echoWord > 0);
  int echoWord = (skip_echoWord < 0) && echoing;
  int echoDone = !echoWord;

  if(nKeys == 1) {
    save_keyWord_1 = *keyWord;
    save_keyWord = &save_keyWord_1;
  } else {
    save_keyWord = keyWord;
  }
  save_nKeyWords = nKeys;
  keyRead = -1;

  for(k=0; k<nKeys; k++) {
    if(keyWord[k][0] == '\0') { keyRead = k; goto exit1; }
  }
  possibleKey = new int[nKeys];
  if(possibleKey == NULL) { error = ERR_MEMORY; goto exit; }

  do {
    error = skipWS(0,0);
    if(error) {
      if(error == ERR_EOF || error == ERR_SECTION_READ
	 || error == ERR_SUBSECTION_READ
	 || error == ERR_SUBSUBSECTION_READ) error = ERR_KEY_NOT_READ;
      echoWord = 0;
      goto exit;
    }

    for(k=0; k<nKeys; k++) possibleKey[k] = 1;
    next = -1;

    do {
      c = getc();
      if(c == TAB || c == NL || c == CR || instrm.eof()) {
	ungetc(c);
        echoDone = 1;
	break;
      }
      if(echoWord && c != SPACE) *echoFile << c;

      next++;
      stillPossible = 0;

      for(k=0; k<nKeys; k++) {
	if(possibleKey[k]) {
	  if(c == keyWord[k][next]) {
	    if(keyWord[k][next+1] == '\0') {
	      cc = getc();
	      ungetc(cc);
	      if(cc==SPACE || cc==TAB || cc==NL || cc == CR || instrm.eof()) {
		if(echoing && !echoWord) *echoFile << keyWord[k] << endl;
		keyRead = k;
		error = 0;
		goto exit;
	      }
	      possibleKey[k] = 0;
	    }
	  } else {
	    possibleKey[k] = 0;
	  }
	}
	stillPossible += possibleKey[k];
      }

      if(echoWord && c == SPACE) {
	if(stillPossible) {
	  *echoFile << c;
	} else {
	  echoDone = 1;
	}
      }
    } while(stillPossible);
  } while(skip);

  error = ERR_KEY_NOT_READ;
  while(!echoDone) {
    c = getc();
    if(c == SPACE || c == TAB || c == NL || c == CR || instrm.eof()) {
      ungetc(c);
      goto exit;
    }
    *echoFile << c;
  }

 exit:
  if(echoWord) *echoFile << endl;
  delete[] possibleKey;

 exit1:
  if(nKeys <= 1) keyRead = 0;
  return error;
}

/******************************************************************************
Title: readKey
Description: Skip or read until a key word is read.
             Processes commands along the way.
             Returns the C-number of the key read in "keyRead".
             If reading, not skipping, the key must be read next from
             the input, except commands, WS, and explicit comments.
Input:  keyWord         The key word to search for.
                        This key word should be statically defined,
			so that error conditions can be processed properly.
        skip_echoWord   Flag:  1 = skip to the key (otherwise read to the key)
                              -1 = echo the word read (if echoing is enabled).
Return: Error flag (see reader.h).
*/
int readKey(char* keyWord, int skip_echoWord) {
  int tmp;
  return readKeys(&keyWord,1,tmp,skip_echoWord);
}

/******************************************************************************
Title: readSection
Description: Skip or read until a SECTION command is read.
             Processes other commands along the way.
Input: skip     Flag: skip to the sectional command (otherwise read to it)
Return: Error flag (see reader.h).
*/
int readSection(int skip) {
  int error = 0;

  while(1) {
    error = skipWS(0,0);
    switch(error) {
    case ERR_SECTION_READ:
      return 0;
    case ERR_EOF:
      return ERR_SECTION_NOT_READ;
    case ERR_SUBSECTION_READ:
    case ERR_SUBSUBSECTION_READ:
      readPWS(0);
      skipWord(1,0,1);
      skip = 1;
      break;
    case 0:
      if(!skip) return ERR_SECTION_NOT_READ;
      skipWord(1,0,0);
      break;
    default:
      return error;
    }
  }
}

/******************************************************************************
Title: readSubSection
Description: Skip or read until a SUBSECTION command is read.
             Processes other commands along the way.
Input: skip     Flag: skip to the sectional command (otherwise read to it)
Return: Error flag (see reader.h).
*/
int readSubSection(int skip) {
  int error = 0;

  while(1) {
    error = skipWS(0,0);
    switch(error) {
    case ERR_SUBSECTION_READ:
      return 0;
    case ERR_EOF:
    case ERR_SECTION_READ:
      return ERR_SUBSECTION_NOT_READ;
    case ERR_SUBSUBSECTION_READ:
      readPWS(0);
      skipWord(1,0,1);
      skip = 1;
      break;
    case 0:
      if(!skip) return ERR_SUBSECTION_NOT_READ;
      skipWord(1,0,0);
      break;
    default:
      return error;
    }
  }
}

/******************************************************************************
Title: readSubSubSection
Description: Skip or read until a SUBSUBSECTION command is read.
             Processes other commands along the way.
Input: skip     Flag: skip to the sectional command (otherwise read to it)
Return: Error flag (see reader.h).
*/
int readSubSubSection(int skip) {
  int error = 0;

  while(1) {
    error = skipWS(0,0);
    switch(error) {
    case ERR_SUBSUBSECTION_READ:
      return 0;
    case ERR_EOF:
    case ERR_SECTION_READ:
    case ERR_SUBSECTION_READ:
      return ERR_SUBSUBSECTION_NOT_READ;
    case 0:
      if(!skip) return ERR_SUBSUBSECTION_NOT_READ;
      skipWord(1,0,0);
      break;
    default:
      return error;
    }
  }
}

/******************************************************************************
Title: skipToSections
Description: Skip to one of several new sections.
             We assume we have just read a section command.
*/
int skipToSections(char** keys, int nKeys, int& keyRead) {
  int error = 0;

  beginSkip();

  while(1) {
    error = readKeys(keys,nKeys,keyRead,-1);
    if(error != ERR_KEY_NOT_READ) goto exit;
    error = readSection(1);
    if(error == ERR_SECTION_NOT_READ) error = ERR_SECTION_KEY_NOT_READ;
    if(error) goto exit;
  }

exit:
  if(!error) endSkip();
  return error;
}

/******************************************************************************
Title: skipToSubSections
Description: Skip to one of several new subsections.
             We assume we have just read a subsection command.
*/
int skipToSubSections(char** keys, int nKeys, int& keyRead) {
  int error = 0;

  beginSkip();

  while(1) {
    error = readKeys(keys,nKeys,keyRead,-1);
    if(error != ERR_KEY_NOT_READ) goto exit;
    error = readSubSection(1);
    if(error == ERR_SUBSECTION_NOT_READ) error = ERR_SUBSECTION_KEY_NOT_READ;
    if(error) goto exit;
  }

exit:
  if(!error) endSkip();
  return error;
}

/******************************************************************************
Title: skipToSubSubSections
Description: Skip to one of several new subsubsections.
             We assume we have just read a subsubsection command.
*/
int skipToSubSubSections(char** keys, int nKeys, int& keyRead) {
  int error = 0;

  beginSkip();

  while(1) {
    error = readKeys(keys,nKeys,keyRead,-1);
    if(error != ERR_KEY_NOT_READ) goto exit;
    error = readSubSubSection(1);
    if(error == ERR_SUBSUBSECTION_NOT_READ)
      error = ERR_SUBSUBSECTION_KEY_NOT_READ;
    if(error) goto exit;
  }

exit:
  if(!error) endSkip();
  return error;
}

/******************************************************************************
Title: skipToSection
Description: Skip to a new section.
             We assume we have just read a section command.
*/
int skipToSection(char* keyWord) {
  int tmp;
  return skipToSections(&keyWord,1,tmp);
}

/******************************************************************************
Title: skipToSubSection
Description: Skip to a new subsection.
             We assume we have just read a subsection command.
*/
int skipToSubSection(char* keyWord) {
  int tmp;
  return skipToSubSections(&keyWord,1,tmp);
}

/******************************************************************************
Title: skipToSubSubSection
Description: Skip to a new subsubsection.
             We assume we have just read a subsubsection command.
*/
int skipToSubSubSection(char* keyWord) {
  int tmp;
  return skipToSubSubSections(&keyWord,1,tmp);
}

/******************************************************************************
Title: beginList
Description: Skip WS and then read a BEGIN_LIST symbol.
             Processes commands along the way.
Return: Error flag (see reader.h).
*/
int beginList() {
  int error = 0;
  char c;

  error = skipWS(0,0);
  if(error) {
    if(error == ERR_EOF) error = ERR_BEGIN_LIST_NOT_READ;
    goto exit;
  }
  c = getc();
  if(c != BEGIN_LIST) {
    error = ERR_BEGIN_LIST_NOT_READ;
    ungetc(c);
    goto exit;
  }
  if(echoing) *echoFile << c << endl;

 exit:
  return error;
}

/******************************************************************************
Title: endList
Description: Skip WS and then read an END_LIST symbol.
             Processes commands along the way.
Return: Error flag (see reader.h).
*/
int endList() {
  int error = 0;
  char c;

  error = skipWS(0,0);
  if(error) {
    if(error == ERR_EOF || error == ERR_SECTION_READ
       || error == ERR_SUBSECTION_READ
       || error == ERR_SUBSUBSECTION_READ) error = ERR_END_LIST_NOT_READ;
    goto exit;
  }
  c = getc();
  if(c != END_LIST) {
    error = ERR_END_LIST_NOT_READ;
    ungetc(c);
    goto exit;
  }
  if(echoing) *echoFile << c << endl;

 exit:
  return error;
}

/******************************************************************************
Title: skipToEndList
Description: Skip to just before an unmatched END_LIST symbol.
             Processes commands along the way.
Return: Error flag (see reader.h).
*/
int skipToEndList() {
  int error = 0;
  char c;
  int nBeginLists = 1;

  while(1) {
    error = skipWS(0,0);
    if(error) {
      if(error == ERR_EOF || error == ERR_SECTION_READ
	 || error == ERR_SUBSECTION_READ
	 || error == ERR_SUBSUBSECTION_READ)
	error = ERR_END_LIST_NOT_READ;
      goto exit;
    }
    c = getc();
    switch(c) {
    case BEGIN_LIST:
      nBeginLists++;
      if(echoing) *echoFile << c;
      break;
    case END_LIST:
      nBeginLists--;
      if(nBeginLists == 0) {
	ungetc(c);
	goto exit;
      }
      if(echoing) *echoFile << c;
      break;
    default:
      skipWord(1,1,0);
    }
  }

 exit:
  if(echoing) *echoFile << endl;
  return error;
}

/******************************************************************************
Title: readUnits
Description: Read the units expression of a real number.
Inputs: dimen    The required dimensions of the units expression.
Output: value    The value of the units expression.
Return: Error flag (see reader.h)
*/
int readUnits(double& value, int* dimen) {
  int error = 0;
  int i;
  char c,str[MAX_LENGTH_UNITS_STR];
  int mode = 1;
  int dimenCheck;

  if(readPWS(screenIn())) {
    value = 1;
    return NO_UNITS;
  }
  c = getc();
  if(c != BEGIN_UNIT) {
    ungetc(c);
    value = 1;
    return NO_UNITS;
  }
  c = getc();
  if(c == BEGIN_UNIT) {
    mode = 2;
  } else {
    ungetc(c);
  }

  for(i = 0; i < MAX_LENGTH_UNITS_STR-1; i++) {
    readPWS(0);
    c = getc();
    if(c == END_UNIT) break;
    str[i] = c;
  }
  if(i >= MAX_LENGTH_UNITS_STR-1) return ERR_UNITS_TOO_LONG;
  str[i] = '\0';
  if(mode == 2) {
    c = getc();
    if(c != END_UNIT) return ERR_ON_LINE;
  }

  if(echoing) {
    if(mode == 2) {
      *echoFile << " [[" << str << "]]";
    } else {
      *echoFile << " [" << str << "]";
    }
  }

  dimenCheck = (mode == 1);
  error = expr_eval(str,value,mode,dimenCheck,dimen);
  switch(error) {
  case 0:                           return 0;
  case EXPR_NEGATIVE:               return ERR_UNITS_NEGATIVE;
  case EXPR_ADDITION_FOUND:         return ERR_UNITS_ADDITION_FOUND;
  case EXPR_DIMENS_FAILED:          return ERR_UNITS_WRONG_DIMEN;
  case EXPR_TOO_BIG:                return ERR_UNITS_TOO_LONG;
  case EXPR_NONNESTED_PARENTHESES:  return ERR_NONNESTED_PARENTHESES;
  case EXPR_ARITHMETIC_OP:          return ERR_ARITHMETIC_OP;
  case EXPR_IMPROPER:               return ERR_UNITS_IMPROPER_EXPR;
  case EXPR_UNKNOWN_UNIT:           return ERR_UNKNOWN_UNIT;
  default:                          return ERR_UNITS;
  }
}

/******************************************************************************
Title: readTempUnits
Description: Read and convert to the absolute temperature expression.
Inputs:   temp            The absolute temperature in the given units.
          echoConversion  Flag to echo the conversion comment.
Outputs:  temp            The absolute temperature in the internal units.
Return: Error flag (see reader.h)
*/
int readTempUnits(double& temp, int echoConversion) {
  int error = 0;
  int i;
  char c;
  double temp_in;
  char degrees[3][MAX_DEGREE_NAME_LENGTH] = {DEG_K, DEG_C, DEG_F};
  char str[MAX_LENGTH_UNITS_STR];
  int base = expr_baseTempIndex();

  if(base < 0) return ERR_IN_CODE;

  readPWS(screenIn());
  c = getc();
  if(c != BEGIN_UNIT) {
    ungetc(c);
    return NO_UNITS;
  }
  temp_in = temp;

  for(i = 0; i < MAX_LENGTH_UNITS_STR-1; i++) {
    readPWS(0);
    c = getc();
    if(c == END_UNIT) break;
    str[i] = c;
  }
  if(i >= MAX_LENGTH_UNITS_STR-1) return ERR_UNITS_TOO_LONG;
  str[i] = '\0';

  error = expr_convertTempBase(temp,str);
  switch(error) {
  case 0:                  break;
  case EXPR_UNKNOWN_UNIT:  return ERR_UNKNOWN_UNIT;
  default:                 return ERR_UNITS;
  }

  if(echoing) {
    *echoFile << " [" << str << "]";
    if(echoConversion) *echoFile << "   # " << temp_in << " [" << str << "] = "
				 << temp << " [" << degrees[base] << "]";
  }

  return 0;
}

/******************************************************************************
Title: readScalar
Description: Read a scalar.
*/
int readScalar(int& i) {
  int error = skipWS(0,0); if(error) return error;
  instrm >> i;
  if(!instrm.good()) return ERR_ON_LINE;
  char c = getc();
  ungetc(c);
  if(c == '.') return ERR_ON_LINE;
  if(echoing) *echoFile << i << endl;
  return 0;
}
int readScalar(double& x) {
  int error = skipWS(0,0); if(error) return error;
  instrm >> x;
  if(!instrm.good()) return ERR_ON_LINE;
  if(echoing) *echoFile << x << endl;
  return 0;
}
int readScalar(float& x) {
  int error = skipWS(0,0); if(error) return error;
  instrm >> x;
  if(!instrm.good()) return ERR_ON_LINE;
  if(echoing) *echoFile << x << endl;
  return 0;
}

int readScalar(double& var, int* dimen) {
  double unit,value;

  int error = skipWS(0,0); if(error) return error;
  instrm >> value;
  if(!instrm.good()) return ERR_ON_LINE;
  if(echoing) *echoFile << value;

  error = readUnits(unit,dimen);
  if(error && error != NO_UNITS) {
    if(echoing) *echoFile << endl;
    return error;
  }

  var = value*unit;
  if(echoing) {
    if(error == NO_UNITS) {
      *echoFile << endl;
    } else {
      *echoFile << "   # " << value << " * " << unit << " = " << var << endl;
    }
  }

  return 0;
}

int readScalar(double& var, int* dimen, double& unitVal) {
  double value;

  int error = skipWS(0,0); if(error) return error;
  instrm >> value;
  if(!instrm.good()) return ERR_ON_LINE;
  if(echoing) *echoFile << value;

  error = readUnits(unitVal,dimen);
  if(error && error != NO_UNITS) {
    if(echoing) *echoFile << endl;
    return error;
  }

  var = value*unitVal;
  if(echoing) {
    if(error == NO_UNITS) {
      *echoFile << endl;
    } else {
      *echoFile << "   # " << value << " * " << unitVal <<" = " << var << endl;
    }
  }

  return 0;
}

int readScalar(float& var, int* dimen) {
  double unit;
  float value;

  int error = skipWS(0,0); if(error) return error;
  instrm >> value;
  if(!instrm.good()) return ERR_ON_LINE;
  if(echoing) *echoFile << value;

  error = readUnits(unit,dimen);
  if(error && error != NO_UNITS) {
    if(echoing) *echoFile << endl;
    return error;
  }

  var = value*unit;
  if(echoing) {
    if(error == NO_UNITS) {
      *echoFile << endl;
    } else {
      *echoFile << "   # " << value << " * " << unit << " = " << var << endl;
    }
  }

  return 0;
}

int readScalar(float& var, int* dimen, double& unitVal) {
  float value;

  int error = skipWS(0,0); if(error) return error;
  instrm >> value;
  if(!instrm.good()) return ERR_ON_LINE;
  if(echoing) *echoFile << value;

  error = readUnits(unitVal,dimen);
  if(error && error != NO_UNITS) {
    if(echoing) *echoFile << endl;
    return error;
  }

  var = value*unitVal;
  if(echoing) {
    if(error == NO_UNITS) {
      *echoFile << endl;
    } else {
      *echoFile << "   # " << value << " * " << unitVal <<" = " << var << endl;
    }
  }

  return 0;
}

/******************************************************************************
Title: readTemp
Description: Read a temperature.
*/
int readTemp(double& var) {
  int error = skipWS(0,0); if(error) return error;
  instrm >> var;
  if(!instrm.good()) return ERR_ON_LINE;
  if(echoing) *echoFile << var;

  error = readTempUnits(var,1);
  if(error && error != NO_UNITS) return error;

  if(echoing) *echoFile << endl;

  return 0;
}

int readTemp(float& var) {
  int error = skipWS(0,0); if(error) return error;
  instrm >> var;
  if(!instrm.good()) return ERR_ON_LINE;
  if(echoing) *echoFile << var;

  double dvar = var;
  error = readTempUnits(dvar,1);
  if(error && error != NO_UNITS) return error;
  var = dvar;

  if(echoing) *echoFile << endl;

  return 0;
}

/******************************************************************************
Title: readString
Description: Read a std::string, up to non-punctuation white space.
Output: string  The string read.
Return: Error flag (see reader.h).
*/
int readString(std::string& str) {
  char c;
  int i;
  int error = 0;

  error = skipWS(0,0); if(error) return error;
  for(i=0; 1; i++) {
    c = getc();
    //std::cout << "[" << c << "]  " << (int)c << "\n";
    if(c == NL || c == CR || instrm.eof() || c == SPACE || c == TAB ) {
      ungetc(c);
      break;
    }
    str += c;
  }
  if(i==0) return ERR_BAD_DATA;
  if(echoing) *echoFile << str << endl;

  return 0;
}

/******************************************************************************
Title: readString
Description: Read a string, up to non-punctuation white space or one less
             than the length of the string, and add the string
             terminating character.
Input:  len     The length of the string.
Output: string  The string read.
Return: Error flag (see reader.h).
*/
int readString(char string[], int len) {
  char c;
  int i;
  int error = 0;

  error = skipWS(0,0); if(error) return error;
  for(i=0; i<len-1; i++) {
    c = getc();
    if(c == NL || c == CR || instrm.eof() || c == SPACE || c == TAB ) {
      ungetc(c);
      break;
    }
    string[i] = c;
  }
  string[i] = '\0';
  if(echoing) *echoFile << string << endl;

  return 0;
}

/******************************************************************************
Title: readLine
Description: Read a line, up to EOL, EOF, or one less than the length of
             the string, and add the string terminating character.
Input:  len     The length of the string.
Output: string  The string read.
Return: Error flag (see reader.h).
*/
int readLine(char string[], int len) {
  char c;
  int i;
  int error = 0;

  error = skipWS(0,0); if(error) return error;
  for(i=0; i<len-1; i++) {
    c = getc();
    if(c == NL || c == CR || instrm.eof()) {
      ungetc(c);
      break;
    }
    string[i] = c;
  }
  string[i] = '\0';
  if(echoing) *echoFile << string << endl;

  return 0;
}

/******************************************************************************
Title: readArray and readVector
Description: Read a <TYPE> array up to the size of the given array. 
             Handles the repetition symbol in the input.
             Possible options are:
             1. constant <value>
             2. linear <value1> <value2> (except int's)
             3. uniform <value1> <value2> (same as linear)
             4. { <value> ... } (no units allowed)
             5. { [<units>] <value> ... } (units apply to all numbers)
             6. ??? (return ERR_UNSUPPORTED_OPTION for possible further
                     processing)
Input:  len      The length of the array (omit for std::vector).
        dimen    The required dimensions of the units expression.
Output: array    The array read.
        unitVal  The value of the units (for value2 in linear/uniform cases)
Return: Error flag (see reader.h).
*/
int readArray(int* array, int len) {
  int error = skipWS(0,0); if(error) return error;
  char c = getc();
  ungetc(c);

  switch(c) {
  case 'c': {
    char str[SIZE_WORD];
    readString(str,SIZE_WORD);
    if(!strcmp(str,"constant")) {
      int value;
      error = readScalar(value); if(error) return error;
      for(int i=0; i<len; i++) array[i] = value;
    } else return ERR_UNSUPPORTED_OPTION;
    return 0;
  }
  case BEGIN_LIST: {
    error = beginList(); if(error) return error;
    int nRead = 0;
    error = readBlock(array,len,nRead,len); if(error) return error;
    error = endList(); if(error) return error;
    return 0;
  }
  default:
    return ERR_UNSUPPORTED_OPTION;
  }
}

int readArray(double* array, int len) {
  int error = skipWS(0,0); if(error) return error;
  char c = getc();
  ungetc(c);

  switch(c) {
  case 'c':
  case 'l':
  case 'u': {
    char str[SIZE_WORD];
    readString(str,SIZE_WORD);
    if(!strcmp(str,"constant")) {
      double value;
      error = readScalar(value); if(error) return error;
      for(int i=0; i<len; i++) array[i] = value;
    } else if(!strcmp(str,"uniform") || !strcmp(str,"linear")) {
      double value1,value2;
      error = readScalar(value1); if(error) return error;
      error = readScalar(value2); if(error) return error;
      array[0] = value1;
      if(len==1) return 0;
      double da = (value2 - value1)/(len-1);
      for(int i=1; i<len-1; i++) array[i] = array[i-1] + da;
      array[len-1] = value2;
    } else return ERR_UNSUPPORTED_OPTION;
    return 0;
  }
  case BEGIN_LIST: {
    error = beginList(); if(error) return error;
    int nRead = 0;
    error = readBlock(array,len,nRead,len); if(error) return error;
    error = endList(); if(error) return error;
    return 0;
  }
  default:
    return ERR_UNSUPPORTED_OPTION;
  }
}

int readArray(float* array, int len) {
  int error = skipWS(0,0); if(error) return error;
  char c = getc();
  ungetc(c);

  switch(c) {
  case 'c':
  case 'l':
  case 'u': {
    char str[SIZE_WORD];
    readString(str,SIZE_WORD);
    if(!strcmp(str,"constant")) {
      float value;
      error = readScalar(value); if(error) return error;
      for(int i=0; i<len; i++) array[i] = value;
    } else if(!strcmp(str,"uniform") || !strcmp(str,"linear")) {
      double value1,value2;
      error = readScalar(value1); if(error) return error;
      error = readScalar(value2); if(error) return error;
      array[0] = value1;
      if(len==1) return 0;
      float da = (value2 - value1)/(len-1);
      for(int i=1; i<len-1; i++) array[i] = array[i-1] + da;
      array[len-1] = value2;
    } else return ERR_UNSUPPORTED_OPTION;
    return 0;
  }
  case BEGIN_LIST: {
    error = beginList(); if(error) return error;
    int nRead = 0;
    error = readBlock(array,len,nRead,len); if(error) return error;
    error = endList(); if(error) return error;
    return 0;
  }
  default:
    return ERR_UNSUPPORTED_OPTION;
  }
}

static int getUnits(double& unitVal, int* dimen) {
  int error;

  error = skipWS(0,0); if(error) return error;
  error = readUnits(unitVal,dimen);
  if(error == NO_UNITS) {
    return 0;
  } else if(error) return error;
  //char str[34];
  reader_echo("   # * ");
  reader_echo(unitVal);
  reader_echo("\n");

  return 0;
}

int readArray(double* array, int len, int* dimen, double& unitVal) {
  int error = skipWS(0,0); if(error) return error;
  char c = getc();
  ungetc(c);

  switch(c) {
  case 'c':
  case 'l':
  case 'u': {
    char str[SIZE_WORD];
    readString(str,SIZE_WORD);
    if(!strcmp(str,"constant")) {
      double value;
      error = readScalar(value,dimen,unitVal); if(error) return error;
      for(int i=0; i<len; i++) array[i] = value;
    } else if(!strcmp(str,"uniform") || !strcmp(str,"linear")) {
      double value1,value2;
      error = readScalar(value1,dimen,unitVal); if(error) return error;
      error = readScalar(value2,dimen,unitVal); if(error) return error;
      array[0] = value1;
      if(len==1) return 0;
      double da = (value2 - value1)/(len-1);
      for(int i=1; i<len-1; i++) array[i] = array[i-1] + da;
      array[len-1] = value2;
    } else return ERR_UNSUPPORTED_OPTION;
    return 0;
  }
  case BEGIN_LIST: {
    error = beginList(); if(error) return error;
    error = getUnits(unitVal,dimen); if(error) return error;
    int nRead = 0;
    error = readBlock(array,len,nRead,len,unitVal); if(error) return error;
    error = endList(); if(error) return error;
    return 0;
  }
  default:
    return ERR_UNSUPPORTED_OPTION;
  }
}

int readArray(float* array, int len, int* dimen, double& unitVal) {
  int error = skipWS(0,0); if(error) return error;
  char c = getc();
  ungetc(c);

  switch(c) {
  case 'c':
  case 'l':
  case 'u': {
    char str[SIZE_WORD];
    readString(str,SIZE_WORD);
    if(!strcmp(str,"constant")) {
      float value;
      error = readScalar(value,dimen,unitVal); if(error) return error;
      for(int i=0; i<len; i++) array[i] = value;
    } else if(!strcmp(str,"uniform") || !strcmp(str,"linear")) {
      double value1,value2;
      error = readScalar(value1,dimen,unitVal); if(error) return error;
      error = readScalar(value2,dimen,unitVal); if(error) return error;
      array[0] = value1;
      if(len==1) return 0;
      float da = (value2 - value1)/(len-1);
      for(int i=1; i<len-1; i++) array[i] = array[i-1] + da;
      array[len-1] = value2;
    } else return ERR_UNSUPPORTED_OPTION;
    return 0;
  }
  case BEGIN_LIST: {
    error = beginList(); if(error) return error;
    error = getUnits(unitVal,dimen); if(error) return error;
    int nRead = 0;
    error = readBlock(array,len,nRead,len,unitVal); if(error) return error;
    error = endList(); if(error) return error;
    return 0;
  }
  default:
    return ERR_UNSUPPORTED_OPTION;
  }
}

int readVector(std::vector<int>& array) {
  int error = skipWS(0,0); if(error) return error;
  char c = getc();
  ungetc(c);

  switch(c) {
  case 'c': {
    char str[SIZE_WORD];
    readString(str,SIZE_WORD);
    if(!strcmp(str,"constant")) {
      int value, len;
      error = readScalar(value); if(error) return error;
      error = readScalar(len); if(error) return error;
      array.assign(len,value);
    } else return ERR_UNSUPPORTED_OPTION;
    return 0;
  }
  case BEGIN_LIST: {
    error = beginList(); if(error) return error;
    int nRead = 0;
    error = readBlock(array,nRead); if(error) return error;
    error = endList(); if(error) return error;
    return 0;
  }
  default:
    return ERR_UNSUPPORTED_OPTION;
  }
}

int readVector(std::vector<double>& array) {
  int error = skipWS(0,0); if(error) return error;
  char c = getc();
  ungetc(c);

  switch(c) {
  case 'c': {
    char str[SIZE_WORD];
    readString(str,SIZE_WORD);
    if(!strcmp(str,"constant")) {
      double value;
      int len;
      error = readScalar(value); if(error) return error;
      error = readScalar(len); if(error) return error;
      array.assign(len,value);
    } else return ERR_UNSUPPORTED_OPTION;
    return 0;
  }
  case BEGIN_LIST: {
    error = beginList(); if(error) return error;
    int nRead = 0;
    error = readBlock(array,nRead); if(error) return error;
    error = endList(); if(error) return error;
    return 0;
  }
  default:
    return ERR_UNSUPPORTED_OPTION;
  }
}

int readVector(std::vector<float>& array) {
  int error = skipWS(0,0); if(error) return error;
  char c = getc();
  ungetc(c);

  switch(c) {
  case 'c': {
    char str[SIZE_WORD];
    readString(str,SIZE_WORD);
    if(!strcmp(str,"constant")) {
      float value;
      int len;
      error = readScalar(value); if(error) return error;
      error = readScalar(len); if(error) return error;
      array.assign(len,value);
    } else return ERR_UNSUPPORTED_OPTION;
    return 0;
  }
  case BEGIN_LIST: {
    error = beginList(); if(error) return error;
    int nRead = 0;
    error = readBlock(array,nRead); if(error) return error;
    error = endList(); if(error) return error;
    return 0;
  }
  default:
    return ERR_UNSUPPORTED_OPTION;
  }
}

int readVector(std::vector<double>& array, int* dimen, double& unitVal) {
  int error = skipWS(0,0); if(error) return error;
  char c = getc();
  ungetc(c);

  switch(c) {
  case 'c': {
    char str[SIZE_WORD];
    readString(str,SIZE_WORD);
    if(!strcmp(str,"constant")) {
      double value;
      int len;
      error = readScalar(value,dimen,unitVal); if(error) return error;
      error = readScalar(len); if(error) return error;
      array.assign(len,value);
    } else return ERR_UNSUPPORTED_OPTION;
    return 0;
  }
  case BEGIN_LIST: {
    error = beginList(); if(error) return error;
    error = getUnits(unitVal,dimen); if(error) return error;
    int nRead = 0;
    error = readBlock(array,nRead,unitVal); if(error) return error;
    error = endList(); if(error) return error;
    return 0;
  }
  default:
    return ERR_UNSUPPORTED_OPTION;
  }
}

int readVector(std::vector<float>& array, int* dimen, double& unitVal) {
  int error = skipWS(0,0); if(error) return error;
  char c = getc();
  ungetc(c);

  switch(c) {
  case 'c': {
    char str[SIZE_WORD];
    readString(str,SIZE_WORD);
    if(!strcmp(str,"constant")) {
      float value;
      int len;
      error = readScalar(value,dimen,unitVal); if(error) return error;
      error = readScalar(len); if(error) return error;
      array.assign(len,value);
    } else return ERR_UNSUPPORTED_OPTION;
    return 0;
  }
  case BEGIN_LIST: {
    error = beginList(); if(error) return error;
    error = getUnits(unitVal,dimen); if(error) return error;
    int nRead = 0;
    error = readBlock(array,nRead,unitVal); if(error) return error;
    error = endList(); if(error) return error;
    return 0;
  }
  default:
    return ERR_UNSUPPORTED_OPTION;
  }
}

/******************************************************************************
Title: readBlock
Description: Read a <TYPE> array up to the size of the given array. Supports
             multiple calls.  Handles the repetition symbol in the input.
Input:  len            The length of the array.
        nToReadTotal   The total number of items to read through
	               sucessive calls.
Output: array          The array read.
        nRead          The number of items read so far.
Return: Error flag (see reader.h).
*/
/* <TYPE == int */
int readBlock(int* array, int len, int& nRead, int nToReadTotal) {
  int error = 0;
  char c;
  static int nReadTotal = 0;
  static int nRepeated = 0;
  static int valueRepeated;
  int i=0;

  nRead = (len >= nToReadTotal-nReadTotal) ? nToReadTotal-nReadTotal : len;
  while(i<nRead) {
    if(nRepeated <= 0) {
      error = skipWS(0,0); if(error) return error;
      instrm >> array[i]; if(!instrm.good()) return ERR_IN_DATA_BLOCK;
      do { c = getc(); } while(c == SPACE);
      if(c != REPETITION) {
	ungetc(c);
	i++;
      } else {
	nRepeated = array[i];
	error = skipWS(0,0); if(error) return error;
	instrm >> valueRepeated; if(!instrm.good()) return ERR_IN_DATA_BLOCK;
      }
    } else {
      do {
	array[i] = valueRepeated;
	i++; nRepeated--;
      } while(i < nRead && nRepeated > 0);
    }
  }
  nReadTotal += nRead;
  if(nReadTotal >= nToReadTotal) {
    nReadTotal = 0;
    if(nRepeated > 0) return ERR_IN_DATA_BLOCK;
  }

  if(echoing) {
    static int j = 0;
    for(i=0; i<nRead; i++) {
      *echoFile << " " << array[i];
      j++;
      if(j == 5) { *echoFile << endl; j = 0; }
    }
    if(nReadTotal == 0 && j != 5) { *echoFile << endl; j = 0; }
  }
  return 0;
}

/* <TYPE == double */
static void multByUnits(double* array, int n, double units) {
  int i;
  if(units != 1) for(i=0; i<n; i++) array[i] *= units;
}

int readBlock(double* array, int len, int& nRead,
		  int nToReadTotal, double units) {
  int error = 0;
  char c;
  static int nReadTotal = 0;
  static int nRepeated = 0;
  static double valueRepeated;
  int i=0;

  nRead = (len >= nToReadTotal-nReadTotal) ? nToReadTotal-nReadTotal : len;
  while(i<nRead) {
    if(nRepeated <= 0) {
      error = skipWS(0,0); if(error) return error;
      instrm >> array[i]; if(!instrm.good()) return ERR_IN_DATA_BLOCK;
      do { c = getc(); } while(c == SPACE);
      if(c != REPETITION) {
	ungetc(c);
	i++;
      } else {
	nRepeated = (int)array[i];
	error = skipWS(0,0); if(error) return error;
	instrm >> valueRepeated; if(!instrm.good()) return ERR_IN_DATA_BLOCK;
      }
    } else {
      do {
	array[i] = valueRepeated;
	i++; nRepeated--;
      } while(i < nRead && nRepeated > 0);
    }
  }
  nReadTotal += nRead;
  if(nReadTotal >= nToReadTotal) {
    nReadTotal = 0;
    if(nRepeated > 0) return ERR_IN_DATA_BLOCK;
  }
  multByUnits(array,nRead,units);

  if(echoing) {
    static int j = 0;
    double unitsInv = (units != 0) ? 1/units : 1;
    for(i=0; i<nRead; i++) {
      *echoFile << " " << array[i]*unitsInv;
      j++;
      if(j == 5) { *echoFile << endl; j = 0; }
    }
    if(nReadTotal == 0 && j != 5) { *echoFile << endl; j = 0; }
  }
  return 0;
}

int readBlock(double* array, int len, int& nRead, int nToReadTotal) {
  double one = 1;
  return readBlock(array, len, nRead, nToReadTotal, one);
}

/* <TYPE == float */
static void multByUnits(float* array, int n, double units) {
  int i;
  if(units != 1) for(i=0; i<n; i++) array[i] *= units;
}

int readBlock(float* array, int len, int& nRead,
		  int nToReadTotal, double units) {
  int error = 0;
  char c;
  static int nReadTotal = 0;
  static int nRepeated = 0;
  static float valueRepeated;
  int i=0;

  nRead = (len >= nToReadTotal-nReadTotal) ? nToReadTotal-nReadTotal : len;
  while(i<nRead) {
    if(nRepeated <= 0) {
      error = skipWS(0,0); if(error) return error;
      instrm >> array[i]; if(!instrm.good()) return ERR_IN_DATA_BLOCK;
      do { c = getc(); } while(c == SPACE);
      if(c != REPETITION) {
	ungetc(c);
	i++;
      } else {
	nRepeated = (int)array[i];
	error = skipWS(0,0); if(error) return error;
	instrm >> valueRepeated; if(!instrm.good()) return ERR_IN_DATA_BLOCK;
      }
    } else {
      do {
	array[i] = valueRepeated;
	i++; nRepeated--;
      } while(i < nRead && nRepeated > 0);
    }
  }
  nReadTotal += nRead;
  if(nReadTotal >= nToReadTotal) {
    nReadTotal = 0;
    if(nRepeated > 0) return ERR_IN_DATA_BLOCK;
  }
  multByUnits(array,nRead,units);

  if(echoing) {
    static int j = 0;
    float unitsInv = (units != 0) ? 1/units : 1;
    for(i=0; i<nRead; i++) {
      *echoFile << " " << array[i]*unitsInv;
      j++;
      if(j == 5) { *echoFile << endl; j = 0; }
    }
    if(nReadTotal == 0 && j != 5) { *echoFile << endl; j = 0; }
  }
  return 0;
}

int readBlock(float* array, int len, int& nRead, int nToReadTotal) {
  double one = 1;
  return readBlock(array, len, nRead, nToReadTotal, one);
}

/******************************************************************************
Title: readBlock
Description: Read a <TYPE> vector. Handles the repetition symbol in the input.
Output: array          The vector read.
        nRead          The number of items read.
Return: Error flag (see reader.h).
*/
/* <TYPE == int */
int readBlock(std::vector<int>& array, int& nRead) {
  int error = 0;
  char c;
  int val;
  int startIndex = array.size();

  nRead=0;
  while(1) {
    error = skipWS(0,0); if(error) return error;
    c = getc();
    ungetc(c);
    if(c == END_LIST) break;

    instrm >> val; if(!instrm.good()) return ERR_IN_DATA_BLOCK;
    error = skipWS(0,0); if(error) return error;
    c = getc();
    if(c != REPETITION) {
      ungetc(c);
      array.push_back(val); 
      nRead++;
    } else {
      int nRepeated = val;
      error = skipWS(0,0); if(error) return error;
      instrm >> val; if(!instrm.good()) return ERR_IN_DATA_BLOCK;
      for(int i=0; i<nRepeated; i++) {
	array.push_back(val); 
	nRead++;
      }
    }
  }

  if(echoing) {
    int j = 0;
    for(int i=0; i<nRead; i++) {
      *echoFile << " " << array[startIndex+i];
      j++;
      if(j == 5) { *echoFile << endl; j = 0; }
    }
    if(j != 5) { *echoFile << endl; }
  }
  return 0;
}

/* <TYPE == double */
static void multByUnits(std::vector<double>& array, int startIndex,
			int nRead, double units) {
  if(units != 1) for(int i=0; i<nRead; i++) array[startIndex+i] *= units;
}

int readBlock(std::vector<double>& array, int& nRead, double units) {
  int error = 0;
  char c;
  double val;
  int startIndex = array.size();

  nRead=0;
  while(1) {
    error = skipWS(0,0); if(error) return error;
    c = getc();
    ungetc(c);
    if(c == END_LIST) break;

    instrm >> val; if(!instrm.good()) return ERR_IN_DATA_BLOCK;
    error = skipWS(0,0); if(error) return error;
    c = getc();
    if(c != REPETITION) {
      ungetc(c);
      array.push_back(val); 
      nRead++;
    } else {
      int nRepeated = val;
      error = skipWS(0,0); if(error) return error;
      instrm >> val; if(!instrm.good()) return ERR_IN_DATA_BLOCK;
      for(int i=0; i<nRepeated; i++) {
	array.push_back(val); 
	nRead++;
      }
    }
  }
  multByUnits(array,startIndex,nRead,units);

  if(echoing) {
    int j = 0;
    for(int i=0; i<nRead; i++) {
      *echoFile << " " << array[startIndex+i];
      j++;
      if(j == 5) { *echoFile << endl; j = 0; }
    }
    if(j != 5) { *echoFile << endl; }
  }
  return 0;
}

int readBlock(std::vector<double>& array, int& nRead) {
  double one = 1;
  return readBlock(array, nRead, one);
}

/* <TYPE == float */
static void multByUnits(std::vector<float>& array, int startIndex,
			int nRead, double units) {
  if(units != 1) for(int i=0; i<nRead; i++) array[startIndex+i] *= units;
}

int readBlock(std::vector<float>& array, int& nRead, double units) {
  int error = 0;
  char c;
  float val;
  int startIndex = array.size();

  nRead=0;
  while(1) {
    error = skipWS(0,0); if(error) return error;
    c = getc();
    ungetc(c);
    if(c == END_LIST) break;

    instrm >> val; if(!instrm.good()) return ERR_IN_DATA_BLOCK;
    error = skipWS(0,0); if(error) return error;
    c = getc();
    if(c != REPETITION) {
      ungetc(c);
      array.push_back(val); 
      nRead++;
    } else {
      int nRepeated = val;
      error = skipWS(0,0); if(error) return error;
      instrm >> val; if(!instrm.good()) return ERR_IN_DATA_BLOCK;
      for(int i=0; i<nRepeated; i++) {
	array.push_back(val); 
	nRead++;
      }
    }
  }
  multByUnits(array,startIndex,nRead,units);

  if(echoing) {
    int j = 0;
    for(int i=0; i<nRead; i++) {
      *echoFile << " " << array[startIndex+i];
      j++;
      if(j == 5) { *echoFile << endl; j = 0; }
    }
    if(j != 5) { *echoFile << endl; }
  }
  return 0;
}

int readBlock(std::vector<float>& array, int& nRead) {
  double one = 1;
  return readBlock(array, nRead, one);
}

/******************************************************************************
Title: openReader
Description: Initialize the routines for reading input files.
Inputs: inFileName0      The initial input file name.
        echoFileName     The echo file name.
	length_fileName  The length of the two file names (excluding '\0').
	batch0           Boolean for batch input.
	echo0            Boolean for echoing the input.
Return: Error flag (see reader.h).
*/
int openReader(const std::string& inFileName0, const std::string& echoFileName,
	       int length_fileName, int batch0, int echo0) {
  if(length_fileName < (int)inFileName0.size()) length_fileName = inFileName0.size();
  if(length_fileName < (int)echoFileName.size()) length_fileName = echoFileName.size();

  char nullEcho[1]; nullEcho[0] = '\0';
  return (echoFileName[0]!='\0') ?
	   openReader(inFileName0.c_str(),echoFileName.c_str(),length_fileName,batch0,echo0)
	   : openReader(inFileName0.c_str(),nullEcho,length_fileName,batch0,echo0);
}
int openReader(const char* inFileName0, const char* echoFileName,
	       int length_fileName, int batch0, int echo0) {
  int error = 0;
  int i;

  /* INITIALIZE GENERAL INPUT AND ECHO FILES */
  SIZE_FILENAME = length_fileName + 1;
  for(i=0; i<DEPTH_INCLUDES; i++) {
    inFileName[i] = new char[SIZE_FILENAME];
    if(inFileName[i] == NULL) return ERR_MEMORY;
  }
  inFileNo = -1;
  strcpy(inFileName[0],inFileName0);
  error = beginFile(); if(error) return error;

  batch = batch0;
  echoing = echo0;
  if(echoing) {
    if(echoFileName[0] != '\0') {
      echoFile_ofstream.open(echoFileName,ios::out);
      echoFile = &echoFile_ofstream;
      if(!echoFile_ofstream || !echoFile_ofstream.good()) {
	cerr << "Error opening echo file " << echoFileName << endl;
	echoing = 0;
      }
    } else {
      echoFile = &cout;
    }
  }

  return 0;
}

/******************************************************************************
Title: closeReader
Description: Close the routines for reading input files and process error flag.
Input:  error  The error flag (see reader.h).  Call with error = 0 if no
               error information is desired.
*/
void closeReader(int error) {
  int i;
  
  processReaderError(error);

  for(i=inFileNo; i>=0; i--) {
    error = finishFile();
    if(error && error != ERR_EOF) processReaderError(error);
  }

  if(echoing) {
    *echoFile << endl;
    if(echoFile != &cout) echoFile_ofstream.close();
  }

  for(i=0; i<DEPTH_INCLUDES; i++) delete[] inFileName[i];
}

/******************************************************************************
Title: screenIn
Description: Determines whether input comes from the screen (stdin).
Return: screenIn   Boolean
*/
int screenIn() {
  return (inFile[inFileNo]==&cin);
}

/******************************************************************************
Title: processReaderError
Description: Process error flag by writing the appropriate message.
Input: error  The error flag (see reader.h).
*/
int processReaderError(const int error) {
  if(error) {
    switch(error) {
    case ERR_UNMATCHED_END_COMMENT_CMD:
    case ERR_UNKNOWN_CMD:
    case ERR_IN_DATA_BLOCK:
    case ERR_IN_GRID_ARRAY:
    case ERR_UNSUPPORTED_OPTION:
    case ERR_BEGIN_LIST_NOT_READ:
    case ERR_SECTION_READ:
    case ERR_SUBSECTION_READ:
    case ERR_SUBSUBSECTION_READ:
    case ERR_NONNESTED_PARENTHESES:
    case ERR_ARITHMETIC_OP:
    case ERR_UNITS_NEGATIVE:
    case ERR_UNITS_ADDITION_FOUND:
    case ERR_UNITS_WRONG_DIMEN:
    case ERR_UNITS_TOO_LONG:
    case ERR_UNITS_IMPROPER_EXPR:
    case ERR_UNKNOWN_UNIT:
    case ERR_UNITS:
    case ERR_ON_LINE:
      switch(error) {
      case ERR_UNMATCHED_END_COMMENT_CMD:
	cerr << "Unmatched end comment command" << endl;
	break;
      case ERR_UNKNOWN_CMD:
	cerr << "Unknown " << COMMAND << " command" << endl;
	break;
      case ERR_IN_DATA_BLOCK:
	cerr << "Error in data block" << endl;
	break;
      case ERR_IN_GRID_ARRAY:
	cerr << "Grid array expected/bad grid" << endl;
	break;
      case ERR_UNSUPPORTED_OPTION:
	cerr << "Unsupported input option" << endl;
	break;
      case ERR_BEGIN_LIST_NOT_READ:
	cerr << "Beginning of list symbol "<< BEGIN_LIST <<" not read"<<endl;
	break;
      case ERR_SECTION_READ:
	cerr << "Premature beginning of SECTION command found" << endl;
	break;
      case ERR_SUBSECTION_READ:
	cerr << "Premature beginning of SUBSECTION command found" << endl;
	break;
      case ERR_SUBSUBSECTION_READ:
	cerr << "Premature beginning of SUBSUBSECTION command found" << endl;
	break;
      case ERR_NONNESTED_PARENTHESES:
	cerr << "Expression has nonnested parentheses" << endl;
	break;
      case ERR_ARITHMETIC_OP:
	cerr << "Arithmetic operation error in evaluating expression"<< endl;
	break;
      case ERR_UNITS_NEGATIVE:
	cerr << "Units expression is not positive" << endl;
	break;
      case ERR_UNITS_ADDITION_FOUND:
	cerr << "Binary addition or subtraction in units expression" << endl;
	break;
      case ERR_UNITS_WRONG_DIMEN:
	cerr << "Unit expression has the wrong dimensions" << endl;
	break;
      case ERR_UNITS_TOO_LONG:
	cerr << "Units expression is too long" << endl;
	break;
      case ERR_UNITS_IMPROPER_EXPR:
	cerr << "Units expression improperly formed or unknown unit" << endl;
	break;
      case ERR_UNKNOWN_UNIT:
	cerr << "Unknown unit" << endl;
	break;
      case ERR_UNITS:
	cerr << "Unit expression error" << endl;
	break;
      }
      if(inFileNo >= 0) {
	if(screenIn()) {
	  cerr << "Error in input on line number "
	       << inFileLineNo[inFileNo] << endl;
	} else {
	  cerr << "Error in input file " << inFileName[inFileNo]
	       << " on line number " << inFileLineNo[inFileNo] << endl;
	}
      }
      break;
    case ERR_BEFORE_LINE:
    case ERR_BAD_DATA:
    case ERR_OUT_OF_RANGE:
      switch(error) {
      case ERR_BAD_DATA:
	cerr << "Bad data detected" << endl;
	break;
      case ERR_OUT_OF_RANGE:
	cerr << "Data out of acceptable range" << endl;
	break;
      }
      if(screenIn()) {
	cerr << "Error in input on or before line number "
	     << inFileLineNo[inFileNo] << endl;
      } else {
	cerr << "Error in input file " << inFileName[inFileNo]
	     << " on or before line number " << inFileLineNo[inFileNo] << endl;
      }
      break;
    case ERR_MEMORY:
      cerr << "Error allocating memory" << endl;
      break;
    case ERR_OPENING_FILE:
      cerr << "Error opening input file " << inFileName[inFileNo+1] << endl;
      break;
    case ERR_FILESYSTEM:
      cerr << "Error opening file or directory" << endl;
      break;
    case ERR_NESTING_FILE:
      cerr << "Too many nested include files: maximum " << DEPTH_INCLUDES
	   << " allowed" << endl;
      break;
    case ERR_EOF:
      break;
    case ERR_IN_CODE:
      cerr << "Sorry, there is an error in the computer code" << endl;
      break;
    case ERR_SECTION_NOT_READ:
    case ERR_SUBSECTION_NOT_READ:
    case ERR_SUBSUBSECTION_NOT_READ:
    case ERR_SECTION_KEY_NOT_READ:
    case ERR_SUBSECTION_KEY_NOT_READ:
    case ERR_SUBSUBSECTION_KEY_NOT_READ:
    case ERR_KEY_NOT_READ:
    case ERR_END_LIST_NOT_READ:
      switch(error) {
      case ERR_SECTION_NOT_READ:
	cerr << "Beginning of SECTION command not found" << endl;
	break;
      case ERR_SUBSECTION_NOT_READ:
	cerr << "Beginning of SUBSECTION command not found" << endl;
	break;
      case ERR_SUBSUBSECTION_NOT_READ:
	cerr << "Beginning of SUBSUBSECTION command not found" << endl;
	break;
      case ERR_SECTION_KEY_NOT_READ:
      case ERR_SUBSECTION_KEY_NOT_READ:
      case ERR_SUBSUBSECTION_KEY_NOT_READ:
      case ERR_KEY_NOT_READ:
	switch(error) {
	case ERR_SECTION_KEY_NOT_READ:
	  cerr << "Section k";
	  break;
	case ERR_SUBSECTION_KEY_NOT_READ:
	  cerr << "Subsection k";
	  break;
	case ERR_SUBSUBSECTION_KEY_NOT_READ:
	  cerr << "Subsubsection k";
	  break;
	case ERR_KEY_NOT_READ:
	  cerr << "K";
	  break;
	}
	{
	  int k;
	  cerr << "ey not found: " << save_keyWord[0] << endl;
	  for(k=1; k<save_nKeyWords; k++)
	    cerr << "               " << save_keyWord[k] << endl;
	}
	break;
      case ERR_END_LIST_NOT_READ:
	cerr << "End of list symbol " << END_LIST << " not found" << endl;
	break;
      }
      if(skipping) {
	if(inFileNo >= 0) {
	  cerr << "Error between input file " << beginSkipFileName
	       << " line number " << beginSkipFileLineNo << endl;
	  cerr << "          and input file " << inFileName[inFileNo]
	       << " line number " << inFileLineNo[inFileNo] << endl;
	} else {
	  cerr << "Error after input file " << beginSkipFileName
	       << " line number " << beginSkipFileLineNo << endl;
	}
      } else {
	if(inFileNo >= 0) {
	  if(screenIn()) {
	    cerr << "Error in input on line number "
		 << inFileLineNo[inFileNo] << endl;
	  } else {
	    cerr << "Error in input file " << inFileName[inFileNo]
		 << " on line number " << inFileLineNo[inFileNo] << endl;
	  }
	}
      }
      break;
    case READ_LITERAL_CMD:
      break;
    default:
      cerr << "Unspecified error" << endl;
    }
    if(error == ERR_EOF || inFileNo < 0)
      cerr << "Premature end of data" << endl;
  }

  return error;
}
