// This class provide a simple way to chain a list of files
// on the grid. It takes as input a string where 
// all the files are separated by coma (,)
//
#ifndef ChainMaker_h
#define ChainMaker_h

#include <TChain.h>
#include <TString.h>
#include <iostream>
#include <stdlib.h>

/* using namespace std; */

class ChainMaker {

 public:
  ChainMaker();
  ChainMaker( TString );
  ~ChainMaker();

  TChain* GetChain() {return _theChain;}
  void    SetTreeFiles(std::string list);
  
 private:
  TString m_ChainName;

 protected:
  TChain *_theChain;   

};

#endif 
