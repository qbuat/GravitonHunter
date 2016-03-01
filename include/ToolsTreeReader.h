/* $Id: TreeReader.h 246 2008-06-24 15:31:47Z /O=GRID-FR/C=FR/O=CNRS/OU=LPSC/CN=Benoit Clement $ */
#ifndef TREEREADER_H
#define TREEREADER_H

//////////////////////////////////////////////////
//
//		Class TreeReader
//		TreeReader.h
//
// Class for Tree reading through TFomula
//////////////////////////////////////////////////

#include "TTree.h"
#include "TTreeFormula.h"
#include "TTreeFormulaManager.h"
#include "TString.h"
#include <map>
#include <vector>
/*
class std::vector<double>;
class std::vector<int>;
class std::vector<std::vector<double> >;
class std::vector<std::vector<int> >;*/

//////////////////////////////////////////////////
class TreeReader //: public TTreeFormulaManager
{
 private:

  TTree* fTree;
  int              fCurrentEntry; 			// current ntuple entry stored in buffer
  int              fEntries;      			// total number of entries
  TString          fName;
  bool             fIsChain;
  int              fCurrentTree;
  std::map<std::string, TTreeFormula*>     fFormulae;	// known formulae
  TTreeFormulaManager* fManager;

 public:

  TreeReader();               // Default ctor
  TreeReader(TTree* n);       // ctor with ntuple
  virtual ~TreeReader();      // dtor

  void SetTree(TTree* n);       //
  double GetVariable(const char* c, int entry=-2); // return variable s for a given entry (<0 -> current entry)
  Bool_t GetEntry(int entry=-1);     // Read a given entry in the buffer (-1 -> next entry);
  int      GetEntries()             { return fEntries ; }
  TTree*   GetTree()                { return fTree    ; }
  TString  GetName()                { return fName    ; }
  void     SetName(TString name)    { fName = name    ; }              

  ClassDef(TreeReader,0);  // Integrate this class into ROOT (must be the last member)
};

#endif
