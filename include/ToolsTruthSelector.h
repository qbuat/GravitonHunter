
#ifndef TruthSelector_h
#define TruthSelector_h

#include <iostream>
#include <vector>
#include <TTree.h>

class TruthSelector
{
 private:
  TTree* m_tree;
  std::vector<float>              *mc_pt;
  std::vector<float>              *mc_eta;
  std::vector<float>              *mc_phi;
  std::vector<float>              *mc_m;
  std::vector<int>                *mc_pdgId;
  std::vector<int>                *mc_status;
  std::vector< std::vector<int> > *mc_child_index;
  std::vector< std::vector<int> > *mc_parent_index;

  TBranch        *b_mc_pt;   //!
  TBranch        *b_mc_eta;   //!
  TBranch        *b_mc_phi;   //!
  TBranch        *b_mc_m;   //!
  TBranch        *b_mc_pdgId;   //!
  TBranch        *b_mc_status;   //!
  TBranch        *b_mc_child_index;   //!
  TBranch        *b_mc_parent_index;   //!

  void InitBranches();
  std::vector<int> FindChildIndex(int index, 
				  bool VERBOSE= false);
 public :

  TruthSelector( TTree* tree );
  bool GetEntry(int entry);
  virtual ~TruthSelector();
  std::vector<int> GetRSGravChildIndex(bool VERBOSE= false);
  std::vector<int> GetDirectGamGamChildIndex(bool VERBOSE = false);
  std::vector<int> GetGamJetChildIndex(bool VERBOSE = false);
  std::vector<int> GetFinalPhotonIndices(bool VERBOSE = false);
  void PrintParentChain(int index, bool VERBOSE = false);
  void Show(int entry,
	    bool sm_1_rs_0=true);
  float GetPartonETIso(int index, float cone_size, float parton_minpt_threshold);
  bool IsPartonIsolated(int index, float cone_size, float parton_minpt_threshold, float E_iso_cut);

  
};

#endif

