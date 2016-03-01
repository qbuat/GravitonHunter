#include "ToolsTruthSelector.h"
#include <TError.h>
#include <TLorentzVector.h>
/////////////////////////////////////////
TruthSelector::TruthSelector(TTree* tree) 
////////////////////////////////////////
{ 
  m_tree = (TTree*)tree->Clone( "clone") ;
  InitBranches();

}
///////////////////////////////
TruthSelector::~TruthSelector() 
///////////////////////////////
{ 

  delete mc_pt;
  delete mc_eta;
  delete mc_phi;
  delete mc_pdgId;
  delete mc_status;
  delete mc_child_index;
  delete mc_parent_index;
  m_tree = 0;
  // delete b_mc_pt;
  // delete b_mc_eta; 
  // delete b_mc_phi;
  // delete b_mc_pdgId;
  // delete b_mc_status;
  // delete b_mc_child_index;
  // delete b_mc_parent_index;
  // delete m_tree;

}


//////////////////////////////////
void TruthSelector::InitBranches() 
//////////////////////////////////
{
  mc_pt  = 0;
  mc_eta = 0;
  mc_phi = 0;

  mc_pdgId = 0;
  mc_status = 0;
  mc_child_index = 0;
  mc_parent_index = 0;

  m_tree->SetBranchStatus("*",0);
  m_tree->SetBranchStatus("mc_pt",1);
  m_tree->SetBranchStatus("mc_eta",1);
  m_tree->SetBranchStatus("mc_phi",1);
  m_tree->SetBranchStatus("mc_pdgId",1);
  m_tree->SetBranchStatus("mc_status",1);
  m_tree->SetBranchStatus("mc_child_index",1);
  m_tree->SetBranchStatus("mc_parent_index",1);

  m_tree->SetBranchAddress("mc_pt", &mc_pt, &b_mc_pt);
  m_tree->SetBranchAddress("mc_eta", &mc_eta, &b_mc_eta);
  m_tree->SetBranchAddress("mc_phi", &mc_phi, &b_mc_phi);
  m_tree->SetBranchAddress("mc_pdgId", &mc_pdgId, &b_mc_pdgId);
  m_tree->SetBranchAddress("mc_status", &mc_status, &b_mc_status);
  m_tree->SetBranchAddress("mc_child_index", &mc_child_index, &b_mc_child_index);
  m_tree->SetBranchAddress("mc_parent_index", &mc_parent_index, &b_mc_parent_index);


}
//////////////////////////////////////
bool TruthSelector::GetEntry(int entry) 
///////////////////////////////////////
{
  bool temp;
  int test = m_tree->GetEntry(entry);
  if(test>0) temp = true; 
  else{
    temp = false;
    Fatal("TruthSelector::GetEntry","CANNOT GET THE ENTRY !");
  }

  return temp;

}
/////////////////////////////////////////////////////////////////
std::vector<int> TruthSelector::GetRSGravChildIndex(bool VERBOSE)
/////////////////////////////////////////////////////////////////
{
  //--------------- Select the graviton ---------------------------
  std::vector<int> grav_index;
  for( int imc=0 ; imc< (int)(*mc_pdgId).size();imc++){
   if( (*mc_pdgId)[imc] != 5100039) continue;
   if( (*mc_status)[imc] != 62 ) continue;
   if(VERBOSE)
     std::cout << imc << " " << (*mc_pdgId)[imc] << std::endl;
    grav_index.push_back(imc);
  }
  //-----------------------------------------------------------------
  if( grav_index.size() !=1 )
    Fatal("TruthSelector::GetRSGravChildIndex", "MORE THAN ONE GRAVITON !!");


  //--------------- Select the children as photon of interest---------------------------
  std::vector<int> child_index;
  if(VERBOSE)
    std::cout << "Grav Index " << grav_index[0] << std::endl;
  for( int ichild=0; ichild< (int)(*mc_child_index)[grav_index[0]].size();ichild++){
    if(VERBOSE)
      std::cout  << "ichild : "<< ichild << std::endl;
    child_index.push_back( (*mc_child_index)[grav_index[0]][ichild] );
   }
  //-----------------------------------------------------------------------------------
  return child_index;
}
///////////////////////////////////////////////////////////////////////
std::vector<int> TruthSelector::GetDirectGamGamChildIndex(bool VERBOSE)
///////////////////////////////////////////////////////////////////////
{
  std::vector<int> DirectGamGamChildIndex;

  //-------------- Some Print-out ----------------------------------------------------
  if(VERBOSE) {
    for( int imc=0 ; imc< (int)(*mc_pdgId).size();imc++){
      if( (*mc_status)[imc] != 21 ) continue;
      std::cout << "-----------------------------" << std::endl;
      std::cout << "parton "<< imc << " --> Id = "<< (*mc_pdgId)[imc] << std::endl;
      for( int ichild = 0 ; ichild< (int)(*mc_child_index)[imc].size();ichild++)
	std::cout << "child "<< ichild  
		  << " : --> Index = " <<  (*mc_child_index)[imc][ichild] 
		  << " Id = " << (*mc_pdgId)[ (*mc_child_index)[imc][ichild] ] 
		  << " status = " << (*mc_status)[ (*mc_child_index)[imc][ichild] ] 
		  << std::endl;
    }
      std::cout << "-----------------------------" << std::endl;
  } 
  //-------------------------------------------------------------------------------

  //--------------- Find the two photon from hardproc direct children from hp glu or qua ------
  std::vector<int> firstchild;
  for( int imc=0 ; imc< (int)(*mc_pdgId).size();imc++){
    if( (*mc_status)[imc] != 21 ) continue;
    for( int ichild = 0 ; ichild< (int)(*mc_child_index)[imc].size();ichild++){
      bool is_already_in_vec = false;
      for(int ic=0; ic<(int)firstchild.size();ic++){
	if( firstchild[ic]== (*mc_child_index)[imc][ichild]){
	  is_already_in_vec=true;
	  break;
	}
      }
      if(!is_already_in_vec)
	firstchild.push_back( (*mc_child_index)[imc][ichild] );
    }
  }
  if( firstchild.size() !=2 ) 
    Fatal( "TruthSelector::GetDirectGamGamChildIndex()",
	   "WRONG NUMBER OF FIRST CHILDREN" );
  //----------------------------------------------------------------------------------------    

  //--------- Find the second child of each photon (start of the virtual chain --------------
  std::vector<int> secondchild0 = FindChildIndex( firstchild[0],VERBOSE );
  std::vector<int> secondchild1 = FindChildIndex( firstchild[1],VERBOSE );
  //-------------------------------------------------------------------------------
  if( VERBOSE ){
    std::cout << "First child chain " << std::endl;
    std::cout << "Index = "  << firstchild[0]
	      << " pt = "    << (*mc_pt)[firstchild[0]]
	      << " eta = "   << (*mc_eta)[firstchild[0]]
	      << " phi = "   << (*mc_phi)[firstchild[0]] 
	      << " pdgId = " << (*mc_pdgId)[firstchild[0]] 
	      << " status = "<< (*mc_status)[firstchild[0]] 
	      << std::endl;
    std::cout << "Index = "  << secondchild0[0]
	      << " pt = "    << (*mc_pt)[secondchild0[0]]
	      << " eta = "   << (*mc_eta)[secondchild0[0]]
	      << " phi = "   << (*mc_phi)[secondchild0[0]] 
	      << " pdgId = " << (*mc_pdgId)[secondchild0[0]] 
	      << " status = "<< (*mc_status)[secondchild0[0]] 
	      << std::endl;
  }


  //------------- Iteration over the children of photon 0 up to the last one ----------------------
  while( secondchild0.size() != 0) {
    if( secondchild0.size() == 0 ) break;
    std::vector<int> temp = FindChildIndex(secondchild0[0]);
    if( temp.size() != 1 ) break;
    else secondchild0 = temp;
  
    if(VERBOSE) 
      std::cout << "Index = "  << secondchild0[0]
		<< " pt = "    << (*mc_pt)[secondchild0[0]]
		<< " eta = "   << (*mc_eta)[secondchild0[0]]
		<< " phi = "   << (*mc_phi)[secondchild0[0]] 
		<< " pdgId = " << (*mc_pdgId)[secondchild0[0]] 
		<< " status = "<< (*mc_status)[secondchild0[0]]  
		<< std::endl;
  }
  DirectGamGamChildIndex.push_back( secondchild0[0] );
  //-------------------------------------------------------------------------------

  if(VERBOSE){
    std::cout << "Second child chain " << std::endl;
    std::cout << "Index = "  << firstchild[1]
	      << " pt = "    << (*mc_pt)[firstchild[1]]
	      << " eta = "   << (*mc_eta)[firstchild[1]]
	      << " phi = "   << (*mc_phi)[firstchild[1]] 
	      << " pdgId = " << (*mc_pdgId)[firstchild[1]] 
	      << " status = "<< (*mc_status)[firstchild[1]] 
	      << std::endl;
    std::cout << "Index = "  << secondchild1[0]
	      << " pt = "    << (*mc_pt)[secondchild1[0]]
	      << " eta = "   << (*mc_eta)[secondchild1[0]]
	      << " phi = "   << (*mc_phi)[secondchild1[0]] 
	      << " pdgId = " << (*mc_pdgId)[secondchild1[0]] 
	      << " status = "<< (*mc_status)[secondchild1[0]] 
	      << std::endl;
  }

  //------------- Iteration over the children of photon 1 up to the last one ----------------------
  while( secondchild1.size() != 0) {
    if( secondchild1.size() == 0 ) break;
    std::vector<int> temp= FindChildIndex(secondchild1[0]);
    if( temp.size() != 1) break ;
    else secondchild1 = temp;

    if(VERBOSE)
      std::cout << "Index = "  << secondchild1[0]
		<< " pt = "    << (*mc_pt)[secondchild1[0]]
		<< " eta = "   << (*mc_eta)[secondchild1[0]]
		<< " phi = "   << (*mc_phi)[secondchild1[0]] 
		<< " pdgId = " << (*mc_pdgId)[secondchild1[0]] 
		<< " status = "<< (*mc_status)[secondchild1[0]] 
		<< std::endl;
  }
  DirectGamGamChildIndex.push_back( secondchild1[0] );
  //-------------------------------------------------------------------------------

  return DirectGamGamChildIndex;
}

///////////////////////////////////////////////////////////////////////
std::vector<int> TruthSelector::FindChildIndex(int index,bool VERBOSE)
//////////////////////////////////////////////////////////////////////
{
  std::vector<int> children;
  for( int ichild = 0 ; ichild< (int)(*mc_child_index)[index].size();ichild++){
    if( VERBOSE ){
      std::cout << "child of part with index " << index   
		<< " : --> Index = " <<  (*mc_child_index)[index][ichild] 
		<< " Id = " << (*mc_pdgId)[ (*mc_child_index)[index][ichild] ] 
		<< " status = " << (*mc_status)[ (*mc_child_index)[index][ichild] ] 
		<< std::endl;
    }
    children.push_back( (*mc_child_index)[index][ichild] );
  }
  return children;

}

//////////////////////////////////////////////////
void TruthSelector::Show(int entry,bool sm_1_rs_0)
//////////////////////////////////////////////////
{
  m_tree->GetEntry(entry);
  std::vector<int> indices;
  std::cout << "mc_pt : mc_eta : mc_phi : mc_status " << std::endl;
  for(int imc=0;imc<(int)(*mc_pt).size();imc++){
    if( (*mc_pdgId)[imc] != 22 ) continue;
    if( (*mc_status)[imc] != 1 ) continue;
    if( (*mc_pt)[imc] < 15000. ) continue;
    std::cout << (*mc_pt)[imc] 
	      << " "
	      << (*mc_eta)[imc] 
	      << " "
	      << (*mc_phi)[imc] 
	      << " "
	      << (*mc_status)[imc]
	      << std::endl;
    indices.push_back(imc);
  }
  if(indices.size()>1){
    TLorentzVector p1_old, p2_old;
    p1_old.SetPtEtaPhiM((*mc_pt)[indices[0]],(*mc_eta)[indices[0]],(*mc_phi)[indices[0]],0);
    p2_old.SetPtEtaPhiM((*mc_pt)[indices[1]],(*mc_eta)[indices[1]],(*mc_phi)[indices[1]],0);
    TLorentzVector sys_old = p1_old+p2_old;
    std::cout << "Mass(GeV) = " << sys_old.M()/1000. << std::endl;
  }

  std::vector<int> vec;
  if(sm_1_rs_0)
    vec = GetDirectGamGamChildIndex();
  else 
    vec = GetRSGravChildIndex();

  if( vec.size() >1 ){
    std::cout << "************************************" << std::endl;
    for(int i=0 ; i<(int)vec.size();i++)
      std::cout << (*mc_pt)[vec[i]] 
	      << " "
		<< (*mc_eta)[vec[i]] 
		<< " "
		<< (*mc_phi)[vec[i]] 
		<< " "
		<< (*mc_status)[vec[i]]
		<< std::endl;
    TLorentzVector p1, p2;
    p1.SetPtEtaPhiM((*mc_pt)[vec[0]],(*mc_eta)[vec[0]],(*mc_phi)[vec[0]],0);
    p2.SetPtEtaPhiM((*mc_pt)[vec[1]],(*mc_eta)[vec[1]],(*mc_phi)[vec[1]],0);
    TLorentzVector sys = p1+p2;
    std::cout << "Mass (GeV) = " << sys.M()/1000. << std::endl;
  }

}




///////////////////////////////////////////////////////////////////////
std::vector<int> TruthSelector::GetFinalPhotonIndices(bool VERBOSE)
//////////////////////////////////////////////////////////////////////
{
  std::vector<int> v_indices;
  for(unsigned int imc = 0; imc < (*mc_pdgId).size(); imc++)
    {
      if((*mc_pdgId)[imc] != 22) continue;
      if((*mc_status)[imc] != 1) continue;
      v_indices.push_back(imc);
    }
  return v_indices;

}

///////////////////////////////////////////////////////////////////////
void TruthSelector::PrintParentChain(int index, bool VERBOSE)
///////////////////////////////////////////////////////////////////////
{

  //print current particle
  if(VERBOSE)
    std::cout<<"particle " << index << " : "
	     <<" status " << (*mc_status)[index]
	     <<" pdgId "  << (*mc_pdgId)[index]
	     <<" pt "     << (*mc_pt)[index]
	     <<" eta "    << (*mc_eta)[index]
	     <<" phi "    << (*mc_phi)[index]
	     << std::endl;

  int tmp_index = index;
  bool isfirstchild = false;

  while(!isfirstchild) {
    int n_parents = (*mc_parent_index)[tmp_index].size();
    if(n_parents == 1) {
      int parentindex = (*mc_parent_index)[tmp_index][0];
        if(VERBOSE)
	  std::cout<<"particle " << parentindex << " : " 
		   <<" status " << (*mc_status)[parentindex]
		   <<" pdgId "  << (*mc_pdgId)[parentindex]
		   <<" pt "     << (*mc_pt)[parentindex]
		   <<" eta "    << (*mc_eta)[parentindex]
		   <<" phi "    << (*mc_phi)[parentindex]
		   <<std::endl;	
      tmp_index = parentindex;
    }
    
    else if(n_parents > 1) {
      if(VERBOSE) {
	std::cout<<"found " << n_parents << " incoming particles with status "<<std::endl;
	for(int p = 0; p < n_parents; p++)
	  std::cout<<(*mc_pdgId)[p]<<", ";
	std::cout<<std::endl;
      }
      isfirstchild = true;
    }

    else if(n_parents== 0) {
      if(VERBOSE) std::cout<<"this photon has no parents!!!!"<<std::endl;
      break; 
    }
    
  }
  
}


///////////////////////////////////////////////////////////////////////
std::vector<int> TruthSelector::GetGamJetChildIndex(bool VERBOSE)
///////////////////////////////////////////////////////////////////////
{

  std::vector<int> indices;
  for(unsigned int imc = 0; imc < (*mc_pdgId).size(); imc++)
    {
      if((*mc_pdgId)[imc] != 22) continue;
      if((*mc_status)[imc] != 1) continue;
      //      indices.push_back(imc);
      
      //print current particle
      if(VERBOSE)
	std::cout<<"particle " << imc << " : "
		 <<" status " << (*mc_status)[imc]
		 <<" pdgId "  << (*mc_pdgId)[imc]
		 <<" pt "     << (*mc_pt)[imc]
		 <<" eta "    << (*mc_eta)[imc]
		 <<" phi "    << (*mc_phi)[imc]
		 <<std::endl;
      
      int tmp_index = imc;
      bool isfirstchild = false;
      
      while(!isfirstchild) {
	int n_parents = (*mc_parent_index)[tmp_index].size();
	if(n_parents == 1) {
	  int parentindex = (*mc_parent_index)[tmp_index][0];
	  if(VERBOSE)
	    std::cout<<"particle " << parentindex << " : " 
		     <<" status " << (*mc_status)[parentindex]
		     <<" pdgId "  << (*mc_pdgId)[parentindex]
		     <<" pt "     << (*mc_pt)[parentindex]
		     <<" eta "    << (*mc_eta)[parentindex]
		     <<" phi "    << (*mc_phi)[parentindex]
		     <<std::endl;	
	  tmp_index = parentindex;
	}
	
	else if(n_parents > 1) {
	  if(VERBOSE) {
	    std::cout<<"found " << n_parents << " incoming particles with status "<<std::endl;
	    for(int p = 0; p < n_parents; p++)
	      std::cout<<(*mc_pdgId)[p]<<", ";
	    std::cout<<std::endl;
	  }
	  isfirstchild = true;
	  indices.push_back(imc);
	}
	
	else if(n_parents== 0) {
	  if(VERBOSE) std::cout<<"this photon has no parents!!!!"<<std::endl;
	  break; 
	}
	
      }
      
    }
  
  return indices;
  
}







//////////////////////////////////////
bool TruthSelector::IsPartonIsolated(int index, float cone_size, float parton_minpt_threshold, float E_iso_cut)
///////////////////////////////////////
{

  return ( GetPartonETIso( index, cone_size, parton_minpt_threshold ) < E_iso_cut );

}





//////////////////////////////////////
float TruthSelector::GetPartonETIso(int index, float cone_size, float parton_minpt_threshold)
///////////////////////////////////////
{

  std::cout<<"AAAAAAAAAAA"<<std::endl;
  
  float Eiso = 0;
  TLorentzVector g;
  std::cout<<"BBBBBBBBBBB"<<std::endl;
  g.SetPtEtaPhiM((*mc_pt)[index],(*mc_eta)[index],(*mc_phi)[index],(*mc_m)[index]);

  for(int i = 0; i < (int)(*mc_pdgId).size(); i++) {
    if((*mc_status)[i] != 1) continue;
    if(abs((*mc_pdgId)[i]) == 12 ||
       abs((*mc_pdgId)[i]) == 14 ||
       abs((*mc_pdgId)[i]) == 16    ) continue;  // exclude neutrino
    if(i == index) continue;
    
    TLorentzVector p;
    p.SetPtEtaPhiM((*mc_pt)[i],(*mc_eta)[i],(*mc_phi)[i],(*mc_m)[i]);
    
    if(p.Pt() < parton_minpt_threshold) continue;
    
    if(g.DeltaR(p) < cone_size) {
      float k = 1;
      if( abs((*mc_pdgId)[i]) == 15) k=0.1783 * (1+0.1);
      else if( abs((*mc_pdgId)[i]) == 13) k=0.1;
      Eiso += k*p.Pt();
    }
    
  }
  
  std::cout<<"before return in GetPartonETIso"<<std::endl;
  return Eiso;
  
  
}
