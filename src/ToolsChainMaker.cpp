#define ChainMaker_cxx
#include "ToolsChainMaker.h" 

/////////////////////////
ChainMaker::ChainMaker()
////////////////////////
{
  //default Constructor
  m_ChainName="DUMMY";
}
/////////////////////////////////////////
ChainMaker::ChainMaker(TString ChainName)
/////////////////////////////////////////
{
  m_ChainName=ChainName;
}
/////////////////////////
ChainMaker::~ChainMaker()
/////////////////////////
{ 
  //destructor
  delete _theChain;
}
///////////////////////////////////////////
void ChainMaker::SetTreeFiles(std::string list)
///////////////////////////////////////////
{
  std::string MyList = list;
  TChain* chain = new TChain(m_ChainName);
  std::string runnumber;
  std::vector<std::string> fileList;
  for (size_t i=0,nchar; i <= MyList.length(); i=nchar+1){
    nchar = MyList.find_first_of(',',i);
    if (nchar == std::string::npos)
      nchar = MyList.length();
    std::string tmp = MyList.substr(i,nchar-i);
    fileList.push_back(tmp);
  }

  // open input files
  for (unsigned int iFile=0; iFile<fileList.size(); ++iFile) {
    chain->Add(fileList[iFile].c_str());
  }
  _theChain = chain;
  
  
}//void ChainMaker::SetTreeFiles(string list){
