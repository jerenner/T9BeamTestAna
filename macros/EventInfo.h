//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jul 24 20:45:44 2023 by ROOT version 6.24/06
// from TTree EventInfo/
// found on file: ntuple_000373.root
//////////////////////////////////////////////////////////

#ifndef EventInfo_h
#define EventInfo_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class EventInfo {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Long64_t        RunNumber;
   Int_t           EventNumber;
   Int_t           SpillNumber;

   // List of branches
   TBranch        *b_RunNumber;   //!
   TBranch        *b_EventNumber;   //!
   TBranch        *b_SpillNumber;   //!

   EventInfo(TFile *infile, TString treeName);
   virtual ~EventInfo();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};


EventInfo::EventInfo(TFile *infile, TString treeName) : fChain(0) 
{
 if (!infile || infile -> IsZombie()) {
    cout << "ERROR, got null TFile pointer!" << endl;
      return;
  }
  TTree* tree = (TTree*) infile -> Get(treeName);
   if (tree == 0) {
     cout << "ERROR getting the tree named " << treeName.Data() << " !" << endl;
     return;
   }
   Init(tree);
}

EventInfo::~EventInfo()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t EventInfo::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t EventInfo::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void EventInfo::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("SpillNumber", &SpillNumber, &b_SpillNumber);
   Notify();
}

Bool_t EventInfo::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void EventInfo::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t EventInfo::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
void EventInfo::Loop()
{
}



#endif // #ifdef EventInfo_cxx
