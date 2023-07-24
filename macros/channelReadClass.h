//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Jul 23 20:18:36 2023 by ROOT version 6.24/06
// from TTree ACT0L/
// found on file: ntuple_000250.root
//////////////////////////////////////////////////////////
// according to new ntuple multi tree format by Nick Prouse
// modified by Jiri Kvita
// modified on 24/07/2023


#ifndef channelReadClass_h
#define channelReadClass_h

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"


// ______________________________________________________________

const int maxnPeaks = 30; // orig: 17

// Header file for the classes stored in the TTree if any.

// ______________________________________________________________
// ______________________________________________________________
// ______________________________________________________________


class channelReadClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        Pedestal;
   Double_t        PedestalSigma;
   Long64_t        nPeaks;
   Double_t        PeakVoltage[maxnPeaks];   //[nPeakVoltage]
   Double_t        PeakTime[maxnPeaks];   //[nPeakTime]
   Double_t        SignalTime[maxnPeaks];   //[nSignalTime]
   Double_t        IntCharge[maxnPeaks];   //[nIntCharge]

   // List of branches
   TBranch        *b_Pedestal;   //!
   TBranch        *b_PedestalSigma;   //!
   TBranch        *b_nPeaks;   //!
   TBranch        *b_PeakVoltage;   //!
   TBranch        *b_PeakTime;   //!
   TBranch        *b_SignalTime;   //!
   TBranch        *b_IntCharge;   //!

   channelReadClass(TFile *infile, TString treeName);
   virtual ~channelReadClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

// ______________________________________________________________
channelReadClass::channelReadClass(TFile *infile, TString treeName) : fChain(0) 
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

// ______________________________________________________________

channelReadClass::~channelReadClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

// ______________________________________________________________
Int_t channelReadClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

// ______________________________________________________________
Long64_t channelReadClass::LoadTree(Long64_t entry)
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

// ______________________________________________________________
void channelReadClass::Init(TTree *tree)
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

   fChain->SetBranchAddress("Pedestal", &Pedestal, &b_Pedestal);
   fChain->SetBranchAddress("PedestalSigma", &PedestalSigma, &b_PedestalSigma);
   fChain->SetBranchAddress("nPeaks", &nPeaks, &b_nPeaks);
   fChain->SetBranchAddress("PeakVoltage", PeakVoltage, &b_PeakVoltage);
   fChain->SetBranchAddress("PeakTime", PeakTime, &b_PeakTime);
   fChain->SetBranchAddress("SignalTime", SignalTime, &b_SignalTime);
   fChain->SetBranchAddress("IntCharge", IntCharge, &b_IntCharge);
   Notify();
}

// ______________________________________________________________
Bool_t channelReadClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

// ______________________________________________________________
void channelReadClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

// ______________________________________________________________
Int_t channelReadClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void  channelReadClass::Loop() {

}

// ______________________________________________________________
// ______________________________________________________________
// ______________________________________________________________


#endif // #ifdef channelReadClass_h
