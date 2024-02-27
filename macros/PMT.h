//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Nov 26 18:20:05 2023 by ROOT version 6.28/04
// from TTree PMT/
// found on file: ../prodscripts/ntuple_000435.root
//////////////////////////////////////////////////////////

#ifndef PMT_h
#define PMT_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class PMT {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        Pedestal;
   Double_t        PedestalSigma;
   Int_t           nPeaks;
   Double_t        PeakVoltage[100];   //[nPeaks]
   Double_t        PeakTime[100];   //[nPeaks]
   Double_t        SignalTime[100];   //[nPeaks]
   Double_t        IntCharge[100];   //[nPeaks]
   UInt_t          timeStamp;
   UInt_t          triggerTime;

   // List of branches
   TBranch        *b_Pedestal;   //!
   TBranch        *b_PedestalSigma;   //!
   TBranch        *b_nPeaks;   //!
   TBranch        *b_PeakVoltage;   //!
   TBranch        *b_PeakTime;   //!
   TBranch        *b_SignalTime;   //!
   TBranch        *b_IntCharge;   //!
   TBranch        *b_timeStamp;   //!
   TBranch        *b_triggerTime;   //!

   PMT(TTree *tree=0);
   virtual ~PMT();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef PMT_cxx
PMT::PMT(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../prodscripts/ntuple_000435.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../prodscripts/ntuple_000435.root");
      }
      f->GetObject("PMT",tree);

   }
   Init(tree);
}

PMT::~PMT()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t PMT::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t PMT::LoadTree(Long64_t entry)
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

void PMT::Init(TTree *tree)
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
   fChain->SetBranchAddress("timeStamp", &timeStamp, &b_timeStamp);
   fChain->SetBranchAddress("triggerTime", &triggerTime, &b_triggerTime);
   Notify();
}

Bool_t PMT::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void PMT::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t PMT::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef PMT_cxx
