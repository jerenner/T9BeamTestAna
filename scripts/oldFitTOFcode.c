
// ______________________________________________________
// old , obsolete
void SetParsSpecialMomenta(int p, double A, TF1 *func, double tof_el, double tof_mu, double tof_pi)
{
  // electrons:  norm [0], mean [2], sigma [3]
  // muons:      norm [1], mean [5], sigma [7]
  // pions:      norm [4], mean [6], sigma [7] -- same!
  
    
     // 200n, 220n, jk 22.11.2022, 24.11.2022
    if (fabs(p + 200) < 1.e-3 || fabs(p + 220) < 1.e-3 ) {
	func->SetParLimits(7, 0.1, 0.4);
      }
     
    // 200p, jk 22.11.2022
    if (fabs(p - 200) < 1.e-3) {
	func->SetParameter(0, A); 
	func->SetParLimits(0, 0.1*A, 2*A);
	func->SetParameter(1, A);
	func->SetParLimits(1, 0., 0.2*A);
	func->SetParLimits(4, 0., 0.1*A);
	func->SetParameter(5, 0.5*( tof_mu + tof_pi) );
	func->SetParameter(7, 0.9);
      }
    
    // 300n
    if (fabs(p + 300) < 1.e-3) {
	func->SetParameter(0, A); 
	func->SetParLimits(0, 0.1*A, 2*A);
	func->SetParameter(1, A);
	func->SetParLimits(1, 0.*A, 2*A);
	//func->SetParameter(5, 0.5*( tof_mu + tof_pi) );
	//func->SetParameter(5, 0.5*( tof_mu + tof_pi) );
	func->SetParameter(5, tof_mu + 0.2);
	cout << "tof_el=" << tof_el << " tof_mu=" << tof_mu << endl;
	func->SetParameter(7, 0.4);
	func->SetParLimits(7, 0.1, 0.5);
	// pi:
	func->SetParameter(6, tof_pi + 0.2);	
      }

    // 340n
    if (fabs(p + 340) < 1.e-3) {
	func->SetParameter(0, 1.5*A); 
	func->SetParLimits(0, 0.1*A, 1.5*A);
	func->SetParameter(1, A/4);
	func->SetParLimits(1, 0.*A, 1.3*A);
	func->SetParameter(6, 1.1*tof_mu);
      }

    // 340p
    if (fabs(p - 340) < 1.e-3) {
	func->SetParameter(0, 1.*A); 
	func->SetParLimits(0, 0.1*A, 1.5*A);
	func->SetParameter(1, A/4);
	func->SetParLimits(1, 0.*A, 1.*A);
	func->SetParameter(6, 1.1*tof_mu);
      }

    // 360n
    if (fabs(p + 360) < 1.e-3) {
	func->SetParameter(0, 1.*A); 
	func->SetParLimits(0, 0.1*A, 1.5*A);
	func->SetParameter(1, A/4);
	func->SetParLimits(1, 0.*A, 1.*A);
	func->SetParameter(6, 1.1*tof_mu);
      }

    // 360p
    if (fabs(p - 360) < 1.e-3) {
	func->SetParameter(0, 1.*A); 
	func->SetParLimits(0, 0.1*A, 1.5*A);
	func->SetParameter(1, A/4);
	func->SetParLimits(1, 0.*A, 1.*A);
	func->SetParameter(6, 1.1*tof_mu);
      }
}




  /*

  // OLD way/attempt by Jiri
	// integrate the ASC2 distribution below non-ele thr in the TOF time widow to get the ASC muon and pion ID efficiency
	// mu

	double nACT2_act2cut_mu = IntegrateAlongYInGivenXWindow(hACT2_act2cut, func->GetParameter(5), tof_reso);
	double nACT2_all_mu    = IntegrateAlongYInGivenXWindow(hACT2_all, func->GetParameter(5), tof_reso);
	double eff_mu          = nACT2_act2cut_mu / nACT2_all_mu;
  cout << " nACT2_act2cut_mu=" << nACT2_act2cut_mu
       << " nACT2_all_mu=" << nACT2_all_mu
       << endl;

	// pi

	double nACT3_act3cut_pi = IntegrateAlongYInGivenXWindow(hACT3_act3cut, func->GetParameter(6) + pioff, tof_reso);
	double nACT3_all_pi    = IntegrateAlongYInGivenXWindow(hACT3_all, func->GetParameter(6) + pioff, tof_reso);
	double eff_pi          = nACT3_act3cut_pi / nACT3_all_pi;
cout << " nACT3_act3cut_pi=" << nACT3_act3cut_pi
     << " nACT3_all_pi=" << nACT3_all_pi
     << endl;
*/


	// lines:
	/*
	x1 = func->GetParameter(5) - tof_reso;
 	x2 = func->GetParameter(5) + tof_reso;
	y1 = 	hACT3_act3cut -> GetYaxis() -> GetXmax();
	y2 = 	hACT3_act3cut -> GetYaxis() -> GetXmin();

	h2can -> cd(2);
	
	TLine *xlmu1 = new TLine(x1, y1, x1, y2);
	xlmu1 -> SetLineColor(kBlue+1);
	xlmu1 -> SetLineStyle(ls);
	xlmu1 -> SetLineWidth(lw);
	xlmu1 -> Draw();

	TLine *xlmu2 = new TLine(x2, y1, x2, y2);
	xlmu2 -> SetLineColor(kBlue+1);
	xlmu2 -> SetLineStyle(ls);
	xlmu2 -> SetLineWidth(lw);
	xlmu2 -> Draw();

	
	x1 = pioff+ func->GetParameter(6) - tof_reso;
 	x2 = pioff + func->GetParameter(6) + tof_reso;
	
	TLine *xlpi1 = new TLine(x1, y1, x1, y2);
	xlpi1 -> SetLineColor(kGreen+1);
	xlpi1 -> SetLineStyle(ls);
	xlpi1 -> SetLineWidth(lw);
	xlpi1 -> Draw();

	TLine *xlpi2 = new TLine(x2, y1, x2, y2);
	xlpi2 -> SetLineColor(kGreen+1);
	xlpi2 -> SetLineStyle(ls);
	xlpi2 -> SetLineWidth(lw);
	xlpi2 -> Draw();
	*/
	


  /*
      	// mu+pi
	double nACT2_act2cut_mupi = IntegrateAlongYInGivenXWindow(hACT2_act2cut, func->GetParameter(5), tof_reso);
	double nACT2_all_mupi    = IntegrateAlongYInGivenXWindow(hACT2_all, func->GetParameter(5), tof_reso);
	double eff_mupi          = nACT2_act2cut_mupi / nACT2_all_mupi;
      */
