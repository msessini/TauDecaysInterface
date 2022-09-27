#include "HiggsCPinTauDecays/TauDecaysInterface/interface/SCalculator.h"
#include <iostream>

SCalculator::SCalculator(){
}

SCalculator::SCalculator(string type):
  type_(type)
{
}

void 
SCalculator::Configure(vector<TLorentzVector> TauAndProd, TLorentzVector Frame, int charge){
  for(unsigned int i=0; i<TauAndProd.size(); i++){
    TauAndProd_HRF.push_back(Boost(TauAndProd.at(i), Frame));

  }
  charge_=charge;
}

bool
SCalculator::isConfigured(){
  if(TauAndProd_LF.size()!=2){ std::cout<<"Error:   SCalculator is not Configured! Check  the size of input vector!  Size =  "<< TauAndProd_LF.size() <<std::endl; return false;} return true;
}

SCalculator::~SCalculator(){
}

TLorentzVector 
SCalculator::Boost(TLorentzVector pB, TLorentzVector frame){
  TMatrixT<double> transform(4,4);
  TMatrixT<double> result(4,1);
  TVectorT<double> vec(4); 
  TVector3 b;
  
  if(frame.Vect().Mag()==0){ std::cout<<"RH Boost is not set, perfrom calculation in the Lab Frame   "<<std::endl; return pB;}
  if(frame.E()==0){ std::cout<<" Caution: Please check that you perform boost correctly!  " <<std::endl; return pB;} 
  else   b=frame.Vect()*(1/frame.E());
  
  
  vec(0)  = pB.E();    vec(1)  = pB.Px();
  vec(2)  = pB.Py();  vec(3)  = pB.Pz();

  double gamma  = 1/sqrt( 1 - b.Mag2());
 
  transform(0,0)=gamma; transform(0,1) =- gamma*b.X() ;  transform(0,2) =  - gamma*b.Y();  transform(0,3) = - gamma*b.Z(); 
 
  transform(1,0)=-gamma*b.X(); transform(1,1) =(1+ (gamma-1)*b.X()*b.X()/b.Mag2()) ;  transform(1,2) = ((gamma-1)*b.X()*b.Y()/b.Mag2());  transform(1,3) = ((gamma-1)*b.X()*b.Z()/b.Mag2());
  transform(2,0)=-gamma*b.Y(); transform(2,1) = ((gamma-1)*b.Y()*b.X()/b.Mag2());  transform(2,2) = (1 + (gamma-1)*b.Y()*b.Y()/b.Mag2());  transform(2,3) =  ((gamma-1)*b.Y()*b.Z()/b.Mag2()); 
  transform(3,0)=-gamma*b.Z(); transform(3,1) =((gamma-1)*b.Z()*b.X()/b.Mag2()) ;  transform(3,2) = ((gamma-1)*b.Z()*b.Y()/b.Mag2());  transform(3,3) = (1 + (gamma-1)*b.Z()*b.Z()/b.Mag2()); 
  result=transform*convertToMatrix(vec);
  return TLorentzVector(result(1,0), result(2,0) ,result(3,0), result(0,0));
}
TMatrixT<double> SCalculator::convertToMatrix(TVectorT<double> V){
  TMatrixT<double> M(V.GetNrows(),1);
  for(int i=0; i < M.GetNrows(); i++){
    M(i,0)=V(i);
  } return M;
}

TVector3
SCalculator::Rotate(TVector3 LVec, TVector3 Rot){
  TVector3 vec = LVec;
  vec.RotateZ(0.5*TMath::Pi() - Rot.Phi());  // not 0.5, to avoid warnings about 0 pT
  vec.RotateX(Rot.Theta());
  return vec;
}

TVector3
SCalculator::pv(){
  TVector3 out(0,0,0);
  if(type_=="pion") out = TauAndProd_HRF.at(1).Vect();
  if(type_=="rho"){
    TLorentzVector pi  = TauAndProd_HRF.at(1);
    TLorentzVector pi0 = TauAndProd_HRF.at(2);
    TLorentzVector Tau = TauAndProd_HRF.at(0);
    TLorentzVector q= pi  - pi0;
    TLorentzVector P= Tau;
    TLorentzVector N= Tau - pi - pi0;
    out = P.M()*(2*(q*N)*q.Vect() - q.Mag2()*N.Vect()) * (1/ (2*(q*N)*(q*P) - q.Mag2()*(N*P)));
  }
  if(type_=="a1"){
    PolarimetricA1  a1pol;
    a1pol.Configure(TauAndProd_HRF,charge_);
    out = -a1pol.PVC().Vect();
  }
  return out;
}

void SCalculator::SortPions(std::vector<TLorentzVector>& pionsvec, std::vector<double>& charges)
{
  
  int npim(0),npip(0), npin(0);
  int OSMCPionIndex(0);
  int SSMCPion1Index(0);
  int SSMCPion2Index(0);
  int OSCharge(0);
  int SS1Charge(0);
  int SS2Charge(0);
 
  TLorentzVector os;
  TLorentzVector ss1;
  TLorentzVector ss2;

  int MCNeutralPionIndex(0);
  int MCChargedPionIndex(0);
  int NeutralPionCharge(0);
  int ChargedPionCharge(0);

  TLorentzVector ChargedPion;
  TLorentzVector NeutralPion;
  //cout<<"test 0"<<endl;
  if(charges.size()==3) //A1
    {
      //cout<<"test 1"<<endl;
      for(unsigned int i=0; i<charges.size(); i++){
	if( charges.at(i)== 1) npip++;
	if( charges.at(i)==-1) npim++;
      }
      if(npip == 1 && npim == 2){
	//cout<<"test 2"<<endl;
	int nss=0;
	for(unsigned int i=0; i<charges.size(); i++){
	  if(charges.at(i)== 1){
	    OSCharge=1;
	    OSMCPionIndex=i;
	  }
	  if(charges.at(i)== -1 && nss == 0){
	    nss++;
	    SS1Charge=-1;
	    SSMCPion1Index=i;
	  }
	  if(charges.at(i)== -1 && nss == 1){
	    SS2Charge=-1;
	    SSMCPion2Index=i;
	  }
	}
      }
      if( npip== 2 && npim==1){
	//cout<<"test 3"<<endl;
	int nss=0;
	for(unsigned int i=0; i<charges.size(); i++){
	  if(charges.at(i) == -1){
	    OSCharge=-1;
	    OSMCPionIndex=i;
	  }
	  if(charges.at(i) == 1 && nss ==0){
	    nss++;
	    SS1Charge=1;
	    SSMCPion1Index=i;
	  }
	  if(charges.at(i) == 1 && nss == 1){
	    SS2Charge=1;
	    SSMCPion2Index=i;
	  }
	}
      }
      os=pionsvec.at(OSMCPionIndex);
      ss1=pionsvec.at(SSMCPion1Index);
      ss2=pionsvec.at(SSMCPion2Index);
      
      charges.clear();
      charges.push_back(OSCharge);
      charges.push_back(SS1Charge);
      charges.push_back(SS2Charge);
      
      pionsvec.clear();
      pionsvec.push_back(os);
      pionsvec.push_back(ss1);
      pionsvec.push_back(ss2);
    }
  
  if(charges.size()==2) //Rho
    {
      for(unsigned int i=0; i<charges.size(); i++){
	if( charges.at(i)== 1) npip++;
	if( charges.at(i)==-1) npim++;
	if( charges.at(i)==0) npin++;
      }
      if(npip == 1 && npin == 1){
	for(unsigned int i=0; i<charges.size(); i++){
	  if(charges.at(i)== 1){
	    ChargedPionCharge=1;
	    MCChargedPionIndex=i;
	  }
	  if(charges.at(i)== 0){
	    NeutralPionCharge=0;
	    MCNeutralPionIndex=i;
	  }
	}
      }
      if( npim== 1 && npin==1){
	for(unsigned int i=0; i<charges.size(); i++){
	  if(charges.at(i) == -1){
	    ChargedPionCharge=-1;
	    MCChargedPionIndex=i;
	  }
	  if(charges.at(i) == 0){
	    NeutralPionCharge=0;
	    MCNeutralPionIndex=i;
	  }
	}
      }
      
      ChargedPion=pionsvec.at(MCChargedPionIndex);
      NeutralPion=pionsvec.at(MCNeutralPionIndex);
      
      charges.clear();
      charges.push_back(ChargedPionCharge);
      charges.push_back(NeutralPionCharge);
      
      pionsvec.clear();
      pionsvec.push_back(ChargedPion);
      pionsvec.push_back(NeutralPion);
    }
  //return pionsvec;
}

bool SCalculator::isOk(TString type1, TString type2, TLorentzVector tauMinus, std::vector<TLorentzVector> sumPionsMinus, std::vector<double> sumPionsChargeMinus, TLorentzVector tauPlus, std::vector<TLorentzVector> sumPionsPlus, std::vector<double> sumPionsChargePlus)
{
  SCalculator Scalc1(type1.Data());
  SCalculator Scalc2(type2.Data());
  
  TLorentzVector zeroLV(0,0,0,0);
  TLorentzVector HadLVMinus(0,0,0,0);
  TLorentzVector HadLVPlus(0,0,0,0);
  
  if(type1!="pion") Scalc1.SortPions(sumPionsMinus, sumPionsChargeMinus);
  if(type2!="pion") Scalc2.SortPions(sumPionsPlus, sumPionsChargePlus);
  //cout<<"----------"<<endl;
  //if(type1=="pion" && type2!="a1")
  vector<TLorentzVector> tauandprodminus;
  vector<TLorentzVector> tauandprodplus;
  bool pionszero=false;
  tauandprodminus.push_back(tauMinus);
  for(unsigned int i=0; i<sumPionsMinus.size();i++) {tauandprodminus.push_back(sumPionsMinus.at(i))/*;sumPionsMinus.at(i).Print()*/;if (sumPionsMinus.at(i)==zeroLV)pionszero=true;}
  
  tauandprodplus.push_back(tauPlus); 
  for(unsigned int i=0; i<sumPionsPlus.size();i++) {tauandprodplus.push_back(sumPionsPlus.at(i))/*;sumPionsPlus.at(i).Print()*/;if (sumPionsPlus.at(i)==zeroLV)pionszero=true;}   
  
  Scalc1.Configure(tauandprodminus,tauandprodminus.at(0)+tauandprodplus.at(0), -1);
  TVector3 h1=Scalc1.pv();
  
  Scalc2.Configure(tauandprodplus,tauandprodminus.at(0)+tauandprodplus.at(0), +1);
  TVector3 h2=Scalc2.pv();
  
  if(std::isnan(h1.Mag())==true || std::isnan(h2.Mag())==true || tauMinus==zeroLV || tauPlus==zeroLV || tauMinus==tauPlus || pionszero ||sumPionsPlus==sumPionsMinus) return false;
  else return true;
}
double SCalculator::AcopAngle(TString type1, TString type2, TLorentzVector tauMinus, std::vector<TLorentzVector> sumPionsMinus, std::vector<double> sumPionsChargeMinus, TLorentzVector tauPlus, std::vector<TLorentzVector> sumPionsPlus, std::vector<double> sumPionsChargePlus)
{  
  
  SCalculator Scalc1(type1.Data());
  SCalculator Scalc2(type2.Data());

  TLorentzVector zeroLV(0,0,0,0);
  TLorentzVector HadLVMinus(0,0,0,0);
  TLorentzVector HadLVPlus(0,0,0,0);

  if(type1!="pion") Scalc1.SortPions(sumPionsMinus, sumPionsChargeMinus);
  if(type2!="pion") Scalc2.SortPions(sumPionsPlus, sumPionsChargePlus);
  
  //if(type1=="pion" && type2!="a1")
  vector<TLorentzVector> tauandprodminus;
  vector<TLorentzVector> tauandprodplus;
  
  tauandprodminus.push_back(tauMinus);
  for(unsigned int i=0; i<sumPionsMinus.size();i++) {tauandprodminus.push_back(sumPionsMinus.at(i));}
  
  
  tauandprodplus.push_back(tauPlus); 
  for(unsigned int i=0; i<sumPionsPlus.size();i++) {tauandprodplus.push_back(sumPionsPlus.at(i));}  
  
  TLorentzVector Frame;
  if(type1=="pion" && type2=="pion"){
    Frame = tauandprodminus.at(1)+tauandprodplus.at(1);
  }
  else Frame = tauandprodminus.at(0)+tauandprodplus.at(0);
  //TLorentzVector Frame = tauandprodminus.at(1)+tauandprodplus.at(1);
  
  Scalc1.Configure(tauandprodminus,Frame, -1);
  TVector3 h1=Scalc1.pv();

  Scalc2.Configure(tauandprodplus,Frame, +1);
  TVector3 h2=Scalc2.pv();

  if(tauandprodminus.at(0)==zeroLV){ cout<<endl;tauandprodminus.at(0).Print();}
  if(tauandprodplus.at(0)==zeroLV) {cout<<endl;tauandprodplus.at(0).Print();}
  
  TLorentzVector tauminus_HRF = Scalc1.Boost(tauandprodminus.at(0),Frame);

  TLorentzVector tauplus_HRF  = Scalc2.Boost(tauandprodplus.at(0),Frame);
  
  // tauandprodminus.at(0).Print();
  // tauandprodplus.at(0).Print();
  // h1.Print();
  // h2.Print();
  
  double h1Norm=1./h1.Mag();
  double h2Norm=1./h2.Mag();
  h1=h1*h1Norm;
  h2=h2*h2Norm;
  double k1Norm=1./((h1.Cross(tauminus_HRF.Vect().Unit())).Mag());
  double k2Norm=1./((h2.Cross(tauplus_HRF.Vect().Unit())).Mag());
  TVector3 k1 = (h1.Cross(tauminus_HRF.Vect().Unit()))*k1Norm;
  TVector3 k2 = (h2.Cross(tauplus_HRF.Vect().Unit()))*k2Norm;
  
  for (unsigned int i=0;i<sumPionsMinus.size();i++){HadLVMinus+=sumPionsMinus.at(i);}
  for (unsigned int i=0;i<sumPionsPlus.size();i++){HadLVPlus+=sumPionsPlus.at(i);}
 
  // tauminus_HRF.Vect().Unit().Print();
  // tauplus_HRF.Vect().Unit().Print();
  // cout<<"k1norm: "<<k1Norm<<endl;
  // cout<<"k2norm: "<<k2Norm<<endl;
  // k1.Print();
  // k2.Print();
  // cout<<"atan2: "<<TMath::ATan2((k1.Cross(k2)).Mag(),k1*k2)<<endl;
  //if(((h1.Cross(h2))*(tauminus_HRF.Vect().Unit()))<=0){cout<<"acos: "<<acos(k1*k2)<<endl; return acos(k1*k2);}
  //else{cout<<"acos: "<<2.*TMath::Pi()-acos(k1*k2)<<endl; return (2.*TMath::Pi()-acos(k1*k2));}
  if(((h1.Cross(h2))*(tauminus_HRF.Vect().Unit()))<=0){// if( TMath::ATan2((k1.Cross(k2)).Mag(),k1*k2)>2.1) {tauandprodminus.at(0).Print();HadLVMinus.Print();
    // tauandprodplus.at(0).Print();HadLVPlus.Print();
    // h1.Print();
    // h2.Print();tauminus_HRF.Vect().Unit().Print();
    // tauplus_HRF.Vect().Unit().Print();
    // cout<<"k1norm: "<<k1Norm<<endl;
    // cout<<"k2norm: "<<k2Norm<<endl;
    // k1.Print();
    // k2.Print();cout<<"ATan2: "<<TMath::ATan2((k1.Cross(k2)).Mag(),k1*k2)<<endl;} 
    return TMath::ATan2((k1.Cross(k2)).Mag(),k1*k2);}
  else{// if((2.*TMath::Pi()- TMath::ATan2((k1.Cross(k2)).Mag(),k1*k2))<4.1){tauandprodminus.at(0).Print();HadLVMinus.Print();
    // tauandprodplus.at(0).Print();HadLVPlus.Print();
    // h1.Print();
    // h2.Print();tauminus_HRF.Vect().Unit().Print();
    // tauplus_HRF.Vect().Unit().Print();
    // cout<<"k1norm: "<<k1Norm<<endl;
    // cout<<"k2norm: "<<k2Norm<<endl;
    // k1.Print();
    // k2.Print();cout<<"ATan2: "<<2.*TMath::Pi()-TMath::ATan2((k1.Cross(k2)).Mag(),k1*k2)<<endl;} 
    return (2.*TMath::Pi()-TMath::ATan2((k1.Cross(k2)).Mag(),k1*k2));}
}

double SCalculator::AcopAngle_DP(TString type1, TString type2, std::vector<TLorentzVector> sumPionsMinus, std::vector<TLorentzVector> sumPionsPlus)
{

  SCalculator Scalc1(type1.Data());
  SCalculator Scalc2(type2.Data());

  TLorentzVector PiMinus, PiPlus, PiZeroMinus, PiZeroPlus;
  double Y1 = 1., Y2 = 1., Y = 1.;

  if(type1 == "a1"){
    double Minv1 = (sumPionsMinus.at(0)+sumPionsMinus.at(1)).Mag();
    double Minv2 = (sumPionsMinus.at(0)+sumPionsMinus.at(2)).Mag();
	
    unsigned int idx_piminus;
    //Compare to rho zero mass
    if(fabs(0.77526-Minv1)<fabs(0.77526-Minv2)) {idx_piminus = 1;}
    else {idx_piminus = 2;}

    PiMinus = sumPionsMinus.at(idx_piminus);
    PiZeroMinus = sumPionsMinus.at(0);
  }
  else if(type1 == "rho"){
    PiMinus = sumPionsMinus.at(0);
    PiZeroMinus = sumPionsMinus.at(1);
  }
  else{cout<<"No Decay Plane for this decay mode"<<endl;}

  if(type2 == "a1"){
    double Minv1 = (sumPionsPlus.at(0)+sumPionsPlus.at(1)).Mag();
    double Minv2 = (sumPionsPlus.at(0)+sumPionsPlus.at(2)).Mag();

    unsigned int idx_piplus;
    //Compare to rho zero mass
    if(fabs(0.77526-Minv1)<fabs(0.77526-Minv2)) {idx_piplus = 1;}
    else {idx_piplus = 2;}

    PiPlus = sumPionsPlus.at(idx_piplus);
    PiZeroPlus = sumPionsPlus.at(0);
  }
  else if(type2 == "rho"){
    PiPlus = sumPionsPlus.at(0);
    PiZeroPlus = sumPionsPlus.at(1);
  }
  else{cout<<"No Decay Plane for this decay mode"<<endl;}

  TLorentzVector ZMF = PiPlus + PiMinus;
  //Minus side
  TLorentzVector PiMinus_ZMF = Scalc1.Boost(PiMinus,ZMF);
  TLorentzVector PiZeroMinus_ZMF = Scalc1.Boost(PiZeroMinus,ZMF);
  TVector3 vecPiMinus = PiMinus_ZMF.Vect();
  TVector3 vecPiZeroMinus = PiZeroMinus_ZMF.Vect();
  vecPiMinus *= 1/vecPiMinus.Mag();
  vecPiZeroMinus *= 1/vecPiZeroMinus.Mag();
  TVector3 vecPiZeroMinustransv = vecPiZeroMinus - vecPiMinus*(vecPiMinus*vecPiZeroMinus);
  vecPiZeroMinustransv *= 1/vecPiZeroMinustransv.Mag();
  //Plus side
  TLorentzVector PiPlus_ZMF = Scalc2.Boost(PiPlus,ZMF);
  TLorentzVector PiZeroPlus_ZMF = Scalc2.Boost(PiZeroPlus,ZMF);
  TVector3 vecPiPlus = PiPlus_ZMF.Vect();
  TVector3 vecPiZeroPlus = PiZeroPlus_ZMF.Vect();
  vecPiPlus *= 1/vecPiPlus.Mag();
  vecPiZeroPlus *= 1/vecPiZeroPlus.Mag();
  TVector3 vecPiZeroPlustransv = vecPiZeroPlus - vecPiPlus*(vecPiPlus*vecPiZeroPlus);
  vecPiZeroPlustransv *= 1/vecPiZeroPlustransv.Mag();
  //Y variable
  Y1 = (PiMinus.E() - PiZeroMinus.E())/(PiMinus.E() + PiZeroMinus.E());
  Y2 = (PiPlus.E() - PiZeroPlus.E())/(PiPlus.E() + PiZeroPlus.E());
  Y = Y1*Y2;
  //angle
  double acop_DP = TMath::ACos(vecPiZeroPlustransv*vecPiZeroMinustransv);
  double sign_DP = vecPiMinus * (vecPiZeroPlustransv.Cross(vecPiZeroMinustransv));
  if (sign_DP<0) acop_DP = 2.0*TMath::Pi() - acop_DP;
  if(Y<0){
    acop_DP = acop_DP + TMath::Pi();
    if (acop_DP>2*TMath::Pi()) acop_DP = acop_DP - 2*TMath::Pi();
  }

  return acop_DP;
}

TVector3 SCalculator::GetIP(TLorentzVector p4, TVector3 ref)
{
  TVector3 dir = p4.Vect();
  double proj = ref*dir/dir.Mag2();
  return ref-dir*proj;
}


double SCalculator::AcopAngle_IP(TLorentzVector pion1, TVector3 IP1, TLorentzVector pion2, TVector3 IP2)
{
  SCalculator Scalc1("pion");
  SCalculator Scalc2("pion");
  //ZMF
  TLorentzVector ZMF = pion1 + pion2;
  //IP Minus
  TLorentzVector eta1(IP1,0.);
  //IP Plus
  TLorentzVector eta2(IP2,0.);
  //Minus side
  TLorentzVector eta1_ZMF = Scalc1.Boost(eta1,ZMF);
  TLorentzVector pion1_ZMF = Scalc1.Boost(pion1,ZMF);
  TVector3 eta1Vec = eta1_ZMF.Vect();
  TVector3 pion1Vec = pion1_ZMF.Vect();
  TVector3 eta1Vectransv = eta1Vec - pion1Vec*(eta1Vec*pion1Vec/pion1Vec.Mag2());
  eta1Vectransv *= 1/eta1Vectransv.Mag();
  //Plus side
  TLorentzVector eta2_ZMF = Scalc2.Boost(eta2,ZMF);
  TLorentzVector pion2_ZMF = Scalc2.Boost(pion2,ZMF);
  TVector3 eta2Vec = eta2_ZMF.Vect();
  TVector3 pion2Vec = pion2_ZMF.Vect();
  TVector3 eta2Vectransv = eta2Vec - pion2Vec*(eta2Vec*pion2Vec/pion2Vec.Mag2());
  eta2Vectransv *= 1/eta2Vectransv.Mag();
  //Angle 
  double acop_IP = TMath::ACos(eta2Vectransv*eta1Vectransv);
  double sign_IP = pion1Vec * (eta2Vectransv.Cross(eta1Vectransv));
  if (sign_IP<0) acop_IP = 2.0*TMath::Pi() - acop_IP;
  return acop_IP;
}

double SCalculator::AcopAngle_PVIP(TString type1, TString type2, TLorentzVector tau1, std::vector<TLorentzVector> sumPions, std::vector<double> sumPionsCharge, TLorentzVector pion, TVector3 pion_IP)
{
  if(!(type2 == "pion" || type2 == "muon") || type1 == "none") return -99;

  SCalculator Scalc1(type1.Data());
  SCalculator Scalc2;
  //ZMF
  TLorentzVector ZMF = tau1 + pion;
  //PV
  Scalc1.SortPions(sumPions, sumPionsCharge);
  //
  vector<TLorentzVector> tauandprod;
  int taucharge = 0;
  tauandprod.push_back(tau1);
  for(unsigned int i=0; i<sumPions.size();i++) {
    tauandprod.push_back(sumPions.at(i));
    taucharge+=sumPionsCharge.at(i);
  }
  //
  Scalc1.Configure(tauandprod,ZMF,taucharge);
  //
  TVector3 h1=Scalc1.pv();
  TLorentzVector tau_HRF = Scalc1.Boost(tau1,ZMF);
  h1*=1./h1.Mag();
  TVector3 k1 = (h1.Cross(tau_HRF.Vect().Unit())).Unit();
  //IP
  TLorentzVector eta(pion_IP,0.);
  //
  TLorentzVector eta_ZMF = Scalc1.Boost(eta,ZMF);
  TLorentzVector pion_ZMF = Scalc1.Boost(pion,ZMF);
  TVector3 etaVec = eta_ZMF.Vect();
  TVector3 pionVec = pion_ZMF.Vect();
  TVector3 etaVectransv = etaVec - pionVec*(etaVec*pionVec/pionVec.Mag2());
  etaVectransv *= 1/etaVectransv.Mag();
  //Angle
  double acop = TMath::ACos(k1*etaVectransv);
  double sign = pionVec * (k1.Cross(etaVectransv));
  if (sign<0) acop = 2.0*TMath::Pi() - acop;
  acop = acop + 0.5*TMath::Pi();
  if(acop>2.0*TMath::Pi()) acop = acop - 2.0*TMath::Pi();
  return acop;
}

double SCalculator::AcopAngle_DPIP(TString type1, TString type2, std::vector<TLorentzVector> sumPions, TLorentzVector pion, TVector3 pion_IP)
{
  if(!(type2 == "pion" || type2 == "muon") || type1=="none") return -99;

  SCalculator Scalc1(type1.Data());

  unsigned int idx_pi = -1;
  unsigned int idx_a1 = -1;

  TLorentzVector Pi, PiZero;
  double Y = 1.;

  if(type1 == "a1"){
    double Minv1 = (sumPions.at(0)+sumPions.at(1)).Mag();
    double Minv2 = (sumPions.at(0)+sumPions.at(2)).Mag();

    //Compare to rho zero mass
    if(fabs(0.77526-Minv1)<fabs(0.77526-Minv2)) {idx_pi = 1; idx_a1 = 2;}
    else {idx_pi = 2; idx_a1 = 1;}

    Pi = sumPions.at(idx_pi);
    PiZero = sumPions.at(0);
  }
  else if(type1 == "rho"){
    Pi = sumPions.at(0);
    PiZero = sumPions.at(1);
  }
  else{cout<<"No Decay Plane for this decay mode"<<endl;}

  Y = (Pi.E() - PiZero.E())/(Pi.E() + PiZero.E());

  TLorentzVector ZMF;

  if(type1 == "a1"){
    if(sumPions.at(idx_pi).E()>sumPions.at(idx_a1).E()){
      ZMF = sumPions.at(idx_pi) + pion;
    }
    if(sumPions.at(idx_pi).E()<sumPions.at(idx_a1).E()){
      ZMF = sumPions.at(idx_a1) + pion;
    }
  }
  if(type1 == "rho"){
    ZMF = Pi + pion;
  }

  TLorentzVector eta1_ZMF = Scalc1.Boost(PiZero,ZMF);
  TLorentzVector pi_ZMF = Scalc1.Boost(Pi,ZMF);
  TVector3 eta1Transv = eta1_ZMF.Vect() - pi_ZMF.Vect()*(pi_ZMF.Vect()*eta1_ZMF.Vect()/pi_ZMF.Vect().Mag2());
  eta1Transv *= 1/eta1Transv.Mag();

  SCalculator Scalc2;
  TLorentzVector eta2(pion_IP,0.);
  TLorentzVector eta2_ZMF = Scalc2.Boost(eta2,ZMF);
  TLorentzVector pion_ZMF = Scalc2.Boost(pion,ZMF);
  TVector3 eta2Transv = eta2_ZMF.Vect() - pion_ZMF.Vect()*(pion_ZMF.Vect()*eta2_ZMF.Vect()/pion_ZMF.Vect().Mag2());
  eta2Transv *= 1/eta2Transv.Mag();

  double acop = TMath::ACos(eta1Transv*eta2Transv);
  double sign = pion_ZMF.Vect() * (eta1Transv.Cross(eta2Transv));
  if (sign<0) acop = 2.0*TMath::Pi() - acop;
  if (Y<0) {
    acop = acop + TMath::Pi();
    if (acop>2*TMath::Pi()) acop = acop - 2*TMath::Pi();
  }
  return acop;
}

double SCalculator::M(TLorentzVector LV) 
{
  Double_t mm = LV.T()*LV.T()-LV.X()*LV.X()-LV.Y()*LV.Y()-LV.Z()*LV.Z();
  return mm < 0.0 ? -TMath::Sqrt(-mm) : TMath::Sqrt(mm);
}
