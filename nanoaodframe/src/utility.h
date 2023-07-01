/*
 * utility.h
 *
 *  Created on: Dec 4, 2018
 *      Author: suyong
 */

#ifndef UTILITY_H_
#define UTILITY_H_

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "Math/Vector4D.h"
#include <string>

using floats =  ROOT::VecOps::RVec<float>;
using floatsVec =  ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>>;
using doubles =  ROOT::VecOps::RVec<double>;
using doublesVec =  ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>>;
using ints =  ROOT::VecOps::RVec<int>;
using bools = ROOT::VecOps::RVec<bool>;
using uchars = ROOT::VecOps::RVec<unsigned char>;
using strings = ROOT::VecOps::RVec<std::string>;

using FourVector = ROOT::Math::PtEtaPhiMVector;
using FourVectorVec = std::vector<FourVector>;
using FourVectorRVec = ROOT::VecOps::RVec<FourVector>;


struct hist1dinfo
{
  ROOT::RDF::TH1DModel hmodel;
  std::string varname;
  std::string weightname;
  std::string systname;
  std::string mincutstep;
  std::string maxcutstep;
} ;


struct varinfo
{
  std::string varname;
  std::string vardefinition;
  std::string mincutstep;
};

struct cutinfo
{
  std::string cutdefinition;
  std::string idx;
};


// generates vectors of 4 vectors given vectors of pt, eta, phi, mass
FourVectorVec gen4vec(floats &pt, floats &eta, floats &phi, floats &mass);

FourVectorVec genmet4vec(float met_pt, float met_phi);

// return a vector size equal to length of x all filled with evWeight value
floats weightv(floats &x, float evWeight);

floats sphericity(FourVectorVec &p);

double foxwolframmoment(int l, FourVectorVec &p, int minj=0, int maxj=-1);

ints good_idx(ints good);

floats lqtop_reconstruction( FourVectorVec &cjet, FourVectorVec &mu, FourVectorVec &tau);

floats top_reconstruction_STLFV(FourVectorVec &jets, FourVectorVec &bjets, FourVectorVec &muons, FourVectorVec &taus);

floats top_reconstruction_TTLFV(FourVectorVec &jets, FourVectorVec &bjets, FourVectorVec &muons, FourVectorVec &taus);

floats top_reco_products_STLFV(FourVectorVec &jets, FourVectorVec &muons, FourVectorVec &taus, floats topreco);

floats top_reco_products_TTLFV(FourVectorVec &jets, FourVectorVec &muons, FourVectorVec &taus, floats topreco);

float calculate_deltaEta( FourVector &p1, FourVector &p2);

float calculate_deltaPhi( FourVector &p1, FourVector &p2);

float calculate_deltaR( FourVector &p1, FourVector &p2);

float calculate_invMass( FourVector &p1, FourVector &p2);

float calculate_MT( FourVectorVec &muons, float met, float metphi);

FourVector sum_4vec( FourVector &p1, FourVector &p2);

floats sort_discriminant( floats discr, floats obj );

ints find_element(ints vec, int a);

ints find_element_binary( ints vec, int a);

ints LastGenPart_idx( int target_id, ints GenPart_pdgId, ints GenPart_genPartIdxMother);

int lastgenpart_idx(int target_i, ints GenPart_pdgId, ints GenPart_genPartIdxMother);

ints FinalGenPart_idx( ints GenPart_pdgId, ints GenPart_genPartIdxMother);

ints dRmatching_binary( int origin_i,float maxdR,  floats origin_pt, floats origin_eta, floats origin_phi, floats origin_mass, floats target_pt, floats target_eta, floats target_phi, floats target_mass);

float SumMass( FourVectorVec object1, FourVectorVec object2, FourVectorVec object3);

FourVector select_leadingvec( FourVectorVec &v );

floats addMuonUnc( floats &input );
#endif /* UTILITY_H_ */
