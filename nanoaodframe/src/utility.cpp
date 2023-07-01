/*
 * utility.cpp
 *
 *  Created on: Dec 4, 2018
 *      Author: suyong
 */
#include "utility.h"
#include "TMatrixDSym.h"
#include "TVectorT.h"
#include "Math/SpecFuncMathMore.h"
#include "Math/GenVector/VectorUtil.h"
#include "Math/GenVector/Rotation3D.h"
#include "Math/Math.h"

// Utility function to generate fourvector objects for thigs that pass selections

FourVectorVec gen4vec(floats &pt, floats &eta, floats &phi, floats &mass)
{
	const int nsize = pt.size();
	FourVectorVec fourvecs;
	fourvecs.reserve(nsize);
	for (auto i=0; i<nsize; i++)
	{
		fourvecs.emplace_back(pt[i], eta[i], phi[i], fabs(mass[i]));
	}

	return fourvecs;
}

FourVectorVec genmet4vec(float met_pt, float met_phi)
{
        FourVectorVec vecs;
        float met_px = met_pt*cos(met_phi);
        float met_py = met_pt*sin(met_phi);
        for(int i = -500; i <= 500; i+=50){
                FourVector metfourvec;
                metfourvec.SetPxPyPzE(met_px, met_py, i, met_pt);
                vecs.emplace_back(metfourvec);
        }
        return vecs;
}

floats weightv(floats &x, float evWeight)
{
	const int nsize = x.size();
	floats weightvector(nsize, evWeight);
	return weightvector;
}

floats sphericity(FourVectorVec &p)
{
	TMatrixDSym NormMomTensor(3);

	NormMomTensor = 0.0;
	double p2sum = 0.0;
	for (auto x: p)
	{
		p2sum += x.P2();
		double mom[3] = {x.Px(), x.Py(), x.Pz()};
		for (int irow=0; irow<3; irow++)
		{
			for (int icol=irow; icol<3; icol++)
			{
				NormMomTensor(irow, icol) += mom[irow] * mom[icol];
			}
		}
	}
	NormMomTensor *= (1.0/p2sum);
	TVectorT<double> Qrev;
	NormMomTensor.EigenVectors(Qrev);
	floats Q(3);
	for (auto i=0; i<3; i++) Q[i] = Qrev[2-i];

	return Q;
}


double foxwolframmoment(int l, FourVectorVec &p, int minj, int maxj)
{   // PRD 87, 073014 (2013)
	double answer = 0.0;

	double ptsum=0.0;

	if (maxj==-1) // process everything
	{
		maxj = p.size();
	}
	//for (auto x: p)
	for (auto i=minj; i<maxj; i++)
	{
		auto x = p[i];
		ptsum += x.Pt();
		//for (auto y: p)
		for (auto j=minj; j<maxj; j++)
		{
			auto y = p[j];
			double wij = x.Pt() * y.Pt();
			double cosdOmega = x.Vect().Dot(y.Vect()) / (x.P() * y.P());
			if (cosdOmega>1.0) cosdOmega=1.0;
			if (cosdOmega<-1.0) cosdOmega=-1.0;
			answer += wij * ROOT::Math::legendre(l, cosdOmega);
		}
	}
	answer /= ptsum*ptsum;
	if (fabs(answer)>1.0) std::cout << "FW>1 " << answer << std::endl;
	return answer;
}


ints good_idx(ints good)
{
        ints out;
        for(int i = 0; i < int(good.size()); i++){
                if( good[i] ){
                        out.emplace_back(i);
                }
        }
        return out;
}


floats top_reconstruction_STLFV(FourVectorVec &jets, FourVectorVec &bjets, FourVectorVec &muons, FourVectorVec &taus){
        
        floats out;
        float SMW_mass, SMtop_mass;
        float X_SMW, X_SMtop;
        float X_min=9999999999, X_min_SMW_mass=-1, X_min_SMtop_mass=-1;
        float X_min_SMW=999999999, X_min_SMtop=999999999;
        float wj1_idx=-1, wj2_idx=-1;
        const float MT = 165.2;
        const float MW = 80.8;
        const float WT = 21.3;
        const float WW = 11.71;	
        
        // Jets from W-1
        for(int j1 = 0; j1<int(jets.size()); j1++){
            if(jets[j1].Pt() == bjets[0].Pt()) continue;
            // Jets from W-2
            for(int j2 = 0; j2<int(jets.size()); j2++){
                if(jets[j2].Pt() == jets[j1].Pt() || jets[j2].Pt() == bjets[0].Pt()) continue;
                SMW_mass = (jets[j1]+jets[j2]).M();
                X_SMW = std::pow((MW-SMW_mass)/WW,2);
                SMtop_mass = (bjets[0]+jets[j1]+jets[j2]).M();
                X_SMtop = std::pow((MT-SMtop_mass)/WT,2);
                if (X_SMW + X_SMtop < X_min){
                    X_min = X_SMW + X_SMtop;
                    X_min_SMW = X_SMW;
                    X_min_SMtop = X_SMtop;
                    X_min_SMW_mass = SMW_mass;
                    X_min_SMtop_mass = SMtop_mass;
                    wj1_idx = float(j1);
                    wj2_idx = float(j2);
                }
            }
        }
        out.push_back(X_min);               // 0
        out.push_back(X_min_SMW_mass);      // 1
        out.push_back(X_min_SMtop_mass);    // 2
        out.push_back(wj1_idx);             // 3
        out.push_back(wj2_idx);             // 4
        out.push_back(X_min_SMW);           // 5
        out.push_back(X_min_SMtop);         // 6

        return out;
}

floats top_reconstruction_TTLFV(FourVectorVec &jets, FourVectorVec &bjets, FourVectorVec &muons, FourVectorVec &taus){
        
        floats out;
        float LFVtop_mass, SMW_mass, SMtop_mass;
        float X_LFVtop, X_SMW, X_SMtop;
        float X_min=9999999999, X_min_LFVtop_mass=-1, X_min_SMW_mass=-1, X_min_SMtop_mass=-1;
        float X_min_LFVtop=999999999, X_min_SMW=999999999, X_min_SMtop=999999999;
        float lfvj_idx=-1, wj1_idx=-1, wj2_idx=-1;
        
        // Mass and Width of top and w boson
        const float MT_LFV = 150.5;
        const float MT_SM = 165.2;
        const float MW = 80.8;
        const float WT_LFV = 17.8;
        const float WT_SM = 21.3;
        const float WW = 11.71;	
        
        // lfv jet
        for(int j1 = 0; j1<int(jets.size()); j1++){
            LFVtop_mass = (jets[j1]+taus[0]+muons[0]).M();
            X_LFVtop = std::pow((MT_LFV-LFVtop_mass)/WT_LFV,2);
            // Jet1 from W
            for(int j2 = 0; j2<int(jets.size()); j2++){
                if(jets[j2].Pt() == bjets[0].Pt() || jets[j2].Pt() == jets[j1].Pt()) continue;
                // Jet2 from W
                for(int j3 = 0; j3<int(jets.size()); j3++){
                    if(jets[j3].Pt() == jets[j2].Pt() || jets[j3].Pt() == bjets[0].Pt() || jets[j3].Pt() == jets[j1].Pt()) continue;
                    SMW_mass = (jets[j2]+jets[j3]).M();
                    X_SMW = std::pow((MW-SMW_mass)/WW,2);
                    SMtop_mass = (bjets[0]+jets[j2]+jets[j3]).M();
                    X_SMtop = std::pow((MT_SM-SMtop_mass)/WT_SM,2);
                    if (X_LFVtop + X_SMW + X_SMtop < X_min){
                        X_min = X_LFVtop + X_SMW + X_SMtop;
                        X_min_LFVtop = X_LFVtop;
                        X_min_SMW = X_SMW;
                        X_min_SMtop = X_SMtop;
                        X_min_LFVtop_mass = LFVtop_mass;
                        X_min_SMW_mass = SMW_mass;
                        X_min_SMtop_mass = SMtop_mass;
                        lfvj_idx = float(j1);
                        wj1_idx = float(j2);
                        wj2_idx = float(j3);
                    }
                }
            }
        }
        out.push_back(X_min);               // 0
        out.push_back(X_min_LFVtop_mass);   // 1
        out.push_back(X_min_SMW_mass);      // 2
        out.push_back(X_min_SMtop_mass);    // 3
        out.push_back(lfvj_idx);            // 4
        out.push_back(wj1_idx);             // 5
        out.push_back(wj2_idx);             // 6
        out.push_back(X_min_LFVtop);        // 7
        out.push_back(X_min_SMW);           // 8
        out.push_back(X_min_SMtop);         // 9

        return out;
}

floats top_reco_products_STLFV(FourVectorVec &jets, FourVectorVec &muons, FourVectorVec &taus, floats topreco){
        floats out;
        int wjet1_idx = topreco[3];
        int wjet2_idx = topreco[4];

        FourVector wjet1 = jets[wjet1_idx];
        FourVector wjet2 = jets[wjet2_idx];

        float wqq_dEta = wjet1.Eta() - wjet2.Eta();
        float wqq_dPhi = ROOT::Math::VectorUtil::DeltaPhi(wjet1, wjet2);
        float wqq_dR = ROOT::Math::VectorUtil::DeltaR(wjet1, wjet2);

        out.emplace_back(wqq_dEta);
        out.emplace_back(wqq_dPhi);
        out.emplace_back(wqq_dR);
        return out;
}

floats top_reco_products_TTLFV(FourVectorVec &jets, FourVectorVec &muons, FourVectorVec &taus, floats topreco){
        floats out;
        int j_idx = topreco[4];
        int wjet1_idx = topreco[5];
        int wjet2_idx = topreco[6];

        FourVector lfvjet = jets[j_idx];
        FourVector wjet1 = jets[wjet1_idx];
        FourVector wjet2 = jets[wjet2_idx];
        FourVector tau = taus[0];
        FourVector muon = muons[0];

        float wqq_dEta = wjet1.Eta() - wjet2.Eta();
        float wqq_dPhi = ROOT::Math::VectorUtil::DeltaPhi(wjet1, wjet2);
        float wqq_dR = ROOT::Math::VectorUtil::DeltaR(wjet1, wjet2);

        float lfvjmu_dEta = lfvjet.Eta() - muon.Eta();
        float lfvjmu_dPhi = ROOT::Math::VectorUtil::DeltaPhi(lfvjet, muon);
        float lfvjmu_dR = ROOT::Math::VectorUtil::DeltaR(lfvjet, muon);
        float lfvjmu_mass = ROOT::Math::VectorUtil::InvariantMass(lfvjet, muon);

        float lfvjtau_dEta = lfvjet.Eta() - tau.Eta();
        float lfvjtau_dPhi = ROOT::Math::VectorUtil::DeltaPhi(lfvjet, tau);
        float lfvjtau_dR = ROOT::Math::VectorUtil::DeltaR(lfvjet, tau);
        float lfvjtau_mass = ROOT::Math::VectorUtil::InvariantMass(lfvjet, tau);

        FourVector mutau = muon + tau;
        float lfvjmutau_dEta = lfvjet.Eta() - mutau.Eta();
        float lfvjmutau_dPhi = ROOT::Math::VectorUtil::DeltaPhi(lfvjet, mutau);
        float lfvjmutau_dR = ROOT::Math::VectorUtil::DeltaR(lfvjet, mutau);
        float lfvjmutau_mass = ROOT::Math::VectorUtil::InvariantMass(lfvjet, mutau);

        out.emplace_back(wqq_dEta);         //0
        out.emplace_back(wqq_dPhi);         //1
        out.emplace_back(wqq_dR);           //2
        out.emplace_back(lfvjmu_dEta);      //3
        out.emplace_back(lfvjmu_dPhi);      //4
        out.emplace_back(lfvjmu_dR);        //5
        out.emplace_back(lfvjmu_mass);      //6
        out.emplace_back(lfvjtau_dEta);     //7
        out.emplace_back(lfvjtau_dPhi);     //8
        out.emplace_back(lfvjtau_dR);       //9
        out.emplace_back(lfvjtau_mass);     //10
        out.emplace_back(lfvjmutau_dEta);   //11
        out.emplace_back(lfvjmutau_dPhi);   //12
        out.emplace_back(lfvjmutau_dR);     //13
        out.emplace_back(lfvjmutau_mass);   //14

        return out;
}


float calculate_deltaEta( FourVector &p1, FourVector &p2){
        return p1.Eta() - p2.Eta();
}

float calculate_deltaPhi( FourVector &p1, FourVector &p2){
        return ROOT::Math::VectorUtil::DeltaPhi(p1, p2);
}

float calculate_deltaR( FourVector &p1, FourVector &p2){
        return ROOT::Math::VectorUtil::DeltaR(p1, p2);
}

float calculate_invMass( FourVector &p1, FourVector &p2){
        return ROOT::Math::VectorUtil::InvariantMass(p1, p2);
}

float calculate_MT( FourVectorVec &muons, float met, float metphi){
        FourVector muon;
        muon = muons[0];
        float dphi = muon.Phi() - metphi;
        if ( dphi > M_PI ) {
            dphi -= 2.0*M_PI;
        } else if ( dphi <= -M_PI ) {
            dphi += 2.0*M_PI;
        }
        float out = sqrt(2*muon.Pt()*met*(1-cos(dphi)));
        return out;
}

FourVector sum_4vec( FourVector &p1, FourVector &p2){
        return p1+p2;
}

floats sort_discriminant( floats discr, floats obj ){
        auto sorted_discr = Reverse(Argsort(discr));
        floats out;
        for (auto idx : sorted_discr){
                out.emplace_back(obj[idx]);
        }
        return out;
}

ints find_element(ints vec, int a){
    ints idx;
    for(int i = 0; i < int(vec.size()); i++){
        if( vec[i] == a ) idx.emplace_back(i);
    }
    return idx;
}

ints find_element_binary( ints vec, int a){
    ints bin;
    for( int i = 0; i < int(vec.size()); i++){
        if( vec[i] == a ) bin.emplace_back(1);
        else bin.emplace_back(0);
    }
    return bin;
}

/////Find last genparticle using pdgid/////(i : idx, id : pdgId, t : target, d: daughter)
ints LastGenPart_idx( int target_id, ints GenPart_pdgId, ints GenPart_genPartIdxMother){
    ints out;
    for( int idx = 0 ; idx < int(GenPart_pdgId.size()); idx++){
        int p_id = GenPart_pdgId[idx];
        if(abs(p_id)==target_id){
            bool isLastParticle = true;
            ints daughters_idx = find_element(GenPart_genPartIdxMother, idx);
            for(int d_idx:daughters_idx){
                int d_id = GenPart_pdgId[d_idx];
                if (abs(d_id)==target_id) isLastParticle = false;
            }
            if(isLastParticle) out.emplace_back(idx);
        }
    }
    return out;
}


/////Find last genparticle using idx/////
int lastgenpart_idx(int target_i, ints GenPart_pdgId, ints GenPart_genPartIdxMother){
    int out;
    int target_id = GenPart_pdgId[target_i];

    int ttemp_i = target_i;
    while(true){
        bool LP = true;
        ints daughters_i = find_element(GenPart_genPartIdxMother, ttemp_i);
        if ( daughters_i.size() == 0 ){}

        else for( int dtemp_i: daughters_i){
            int dtemp_id = GenPart_pdgId[dtemp_i];
            if( dtemp_id == target_id ){
                LP = false;
                ttemp_i = dtemp_i;
            }
        }
        if( LP ){
        out=ttemp_i;
        break;
        }
    }
    return out;
}

ints FinalGenPart_idx( ints GenPart_pdgId, ints GenPart_genPartIdxMother ){
    ints out;
    int LFVtop_idx = -1, SMtop_idx=-1;
    int up_idx = -1, muon_idx = -1, tau_idx = -1;
    int b_idx = -1, W_idx = -1;
    int Wq1_idx=-1, Wq2_idx=-1;
    ints Wds_i;
    ints LastTop_idx = LastGenPart_idx(6, GenPart_pdgId, GenPart_genPartIdxMother);
    for( int i : LastTop_idx ){
        ints ds_idx = find_element(GenPart_genPartIdxMother, i);
        for( int d_idx : ds_idx ){
            int d_id = GenPart_pdgId[d_idx];
            if( abs(d_id) == 4 || abs(d_id) == 2){
                up_idx = lastgenpart_idx(d_idx, GenPart_pdgId, GenPart_genPartIdxMother);
                LFVtop_idx=i;
            }
            if(abs(d_id)==13){
                muon_idx = lastgenpart_idx(d_idx, GenPart_pdgId, GenPart_genPartIdxMother);
            }
            if(abs(d_id)==15){
                tau_idx = lastgenpart_idx(d_idx, GenPart_pdgId, GenPart_genPartIdxMother);
            }
            if(abs(d_id)==5){
                b_idx = lastgenpart_idx(d_idx, GenPart_pdgId, GenPart_genPartIdxMother);
                SMtop_idx=i;
            }
            if(abs(d_id)==24){
                W_idx = lastgenpart_idx(d_idx, GenPart_pdgId, GenPart_genPartIdxMother);
                Wds_i = find_element(GenPart_genPartIdxMother, W_idx);
                if( Wds_i.size() !=2 ) break;
                Wq1_idx = lastgenpart_idx(Wds_i[0], GenPart_pdgId, GenPart_genPartIdxMother);
                Wq2_idx = lastgenpart_idx(Wds_i[1], GenPart_pdgId, GenPart_genPartIdxMother);
            }
        }
    }
    out.emplace_back(up_idx);
    out.emplace_back(muon_idx);
    out.emplace_back(tau_idx);
    out.emplace_back(b_idx);
    out.emplace_back(Wq1_idx);
    out.emplace_back(Wq2_idx);
    out.emplace_back(LFVtop_idx);
    out.emplace_back(SMtop_idx);

    return out;
}


ints dRmatching_binary( int origin_i,float maxdR,  floats origin_pt, floats origin_eta, floats origin_phi, floats origin_mass, floats target_pt, floats target_eta, floats target_phi, floats target_mass){
    FourVectorVec origin, target;
    ints target_binary;
    int target_i = -1;
    float tempdR = 1, dR = maxdR;

    origin = gen4vec(origin_pt, origin_eta, origin_phi, origin_mass);
    target = gen4vec(target_pt, target_eta, target_phi, target_mass);

    for(int i=0; i<int(target.size()); i++){
        tempdR = ROOT::Math::VectorUtil::DeltaR(origin[origin_i],target[i]); 
        if( tempdR < dR ){
            target_i = i;
            dR = tempdR;
        }
    }
    for(int j=0; j<int(target.size()); j++){
        if( j == int(target_i) )target_binary.emplace_back(1);
        else target_binary.emplace_back(0);
    }
    return target_binary;
}


float SumMass( FourVectorVec object1, FourVectorVec object2, FourVectorVec object3){
    FourVector SumV = object1[0] + object2[0] + object3[0];
    float SumM = SumV.M();
    return SumM;
}

FourVector select_leadingvec( FourVectorVec &v ){
    FourVector vout;
    if(v.size() > 0) return v[0];
    else return vout;
}

//Considering only one muon!!!
floats addMuonUnc( floats &input ) {
    floats out;
    out.emplace_back(input[0]);
    out.emplace_back(input[0] + sqrt(pow(input[1] - input[0], 2) + pow(0.005, 2)));
    out.emplace_back(input[0] - sqrt(pow(input[2] - input[0], 2) + pow(0.005, 2)));
    return out;
}
