#include "TObject.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include <vector>

int ican2 = 0;
void makeCanvas()  {
    TCanvas * can = new TCanvas( TString::Format( "can%d", ican2++ ), "", 900, 900);
    can->SetBottomMargin(0.1);
    can->SetTopMargin(0.1);
    can->SetLeftMargin(0.1);
    can->SetRightMargin(0.1);
}

void plotting() {
    TFile * fo = new TFile( "final_plots.root", "RECREATE" );

    TFile * h5gev = TFile::Open("backhcal_5gev_batch1.root");
    TFile * h4gev = TFile::Open("backhcal_4gev_batch1.root");
    TFile * h3gev = TFile::Open("backhcal_3gev_batch1.root");
    TFile * h2gev = TFile::Open("backhcal_2gev_batch1.root");
    TFile * h1_5gev = TFile::Open("backhcal_1.5gev_batch1.root");
    TFile * h1gev = TFile::Open("backhcal_1gev_batch1.root");
    TFile * h0_8gev = TFile::Open("backhcal_0.8gev_batch1.root");
    TFile * h0_65gev = TFile::Open("backhcal_0.65gev_batch1.root");
    TFile * h0_5gev = TFile::Open("backhcal_0.5gev_batch1.root");
    TFile * h0_4gev = TFile::Open("backhcal_0.4gev_batch1.root");
    TFile * h0_3gev = TFile::Open("backhcal_0.3gev_batch1.root");
    TFile * h0_2gev = TFile::Open("backhcal_0.2gev_batch1.root");
    TFile * h0_1gev = TFile::Open("backhcal_0.1gev_batch1.root");

    h5gev->ls();
    auto h_5gev_hitcont_energy = (TH1D*)h5gev->Get("h_nHCal_hit_contrib_energy");
    auto h_5gev_hitcont_time = (TH1D*)h5gev->Get("h_nHCal_hit_contrib_time");
    auto h_5gev_hitcont_energy_vs_time = (TH2D*)h5gev->Get("h_nHCal_hit_contrib_energy_vs_time");
    auto h_5gev_hitcont_energy_vs_telap = (TH2D*)h5gev->Get("h_nHCal_hit_contrib_energy_vs_telap");
    auto h_5gev_hitcont_energy_vs_time_total = (TH2D*)h5gev->Get("h_nHCal_hit_contrib_energy_vs_time_total");

    h4gev->ls();
    auto h_4gev_hitcont_energy = (TH1D*)h4gev->Get("h_nHCal_hit_contrib_energy");
    auto h_4gev_hitcont_time = (TH1D*)h4gev->Get("h_nHCal_hit_contrib_time");
    auto h_4gev_hitcont_energy_vs_time = (TH2D*)h4gev->Get("h_nHCal_hit_contrib_energy_vs_time");
    auto h_4gev_hitcont_energy_vs_telap = (TH2D*)h4gev->Get("h_nHCal_hit_contrib_energy_vs_telap");
    auto h_4gev_hitcont_energy_vs_time_total = (TH2D*)h4gev->Get("h_nHCal_hit_contrib_energy_vs_time_total");

    h3gev->ls();
    auto h_3gev_hitcont_energy = (TH1D*)h3gev->Get("h_nHCal_hit_contrib_energy");
    auto h_3gev_hitcont_time = (TH1D*)h3gev->Get("h_nHCal_hit_contrib_time");
    auto h_3gev_hitcont_energy_vs_time = (TH2D*)h3gev->Get("h_nHCal_hit_contrib_energy_vs_time");
    auto h_3gev_hitcont_energy_vs_telap = (TH2D*)h3gev->Get("h_nHCal_hit_contrib_energy_vs_telap");
    auto h_3gev_hitcont_energy_vs_time_total = (TH2D*)h3gev->Get("h_nHCal_hit_contrib_energy_vs_time_total");

    h2gev->ls();
    auto h_2gev_hitcont_energy = (TH1D*)h2gev->Get("h_nHCal_hit_contrib_energy");
    auto h_2gev_hitcont_time = (TH1D*)h2gev->Get("h_nHCal_hit_contrib_time");
    auto h_2gev_hitcont_energy_vs_time = (TH2D*)h2gev->Get("h_nHCal_hit_contrib_energy_vs_time");
    auto h_2gev_hitcont_energy_vs_telap = (TH2D*)h2gev->Get("h_nHCal_hit_contrib_energy_vs_telap");
    auto h_2gev_hitcont_energy_vs_time_total = (TH2D*)h2gev->Get("h_nHCal_hit_contrib_energy_vs_time_total");

    h1_5gev->ls();
    auto h_1_5gev_hitcont_energy = (TH1D*)h1_5gev->Get("h_nHCal_hit_contrib_energy");
    auto h_1_5gev_hitcont_time = (TH1D*)h1_5gev->Get("h_nHCal_hit_contrib_time");
    auto h_1_5gev_hitcont_energy_vs_time = (TH2D*)h1_5gev->Get("h_nHCal_hit_contrib_energy_vs_time");
    auto h_1_5gev_hitcont_energy_vs_telap = (TH2D*)h1_5gev->Get("h_nHCal_hit_contrib_energy_vs_telap");
    auto h_1_5gev_hitcont_energy_vs_time_total = (TH2D*)h1_5gev->Get("h_nHCal_hit_contrib_energy_vs_time_total");

    h1gev->ls();
    auto h_1gev_hitcont_energy = (TH1D*)h1gev->Get("h_nHCal_hit_contrib_energy");
    auto h_1gev_hitcont_time = (TH1D*)h1gev->Get("h_nHCal_hit_contrib_time");
    auto h_1gev_hitcont_energy_vs_time = (TH2D*)h1gev->Get("h_nHCal_hit_contrib_energy_vs_time");
    auto h_1gev_hitcont_energy_vs_telap = (TH2D*)h1gev->Get("h_nHCal_hit_contrib_energy_vs_telap");
    auto h_1gev_hitcont_energy_vs_time_total = (TH2D*)h1gev->Get("h_nHCal_hit_contrib_energy_vs_time_total");

    h0_8gev->ls();
    auto h_0_8gev_hitcont_energy = (TH1D*)h0_8gev->Get("h_nHCal_hit_contrib_energy");
    auto h_0_8gev_hitcont_time = (TH1D*)h0_8gev->Get("h_nHCal_hit_contrib_time");
    auto h_0_8gev_hitcont_energy_vs_time = (TH2D*)h0_8gev->Get("h_nHCal_hit_contrib_energy_vs_time");
    auto h_0_8gev_hitcont_energy_vs_telap = (TH2D*)h0_8gev->Get("h_nHCal_hit_contrib_energy_vs_telap");
    auto h_0_8gev_hitcont_energy_vs_time_total = (TH2D*)h0_8gev->Get("h_nHCal_hit_contrib_energy_vs_time_total");

    h0_65gev->ls();
    auto h_0_65gev_hitcont_energy = (TH1D*)h0_65gev->Get("h_nHCal_hit_contrib_energy");
    auto h_0_65gev_hitcont_time = (TH1D*)h0_65gev->Get("h_nHCal_hit_contrib_time");
    auto h_0_65gev_hitcont_energy_vs_time = (TH2D*)h0_65gev->Get("h_nHCal_hit_contrib_energy_vs_time");
    auto h_0_65gev_hitcont_energy_vs_telap = (TH2D*)h0_65gev->Get("h_nHCal_hit_contrib_energy_vs_telap");
    auto h_0_65gev_hitcont_energy_vs_time_total = (TH2D*)h0_65gev->Get("h_nHCal_hit_contrib_energy_vs_time_total");

    h0_5gev->ls();
    auto h_0_5gev_hitcont_energy = (TH1D*)h0_5gev->Get("h_nHCal_hit_contrib_energy");
    auto h_0_5gev_hitcont_time = (TH1D*)h0_5gev->Get("h_nHCal_hit_contrib_time");
    auto h_0_5gev_hitcont_energy_vs_time = (TH2D*)h0_5gev->Get("h_nHCal_hit_contrib_energy_vs_time");
    auto h_0_5gev_hitcont_energy_vs_telap = (TH2D*)h0_5gev->Get("h_nHCal_hit_contrib_energy_vs_telap");
    auto h_0_5gev_hitcont_energy_vs_time_total = (TH2D*)h0_5gev->Get("h_nHCal_hit_contrib_energy_vs_time_total");

    h0_4gev->ls();
    auto h_0_4gev_hitcont_energy = (TH1D*)h0_4gev->Get("h_nHCal_hit_contrib_energy");
    auto h_0_4gev_hitcont_time = (TH1D*)h0_4gev->Get("h_nHCal_hit_contrib_time");
    auto h_0_4gev_hitcont_energy_vs_time = (TH2D*)h0_4gev->Get("h_nHCal_hit_contrib_energy_vs_time");
    auto h_0_4gev_hitcont_energy_vs_telap = (TH2D*)h0_4gev->Get("h_nHCal_hit_contrib_energy_vs_telap");
    auto h_0_4gev_hitcont_energy_vs_time_total = (TH2D*)h0_4gev->Get("h_nHCal_hit_contrib_energy_vs_time_total");

    h0_3gev->ls();
    auto h_0_3gev_hitcont_energy = (TH1D*)h0_3gev->Get("h_nHCal_hit_contrib_energy");
    auto h_0_3gev_hitcont_time = (TH1D*)h0_3gev->Get("h_nHCal_hit_contrib_time");
    auto h_0_3gev_hitcont_energy_vs_time = (TH2D*)h0_3gev->Get("h_nHCal_hit_contrib_energy_vs_time");
    auto h_0_3gev_hitcont_energy_vs_telap = (TH2D*)h0_3gev->Get("h_nHCal_hit_contrib_energy_vs_telap");
    auto h_0_3gev_hitcont_energy_vs_time_total = (TH2D*)h0_3gev->Get("h_nHCal_hit_contrib_energy_vs_time_total");

    h0_2gev->ls();
    auto h_0_2gev_hitcont_energy = (TH1D*)h0_2gev->Get("h_nHCal_hit_contrib_energy");
    auto h_0_2gev_hitcont_time = (TH1D*)h0_2gev->Get("h_nHCal_hit_contrib_time");
    auto h_0_2gev_hitcont_energy_vs_time = (TH2D*)h0_2gev->Get("h_nHCal_hit_contrib_energy_vs_time");
    auto h_0_2gev_hitcont_energy_vs_telap = (TH2D*)h0_2gev->Get("h_nHCal_hit_contrib_energy_vs_telap");
    auto h_0_2gev_hitcont_energy_vs_time_total = (TH2D*)h0_2gev->Get("h_nHCal_hit_contrib_energy_vs_time_total");

    h0_1gev->ls();
    auto h_0_1gev_hitcont_energy = (TH1D*)h0_1gev->Get("h_nHCal_hit_contrib_energy");
    auto h_0_1gev_hitcont_time = (TH1D*)h0_1gev->Get("h_nHCal_hit_contrib_time");
    auto h_0_1gev_hitcont_energy_vs_time = (TH2D*)h0_1gev->Get("h_nHCal_hit_contrib_energy_vs_time");
    auto h_0_1gev_hitcont_energy_vs_telap = (TH2D*)h0_1gev->Get("h_nHCal_hit_contrib_energy_vs_telap");
    auto h_0_1gev_hitcont_energy_vs_time_total = (TH2D*)h0_1gev->Get("h_nHCal_hit_contrib_energy_vs_time_total");

    auto * h_E1t1 = new TH1F("h_E1t1", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * h_E1t2 = new TH1F("h_E1t2", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * h_E1t3 = new TH1F("h_E1t3", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * h_E1t4 = new TH1F("h_E1t4", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * h_E1t5 = new TH1F("h_E1t5", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * h_E2t1 = new TH1F("h_E2t1", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * h_E2t2 = new TH1F("h_E2t2", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * h_E2t3 = new TH1F("h_E2t3", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * h_E2t4 = new TH1F("h_E2t4", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * h_E2t5 = new TH1F("h_E2t5", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * h_E3t1 = new TH1F("h_E3t1", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * h_E3t2 = new TH1F("h_E3t2", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * h_E3t3 = new TH1F("h_E3t3", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * h_E3t4 = new TH1F("h_E3t4", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * h_E3t5 = new TH1F("h_E3t5", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * h_E4t1 = new TH1F("h_E4t1", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * h_E4t2 = new TH1F("h_E4t2", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * h_E4t3 = new TH1F("h_E4t3", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * h_E4t4 = new TH1F("h_E4t4", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * h_E4t5 = new TH1F("h_E4t5", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * h_E5t1 = new TH1F("h_E5t1", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * h_E5t2 = new TH1F("h_E5t2", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * h_E5t3 = new TH1F("h_E5t3", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * h_E5t4 = new TH1F("h_E5t4", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * h_E5t5 = new TH1F("h_E5t5", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * htot_E1t1 = new TH1F("htot_E1t1", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * htot_E1t2 = new TH1F("htot_E1t2", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * htot_E1t3 = new TH1F("htot_E1t3", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * htot_E1t4 = new TH1F("htot_E1t4", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * htot_E1t5 = new TH1F("htot_E1t5", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * htot_E2t1 = new TH1F("htot_E2t1", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * htot_E2t2 = new TH1F("htot_E2t2", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * htot_E2t3 = new TH1F("htot_E2t3", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * htot_E2t4 = new TH1F("htot_E2t4", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * htot_E2t5 = new TH1F("htot_E2t5", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * htot_E3t1 = new TH1F("htot_E3t1", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * htot_E3t2 = new TH1F("htot_E3t2", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * htot_E3t3 = new TH1F("htot_E3t3", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * htot_E3t4 = new TH1F("htot_E3t4", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * htot_E3t5 = new TH1F("htot_E3t5", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * htot_E4t1 = new TH1F("htot_E4t1", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * htot_E4t2 = new TH1F("htot_E4t2", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * htot_E4t3 = new TH1F("htot_E4t3", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * htot_E4t4 = new TH1F("htot_E4t4", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * htot_E4t5 = new TH1F("htot_E4t5", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * htot_E5t1 = new TH1F("htot_E5t1", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * htot_E5t2 = new TH1F("htot_E5t2", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * htot_E5t3 = new TH1F("htot_E5t3", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * htot_E5t4 = new TH1F("htot_E5t4", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * htot_E5t5 = new TH1F("htot_E5t5", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 5.5);
    auto * h_tot = new TH1F("h_tot", ";Incident Neutron Energy (GeV); Efficiency", 1000, 0, 2.5);

cout << "test";

double E_mip = 0.00075;
vector<double> E_th{ 0.5*E_mip, 0.25*E_mip, 0.1*E_mip, 0.05*E_mip, 0};
vector<double> t_max{100, 500, 1000, 10000, 25};
vector<double> E_gun{ 5, 4, 3, 2, 1.5, 1, 0.8, 0.65, 0.5, 0.4, 0.3, 0.2, 0.1};

vector<TH2D*> in_hists_total{ h_5gev_hitcont_energy_vs_time_total, h_4gev_hitcont_energy_vs_time_total, h_3gev_hitcont_energy_vs_time_total, h_2gev_hitcont_energy_vs_time_total, h_1_5gev_hitcont_energy_vs_time_total, h_1gev_hitcont_energy_vs_time_total, h_0_8gev_hitcont_energy_vs_time_total, h_0_65gev_hitcont_energy_vs_time_total, h_0_5gev_hitcont_energy_vs_time_total, h_0_4gev_hitcont_energy_vs_time_total, h_0_3gev_hitcont_energy_vs_time_total, h_0_2gev_hitcont_energy_vs_time_total, h_0_1gev_hitcont_energy_vs_time_total };

vector<TH2D*> in_hists_time{ h_5gev_hitcont_energy_vs_time, h_4gev_hitcont_energy_vs_time, h_3gev_hitcont_energy_vs_time, h_2gev_hitcont_energy_vs_time, h_1_5gev_hitcont_energy_vs_time,  h_1gev_hitcont_energy_vs_time, h_0_8gev_hitcont_energy_vs_time, h_0_65gev_hitcont_energy_vs_time, h_0_5gev_hitcont_energy_vs_time, h_0_4gev_hitcont_energy_vs_time, h_0_3gev_hitcont_energy_vs_time, h_0_2gev_hitcont_energy_vs_time, h_0_1gev_hitcont_energy_vs_time };

vector<TH2D*> in_hists_telap{ h_5gev_hitcont_energy_vs_telap, h_4gev_hitcont_energy_vs_telap, h_3gev_hitcont_energy_vs_telap, h_2gev_hitcont_energy_vs_telap, h_1_5gev_hitcont_energy_vs_telap,  h_1gev_hitcont_energy_vs_telap, h_0_8gev_hitcont_energy_vs_telap, h_0_65gev_hitcont_energy_vs_telap, h_0_5gev_hitcont_energy_vs_telap, h_0_4gev_hitcont_energy_vs_telap, h_0_3gev_hitcont_energy_vs_telap, h_0_2gev_hitcont_energy_vs_telap, h_0_1gev_hitcont_energy_vs_telap };

vector<vector<TH1F*>> out_hists_time{ {h_E1t1, h_E1t2, h_E1t3, h_E1t4, h_E1t5 }, {h_E2t1, h_E2t2, h_E2t3, h_E2t4, h_E2t5 }, {h_E3t1, h_E3t2, h_E3t3, h_E3t4, h_E3t5 }, {h_E4t1, h_E4t2, h_E4t3, h_E4t4, h_E4t5 }, {h_E5t1, h_E5t2, h_E5t3, h_E5t4, h_E5t5 } };

vector<vector<TH1F*>> tot_hists_time{ {htot_E1t1, htot_E1t2, htot_E1t3, htot_E1t4, htot_E1t5 }, {htot_E2t1, htot_E2t2, htot_E2t3, htot_E2t4, htot_E2t5 }, {htot_E3t1, htot_E3t2, htot_E3t3, htot_E3t4, htot_E3t5 }, {htot_E4t1, htot_E4t2, htot_E4t3, htot_E4t4, htot_E4t5 }, {htot_E5t1, htot_E5t2, htot_E5t3, htot_E5t4, htot_E5t5 } };

for ( int i = 0; i < E_th.size(); i++ ) {
    for ( int j = 0; j < t_max.size(); j++ ) {
        for ( int k = 0; k < E_gun.size(); k++ ) {
            double out_bin_no = out_hists_time[i][j]->FindBin(E_gun[k]);
            double out_bin_content = in_hists_telap[k]->GetBinContent( in_hists_telap[k]->FindBin( E_th[i], t_max[j] ) );
            out_hists_time[i][j]->SetBinContent( out_bin_no, out_bin_content );
            double tot_bin_no = tot_hists_time[i][j]->FindBin(E_gun[k]);
            double tot_bin_content = in_hists_total[k]->GetBinContent( in_hists_total[k]->FindBin( E_th[i], t_max[j] ) );
            tot_hists_time[i][j]->SetBinContent( tot_bin_no, tot_bin_content );
        }
    }
}

for ( int i = 0; i < E_th.size(); i++ ) {
    for ( int j = 0; j < t_max.size(); j++ ) {
        out_hists_time[i][j]->Divide( out_hists_time[i][j], tot_hists_time[i][j], 1.0, 1.0, "B" );
    }
}

makeCanvas();
h_2gev_hitcont_energy->SetLineColor(kBlack);
h_2gev_hitcont_energy->SetLineColor(kGreen - 3);
h_0_5gev_hitcont_energy->SetLineColor(kMagenta);
h_0_2gev_hitcont_energy->SetLineColor(kGray);
h_0_1gev_hitcont_energy->SetLineColor(kRed - 5);
h_2gev_hitcont_energy->SetLineWidth(2);
h_1gev_hitcont_energy->SetLineWidth(2);
h_0_5gev_hitcont_energy->SetLineWidth(2);
h_0_2gev_hitcont_energy->SetLineWidth(2);
h_0_1gev_hitcont_energy->SetLineWidth(2);
h_2gev_hitcont_energy->SetStats(0);
h_2gev_hitcont_energy->SetMinimum(100);
gPad->SetLogy();
h_2gev_hitcont_energy->Draw();
h_1gev_hitcont_energy->Draw("SAME");
h_0_5gev_hitcont_energy->Draw("SAME");
h_0_2gev_hitcont_energy->Draw("SAME");
h_0_1gev_hitcont_energy->Draw("SAME");
auto leg1 = new TLegend(0.7,0.85,0.9,0.65);
leg1->AddEntry(h_2gev_hitcont_energy, "2 GeV Neutrons");
leg1->AddEntry(h_1gev_hitcont_energy, "1 GeV Neutrons");
leg1->AddEntry(h_0_5gev_hitcont_energy, "0.5 GeV Neutrons");
leg1->AddEntry(h_0_2gev_hitcont_energy, "0.2 GeV Neutrons");
leg1->AddEntry(h_0_1gev_hitcont_energy, "0.1 GeV Neutrons");
leg1->Draw();
gPad->Print("plots/hitcont_energy.pdf");
gPad->Print("plots/hitcont_energy.png");

makeCanvas();
h_2gev_hitcont_time->SetLineColor(kBlack);
h_2gev_hitcont_time->SetLineColor(kGreen - 3);
h_0_5gev_hitcont_time->SetLineColor(kMagenta);
h_0_2gev_hitcont_time->SetLineColor(kGray);
h_0_1gev_hitcont_time->SetLineColor(kRed - 5);
h_2gev_hitcont_time->SetLineWidth(2);
h_1gev_hitcont_time->SetLineWidth(2);
h_0_5gev_hitcont_time->SetLineWidth(2);
h_0_2gev_hitcont_time->SetLineWidth(2);
h_0_1gev_hitcont_time->SetLineWidth(2);
h_2gev_hitcont_time->SetStats(0);
h_2gev_hitcont_time->SetMinimum(100);
gPad->SetLogy();
h_2gev_hitcont_time->Draw();
h_1gev_hitcont_time->Draw("SAME");
h_0_5gev_hitcont_time->Draw("SAME");
h_0_2gev_hitcont_time->Draw("SAME");
h_0_1gev_hitcont_time->Draw("SAME");
auto leg2 = new TLegend(0.7,0.85,0.9,0.65);
leg2->AddEntry(h_2gev_hitcont_time, "2 GeV Neutrons");
leg2->AddEntry(h_1gev_hitcont_time, "1 GeV Neutrons");
leg2->AddEntry(h_0_5gev_hitcont_time, "0.5 GeV Neutrons");
leg2->AddEntry(h_0_2gev_hitcont_time, "0.2 GeV Neutrons");
leg2->AddEntry(h_0_1gev_hitcont_time, "0.1 GeV Neutrons");
leg2->Draw();
gPad->Print("plots/hitcont_time.pdf");
gPad->Print("plots/hitcont_time.png");

makeCanvas();
h_E1t1->SetMarkerColor(kOrange);
h_E1t1->SetMarkerStyle(53);
h_E1t2->SetMarkerColor(kGreen);
h_E1t2->SetMarkerStyle(53);
h_E1t3->SetMarkerColor(kBlue);
h_E1t3->SetMarkerStyle(53);
h_E2t1->SetMarkerColor(kOrange);
h_E2t1->SetMarkerStyle(54);
h_E2t2->SetMarkerColor(kGreen);
h_E2t2->SetMarkerStyle(54);
h_E2t3->SetMarkerColor(kBlue);
h_E2t3->SetMarkerStyle(54);
h_E3t1->SetMarkerColor(kOrange);
h_E3t1->SetMarkerStyle(55);
h_E3t2->SetMarkerColor(kGreen);
h_E3t2->SetMarkerStyle(55);
h_E3t3->SetMarkerColor(kBlue);
h_E3t3->SetMarkerStyle(55);
h_E4t1->SetMarkerColor(kOrange);
h_E4t1->SetMarkerStyle(56);
h_E4t2->SetMarkerColor(kGreen);
h_E4t2->SetMarkerStyle(56);
h_E4t3->SetMarkerColor(kBlue);
h_E4t3->SetMarkerStyle(56);

h_E1t4->SetMarkerColor(kMagenta);
h_E1t4->SetMarkerStyle(53);
h_E2t4->SetMarkerColor(kMagenta);
h_E2t4->SetMarkerStyle(54);
h_E3t4->SetMarkerColor(kMagenta);
h_E3t4->SetMarkerStyle(55);
h_E4t4->SetMarkerColor(kMagenta);
h_E4t4->SetMarkerStyle(56);

h_E1t5->SetMarkerColor(kRed);
h_E1t5->SetMarkerStyle(53);
h_E2t5->SetMarkerColor(kRed);
h_E2t5->SetMarkerStyle(54);
h_E3t5->SetMarkerColor(kRed);
h_E3t5->SetMarkerStyle(55);
h_E4t5->SetMarkerColor(kRed);
h_E4t5->SetMarkerStyle(56);

h_E5t1->SetMarkerColor(kOrange);
h_E5t1->SetMarkerStyle(58);
h_E5t2->SetMarkerColor(kGreen);
h_E5t2->SetMarkerStyle(58);
h_E5t3->SetMarkerColor(kBlue);
h_E5t3->SetMarkerStyle(58);
h_E5t4->SetMarkerColor(kMagenta);
h_E5t4->SetMarkerStyle(58);
h_E5t5->SetMarkerColor(kRed);
h_E5t5->SetMarkerStyle(58);


h_E1t1->SetStats(0);
h_E1t1->SetMaximum(1);
h_E1t1->SetMinimum(0);
h_E1t1->SetTitle("Efficiency of various E_{th} and t_{int} for n^{0}; n^{0} energy; Efficiency");
auto leg3 = new TLegend(0.55,0.85,.9,0.1);
gStyle->SetLegendTextSize(0.025);
for ( int i = 0; i < (E_th.size()); i++ ) {
    for ( int j = 0; j < (t_max.size()); j++ ) {
        if ( i==0 && j==0 ) { out_hists_time[i][j]->Draw("P"); }else{ out_hists_time[i][j]->Draw("P SAME"); }
    }
}
leg3->AddEntry( h_E1t5, "E_{th}=0.5*E_{MIP}, t_{int}=25 ns");
leg3->AddEntry( h_E1t1, "E_{th}=0.5*E_{MIP}, t_{int}=100 ns");
leg3->AddEntry( h_E1t2, "E_{th}=0.5*E_{MIP}, t_{int}=500 ns");
leg3->AddEntry( h_E1t3, "E_{th}=0.5*E_{MIP}, t_{int}=1000 ns");
leg3->AddEntry( h_E1t4, "E_{th}=0.5*E_{MIP}, t_{int}=10000 ns");
leg3->AddEntry( h_E2t5, "E_{th}=0.25*E_{MIP}, t_{int}=25 ns");
leg3->AddEntry( h_E2t1, "E_{th}=0.25*E_{MIP}, t_{int}=100 ns");
leg3->AddEntry( h_E2t2, "E_{th}=0.25*E_{MIP}, t_{int}=500 ns");
leg3->AddEntry( h_E2t3, "E_{th}=0.25*E_{MIP}, t_{int}=1000 ns");
leg3->AddEntry( h_E2t4, "E_{th}=0.25*E_{MIP}, t_{int}=10000 ns");
leg3->AddEntry( h_E3t5, "E_{th}=0.1*E_{MIP}, t_{int}=25 ns");
leg3->AddEntry( h_E3t1, "E_{th}=0.1*E_{MIP}, t_{int}=100 ns");
leg3->AddEntry( h_E3t2, "E_{th}=0.1*E_{MIP}, t_{int}=500 ns");
leg3->AddEntry( h_E3t3, "E_{th}=0.1*E_{MIP}, t_{int}=1000 ns");
leg3->AddEntry( h_E3t4, "E_{th}=0.1*E_{MIP}, t_{int}=10000 ns");
leg3->AddEntry( h_E4t5, "E_{th}=0.05*E_{MIP}, t_{int}=25 ns");
leg3->AddEntry( h_E4t1, "E_{th}=0.05*E_{MIP}, t_{int}=100 ns");
leg3->AddEntry( h_E4t2, "E_{th}=0.05*E_{MIP}, t_{int}=500 ns");
leg3->AddEntry( h_E4t3, "E_{th}=0.05*E_{MIP}, t_{int}=1000 ns");
leg3->AddEntry( h_E4t4, "E_{th}=0.05*E_{MIP}, t_{int}=10000 ns");
leg3->AddEntry( h_E5t5, "E_{th}=0, t_{int}=25 ns");
leg3->AddEntry( h_E5t1, "E_{th}=0, t_{int}=100 ns");
leg3->AddEntry( h_E5t2, "E_{th}=0, t_{int}=500 ns");
leg3->AddEntry( h_E5t3, "E_{th}=0, t_{int}=1000 ns");
leg3->AddEntry( h_E5t4, "E_{th}=0, t_{int}=10000 ns");
leg3->Draw();
gPad->Print("plot_telap_backHCal_neutron_eff.pdf");

makeCanvas();
h_E1t1->SetMarkerColor(kOrange);
h_E1t1->SetMarkerStyle(53);
h_E1t2->SetMarkerColor(kGreen);
h_E1t2->SetMarkerStyle(53);
h_E1t3->SetMarkerColor(kBlue);
h_E1t3->SetMarkerStyle(53);
h_E2t1->SetMarkerColor(kOrange);
h_E2t1->SetMarkerStyle(54);
h_E2t2->SetMarkerColor(kGreen);
h_E2t2->SetMarkerStyle(54);
h_E2t3->SetMarkerColor(kBlue);
h_E2t3->SetMarkerStyle(54);
h_E3t1->SetMarkerColor(kOrange);
h_E3t1->SetMarkerStyle(55);
h_E3t2->SetMarkerColor(kGreen);
h_E3t2->SetMarkerStyle(55);
h_E3t3->SetMarkerColor(kBlue);
h_E3t3->SetMarkerStyle(55);
h_E4t1->SetMarkerColor(kOrange);
h_E4t1->SetMarkerStyle(56);
h_E4t2->SetMarkerColor(kGreen);
h_E4t2->SetMarkerStyle(56);
h_E4t3->SetMarkerColor(kBlue);
h_E4t3->SetMarkerStyle(56);

h_E1t4->SetMarkerColor(kMagenta);
h_E1t4->SetMarkerStyle(53);
h_E2t4->SetMarkerColor(kMagenta);
h_E2t4->SetMarkerStyle(54);
h_E3t4->SetMarkerColor(kMagenta);
h_E3t4->SetMarkerStyle(55);
h_E4t4->SetMarkerColor(kMagenta);
h_E4t4->SetMarkerStyle(56);

h_E1t5->SetMarkerColor(kRed);
h_E1t5->SetMarkerStyle(53);
h_E2t5->SetMarkerColor(kRed);
h_E2t5->SetMarkerStyle(54);
h_E3t5->SetMarkerColor(kRed);
h_E3t5->SetMarkerStyle(55);
h_E4t5->SetMarkerColor(kRed);
h_E4t5->SetMarkerStyle(56);

h_E5t1->SetMarkerColor(kOrange);
h_E5t1->SetMarkerStyle(58);
h_E5t2->SetMarkerColor(kGreen);
h_E5t2->SetMarkerStyle(58);
h_E5t3->SetMarkerColor(kBlue);
h_E5t3->SetMarkerStyle(58);
h_E5t4->SetMarkerColor(kMagenta);
h_E5t4->SetMarkerStyle(58);
h_E5t5->SetMarkerColor(kRed);
h_E5t5->SetMarkerStyle(58);

h_E1t1->SetStats(0);
h_E1t1->SetMaximum(1);
h_E1t1->SetMinimum(0);
h_E1t1->SetTitle("Efficiency of various E_{th} and t_{int} for n^{0}; n^{0} energy; Efficiency");
auto leg35 = new TLegend(0.55,0.85,.9,0.1);
gStyle->SetLegendTextSize(0.025);
for ( int i = 0; i < (E_th.size()); i++ ) {
    for ( int j = 0; j < (t_max.size()); j++ ) {
        if ( i==0 && j==0 ) { out_hists_time[i][j]->Draw("P"); }else{ if (j==3 || j==4 || j==0) { out_hists_time[i][j]->Draw("P SAME"); } }
    }
}

leg35->AddEntry( h_E1t5, "E_{th}=0.5*E_{MIP}, t_{int}=25 ns");
leg35->AddEntry( h_E1t1, "E_{th}=0.5*E_{MIP}, t_{int}=100 ns");
leg35->AddEntry( h_E1t4, "E_{th}=0.5*E_{MIP}, t_{int}=10000 ns");
leg35->AddEntry( h_E2t5, "E_{th}=0.25*E_{MIP}, t_{int}=25 ns");
leg35->AddEntry( h_E2t1, "E_{th}=0.25*E_{MIP}, t_{int}=100 ns");
leg35->AddEntry( h_E2t4, "E_{th}=0.25*E_{MIP}, t_{int}=10000 ns");
leg35->AddEntry( h_E3t5, "E_{th}=0.1*E_{MIP}, t_{int}=25 ns");
leg35->AddEntry( h_E3t1, "E_{th}=0.1*E_{MIP}, t_{int}=100 ns");
leg35->AddEntry( h_E3t4, "E_{th}=0.1*E_{MIP}, t_{int}=10000 ns");
leg35->AddEntry( h_E4t5, "E_{th}=0.05*E_{MIP}, t_{int}=25 ns");
leg35->AddEntry( h_E4t1, "E_{th}=0.05*E_{MIP}, t_{int}=100 ns");
leg35->AddEntry( h_E4t4, "E_{th}=0.05*E_{MIP}, t_{int}=10000 ns");
leg35->AddEntry( h_E5t5, "E_{th}=0, t_{int}=25 ns");
leg35->AddEntry( h_E5t1, "E_{th}=0, t_{int}=100 ns");
leg35->AddEntry( h_E5t4, "E_{th}=0, t_{int}=10000 ns");
leg35->Draw();
gPad->Print("plot_telap_select_tmax_backHCal_neutron_eff.pdf");

makeCanvas();
h_E1t1->SetMarkerColor(kRed);
h_E2t1->SetMarkerColor(kOrange);
h_E3t1->SetMarkerColor(kGreen);
h_E4t1->SetMarkerColor(kBlue);
h_E5t1->SetMarkerColor(kMagenta);
h_E1t1->SetStats(0);
h_E1t1->SetMaximum(1);
h_E1t1->SetMinimum(0);
h_E1t1->SetTitle("Efficiency of various E_{th} and t_{int} for n^{0}; n^{0} energy; Efficiency");
auto leg4 = new TLegend(0.5,0.5,.9,0.1);
gStyle->SetLegendTextSize(0.03);
for ( int i = 0; i < (E_th.size()); i++ ) {
        if ( i==0 ) { out_hists_time[i][0]->Draw("P"); }else{ out_hists_time[i][0]->Draw("P SAME"); }
}
leg4->AddEntry( h_E1t1, "E_{th}=0.5*E_{MIP}, t_{int}=100 ns");
leg4->AddEntry( h_E2t1, "E_{th}=0.25*E_{MIP}, t_{int}=100 ns");
leg4->AddEntry( h_E3t1, "E_{th}=0.1*E_{MIP}, t_{int}=100 ns");
leg4->AddEntry( h_E4t1, "E_{th}=0.05*E_{MIP}, t_{int}=100 ns");
leg4->AddEntry( h_E5t1, "E_{th}=0, t_{int}=100 ns");
leg4->Draw();
gPad->Print("plot_telap_100ns_backHCal_neutron_eff.pdf");

makeCanvas();
h_E3t1->SetMarkerColor(kOrange);
h_E3t1->SetStats(0);
h_E3t1->SetMaximum(1);
h_E3t1->SetMinimum(0);
h_E3t1->SetTitle("Efficiency of various E_{th} and t_{int} for n^{0}; n^{0} energy; Efficiency");
auto leg5 = new TLegend(0.5,0.5,.9,0.1);
gStyle->SetLegendTextSize(0.03);
for ( int i = 0; i < (t_max.size()); i++ ) {
        if ( i==0 ) { out_hists_time[2][i]->Draw("P"); }else{ out_hists_time[2][i]->Draw("P SAME"); }
}
leg5->AddEntry( h_E3t5, "E_{th}=0.1*E_{MIP}, t_{int}=25 ns");
leg5->AddEntry( h_E3t1, "E_{th}=0.1*E_{MIP}, t_{int}=100 ns");
leg5->AddEntry( h_E3t2, "E_{th}=0.1*E_{MIP}, t_{int}=500 ns");
leg5->AddEntry( h_E3t3, "E_{th}=0.1*E_{MIP}, t_{int}=1000 ns");
leg5->AddEntry( h_E3t4, "E_{th}=0.1*E_{MIP}, t_{int}=10000 ns");
leg5->Draw();
gPad->Print("plot_telap_0.1Emip_backHCal_neutron_eff.pdf");

}
