// c++
#include <iostream>
#include <sys/stat.h>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <set>

// ROOT
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TDatabasePDG.h>
#include <TParticle.h>

// Custom headers
#include "config.h"
#include "ana_helper.h"
#include "paths.h"

Config& conf = Config::getInstance();

void analyze(TString path, TString particle){    
    // +---------+
    // | setting |
    // +---------+
    gROOT->GetColor(kBlue)->SetRGB(0.12156862745098039, 0.4666666666666667, 0.7058823529411765);
    gROOT->GetColor(kOrange)->SetRGB(1.0, 0.4980392156862745, 0.054901960784313725);
    gROOT->GetColor(kGreen)->SetRGB(44.0/256, 160.0/256, 44.0/256);
    
    gStyle->SetOptStat(0);
    gStyle->SetLabelSize(0.06, "XY");
    gStyle->SetTitleSize(0.06, "x"); // x軸のタイトルサイズ
    gStyle->SetTitleSize(0.06, "y"); // y軸のタイトルサイズ
    gStyle->SetTitleFontSize(0.06);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gROOT->GetColor(0)->SetAlpha(0.01);


    // +-----------+
    // | load file |
    // +-----------+
    auto *f = new TFile(path.Data());
    if (!f || f->IsZombie()) {
        std::cerr << "Error: Could not open file : " << path << std::endl;
        return;
    }
    TTreeReader reader("hodo", f);
    TTreeReaderValue<unsigned int> run_number(reader, "run_number");
    reader.SetEntry(0);
    Int_t run_num = *run_number;

    // +--------------------------+
    // | prepare output root file |
    // +--------------------------+
    TString output_path = Form("%s/root/run%05d_T0_HDPEM_%s.root", OUTPUT_DIR.Data(), run_num, particle.Data());
    if (std::ifstream(output_path.Data())) std::remove(output_path.Data());
    TFile* fout = new TFile(output_path.Data(), "RECREATE");

    // +-------------------+
    // | prepare histogram |
    // +-------------------+    
    // -- tdc ----------
    TH1D *h_t0_tdc[2][conf.num_of_ch.at("t0")];
    for (Int_t i = 0; i < conf.num_of_ch.at("t0"); i++ ) h_t0_tdc[0][i] = (TH1D*)f->Get(Form("T0_TDC_seg%dU_%s", i, particle.Data()));
    for (Int_t i = 0; i < conf.num_of_ch.at("t0"); i++ ) h_t0_tdc[1][i] = (TH1D*)f->Get(Form("T0_TDC_seg%dD_%s", i, particle.Data()));
    // -- adc ----------
    TH1D *h_t0_adc[2][conf.num_of_ch.at("t0")];
    for (Int_t i = 0; i < conf.num_of_ch.at("t0"); i++ ) h_t0_adc[0][i] = (TH1D*)f->Get(Form("T0_ADC_seg%dU_%s", i, particle.Data()));
    for (Int_t i = 0; i < conf.num_of_ch.at("t0"); i++ ) h_t0_adc[1][i] = (TH1D*)f->Get(Form("T0_ADC_seg%dD_%s", i, particle.Data()));

    // -- set tdc range ----------
    TH1D *h_sum_tdc = (TH1D*)h_t0_tdc[0][0]->Clone("h_sum_tdc");
    h_sum_tdc->Reset(); 
    for (Int_t i = 0; i < conf.num_of_ch.at("t0"); i++) {
        h_sum_tdc->Add(h_t0_tdc[0][i]);
        h_sum_tdc->Add(h_t0_tdc[1][i]);
    }
    ana_helper::set_tdc_search_range(h_sum_tdc);

    // +--------------+
    // | fit and plot |
    // +--------------+
    // -- prepare pdf -----
    Int_t nth_pad = 1;
    Int_t rows = 2, cols = 2;
    Int_t max_pads = rows * cols;
    TString pdf_path = Form("%s/img/run%05d_T0_HDPRM_%s.pdf", OUTPUT_DIR.Data(), run_num, particle.Data());

    // -- container -----
    std::vector<FitResult> adc_up;
    std::vector<FitResult> adc_down;
    std::vector<FitResult> tdc_up;
    std::vector<FitResult> tdc_down;

    auto c_t0 = new TCanvas("t0", "", 1500, 1200);
    c_t0->Divide(cols, rows);
    c_t0->Print(pdf_path + "["); // start
    nth_pad = 1;
    for (Int_t i = 0; i < conf.num_of_ch.at("t0"); i++) {
        if (nth_pad > max_pads) {
            c_t0->Print(pdf_path);
            c_t0->Clear();
            c_t0->Divide(cols, rows);
            nth_pad = 1;
        }

        FitResult result;
        // -- UP -----
        result = ana_helper::tdc_fit(h_t0_tdc[0][i], c_t0, nth_pad);
        tdc_up.push_back(result);
        nth_pad++;

        result = ana_helper::t0_adc_fit(h_t0_adc[0][i], c_t0, nth_pad, conf.t0_ped_mip_distance[0][i]);
        adc_up.push_back(result);
        nth_pad++;

        // -- DOWN -----
        result = ana_helper::tdc_fit(h_t0_tdc[1][i], c_t0, nth_pad);
        tdc_down.push_back(result);
        nth_pad++;

        result = ana_helper::t0_adc_fit(h_t0_adc[1][i], c_t0, nth_pad, conf.t0_ped_mip_distance[1][i]);
        adc_down.push_back(result);
        nth_pad++;        
    }
    c_t0->Print(pdf_path);
    c_t0->Print(pdf_path + "]"); // end
    delete c_t0;

    // +-------+
    // | Write |
    // +-------+
    TTree* tree = new TTree("tree", "");
    Int_t ch;
    std::vector<Double_t> adc_p0_val, adc_p1_val, tdc_p0_val; 
    std::vector<Double_t> adc_p0_err, adc_p1_err, tdc_p0_err; 
    tree->Branch("ch", &ch, "ch/I");
    tree->Branch("adc_p0_val", &adc_p0_val);
    tree->Branch("adc_p1_val", &adc_p1_val);
    tree->Branch("tdc_p0_val", &tdc_p0_val);
    tree->Branch("adc_p0_err", &adc_p0_err);
    tree->Branch("adc_p1_err", &adc_p1_err);
    tree->Branch("tdc_p0_err", &tdc_p0_err);
    
    for (Int_t i = 0; i < conf.num_of_ch.at("t0"); i++) {
        ch = i;
        adc_p0_val.clear(); adc_p1_val.clear(); tdc_p0_val.clear();
        adc_p0_err.clear(); adc_p1_err.clear(); tdc_p0_err.clear();

        // -- pedestal -----
        adc_p0_val.push_back(adc_up[i].par[1]);
        adc_p0_val.push_back(adc_down[i].par[1]);
        adc_p0_err.push_back(adc_up[i].err[1]);
        adc_p0_err.push_back(adc_down[i].err[1]);
    
        // -- mip -----
        adc_p1_val.push_back(adc_up[i].par[4]);
        adc_p1_val.push_back(adc_down[i].par[4]);
        adc_p1_err.push_back(adc_up[i].err[4]);
        adc_p1_err.push_back(adc_down[i].err[4]);
    
        // -- tdc -----
        tdc_p0_val.push_back(tdc_up[i].par[1]);
        tdc_p0_val.push_back(tdc_down[i].par[1]);
        tdc_p0_err.push_back(tdc_up[i].err[1]);
        tdc_p0_err.push_back(tdc_down[i].err[1]);
    
        tree->Fill();
    }
    
    fout->cd();
    tree->Write();
    fout->Close(); 
}

Int_t main(int argc, char** argv) {

    // -- check argments -----
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <root file path> <particle>" << std::endl;
        return 1;
    }
    TString path = argv[1];
    TString particle = argv[2];
    if (particle != "Pi" && particle != "K") {
        std::cerr << "Error: Unexpected particle name: " << particle << std::endl;
        return 1;
    }

    conf.detector = "t0";
    analyze(path, particle);
    return 0;
}