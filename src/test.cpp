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
#include <TTreeReaderArray.h>
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

void adjust_range(TH1D *h)
{
    Double_t mean  = h->GetMean();        
    Double_t stdev = h->GetStdDev();        
    Double_t lower = std::max(mean - 3.0 * stdev, conf.adc_min);
    Double_t upper = std::min(mean + 7.0 * stdev, conf.adc_max);
    h->GetXaxis()->SetRangeUser(lower, upper);
}

void analyze(Int_t run_num){
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
    TString root_file_path = Form("%s/hodo_run%05d.root", DATA_DIR.Data(), run_num);
    auto *f = new TFile( root_file_path.Data() );
    if (!f || f->IsZombie()) {
        std::cerr << "Error: Could not open file : " << root_file_path << std::endl;
        return;
    }
    TTreeReader reader("hodo", f);
    TTreeReaderValue<unsigned int> run_number(reader, "run_number");
    TTreeReaderValue<std::vector<Double_t>> bac_raw_seg(reader, "bac_raw_seg");
    TTreeReaderValue<std::vector<Double_t>> bac_adc(reader, "bac_adc_u");

    TTreeReaderValue<std::vector<Double_t>> kvc2_raw_seg(reader, "kvc2_raw_seg");
    TTreeReaderValue<std::vector<Double_t>> kvc2_adc_a(reader, "kvc2_adc_a");
    TTreeReaderValue<std::vector<Double_t>> kvc2_adc_b(reader, "kvc2_adc_b");
    TTreeReaderValue<std::vector<Double_t>> kvc2_adc_c(reader, "kvc2_adc_c");
    TTreeReaderValue<std::vector<Double_t>> kvc2_adc_d(reader, "kvc2_adc_d");

    TTreeReaderValue<std::vector<Double_t>> kvc1_raw_seg(reader, "kvc1_raw_seg");
    TTreeReaderValue<std::vector<Double_t>> kvc1_adc_u(reader, "kvc1_adc_u");
    TTreeReaderValue<std::vector<Double_t>> kvc1_adc_d(reader, "kvc1_adc_d");

    // +-------------------+
    // | prepare histogram |
    // +-------------------+    
    // -- bac ----------
    TH1D *h_bac_adc[conf.num_of_ch.at("bac")];
    for (Int_t ch = 0; ch < conf.num_of_ch.at("bac"); ch++ ) {
        TString name  = Form("BACa_%d_%d", run_num, ch + 1);
        TString title = Form("run%05d BAC(ADC) ch%d;ADC;", run_num, ch + 1);
        h_bac_adc[ch] = new TH1D(name, title, conf.adjust_adc_bin_num, conf.adc_min, conf.adc_max);
    }   
    // -- kvc2 ----------
    TH1D *h_kvc2_adc[conf.num_of_ch.at("kvc2")][conf.num_of_UorD.at("kvc2")];
    for (Int_t ch = 0; ch < conf.num_of_ch.at("kvc2"); ch++ ) {
        for (Int_t UorD = 0; UorD < conf.num_of_UorD.at("kvc2"); UorD++ ){
            TString name  = Form("KVC2a_%d_%d_%d", run_num, ch + 1, UorD + 1);
            TString title = Form("run%05d KVC2(ADC) ch%d board%d;ADC;", run_num, ch + 1, UorD+1);
            h_kvc2_adc[ch][UorD] = new TH1D(name, title, conf.adjust_adc_bin_num, conf.adc_min, conf.adc_max);
        }
    }
    // -- kvc1 ----------
    TH1D *h_kvc1_adc[conf.num_of_ch.at("kvc1")][conf.num_of_UorD.at("kvc1")];
    for (Int_t ch = 0; ch < conf.num_of_ch.at("kvc1"); ch++ ) {
        for (Int_t UorD = 0; UorD < conf.num_of_UorD.at("kvc1"); UorD++ ){
            TString name  = Form("KVC1a_%d_%d_%d", run_num, ch + 1, UorD + 1);
            TString title = Form("run%05d KVC1(ADC) ch%d board%d;ADC;", run_num, ch + 1, UorD+1);
            h_kvc1_adc[ch][UorD] = new TH1D(name, title, conf.adjust_adc_bin_num, conf.adc_min, conf.adc_max);
        }
    }


    // +------------+
    // | Fill event |
    // +------------+
    reader.Restart();
    while (reader.Next()){
        // -- BAC -----
        for (Int_t i = 0, n = (*bac_raw_seg).size(); i < n; i++) {
            Int_t index = static_cast<Int_t>((*bac_raw_seg)[i]);
            if (0 <= index && index < conf.num_of_ch.at("bac"))
                h_bac_adc[index]->Fill((*bac_adc)[index]);
        }

        // -- KVC2 -----
        for (Int_t i = 0, n = (*kvc2_raw_seg).size(); i < n; i++) {
            Int_t index = static_cast<Int_t>((*kvc2_raw_seg)[i]);
            if (0 <= index && index < conf.num_of_ch.at("kvc2")) {
                h_kvc2_adc[index][0]->Fill((*kvc2_adc_a)[index]);
                h_kvc2_adc[index][1]->Fill((*kvc2_adc_b)[index]);
                h_kvc2_adc[index][2]->Fill((*kvc2_adc_c)[index]);
                h_kvc2_adc[index][3]->Fill((*kvc2_adc_d)[index]);
            }
        }
        
        // -- KVC1 -----
        for (Int_t i = 0, n = (*kvc1_raw_seg).size(); i < n; i++) {
            Int_t index = static_cast<Int_t>((*kvc1_raw_seg)[i]);
            if (0 <= index && index < conf.num_of_ch.at("kvc1")) {
                h_kvc1_adc[index][0]->Fill((*kvc1_adc_u)[index]);
                h_kvc1_adc[index][1]->Fill((*kvc1_adc_d)[index]);
            }
        }
        
    }

    // +------+
    // | plot |
    // +------+
    // メインフレームを作成
    TGMainFrame *main = new TGMainFrame(gClient->GetRoot(), 1000, 800);
    // Handle the window close event to terminate the application
    main->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
    // タブウィジェットを作成
    TGTab *tab = new TGTab(main, 1000, 800);

    // -- BAC -----
    TCanvas *c_bac = ana_helper::add_tab(tab, "bac");
    c_bac->Divide(2, 2);
    for (Int_t ch = 0; ch < conf.num_of_ch.at("bac"); ch++ ) {
        c_bac->cd(ch+1);
        adjust_range(h_bac_adc[ch]);
        h_bac_adc[ch]->Draw();
    }

    // -- KVC -----
    TCanvas *c_kvc2_seg1 = ana_helper::add_tab(tab, "kvc2_seg1");
    c_kvc2_seg1->Divide(2, 2);
    TCanvas *c_kvc2_seg2 = ana_helper::add_tab(tab, "kvc2_seg2");
    c_kvc2_seg2->Divide(2, 2);
    TCanvas *c_kvc2_seg3 = ana_helper::add_tab(tab, "kvc2_seg3");
    c_kvc2_seg3->Divide(2, 2);
    TCanvas *c_kvc2_seg4 = ana_helper::add_tab(tab, "kvc2_seg4");
    c_kvc2_seg4->Divide(2, 2);

    std::vector<TCanvas*> canvas_cantaier{c_kvc2_seg1, c_kvc2_seg2, c_kvc2_seg3, c_kvc2_seg4};

    for (Int_t ch = 0; ch < conf.num_of_ch.at("kvc2"); ch++ ) {
        for (Int_t UorD = 0; UorD < conf.num_of_UorD.at("kvc2"); UorD++ ) {
            canvas_cantaier[ch]->cd(UorD+1);
            adjust_range(h_kvc2_adc[ch][UorD]);
            h_kvc2_adc[ch][UorD]->Draw();
        }
    }

    // -- KVC -----
    TCanvas *c_kvc1_seg1 = ana_helper::add_tab(tab, "kvc1_seg1");
    c_kvc1_seg1->Divide(2, 1);
    TCanvas *c_kvc1_seg2 = ana_helper::add_tab(tab, "kvc1_seg2");
    c_kvc1_seg2->Divide(2, 1);
    TCanvas *c_kvc1_seg3 = ana_helper::add_tab(tab, "kvc1_seg3");
    c_kvc1_seg3->Divide(2, 1);
    TCanvas *c_kvc1_seg4 = ana_helper::add_tab(tab, "kvc1_seg4");
    c_kvc1_seg4->Divide(2, 1);

    std::vector<TCanvas*> canvas_cantaier1{c_kvc1_seg1, c_kvc1_seg2, c_kvc1_seg3, c_kvc1_seg4};

    for (Int_t ch = 0; ch < conf.num_of_ch.at("kvc1"); ch++ ) {
        for (Int_t UorD = 0; UorD < conf.num_of_UorD.at("kvc1"); UorD++ ) {
            canvas_cantaier1[ch]->cd(UorD+1);
            adjust_range(h_kvc1_adc[ch][UorD]);
            h_kvc1_adc[ch][UorD]->Draw();
        }
    }
    
    // メインフレームにタブウィジェットを追加
    main->AddFrame(tab, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
    // ウィンドウを表示
    main->MapSubwindows();
    main->Resize(main->GetDefaultSize());
    main->MapWindow();
}

Int_t main(int argc, char** argv) {

    // -- check argments -----
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <run number>" << std::endl;
        return 1;
    }
    Int_t run_num = std::atoi(argv[1]);

    TApplication *theApp = new TApplication("App", &argc, argv);    
    analyze(run_num);
    theApp->Run();

    return 0;
}