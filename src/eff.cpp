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
#include <TGraph.h>
#include <TPolyLine.h>

// Custom headers
#include "config.h"
#include "ana_helper.h"
#include "paths.h"
#include "progress_bar.h"

Config& conf = Config::getInstance();

static Double_t width = 120.0, height = 26.0;

void adjust_range(TH2D *h)
{
    h->GetXaxis()->SetRangeUser(
        -150.0, 
        150.0
    );
    h->GetYaxis()->SetRangeUser(
        -75.0, 
        75.0
    );
}

Int_t check_inside_box(Double_t x, Double_t y)
{
    // Double_t width = 114.0, height = 20.0;
    Double_t theta = 1.0 * TMath::Pi() / 180.0;
    Double_t x_shift = 20.0, y_shift = 0.0;

    Int_t hit_seg = -1;
    for (Int_t i = 0; i < 4; ++i) {
        Double_t mean_x = 0.0, mean_y = 26.0*(2.0-i) - 13.0;
        // 回転・平行移動された座標系における相対位置
        Double_t dx = x - (mean_x + x_shift);
        Double_t dy = y - (mean_y + y_shift);

        // ローカル座標へ逆回転（global → box local）
        Double_t x_local = dx * cos(theta) + dy * sin(theta);
        Double_t y_local = -dx * sin(theta) + dy * cos(theta);

        // 判定（長方形内かどうか）
        bool inside = (std::abs(x_local) <= width / 2.0) && (std::abs(y_local) <= height / 2.0);

        if (inside) {
            hit_seg = 3 - i;
        }
    }
    return hit_seg;
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
    // -- hodo -----
    TString root_file_path_hodo = Form("%s/hodo/hodo_run%05d.root", DATA_DIR.Data(), run_num);
    auto *f_hodo = new TFile( root_file_path_hodo.Data() );
    if (!f_hodo || f_hodo->IsZombie()) {
        std::cerr << "Error: Could not open file : " << root_file_path_hodo << std::endl;
        return;
    }
    TTreeReader reader_hodo("hodo", f_hodo);
    Int_t total_entry = reader_hodo.GetEntries();
    TTreeReaderValue<unsigned int> run_number(reader_hodo, "run_number");
    TTreeReaderValue<unsigned int> evnum_hodo(reader_hodo, "event_number");

    TTreeReaderValue<std::vector<Double_t>> t0_raw_seg(reader_hodo, "t0_raw_seg");
    TTreeReaderValue<std::vector<std::vector<Double_t>>> t0_tdc(reader_hodo, "t0_tdc_s");

    TTreeReaderValue<std::vector<std::vector<Double_t>>> bac_tdc(reader_hodo, "bac_tdc_u");

    TTreeReaderValue<std::vector<std::vector<Double_t>>> sac_tdc(reader_hodo, "sac_tdc_u");

    TTreeReaderValue<std::vector<Double_t>> kvc2_raw_seg(reader_hodo, "kvc2_raw_seg");
    TTreeReaderValue<std::vector<Double_t>> kvc2_hit_seg(reader_hodo, "kvc2_hit_seg");
    TTreeReaderValue<std::vector<Double_t>> kvc2_adc(reader_hodo, "kvc2_adc_s");
    TTreeReaderValue<std::vector<std::vector<Double_t>>> kvc2_tdc(reader_hodo, "kvc2_tdc_s");

    TTreeReaderValue<std::vector<Double_t>> bh2_raw_seg(reader_hodo, "bh2_raw_seg");
    TTreeReaderValue<std::vector<std::vector<Double_t>>> bh2_tdc(reader_hodo, "bh2_tdc_s");

    TTreeReaderValue<std::vector<Double_t>> htof_raw_seg(reader_hodo, "htof_raw_seg");
    TTreeReaderValue<std::vector<std::vector<Double_t>>> htof_tdc(reader_hodo, "htof_tdc_s");

    TTreeReaderValue<Double_t> btof0(reader_hodo, "btof0");

    // -- bcout-----
    TString root_file_path_bcout = Form("%s/bcout/bcout_run%05d.root", DATA_DIR.Data(), run_num);
    auto *f_bcout = new TFile( root_file_path_bcout.Data() );
    if (!f_bcout || f_bcout->IsZombie()) {
        std::cerr << "Error: Could not open file : " << root_file_path_bcout << std::endl;
        return;
    }
    TTreeReader reader_bcout("bcout", f_bcout);
    TTreeReaderValue<unsigned int> evnum_bcout(reader_bcout, "event_number");

    TTreeReaderValue<std::vector<Double_t>> x0(reader_bcout, "x0");
    TTreeReaderValue<std::vector<Double_t>> u0(reader_bcout, "u0");
    TTreeReaderValue<std::vector<Double_t>> y0(reader_bcout, "y0");
    TTreeReaderValue<std::vector<Double_t>> v0(reader_bcout, "v0");
    

    // +-------------------+
    // | prepare histogram |
    // +-------------------+
    // -- T0 -----    
    std::vector<HistPair> h_t0t;
    for (Int_t ch = 0; ch < conf.num_of_ch.at("t0"); ch++) {
        TString name  = Form("T0t_%d_%d", run_num, ch + 1);
        TString title = Form("run%05d T0(TDC) ch%d;TDC;", run_num, ch + 1);
        h_t0t.emplace_back(name, title, conf.adjust_tdc_bin_num, conf.tdc_min, conf.tdc_max);
    }

    // -- BAC -----    
    std::vector<HistPair> h_bact;
    for (Int_t ch = 0; ch < conf.num_of_ch.at("bacs"); ch++) {
        TString name  = Form("BACt_%d_%d", run_num, ch + 1);
        TString title = Form("run%05d BAC(TDC) ch%d;TDC;", run_num, ch + 1);
        h_bact.emplace_back(name, title, conf.adjust_tdc_bin_num, conf.tdc_min, conf.tdc_max);
    }

    // -- SAC -----    
    std::vector<HistPair> h_sact;
    for (Int_t ch = 0; ch < conf.num_of_ch.at("sacs"); ch++) {
        TString name  = Form("SACt_%d_%d", run_num, ch + 1);
        TString title = Form("run%05d SAC(TDC) ch%d;TDC;", run_num, ch + 1);
        h_sact.emplace_back(name, title, conf.adjust_tdc_bin_num, conf.tdc_min, conf.tdc_max);
    }

    // -- KVC2 -----    
    std::vector<HistPair> h_kvc2a;
    std::vector<HistPair> h_kvc2t;
    for (Int_t ch = 0; ch < conf.num_of_ch.at("kvc2"); ch++) {
        {
            TString name  = Form("KVC2a_%d_%d", run_num, ch + 1);
            TString title = Form("run%05d KVC2(ADC) ch%d;ADC;", run_num, ch + 1);
            h_kvc2a.emplace_back(name, title, conf.adjust_adc_bin_num, conf.adc_min, conf.adc_max);
        }
        {
            TString name  = Form("KVC2t_%d_%d", run_num, ch + 1);
            TString title = Form("run%05d KVC2(TDC) ch%d;TDC;", run_num, ch + 1);
            h_kvc2t.emplace_back(name, title, conf.adjust_tdc_bin_num, conf.tdc_min, conf.tdc_max);
        }
    }

    // -- BH2 -----    
    std::vector<HistPair> h_bh2t;
    for (Int_t ch = 0; ch < conf.num_of_ch.at("bh2"); ch++) {
        TString name  = Form("BH2t_%d_%d", run_num, ch + 1);
        TString title = Form("run%05d BH2(TDC) ch%d;TDC;", run_num, ch + 1);
        h_bh2t.emplace_back(name, title, conf.adjust_tdc_bin_num, conf.tdc_min, conf.tdc_max);
    }

    // -- KVC2 -----    
    std::vector<HistPair2D> h_kvc2_seg_profile;
    for (Int_t ch = 0; ch < conf.num_of_ch.at("kvc2"); ch++) {
        TString name  = Form("KVC2_profile_%d_%d", run_num, ch + 1);
        TString title = Form("run%05d KVC2 seg%d;x;y;", run_num, ch + 1);
        h_kvc2_seg_profile.emplace_back(name, title, 600, -150.0, 150.0, 280, 70.0, 70.0);
    }

    // -- HTOF -----    
    std::vector<HistPair2D> h_htof_seg_profile;
    for (Int_t ch = 0; ch < conf.num_of_ch.at("htof"); ch++) {
        TString name  = Form("HTOF_profile_%d_%d", run_num, ch + 1);
        TString title = Form("run%05d HTOF seg%d;x;y;", run_num, ch + 1);
        h_htof_seg_profile.emplace_back(name, title, 300, -150.0, 150.0, 140, 70.0, 70.0);
    }

    // -- T0 -----    
    std::vector<HistPair2D> h_t0_seg_profile;
    for (Int_t ch = 0; ch < conf.num_of_ch.at("t0")+1; ch++) {
        TString name  = Form("T0_profile_%d_%d", run_num, ch + 1);
        TString title = Form("run%05d T0 seg%d;x;y;", run_num, ch + 1);
        h_t0_seg_profile.emplace_back(name, title, 300, -150.0, 150.0, 140, 70.0, 70.0);
    }
    
    // -- BH2 -----    
    std::vector<HistPair2D> h_bh2_seg_profile;
    for (Int_t ch = 0; ch < conf.num_of_ch.at("bh2"); ch++) {
        TString name  = Form("BH2_profile_%d_%d", run_num, ch + 1);
        TString title = Form("run%05d BH2 seg%d;x;y;", run_num, ch + 1);
        h_bh2_seg_profile.emplace_back(name, title, 600, -150.0, 150.0, 280, 70.0, 70.0);
    }

    // -- To x BH2 -----    
    std::vector<HistPair2D> h_t0_bh2_profile;
    for (Int_t ch_t0 = 0; ch_t0 < conf.num_of_ch.at("t0"); ch_t0++) {
        for (Int_t ch_bh2 = 0; ch_bh2 < conf.num_of_ch.at("bh2"); ch_bh2++) {
            TString name  = Form("T0xBH2_profile_%d_%d", run_num, ch_t0*conf.num_of_ch.at("bh2") + ch_bh2);
            TString title = Form("run%05d T0 seg%d BH2 seg%d;x;y;", run_num, ch_t0+1, ch_bh2 + 1);
            h_t0_bh2_profile.emplace_back(name, title, 600, -150.0, 150.0, 280, 70.0, 70.0);
        }
    }    
    

    // -- tracking -----
    HistPair2D h_kvc2_profile(
        Form("KVC2_profile_%d", run_num),
        Form("run%05d KVC2_profile;x;y", run_num),
        600, -150.0, 150.0,
        280, 70.0, 70.0
    );

    

    // +------------------+
    // | Fill event (1st) |
    // +------------------+
    Int_t evnum = 0;
    reader_hodo.Restart();
    reader_bcout.Restart();
    while (reader_hodo.Next() && reader_bcout.Next()){ displayProgressBar(++evnum, total_entry);
        if (*evnum_hodo != *evnum_bcout) {
            std::cerr << "Error: event numbers are not mutched; hodo : " << *evnum_hodo << ", bcout : " << *evnum_bcout << std::endl;
            return;
        }

        if (*btof0 > conf.btof_threshold) continue; 
    
        // -- T0 -----
        for (Int_t i = 0, n_i = (*t0_raw_seg).size(); i < n_i; i++) {
            Int_t index = static_cast<Int_t>((*t0_raw_seg)[i]);
            if (0 <= index && index < conf.num_of_ch.at("t0"))
                for (Int_t j = 0, n_j = (*t0_tdc)[index].size(); j < n_j; j++) 
                    h_t0t[index].raw->Fill((*t0_tdc)[index][j]);
        }

        // -- BAC -----
        for (Int_t i = 0, n_i = (*bac_tdc).size(); i < n_i; i++) {
            for (Int_t j = 0, n_j = (*bac_tdc)[i].size(); j < n_j; j++) {
                h_bact[0].raw->Fill((*bac_tdc)[i][j]);
            }
        }

        // -- SAC -----
        for (Int_t i = 0, n_i = (*sac_tdc).size(); i < n_i; i++) {
            for (Int_t j = 0, n_j = (*sac_tdc)[i].size(); j < n_j; j++) {
                h_sact[0].raw->Fill((*sac_tdc)[i][j]);
            }
        }

        // -- KVC -----
        for (Int_t i = 0, n_i = (*kvc2_raw_seg).size(); i < n_i; i++) {
            Int_t index = static_cast<Int_t>((*kvc2_raw_seg)[i]);
            if (0 <= index && index < conf.num_of_ch.at("kvc2")) {
                h_kvc2a[index].raw->Fill((*kvc2_adc)[index]);
                for (Int_t j = 0, n_j = (*kvc2_tdc)[index].size(); j < n_j; j++) 
                    h_kvc2t[index].raw->Fill((*kvc2_tdc)[index][j]);
            }
        }
        
        // -- BH2 -----
        for (Int_t i = 0, n_i = (*bh2_raw_seg).size(); i < n_i; i++) {
            Int_t index = static_cast<Int_t>((*bh2_raw_seg)[i]);
            if (0 <= index && index < conf.num_of_ch.at("bh2"))
                for (Int_t j = 0, n_j = (*bh2_tdc)[index].size(); j < n_j; j++) 
                    h_bh2t[index].raw->Fill((*bh2_tdc)[index][j]);
        }

        // -- tracking -----
        for (Int_t i = 0, n = (*x0).size(); i < n; i++) {
            Double_t x = (*x0)[i] + (*u0)[i]*conf.kvc2_pos_z;
            Double_t y = (*y0)[i] + (*v0)[i]*conf.kvc2_pos_z;
            h_kvc2_profile.raw->Fill(x, y);
        }

    }

    // +--------------+
    // | Fit and Plot |
    // +--------------+
    // メインフレームを作成
    TGMainFrame *main = new TGMainFrame(gClient->GetRoot(), 1000, 800);
    // Handle the window close event to terminate the application
    main->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
    // タブウィジェットを作成
    TGTab *tab = new TGTab(main, 1000, 800);

    // -- T0 -----
    TCanvas *c_t0 = ana_helper::add_tab(tab, "t0");
    c_t0->Divide(3, 2);
    std::vector<FitResult> result_t0;
    conf.detector = "t0";
    for (Int_t ch = 0; ch < conf.num_of_ch.at("t0"); ch++) {
        FitResult tmp_result = ana_helper::tdc_fit(h_t0t[ch].raw, c_t0, ch+1);
        result_t0.push_back(tmp_result);
        c_t0->Update();
    }

    // -- BAC -----
    TCanvas *c_bac = ana_helper::add_tab(tab, "bac");
    std::vector<FitResult> result_bac;
    conf.detector = "bacs";
    for (Int_t ch = 0; ch < conf.num_of_ch.at("bacs"); ch++) {
        FitResult tmp_result = ana_helper::tdc_fit(h_bact[ch].raw, c_bac, ch+1);
        result_bac.push_back(tmp_result);
    }
    

    // -- SAC -----
    TCanvas *c_sac = ana_helper::add_tab(tab, "sac");
    std::vector<FitResult> result_sac;
    conf.detector = "sacs";
    for (Int_t ch = 0; ch < conf.num_of_ch.at("sacs"); ch++) {
        FitResult tmp_result = ana_helper::tdc_fit(h_sact[ch].raw, c_sac, ch+1);
        result_sac.push_back(tmp_result);
    }


    // -- KVC2 -----
    TCanvas *c_kvc2 = ana_helper::add_tab(tab, "kvc2");
    c_kvc2->Divide(2, 2);
    std::vector<FitResult> result_kvc2;
    conf.detector = "kvc2";
    for (Int_t ch = 0; ch < conf.num_of_ch.at("kvc2"); ch++) {
        FitResult tmp_result = ana_helper::tdc_fit(h_kvc2t[ch].raw, c_kvc2, ch+1);
        result_kvc2.push_back(tmp_result);
    }

    // -- BH2 -----
    TCanvas *c_bh2 = ana_helper::add_tab(tab, "bh2");
    c_bh2->Divide(4, 3);
    std::vector<FitResult> result_bh2;
    conf.detector = "bh2";
    for (Int_t ch = 0; ch < conf.num_of_ch.at("bh2"); ch++) {
        FitResult tmp_result = ana_helper::tdc_fit(h_bh2t[ch].raw, c_bh2, ch+1);
        result_bh2.push_back(tmp_result);
    }


    for (Int_t index = 0; index<4; index++) {
        Double_t lower  = result_kvc2[index].par[1] - 10.0*result_kvc2[index].par[2];
        Double_t upper = result_kvc2[index].par[1] + 10.0*result_kvc2[index].par[2];
    }
    
    // +------------------+
    // | Fill event (1st) |
    // +------------------+
    evnum = 0;
    std::vector<Int_t> n_kaon(conf.num_of_ch.at("kvc2"), 0);
    std::vector<Int_t> n_trig(conf.num_of_ch.at("kvc2"), 0);
    Double_t n_sigma = 5.0;
    reader_hodo.Restart();
    reader_bcout.Restart();
    while (reader_hodo.Next() && reader_bcout.Next()){ displayProgressBar(++evnum, total_entry);
        // -- T0 -----
        std::vector<Bool_t> t0_hitseg(conf.num_of_ch.at("t0"), false);
        Bool_t flag_t0 = false;
        for (Int_t i = 0, n_i = (*t0_raw_seg).size(); i < n_i; i++) {
            Int_t index = static_cast<Int_t>((*t0_raw_seg)[i]);
            if (0 <= index && index < conf.num_of_ch.at("t0")) {
                for (Int_t j = 0, n_j = (*t0_tdc)[index].size(); j < n_j; j++) {
                    Double_t lower  = result_t0[index].par[1] - 5.0*result_t0[index].par[2];
                    Double_t upper = result_t0[index].par[1] + 5.0*result_t0[index].par[2];
                    if ( lower < (*t0_tdc)[index][j] && (*t0_tdc)[index][j] < upper ) {
                        flag_t0 = true;
                        t0_hitseg[index] = true;
                    }
                }
            }
        }

        // -- BAC -----
        Bool_t flag_bac = false;
        for (Int_t i = 0, n_i = (*bac_tdc).size(); i < n_i; i++) {
            for (Int_t j = 0, n_j = (*bac_tdc)[i].size(); j < n_j; j++) {
                Double_t lower  = result_bac[0].par[1] - 5.0*result_bac[0].par[2];
                Double_t upper = result_bac[0].par[1] + 5.0*result_bac[0].par[2];
                if ( lower < (*bac_tdc)[i][j] && (*bac_tdc)[i][j] < upper ) flag_bac = true;
            }
        }

        // -- SAC -----
        Bool_t flag_sac = false;
        for (Int_t i = 0, n_i = (*sac_tdc).size(); i < n_i; i++) {
            for (Int_t j = 0, n_j = (*sac_tdc)[i].size(); j < n_j; j++) {
                Double_t lower  = result_sac[0].par[1] - 5.0*result_sac[0].par[2];
                Double_t upper = result_sac[0].par[1] + 5.0*result_sac[0].par[2];
                if ( lower < (*sac_tdc)[i][j] && (*sac_tdc)[i][j] < upper ) flag_sac = true;
            }
        }

        // -- KVC -----
        std::vector<Bool_t> flag_kvc2(conf.num_of_ch.at("kvc2"), false);
        for (Int_t i = 0, n_i = (*kvc2_raw_seg).size(); i < n_i; i++) {
            Int_t index = static_cast<Int_t>((*kvc2_raw_seg)[i]);
            if (0 <= index && index < conf.num_of_ch.at("kvc2")) {
                for (Int_t j = 0, n_j = (*kvc2_tdc)[i].size(); j < n_j; j++) {
                    Double_t lower  = result_kvc2[index].par[1] - 10.0*result_kvc2[index].par[2];
                    Double_t upper = result_kvc2[index].par[1] + 10.0*result_kvc2[index].par[2];
                    if ( lower < (*kvc2_tdc)[i][j] && (*kvc2_tdc)[i][j] < upper ) flag_kvc2[index] = true;
                }
            }
        }
        
        // -- BH2 -----
        std::vector<Bool_t> flag_bh2(conf.num_of_ch.at("bh2"), false);
        Bool_t flag_bh2_narrow = false;
        for (Int_t i = 0, n_i = (*bh2_raw_seg).size(); i < n_i; i++) {
            Int_t index = static_cast<Int_t>((*bh2_raw_seg)[i]);
            if (0 <= index && index < conf.num_of_ch.at("bh2")) {
                for (Int_t j = 0, n_j = (*bh2_tdc)[index].size(); j < n_j; j++) {
                    Double_t lower  = result_bh2[index].par[1] - 5.0*result_bh2[index].par[2];
                    Double_t upper = result_bh2[index].par[1] + 5.0*result_bh2[index].par[2];
                    if ( lower < (*bh2_tdc)[i][j] && (*bh2_tdc)[i][j] < upper ) {
                        flag_bh2[index] = true;
                        if (1 < index && index < 9) flag_bh2_narrow = true;
                    }
                } 
            }   
        }

        // -- tracking -----
        std::vector<Bool_t> hitseg(conf.num_of_ch.at("kvc2"), false);
        for (Int_t i = 0, n = (*x0).size(); i < n; i++) {
            Double_t x = (*x0)[i] + (*u0)[i]*conf.kvc2_pos_z;
            Double_t y = (*y0)[i] + (*v0)[i]*conf.kvc2_pos_z;
            Int_t hit_seg_id = check_inside_box(x, y);
            if (hit_seg_id != -1) hitseg[hit_seg_id] = true;
        }

        // -- HTOF -----
        std::vector<Bool_t> flag_htof(conf.num_of_ch.at("htof"), false);
        for (Int_t i = 0, n_i = (*htof_raw_seg).size(); i < n_i; i++) {
            Int_t index = static_cast<Int_t>((*htof_raw_seg)[i]);
            if (0 <= index && index < conf.num_of_ch.at("htof")) {
                for (Int_t j = 0, n_j = (*htof_tdc)[i].size(); j < n_j; j++) {
                    Double_t lower = 690000.0;
                    Double_t upper = 710000.0;
                    if ( lower < (*htof_tdc)[i][j] && (*htof_tdc)[i][j] < upper ) flag_htof[index] = true;
                }
            }
        }
        
        // for (Int_t ch = 0; ch < conf.num_of_ch.at("kvc2"); ch++) {
        //     // if (flag_kvc2[ch]) {
        //     Int_t index = static_cast<Int_t>((*kvc2_hit_seg)[ch]);
        //     if (index == ch) {
        //         for (Int_t i = 0, n = (*x0).size(); i < n; i++) {
        //             Double_t x = (*x0)[i] + (*u0)[i]*conf.kvc2_pos_z;
        //             Double_t y = (*y0)[i] + (*v0)[i]*conf.kvc2_pos_z;
        //             h_kvc2_seg_profile[ch].trig->Fill(x, y);
        //         }
        //     }
        // }

        for (Int_t ch = 0; ch < conf.num_of_ch.at("kvc2"); ch++) {
            if (flag_kvc2[ch]) {
                for (Int_t i = 0, n = (*x0).size(); i < n; i++) {
                    Double_t x = (*x0)[i] + (*u0)[i]*conf.kvc2_pos_z;
                    Double_t y = (*y0)[i] + (*v0)[i]*conf.kvc2_pos_z;
                    h_kvc2_seg_profile[ch].trig->Fill(x, y);
                }
            }
        }

        // for (Int_t ch = 0; ch < conf.num_of_ch.at("htof"); ch++) {
        for (Int_t ch = 0; ch < 6; ch++) {
            if (flag_htof[ch]) {
                for (Int_t i = 0, n = (*x0).size(); i < n; i++) {
                    Double_t x = (*x0)[i] + (*u0)[i]*conf.htof_pos_z;
                    Double_t y = (*y0)[i] + (*v0)[i]*conf.htof_pos_z;
                    h_htof_seg_profile[ch].trig->Fill(x, y);
                    h_htof_seg_profile[6].trig->Fill(x, y);
                }
            }
        }

        for (Int_t ch = 0; ch < conf.num_of_ch.at("bh2"); ch++) {
            if (flag_bh2[ch]) {
                for (Int_t i = 0, n = (*x0).size(); i < n; i++) {
                    Double_t x = (*x0)[i] + (*u0)[i]*conf.bh2_pos_z;
                    Double_t y = (*y0)[i] + (*v0)[i]*conf.bh2_pos_z;
                    h_bh2_seg_profile[ch].trig->Fill(x, y);
                }
            }
        }

        for (Int_t ch = 0; ch < conf.num_of_ch.at("t0"); ch++) {
            if (t0_hitseg[ch]) {
                for (Int_t i = 0, n = (*x0).size(); i < n; i++) {
                    Double_t x = (*x0)[i] + (*u0)[i]*conf.t0_pos_z;
                    Double_t y = (*y0)[i] + (*v0)[i]*conf.t0_pos_z;
                    h_t0_seg_profile[ch].trig->Fill(x, y);
                    h_t0_seg_profile[conf.num_of_ch.at("t0")].trig->Fill(x, y);
                }
            }
        }


        for (Int_t i = 0, n = (*x0).size(); i < n; i++) {
            Double_t x = (*x0)[i] + (*u0)[i]*conf.htof_pos_z;
            Double_t y = (*y0)[i] + (*v0)[i]*conf.htof_pos_z;
            
            for (Int_t ch_t0 = 0; ch_t0 < conf.num_of_ch.at("t0"); ch_t0++) {
                for (Int_t ch_bh2 = 0; ch_bh2 < conf.num_of_ch.at("bh2"); ch_bh2++) {
                    if (t0_hitseg[ch_t0] && flag_bh2[ch_bh2]) {
                        h_t0_bh2_profile[ch_t0*conf.num_of_ch.at("bh2") + ch_bh2].trig->Fill(x, y);
                    }
                }
            } 
        }

        if ( flag_t0 && !flag_bac && !flag_sac && flag_bh2_narrow && *btof0 < conf.btof_threshold) {         
            for (Double_t ch = 0; ch < conf.num_of_ch.at("kvc2"); ch++) {
                if (hitseg[ch]) {
                    n_kaon[ch]++;
                    if (flag_kvc2[ch]) {
                        n_trig[ch]++;
                        h_kvc2a[ch].trig->Fill((*kvc2_adc)[ch]);
                    }
                }
            }
        }
    }


    // -- BH2 -----
    TCanvas *c_kvc2_profile_raw = ana_helper::add_tab(tab, "profile raw");
    c_kvc2_profile_raw->cd(1);
    h_kvc2_profile.raw->Draw("colz");

    TCanvas *c_kvc2_profile_trig = ana_helper::add_tab(tab, "profile trig");
    c_kvc2_profile_trig->cd(1);
    h_kvc2_profile.trig->Draw("colz");

    TCanvas *c_kvc2_profile = ana_helper::add_tab(tab, "profile");
    c_kvc2_profile->Divide(2, 2);
    
    for (Int_t ch = 0; ch < conf.num_of_ch.at("kvc2"); ch++) {
        c_kvc2_profile->cd(ch+1);
        h_kvc2_seg_profile[ch].trig->Draw("colz");

        // 軸に沿った回転＋範囲カット
        // Double_t width = 120.0, height = 26.0;
        
        Double_t theta = 1.0 * TMath::Pi() / 180.0;
        Double_t x_shift = 20.0, y_shift = 0.0;

        for (Int_t i = 0; i < 5; ++i) {
            Double_t mean_x = 0.0, mean_y = 26.0*(2.0-i) - 13.0;

            // 長さ（描画範囲、適宜調整）
            Double_t L = 200.0;

            Double_t x1 = (mean_x - x_shift) - L * cos(theta);
            Double_t y1 = (mean_y - y_shift) - L * sin(theta);
            Double_t x2 = (mean_x - x_shift) + L * cos(theta);
            Double_t y2 = (mean_y - y_shift) + L * sin(theta);

            // 線を描画
            // TLine* line = new TLine(x1, y1, x2, y2);
            // line->SetLineColor(kBlue);
            // line->SetLineWidth(2);
            // line->Draw("SAME");

            if (i< 4) {
                    // 4つの頂点（中心を原点としたローカル座標）
        std::vector<std::pair<Double_t, Double_t>> corners = {
            {-width/2, -height/2},
            { width/2, -height/2},
            { width/2,  height/2},
            {-width/2,  height/2},
            {-width/2, -height/2}  // 閉じるために最初の点をもう一度
        };

        // 回転・平行移動して世界座標へ変換
        std::vector<Double_t> x_coords, y_coords;
        for (auto& [x_local, y_local] : corners) {
            Double_t x_rot = x_local * cos(theta) - y_local * sin(theta);
            Double_t y_rot = x_local * sin(theta) + y_local * cos(theta);
            x_coords.push_back(x_rot + mean_x + x_shift);
            y_coords.push_back(y_rot + mean_y + y_shift);
        }

        // 描画
        TPolyLine* box = new TPolyLine(x_coords.size(), x_coords.data(), y_coords.data());
        box->SetLineColor(kRed);
        box->SetLineWidth(2);
        box->SetLineStyle(1);  // 実線
        box->Draw("L SAME");

            }

        }

        // for (Int_t i = 0; i < n; ++i) {
        //     Double_t dx = x[i] - mean_x;
        //     Double_t dy = y[i] - mean_y;

        //     Double_t x_rot =  dx * cos(theta) + dy * sin(theta);
        //     Double_t y_rot = -dx * sin(theta) + dy * cos(theta);

        //     if (std::abs(x_rot) < width / 2.0 && std::abs(y_rot) < height / 2.0) {
        //         // 通過（長方形内）
        //     }
        // }


        // // 4つの頂点（中心を原点としたローカル座標）
        // std::vector<std::pair<Double_t, Double_t>> corners = {
        //     {-width/2, -height/2},
        //     { width/2, -height/2},
        //     { width/2,  height/2},
        //     {-width/2,  height/2},
        //     {-width/2, -height/2}  // 閉じるために最初の点をもう一度
        // };

        // // 回転・平行移動して世界座標へ変換
        // std::vector<Double_t> x_coords, y_coords;
        // for (auto& [x_local, y_local] : corners) {
        //     Double_t x_rot = x_local * cos(theta) - y_local * sin(theta);
        //     Double_t y_rot = x_local * sin(theta) + y_local * cos(theta);
        //     x_coords.push_back(x_rot + mean_x);
        //     y_coords.push_back(y_rot + mean_y);
        // }

        // // 描画
        // TPolyLine* box = new TPolyLine(x_coords.size(), x_coords.data(), y_coords.data());
        // box->SetLineColor(kRed);
        // box->SetLineWidth(2);
        // box->SetLineStyle(1);  // 実線
        // box->Draw("L SAME");

        // // TBox(x1, y1, x2, y2) → 左下と右上の座標
        // TBox *box = new TBox(-60.0, 26.0*(1.0-ch), 60.0, 26.0*(2.0-ch));        
        // // 枠線の設定（赤・塗りなし）
        // box->SetLineColor(kRed);  // 赤色の枠線
        // box->SetLineWidth(2);     // 線の太さ（任意）
        // box->SetFillStyle(0);     // 塗りつぶしなし
        // // 描画
        // box->Draw("same");
    }


    TCanvas *c_htof_profile = ana_helper::add_tab(tab, "profile");
    c_htof_profile->Divide(2, 3);
    
    for (Int_t ch = 0; ch < 6; ch++) {
        c_htof_profile->cd(ch+1);
        adjust_range(h_htof_seg_profile[ch].trig);
        h_htof_seg_profile[ch].trig->Draw("colz");
    }

    TCanvas *c_bh2_profile = ana_helper::add_tab(tab, "profile");
    c_bh2_profile->Divide(3, 4);
    
    for (Int_t ch = 0; ch < conf.num_of_ch.at("bh2"); ch++) {
        c_bh2_profile->cd(ch+1);
        h_bh2_seg_profile[ch].trig->Draw("colz");
    }

    TCanvas *c_t0_profile = ana_helper::add_tab(tab, "profile");
    c_t0_profile->Divide(3, 2);
    
    for (Int_t ch = 0; ch < conf.num_of_ch.at("t0")+1; ch++) {
        c_t0_profile->cd(ch+1);
        h_t0_seg_profile[ch].trig->Draw("colz");
    }

    TCanvas *c_t0_bh2_profile = ana_helper::add_tab(tab, "profile");
    c_t0_bh2_profile->Divide(7, 8);
    for (Int_t ch_t0 = 0; ch_t0 < conf.num_of_ch.at("t0"); ch_t0++) {
        for (Int_t ch_bh2 = 0; ch_bh2 < conf.num_of_ch.at("bh2"); ch_bh2++) {
            Int_t i = ch_t0*conf.num_of_ch.at("bh2") + ch_bh2;
            c_t0_bh2_profile->cd(i+1);
            h_t0_bh2_profile[i].trig->GetXaxis()->SetRangeUser(-150.0, 150.0);
            h_t0_bh2_profile[i].trig->GetYaxis()->SetRangeUser(-150.0, 150.0);            
            h_t0_bh2_profile[i].trig->Draw("colz");

                    // 軸に沿った回転＋範囲カット
        // Double_t width = 120.0, height = 26.0;
        Double_t width = 116.0, height = 24.0;
        
        Double_t theta = 1.0 * TMath::Pi() / 180.0;
        Double_t x_shift = 20.0, y_shift = 0.0;

        for (Int_t i = 0; i < 5; ++i) {
            Double_t mean_x = 0.0, mean_y = 26.0*(2.0-i) - 13.0;

            // 長さ（描画範囲、適宜調整）
            Double_t L = 200.0;

            Double_t x1 = (mean_x - x_shift) - L * cos(theta);
            Double_t y1 = (mean_y - y_shift) - L * sin(theta);
            Double_t x2 = (mean_x - x_shift) + L * cos(theta);
            Double_t y2 = (mean_y - y_shift) + L * sin(theta);

            // 線を描画
            // TLine* line = new TLine(x1, y1, x2, y2);
            // line->SetLineColor(kBlue);
            // line->SetLineWidth(2);
            // line->Draw("SAME");

            if (i< 4) {
                    // 4つの頂点（中心を原点としたローカル座標）
        std::vector<std::pair<Double_t, Double_t>> corners = {
            {-width/2, -height/2},
            { width/2, -height/2},
            { width/2,  height/2},
            {-width/2,  height/2},
            {-width/2, -height/2}  // 閉じるために最初の点をもう一度
        };

        // 回転・平行移動して世界座標へ変換
        std::vector<Double_t> x_coords, y_coords;
        for (auto& [x_local, y_local] : corners) {
            Double_t x_rot = x_local * cos(theta) - y_local * sin(theta);
            Double_t y_rot = x_local * sin(theta) + y_local * cos(theta);
            x_coords.push_back(x_rot + mean_x + x_shift);
            y_coords.push_back(y_rot + mean_y + y_shift);
        }

        // 描画
        TPolyLine* box = new TPolyLine(x_coords.size(), x_coords.data(), y_coords.data());
        box->SetLineColor(kRed);
        box->SetLineWidth(2);
        box->SetLineStyle(1);  // 実線
        box->Draw("L SAME");

            }

        }
        }
    } 

    for (Double_t ch = 0; ch < conf.num_of_ch.at("kvc2"); ch++) {
        std::cout << n_kaon[ch] << ", " << n_trig[ch] << ", " << (Double_t) n_trig[ch]/n_kaon[ch] << std::endl;
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