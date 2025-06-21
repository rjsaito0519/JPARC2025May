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
static Double_t theta = 0.0 * TMath::Pi() / 180.0;
static Double_t x_shift = 20.0, y_shift = 3.0;


Int_t check_inside_box(Double_t x, Double_t y)
{
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
            hit_seg = i;
        }
    }
    return hit_seg;
}

std::vector<std::vector<Int_t>> analyze(Int_t run_num){
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
    TString root_file_path_hodo = Form("%s/hodo/hodo_run%05d_Pi_t0.root", DATA_DIR.Data(), run_num);
    auto *f_hodo = new TFile( root_file_path_hodo.Data() );
    if (!f_hodo || f_hodo->IsZombie()) {
        std::cerr << "Error: Could not open file : " << root_file_path_hodo << std::endl;
        return std::vector<std::vector<Int_t>>{};
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
        return std::vector<std::vector<Int_t>>{};
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
    
    // +------------------+
    // | Fill event (1st) |
    // +------------------+
    Int_t evnum = 0;
    reader_hodo.Restart();
    reader_bcout.Restart();
    while (reader_hodo.Next() && reader_bcout.Next()){ displayProgressBar(++evnum, total_entry);
        if (*evnum_hodo != *evnum_bcout) {
            std::cerr << "Error: event numbers are not mutched; hodo : " << *evnum_hodo << ", bcout : " << *evnum_bcout << std::endl;
            return std::vector<std::vector<Int_t>>{};
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

    }

    // +--------------+
    // | Fit and Plot |
    // +--------------+
    // メインフレームを作成
    TGMainFrame *main = new TGMainFrame(gClient->GetRoot(), 1000, 800);
    // Handle the window close event to terminate the application
    // main->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
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


    // +--------------------------+
    // | prepare output root file |
    // +--------------------------+
    TString output_path = OUTPUT_DIR + "/root/kvc2_eff.root";
    if (std::ifstream(output_path.Data())) std::remove(output_path.Data());
    TFile fout(output_path.Data(), "create");
    TTree output_tree("tree", ""); 

    // -- prepare root file branch -----
    std::vector<Int_t> n_kaon(conf.num_of_ch.at("kvc2"), 0);
    std::vector<Int_t> n_trig(conf.num_of_ch.at("kvc2"), 0);
    std::vector<Double_t> eff_val, eff_err;

    output_tree.Branch("width", &width, "width/D");
    output_tree.Branch("height", &height, "height/D");
    output_tree.Branch("n_kaon", &n_kaon);
    output_tree.Branch("n_trig", &n_trig);
    output_tree.Branch("eff_val", &eff_val);
    output_tree.Branch("eff_err", &eff_err);
    
    // +------------------+
    // | Fill event (2nd) |
    // +------------------+
    for (Double_t w_i = 0.0; w_i < 10.0; w_i++) {
        for (Double_t h_i = 0.0; h_i < 10.0; h_i++) {
            std::cout << w_i << ", " << h_i << std::endl;
            
            eff_val.clear(); eff_err.clear();
            width  = 120.0 - 120.0/20.0*w_i;
            height =  26.0 -  26.0/20.0*h_i;
            
            evnum = 0;
            for (Int_t ch = 0; ch < conf.num_of_ch.at("kvc2"); ch++) {
                n_kaon[ch] = 0;
                n_trig[ch] = 0;
            }
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
                            Double_t lower = result_t0[index].par[1] - 5.0*result_t0[index].par[2];
                            Double_t upper = result_t0[index].par[1] + 5.0*result_t0[index].par[2];
                            if ( lower < (*t0_tdc)[i][j] && (*t0_tdc)[i][j] < upper ) {
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
                        Double_t lower = result_bac[0].par[1] - 5.0*result_bac[0].par[2];
                        Double_t upper = result_bac[0].par[1] + 5.0*result_bac[0].par[2];
                        if ( lower < (*bac_tdc)[i][j] && (*bac_tdc)[i][j] < upper ) flag_bac = true;
                    }
                }

                // -- SAC -----
                Bool_t flag_sac = false;
                for (Int_t i = 0, n_i = (*sac_tdc).size(); i < n_i; i++) {
                    for (Int_t j = 0, n_j = (*sac_tdc)[i].size(); j < n_j; j++) {
                        Double_t lower = result_sac[0].par[1] - 5.0*result_sac[0].par[2];
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
                            Double_t lower = result_kvc2[index].par[1] - 10.0*result_kvc2[index].par[2];
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
                            Double_t lower = result_bh2[index].par[1] - 5.0*result_bh2[index].par[2];
                            Double_t upper = result_bh2[index].par[1] + 5.0*result_bh2[index].par[2];
                            if ( lower < (*bh2_tdc)[i][j] && (*bh2_tdc)[i][j] < upper ) {
                                flag_bh2[index] = true;
                                if (0 < index && index < 10) flag_bh2_narrow = true;
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

                if ( flag_t0 && !flag_bac && !flag_sac && flag_bh2_narrow && *btof0 < conf.btof_threshold) {         
                    for (Double_t ch = 0; ch < conf.num_of_ch.at("kvc2"); ch++) {
                        if (hitseg[ch]) {
                            n_kaon[ch]++;
                            if (flag_kvc2[ch]) {
                                n_trig[ch]++;
                                // h_kvc2a[ch].trig->Fill((*kvc2_adc)[ch]);
                            }
                        }
                    }
                }
            }


            for (Double_t ch = 0; ch < conf.num_of_ch.at("kvc2"); ch++) {
                Double_t eff = static_cast<Double_t>(n_trig[ch]) / static_cast<Double_t>(n_kaon[ch]);
                Double_t err = TMath::Sqrt(eff*(1.0 - eff) / static_cast<Double_t>(n_kaon[ch]));
                std::cout << n_kaon[ch] << ", " << n_trig[ch] << ", " 
                          << std::fixed << std::setprecision(4) << eff << " +/- " << err << std::endl;
                eff_val.push_back(eff);
                eff_err.push_back(err);
            }
            output_tree.Fill();
        }
    }

    // +------------+
    // | Write data |
    // +------------+
    fout.cd(); // 明示的にカレントディレクトリを設定
    output_tree.Write();
    fout.Close();

    // // メインフレームにタブウィジェットを追加
    // main->AddFrame(tab, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
    // // ウィンドウを表示
    // main->MapSubwindows();
    // main->Resize(main->GetDefaultSize());
    // main->MapWindow();

    return std::vector<std::vector<Int_t>>{};
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