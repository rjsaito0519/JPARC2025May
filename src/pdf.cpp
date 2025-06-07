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
static Int_t nth_pad = 1;
static Int_t rows = 2, cols = 2;
static Int_t max_pads = rows * cols;
static TString pdf_path = Form("%s/img/LED.pdf", OUTPUT_DIR.Data());
    


void adjust_range(TH1D *h)
{
    Double_t mean  = h->GetMean();        
    Double_t stdev = h->GetStdDev();        
    Double_t lower = std::max(mean - 3.0 * stdev, conf.adc_min);
    Double_t upper = std::min(mean + 7.0 * stdev, conf.adc_max);
    h->GetXaxis()->SetRangeUser(lower, upper);
}

void analyze(Int_t run_num, TCanvas* c){
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

        // -- KVC -----
        for (Int_t i = 0, n = (*kvc2_raw_seg).size(); i < n; i++) {
            Int_t index = static_cast<Int_t>((*kvc2_raw_seg)[i]);
            if (0 <= index && index < conf.num_of_ch.at("kvc2")) {
                h_kvc2_adc[index][0]->Fill((*kvc2_adc_a)[index]);
                h_kvc2_adc[index][1]->Fill((*kvc2_adc_b)[index]);
                h_kvc2_adc[index][2]->Fill((*kvc2_adc_c)[index]);
                h_kvc2_adc[index][3]->Fill((*kvc2_adc_d)[index]);
            }
        }
    }

    // +------+
    // | plot |
    // +------+

    // -- BAC -----
    for (Int_t ch = 0; ch < conf.num_of_ch.at("bac"); ch++ ) {
        if (nth_pad > max_pads) {
            c->Print(pdf_path);
            c->Clear();
            c->Divide(cols, rows);
            nth_pad = 1;
        }
        c->cd(nth_pad);
        adjust_range(h_bac_adc[ch]);
        h_bac_adc[ch]->Draw();
        nth_pad++;
    }

    // -- KVC -----
    for (Int_t ch = 0; ch < conf.num_of_ch.at("kvc2"); ch++ ) {
        for (Int_t UorD = 0; UorD < conf.num_of_UorD.at("kvc2"); UorD++ ) {
            if (nth_pad > max_pads) {
                c->Print(pdf_path);
                c->Clear();
                c->Divide(cols, rows);
                nth_pad = 1;
            }
            c->cd(nth_pad);
            adjust_range(h_kvc2_adc[ch][UorD]);
            h_kvc2_adc[ch][UorD]->Draw();
            nth_pad++;
        }
    }
    c->Print(pdf_path);

    delete f;
    
}

Int_t main(int argc, char** argv) {

    // -- prepare pdf -----
    
    auto c = new TCanvas("", "", 1500, 1200);
    c->Divide(cols, rows);
    c->Print(pdf_path + "["); // start

    std::vector<Int_t> runNumbers = {
        626, 627, 628, 629, 630, 631, 632, 633, 634, 635, 636, 637, 642, 643, 644, 645, 646, 647, 648,
        650, 651, 652, 653, 655, 656, 657, 658, 659, 661, 662, 663, 664, 665, 667, 668, 669, 671, 672, 673,
        674, 678, 679, 680, 681, 683, 684, 685, 686, 687, 689, 690, 691, 692, 696, 697, 698, 699, 702, 703,
        704, 705, 706, 708, 709, 710, 711, 716, 717, 718, 719, 720, 722, 723, 724, 725, 726, 728, 729, 730,
        731, 732, 734, 735, 736, 737, 738, 740, 741, 742, 743, 744, 746, 747, 748, 749, 751, 752, 754, 755,
        757, 758, 759, 760, 762, 763, 765, 766, 768, 769, 770, 772, 773, 774, 775
    };

    for (Int_t run_num : runNumbers) {
        analyze(run_num, c);
    }

    c->Print(pdf_path);
    c->Print(pdf_path + "]"); // end

    return 0;
}