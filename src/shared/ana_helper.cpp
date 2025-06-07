#include "ana_helper.h"

namespace ana_helper {

    // ____________________________________________________________________________________________
    TCanvas* add_tab(TGTab *tab, const char* tabName) {
        // タブを作成し、キャンバスを埋め込む
        TGCompositeFrame *tf = tab->AddTab(tabName);
        TRootEmbeddedCanvas *embeddedCanvas = new TRootEmbeddedCanvas(tabName, tf, 1000, 800);
       tf->AddFrame(embeddedCanvas, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
        return embeddedCanvas->GetCanvas();
    }

    // ____________________________________________________________________________________________
    std::pair<Double_t, Double_t> cal_one_photon_gain(std::pair<Double_t, Double_t> mean, std::pair<Double_t, Double_t> pedestal, std::pair<Double_t, Double_t> n_pedestal, Double_t n_total) {
        Double_t mu = -TMath::Log( n_pedestal.first / n_total );
        Double_t gain = ( mean.first - pedestal.first ) / mu;
        
        // -- cal error propagation -----
        Double_t pdv_mean = 1.0/mu;
        Double_t pdv_ped = -1.0/mu;
        Double_t pdv_n_ped = gain / (n_pedestal.first * mu);
        Double_t error = TMath::Sqrt( TMath::Power(pdv_mean*mean.second, 2.0) + TMath::Power(pdv_ped*pedestal.second, 2.0) + TMath::Power(pdv_n_ped*n_pedestal.second, 2.0) );
        
        std::pair<Double_t, Double_t> one_photon_gain{gain, error};
        return one_photon_gain;
    }

}