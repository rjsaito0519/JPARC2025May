// config.h
#ifndef CONFIG_H
#define CONFIG_H

class Config {
public:
    static Config& getInstance() {
        static Config instance;
        return instance;
    }

    TString detector = "none";

    const Int_t max_tdc_hit = 16;

    Double_t btof_threshold = -3.0;

    Double_t kvc2_pos_z = -672.1167;
    Double_t htof_pos_z = -493.363;
    Double_t bh2_pos_z = -554.872;
    Double_t t0_pos_z  = -1086.55;
    Double_t bac_pos_z  = -1041.63;

    Double_t detector_HV = 58.0;

    // ADCのbin設定
    const Int_t adc_bin_num = 4096;
    Int_t adjust_adc_bin_num = 4096; // for adjustment
    const Double_t adc_min = 0.0;
    const Double_t adc_max = 4096.0;

    // TDCのbin設定
    Int_t tdc_bin_num = 65536;
    Int_t adjust_tdc_bin_num = 65536; // for adjustment
    Double_t tdc_min = 0.0;
    Double_t tdc_max = 2097152.0;

    // NPEのbin設定
    Int_t npe_bin_num = 2100;
    Int_t adjust_npe_bin_num = 280;
    Double_t npe_min = -10.;
    Double_t npe_max = 410.;

    const std::unordered_map<std::string, Int_t> num_of_ch{
        { "bht", 63 },
        {  "t0",  5 },
        { "bac",  4 },
        { "bacs",  1 },
        { "sac",  8 },
        { "sacs",  1 },
        { "kvc1",  4 },
        { "kvc2",  4 },
        { "bh2",  11 },
        { "htof",  34 },
    };

    const std::unordered_map<std::string, Int_t> num_of_UorD{
        { "bht", 2 },
        {  "t0", 2 },
        { "bac", 1 },
        { "sac", 1 },
        { "kvc1", 2 },
        { "kvc2", 4 },
        // { "bh2",  1 },
    };


    std::unordered_map<std::string, std::pair<Double_t, Double_t>> tdc_search_range{
        { "bht", {680000.0, 700000.0} },
        {  "t0", {680000.0, 700000.0} },
        { "bacs", {720000.0, 750000.0} },
        { "sacs", {670000.0, 710000.0} },
        { "kvc2", {710000.0, 730000.0} },
        { "bh2", {700000.0, 720000.0} },
        
        // { "bac",  4 },
        // { "SAC",  8 },
        // { "KVC",  4 },
        // { "BH2",  1 },
    };

    const std::vector<std::vector<Double_t>> t0_ped_mip_distance{
        {77.0, 52.0, 75.0, 50.0, 71.0}, // UP
        {55.0, 35.0, 55.0, 53.0, 45.0}  // DOWN
    };

    // KVC one photon gain
    std::unordered_map<Int_t, std::vector<std::pair<Double_t, Double_t>>> kvc2_opg{
        { 56, std::vector<std::pair<Double_t, Double_t>>(16, {TMath::QuietNaN(), TMath::QuietNaN()}) },

        { 57, {
            {7.605, 0.019},  // board0 seg0
            {7.311, 0.023},  // board0 seg1
            {7.195, 0.022},  // board0 seg2
            {TMath::QuietNaN(), TMath::QuietNaN()}, // board0 seg3

            {7.019, 0.024},  // board1 seg0
            {7.158, 0.024},  // board1 seg1
            {5.585, 0.028},  // board1 seg2
            {TMath::QuietNaN(), TMath::QuietNaN()}, // board1 seg3

            {5.582, 0.068},  // board2 seg0
            {7.150, 0.021},  // board2 seg1
            {7.372, 0.023},  // board2 seg2
            {TMath::QuietNaN(), TMath::QuietNaN()}, // board2 seg3

            {7.159, 0.027},  // board3 seg0
            {6.738, 0.020},  // board3 seg1
            {6.477, 0.022},  // board3 seg2
            {TMath::QuietNaN(), TMath::QuietNaN()}  // board3 seg3
        }},

        { 58, {
            {9.392, 0.023},  // board0 seg0
            {9.217, 0.029},  // board0 seg1
            {8.829, 0.027},  // board0 seg2
            {8.631, 0.029},  // board0 seg3

            {8.413, 0.029},  // board1 seg0
            {8.844, 0.029},  // board1 seg1
            {6.756, 0.034},  // board1 seg2
            {8.536, 0.030},  // board1 seg3

            {6.456, 0.079},  // board2 seg0
            {8.928, 0.027},  // board2 seg1
            {8.974, 0.028},  // board2 seg2
            {8.536, 0.035},  // board2 seg3

            {8.629, 0.032},  // board3 seg0
            {8.618, 0.025},  // board3 seg1
            {7.834, 0.026},  // board3 seg2
            {8.101, 0.036}   // board3 seg3
        }}
    };

    std::unordered_map<Int_t, std::vector<Double_t>> kvc2_pedestal{
        { 58, {
            117.852,   // board0 seg0
            98.5445,  // board1 seg0
            51.5,     // board2 seg0
            102.606,   // board3 seg0

            113.605,   // board0 seg1
            96.0246,  // board1 seg1
            120.671,   // board2 seg1
            100.695,   // board3 seg1

            113.523,   // board0 seg2
            116.772,   // board1 seg2
            121.768,   // board2 seg2
            122.725,   // board3 seg2

            102.681,   // board0 seg3
            103.2,     // board1 seg3
            102.065,   // board2 seg3
            111.553    // board3 seg3
        }}
    };


    void kvc2_initialize() {
        adjust_adc_bin_num = 1024;
    }
    
private:
    Config() = default; // コンストラクタをプライベートにして外部からのインスタンス生成を禁止
    Config(const Config&) = delete;
    Config& operator=(const Config&) = delete;
};

#endif
