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
    Double_t t0_pos_z  = -1100.0;

    // ADCのbin設定
    const Int_t adc_bin_num = 4096;
    const Int_t adjust_adc_bin_num = 4096; // for adjustment
    const Double_t adc_min = 0.0;
    const Double_t adc_max = 4096.0;

    // TDCのbin設定
    Int_t tdc_bin_num = 65536;
    Int_t adjust_tdc_bin_num = 65536; // for adjustment
    Double_t tdc_min = 0.0;
    Double_t tdc_max = 2097152.0;

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

    // void beam_initialize() {
    //     beam_generator = -1;
    // }
    
private:
    Config() = default; // コンストラクタをプライベートにして外部からのインスタンス生成を禁止
    Config(const Config&) = delete;
    Config& operator=(const Config&) = delete;
};

#endif
