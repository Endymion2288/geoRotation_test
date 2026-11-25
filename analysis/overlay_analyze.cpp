// Compiled C++17 16-thread OpenMP paired 3-way analyzer for Stage 08 Two-Track Overlay study.
// Supports --mode baseline-only (2-way EOS Baseline-May30 vs Baseline-Jul8) and
// --mode full (3-way adding Latest-GeoRotation) across all 322 matched bin-pairs (6,339,599 events).

#include <Math/QuantFuncMathCore.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TROOT.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <omp.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

static constexpr double CL_1SIGMA = 0.6826894921370859; // erf(1/sqrt(2))
static constexpr double ALPHA_1SIGMA = 1.0 - CL_1SIGMA;

enum RawDatasetIdx { DS_008749 = 0, DS_010698 = 1, DS_CROSS = 2, DS_10004X = 3, N_RAW_DS = 4 };
enum VersionIdx { VER_MAY30 = 0, VER_JUL8 = 1, VER_LATEST = 2, N_MAX_VER = 3 };
enum PairIdx { PAIR_JUL8_MAY30 = 0, PAIR_LATEST_JUL8 = 1, PAIR_LATEST_MAY30 = 2, N_MAX_PAIRS = 3 };
enum MetricIdx { MET_PAWAN_C2R = 0, MET_LT2 = 1, MET_GE2 = 2, MET_MATCHED = 3, N_METRICS = 4 };
enum RegimeIdx { REG_GLOBAL = 0, REG_SMALL_DR1 = 1, REG_PLATEAU = 2, REG_LARGE_DX0 = 3, N_REGIMES = 4 };

struct CurveSpec {
    const char* key;
    const char* xlabel;
    double xmin;
    double xmax;
    int nbins;
};

static constexpr int N_CURVES = 14;
static constexpr CurveSpec CURVE_SPECS[N_CURVES] = {
    {"dR1_wide",     "Track Separation (\\Delta R_1) [mm]",     0.0, 100.0, 200},
    {"dR1_zoom",     "Track Separation (\\Delta R_1) [mm]",     0.0,  20.0,  40},
    {"dR0_wide",     "Track Separation (\\Delta R_0) [mm]",     0.0, 100.0, 200},
    {"dR0_zoom",     "Track Separation (\\Delta R_0) [mm]",     0.0,  20.0,  40},
    {"abs_dY0_wide", "Track Separation (|\\Delta Y_0|) [mm]",   0.0, 100.0, 100},
    {"abs_dY0_zoom", "Track Separation (|\\Delta Y_0|) [mm]",   0.0,  20.0,  40},
    {"abs_dY1_wide", "Track Separation (|\\Delta Y_1|) [mm]",   0.0, 100.0, 100},
    {"abs_dY1_zoom", "Track Separation (|\\Delta Y_1|) [mm]",   0.0,  20.0,  40},
    {"abs_dX0_wide", "Track Separation (|\\Delta X_0|) [mm]",   0.0, 100.0, 100},
    {"dX0_signed",   "Track Separation (\\Delta X_0) [mm]",  -100.0, 100.0, 200},
    {"abs_dX1_wide", "Track Separation (|\\Delta X_1|) [mm]",   0.0, 100.0, 100},
    {"dX1_signed",   "Track Separation (\\Delta X_1) [mm]",  -100.0, 100.0, 200},
    {"dY1_signed",   "Track Separation (\\Delta Y_1) [mm]",  -100.0, 100.0, 200},
    {"dY0_signed",   "Track Separation (\\Delta Y_0) [mm]",  -100.0, 100.0, 200},
};

static constexpr int N_RW_CURVES = 5;
static constexpr int RW_CURVE_IDS[N_RW_CURVES] = {0, 1, 2, 8, 4}; // dR1_wide, dR1_zoom, dR0_wide, abs_dX0_wide, abs_dY0_wide
static constexpr int N_MULT_CURVES = 3;
static constexpr int MULT_CURVE_IDS[N_MULT_CURVES] = {0, 1, 8}; // dR1_wide, dR1_zoom, abs_dX0_wide

static constexpr int N_PZ_BINS = 100;
static constexpr double PZ_LOG_MIN = 0.0;                  // 1 MeV
static const double PZ_LOG_MAX = std::log10(2.0e7);        // 2e7 MeV

static constexpr int N_RES_BINS = 100;
static constexpr double RES_MIN = -0.5;
static constexpr double RES_MAX = 0.5;

static constexpr int N_RW_P0_BINS = 40;
static const double RW_P0_LOG_MIN = std::log10(50.0);      // 50 GeV
static const double RW_P0_LOG_MAX = std::log10(5000.0);    // 5000 GeV

static constexpr int N_PLAT_P0_BINS = 20;
static constexpr int N_PLAT_R0_BINS = 19;                  // [0, 95] mm, step 5
static constexpr int N_PLAT_SEG_BINS = 13;                 // 0..12
static constexpr int N_DX_STRAT_BINS = 20;                 // [0, 100] mm, step 5
static constexpr int N_MAP_BINS = 40;                      // [-100, 100] mm, step 5

inline void clopperPearson1Sigma(long long k, long long n, double& eff, double& err_lo, double& err_hi) {
    if (n <= 0) {
        eff = std::numeric_limits<double>::quiet_NaN();
        err_lo = 0.0;
        err_hi = 0.0;
        return;
    }
    eff = static_cast<double>(k) / static_cast<double>(n);
    double lo = (k == 0) ? 0.0 : ROOT::Math::beta_quantile(0.5 * ALPHA_1SIGMA, static_cast<double>(k), static_cast<double>(n - k + 1));
    double hi = (k == n) ? 1.0 : ROOT::Math::beta_quantile(1.0 - 0.5 * ALPHA_1SIGMA, static_cast<double>(k + 1), static_cast<double>(n - k));
    err_lo = std::max(0.0, eff - lo);
    err_hi = std::max(0.0, hi - eff);
}

inline void mcnemarDelta(long long bp, long long g, long long l, long long bf, double& delta, double& err) {
    long long n = bp + g + l + bf;
    if (n <= 0) {
        delta = std::numeric_limits<double>::quiet_NaN();
        err = std::numeric_limits<double>::quiet_NaN();
        return;
    }
    double dn = static_cast<double>(n);
    double diff = static_cast<double>(g - l);
    delta = diff / dn;
    double var_num = std::max(0.0, static_cast<double>(g + l) - (diff * diff) / dn);
    err = std::sqrt(var_num) / dn;
}

inline std::string fmtDoubleOrNull(double v) {
    if (!std::isfinite(v)) return "null";
    std::ostringstream oss;
    oss << std::setprecision(10) << v;
    return oss.str();
}

struct EffBinSet {
    int nbins = 0;
    double xmin = 0.0;
    double xmax = 1.0;
    bool log_x = false;
    long long underflow_den = 0;
    long long overflow_den = 0;
    long long invalid_den = 0;
    std::vector<long long> den;
    std::array<std::vector<long long>, N_METRICS> pass;

    EffBinSet() = default;
    EffBinSet(int nb, double lo, double hi, bool is_log = false)
        : nbins(nb), xmin(lo), xmax(hi), log_x(is_log), den(nb, 0) {
        for (int m = 0; m < N_METRICS; ++m) pass[m].assign(nb, 0);
    }

    int findBin(double x) {
        if (!std::isfinite(x)) { invalid_den++; return -1; }
        double v = log_x ? (x > 0.0 ? std::log10(x) : -1e30) : x;
        double lo = log_x ? std::log10(xmin) : xmin;
        double hi = log_x ? std::log10(xmax) : xmax;
        if (v < lo) { underflow_den++; return -1; }
        if (v >= hi) { overflow_den++; return -1; }
        int b = static_cast<int>((v - lo) / (hi - lo) * nbins);
        if (b < 0) b = 0;
        if (b >= nbins) b = nbins - 1;
        return b;
    }

    void fill(double x, const bool p[N_METRICS]) {
        int b = findBin(x);
        if (b < 0) return;
        den[b]++;
        for (int m = 0; m < N_METRICS; ++m) {
            if (p[m]) pass[m][b]++;
        }
    }

    void fillSingleMetric(double x, bool p0) {
        int b = findBin(x);
        if (b < 0) return;
        den[b]++;
        if (p0) pass[0][b]++;
    }

    void merge(const EffBinSet& o) {
        underflow_den += o.underflow_den;
        overflow_den += o.overflow_den;
        invalid_den += o.invalid_den;
        for (int i = 0; i < nbins; ++i) {
            den[i] += o.den[i];
            for (int m = 0; m < N_METRICS; ++m) pass[m][i] += o.pass[m][i];
        }
    }
};

struct MultBinSet {
    int nbins = 0;
    double xmin = 0.0;
    double xmax = 1.0;
    std::vector<long long> den, c_lt0, c_lt1, c_lt2_pass, c_lt2_fail, c_ge3;

    MultBinSet() = default;
    MultBinSet(int nb, double lo, double hi)
        : nbins(nb), xmin(lo), xmax(hi),
          den(nb, 0), c_lt0(nb, 0), c_lt1(nb, 0), c_lt2_pass(nb, 0), c_lt2_fail(nb, 0), c_ge3(nb, 0) {}

    void fill(double x, int lt, bool pass_c2r) {
        if (!std::isfinite(x) || x < xmin || x >= xmax) return;
        int b = static_cast<int>((x - xmin) / (xmax - xmin) * nbins);
        if (b < 0) b = 0;
        if (b >= nbins) b = nbins - 1;
        den[b]++;
        if (lt <= 0) c_lt0[b]++;
        else if (lt == 1) c_lt1[b]++;
        else if (lt == 2) {
            if (pass_c2r) c_lt2_pass[b]++;
            else c_lt2_fail[b]++;
        } else {
            c_ge3[b]++;
        }
    }

    void merge(const MultBinSet& o) {
        for (int i = 0; i < nbins; ++i) {
            den[i] += o.den[i];
            c_lt0[i] += o.c_lt0[i];
            c_lt1[i] += o.c_lt1[i];
            c_lt2_pass[i] += o.c_lt2_pass[i];
            c_lt2_fail[i] += o.c_lt2_fail[i];
            c_ge3[i] += o.c_ge3[i];
        }
    }
};

struct PairedCell {
    long long bp = 0, g = 0, l = 0, bf = 0;
    void add(bool a, bool b) {
        if (a && b) bp++;
        else if (!a && b) g++;
        else if (a && !b) l++;
        else bf++;
    }
    void merge(const PairedCell& o) {
        bp += o.bp; g += o.g; l += o.l; bf += o.bf;
    }
};

struct PairedBinSet {
    int nbins = 0;
    double xmin = 0.0;
    double xmax = 1.0;
    // Track MET_PAWAN_C2R (0) and MET_LT2 (1)
    std::array<std::vector<PairedCell>, 2> cells;

    PairedBinSet() = default;
    PairedBinSet(int nb, double lo, double hi) : nbins(nb), xmin(lo), xmax(hi) {
        cells[0].resize(nb);
        cells[1].resize(nb);
    }

    void fill(double x, bool a0, bool b0, bool a1, bool b1) {
        if (!std::isfinite(x) || x < xmin || x >= xmax) return;
        int b = static_cast<int>((x - xmin) / (xmax - xmin) * nbins);
        if (b < 0) b = 0;
        if (b >= nbins) b = nbins - 1;
        cells[0][b].add(a0, b0);
        cells[1][b].add(a1, b1);
    }

    void merge(const PairedBinSet& o) {
        for (int m = 0; m < 2; ++m)
            for (int i = 0; i < nbins; ++i) cells[m][i].merge(o.cells[m][i]);
    }
};

struct StageLossCounts {
    long long n_den = 0;
    long long pass_c2R = 0;
    long long lt2_fail_c2R = 0;
    long long lt_ge3 = 0;
    long long lt0 = 0;
    long long lt1_total = 0;
    long long lt1_miss_seg_st1 = 0;
    long long lt1_miss_seg_st2 = 0;
    long long lt1_miss_seg_st3 = 0;
    long long lt1_miss_seg_multi = 0;
    long long lt1_killed_by_gb = 0;
    long long lt1_all3gb_lost_in_ckf = 0;
    long long lt1_unmatched_kept = 0;

    void merge(const StageLossCounts& o) {
        n_den += o.n_den;
        pass_c2R += o.pass_c2R;
        lt2_fail_c2R += o.lt2_fail_c2R;
        lt_ge3 += o.lt_ge3;
        lt0 += o.lt0;
        lt1_total += o.lt1_total;
        lt1_miss_seg_st1 += o.lt1_miss_seg_st1;
        lt1_miss_seg_st2 += o.lt1_miss_seg_st2;
        lt1_miss_seg_st3 += o.lt1_miss_seg_st3;
        lt1_miss_seg_multi += o.lt1_miss_seg_multi;
        lt1_killed_by_gb += o.lt1_killed_by_gb;
        lt1_all3gb_lost_in_ckf += o.lt1_all3gb_lost_in_ckf;
        lt1_unmatched_kept += o.lt1_unmatched_kept;
    }
};

struct Map2D {
    std::array<std::array<long long, N_MAP_BINS>, N_MAP_BINS> den{};
    std::array<std::array<long long, N_MAP_BINS>, N_MAP_BINS> lost{};

    void fill(double x, double y, bool is_lost) {
        if (!std::isfinite(x) || !std::isfinite(y)) return;
        if (x < -100.0 || x >= 100.0 || y < -100.0 || y >= 100.0) return;
        int ix = static_cast<int>((x + 100.0) / 200.0 * N_MAP_BINS);
        int iy = static_cast<int>((y + 100.0) / 200.0 * N_MAP_BINS);
        if (ix < 0 || ix >= N_MAP_BINS || iy < 0 || iy >= N_MAP_BINS) return;
        den[iy][ix]++;
        if (is_lost) lost[iy][ix]++;
    }

    void merge(const Map2D& o) {
        for (int iy = 0; iy < N_MAP_BINS; ++iy)
            for (int ix = 0; ix < N_MAP_BINS; ++ix) {
                den[iy][ix] += o.den[iy][ix];
                lost[iy][ix] += o.lost[iy][ix];
            }
    }
};

struct VersionDatasetStats {
    int n_files = 0;
    long long n_events = 0;
    long long n_og2 = 0;
    std::array<long long, 6> lt_counts{}; // 0, 1, 2, 3, 4, 5+
    long long lt2_fail_c2R = 0;
    std::array<std::array<long long, N_METRICS>, N_REGIMES> reg_pass{};
    std::array<long long, N_REGIMES> reg_den{};

    std::array<EffBinSet, N_CURVES> curves;
    std::array<MultBinSet, N_MULT_CURVES> mult_curves;

    // Momentum
    std::array<long long, N_PZ_BINS> pz_reco{};
    long long pz_reco_tot = 0, pz_reco_uf = 0, pz_reco_of = 0;
    std::array<long long, N_PZ_BINS> pz_og{};
    long long pz_og_tot = 0, pz_og_uf = 0, pz_og_of = 0;
    std::array<long long, N_RES_BINS> pz_res{};
    long long pz_res_tot = 0, pz_res_uf = 0, pz_res_of = 0;
    double pz_res_sum = 0.0, pz_res_sum2 = 0.0;
    std::array<long long, N_RW_P0_BINS> og_p0_40{};
    long long og_p0_40_tot = 0;

    // Stage-of-loss
    std::array<StageLossCounts, N_REGIMES> stage_reg;
    EffBinSet plat_nSeg{N_PLAT_SEG_BINS, 0.0, 13.0, false};
    EffBinSet plat_nGBSeg{N_PLAT_SEG_BINS, 0.0, 13.0, false};
    EffBinSet plat_min_p0{N_PLAT_P0_BINS, 50.0, 5000.0, true};
    EffBinSet plat_max_R0{N_PLAT_R0_BINS, 0.0, 95.0, false};
    std::array<StageLossCounts, N_DX_STRAT_BINS> abs_dX0_stage{};

    // Large-dX stratified + 2D loss maps
    std::array<EffBinSet, 3> eff_dx0_by_x0;
    std::array<EffBinSet, 3> eff_dx0_by_R0;
    Map2D loss_map_st1;
    Map2D loss_map_st3;

    VersionDatasetStats() {
        for (int c = 0; c < N_CURVES; ++c) {
            curves[c] = EffBinSet(CURVE_SPECS[c].nbins, CURVE_SPECS[c].xmin, CURVE_SPECS[c].xmax, false);
        }
        for (int m = 0; m < N_MULT_CURVES; ++m) {
            int cid = MULT_CURVE_IDS[m];
            mult_curves[m] = MultBinSet(CURVE_SPECS[cid].nbins, CURVE_SPECS[cid].xmin, CURVE_SPECS[cid].xmax);
        }
        for (int s = 0; s < 3; ++s) {
            eff_dx0_by_x0[s] = EffBinSet(N_DX_STRAT_BINS, 0.0, 100.0, false);
            eff_dx0_by_R0[s] = EffBinSet(N_DX_STRAT_BINS, 0.0, 100.0, false);
        }
    }

    void merge(const VersionDatasetStats& o) {
        n_files += o.n_files;
        n_events += o.n_events;
        n_og2 += o.n_og2;
        for (int i = 0; i < 6; ++i) lt_counts[i] += o.lt_counts[i];
        lt2_fail_c2R += o.lt2_fail_c2R;
        for (int r = 0; r < N_REGIMES; ++r) {
            reg_den[r] += o.reg_den[r];
            for (int m = 0; m < N_METRICS; ++m) reg_pass[r][m] += o.reg_pass[r][m];
            stage_reg[r].merge(o.stage_reg[r]);
        }
        for (int c = 0; c < N_CURVES; ++c) curves[c].merge(o.curves[c]);
        for (int m = 0; m < N_MULT_CURVES; ++m) mult_curves[m].merge(o.mult_curves[m]);
        for (int i = 0; i < N_PZ_BINS; ++i) { pz_reco[i] += o.pz_reco[i]; pz_og[i] += o.pz_og[i]; }
        pz_reco_tot += o.pz_reco_tot; pz_reco_uf += o.pz_reco_uf; pz_reco_of += o.pz_reco_of;
        pz_og_tot += o.pz_og_tot; pz_og_uf += o.pz_og_uf; pz_og_of += o.pz_og_of;
        for (int i = 0; i < N_RES_BINS; ++i) pz_res[i] += o.pz_res[i];
        pz_res_tot += o.pz_res_tot; pz_res_uf += o.pz_res_uf; pz_res_of += o.pz_res_of;
        pz_res_sum += o.pz_res_sum; pz_res_sum2 += o.pz_res_sum2;
        for (int i = 0; i < N_RW_P0_BINS; ++i) og_p0_40[i] += o.og_p0_40[i];
        og_p0_40_tot += o.og_p0_40_tot;
        plat_nSeg.merge(o.plat_nSeg);
        plat_nGBSeg.merge(o.plat_nGBSeg);
        plat_min_p0.merge(o.plat_min_p0);
        plat_max_R0.merge(o.plat_max_R0);
        for (int i = 0; i < N_DX_STRAT_BINS; ++i) abs_dX0_stage[i].merge(o.abs_dX0_stage[i]);
        for (int s = 0; s < 3; ++s) {
            eff_dx0_by_x0[s].merge(o.eff_dx0_by_x0[s]);
            eff_dx0_by_R0[s].merge(o.eff_dx0_by_R0[s]);
        }
        loss_map_st1.merge(o.loss_map_st1);
        loss_map_st3.merge(o.loss_map_st3);
    }
};

struct PairedDatasetStats {
    std::array<std::array<PairedCell, N_METRICS>, N_REGIMES> reg_cells{};
    std::array<std::array<std::array<long long, 5>, 5>, N_REGIMES> trans{};
    std::array<PairedBinSet, N_CURVES> curves;

    PairedDatasetStats() {
        for (int c = 0; c < N_CURVES; ++c) {
            curves[c] = PairedBinSet(CURVE_SPECS[c].nbins, CURVE_SPECS[c].xmin, CURVE_SPECS[c].xmax);
        }
    }

    void merge(const PairedDatasetStats& o) {
        for (int r = 0; r < N_REGIMES; ++r) {
            for (int m = 0; m < N_METRICS; ++m) reg_cells[r][m].merge(o.reg_cells[r][m]);
            for (int a = 0; a < 5; ++a)
                for (int b = 0; b < 5; ++b) trans[r][a][b] += o.trans[r][a][b];
        }
        for (int c = 0; c < N_CURVES; ++c) curves[c].merge(o.curves[c]);
    }
};

struct WeightedBinSet {
    int nbins = 0;
    double xmin = 0.0;
    double xmax = 1.0;
    std::vector<double> sumw_den, sumw2_den, sumw_pass, sumw2_pass, sumw2_fail;

    WeightedBinSet() = default;
    WeightedBinSet(int nb, double lo, double hi)
        : nbins(nb), xmin(lo), xmax(hi),
          sumw_den(nb, 0.0), sumw2_den(nb, 0.0), sumw_pass(nb, 0.0), sumw2_pass(nb, 0.0), sumw2_fail(nb, 0.0) {}

    void fill(double x, double w, bool pass) {
        if (!std::isfinite(x) || x < xmin || x >= xmax) return;
        int b = static_cast<int>((x - xmin) / (xmax - xmin) * nbins);
        if (b < 0) b = 0;
        if (b >= nbins) b = nbins - 1;
        double w2 = w * w;
        sumw_den[b] += w;
        sumw2_den[b] += w2;
        if (pass) {
            sumw_pass[b] += w;
            sumw2_pass[b] += w2;
        } else {
            sumw2_fail[b] += w2;
        }
    }

    void merge(const WeightedBinSet& o) {
        for (int i = 0; i < nbins; ++i) {
            sumw_den[i] += o.sumw_den[i];
            sumw2_den[i] += o.sumw2_den[i];
            sumw_pass[i] += o.sumw_pass[i];
            sumw2_pass[i] += o.sumw2_pass[i];
            sumw2_fail[i] += o.sumw2_fail[i];
        }
    }
};

struct ReweightStats {
    long long n_events = 0;
    double sumw = 0.0;
    double sumw2 = 0.0;
    std::array<double, N_RW_P0_BINS> closure_sumw{};
    double closure_tot = 0.0;
    std::array<WeightedBinSet, N_RW_CURVES> curves;

    ReweightStats() {
        for (int i = 0; i < N_RW_CURVES; ++i) {
            int cid = RW_CURVE_IDS[i];
            curves[i] = WeightedBinSet(CURVE_SPECS[cid].nbins, CURVE_SPECS[cid].xmin, CURVE_SPECS[cid].xmax);
        }
    }

    void merge(const ReweightStats& o) {
        n_events += o.n_events;
        sumw += o.sumw;
        sumw2 += o.sumw2;
        closure_tot += o.closure_tot;
        for (int i = 0; i < N_RW_P0_BINS; ++i) closure_sumw[i] += o.closure_sumw[i];
        for (int i = 0; i < N_RW_CURVES; ++i) curves[i].merge(o.curves[i]);
    }
};

struct FileTriple {
    int ds_idx = 0;
    std::string ds_name;
    int pos_bin = 0;
    int neg_bin = 0;
    std::string path_may30;
    std::string path_jul8;
    std::string path_latest;
};

struct EventEval {
    int lt = 0;
    bool pass[N_METRICS] = {false, false, false, false};
    bool lt2_fail_c2r = false;
    int lt1_kept_idx = -1; // 0 or 1 if matched, -1 if unmatched
    int stage_code = 0;    // 0=pass_c2R, 1=lt2_fail, 2=lt>=3, 3=lt0, 4..10=lt1 subtypes
    int nSeg = 0;
    int nGBSeg = 0;
    double matched_res[2] = {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
};

inline void hermitePredict(double sz,
                           double z0, double x0, double y0, double tx0, double ty0,
                           double z1, double x1, double y1, double tx1, double ty1,
                           double& xp, double& yp) {
    double L = z1 - z0;
    if (std::abs(L) < 1e-6) { xp = x0; yp = y0; return; }
    double t = (sz - z0) / L;
    double t2 = t * t;
    double t3 = t2 * t;
    double h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
    double h10 = t3 - 2.0 * t2 + t;
    double h01 = -2.0 * t3 + 3.0 * t2;
    double h11 = t3 - t2;
    xp = h00 * x0 + h10 * L * tx0 + h01 * x1 + h11 * L * tx1;
    yp = h00 * y0 + h10 * L * ty0 + h01 * y1 + h11 * L * ty1;
}

inline int stationOfZ(double z) {
    if (z >= -100.0 && z <= 350.0) return 0;   // St. 1
    if (z >= 950.0 && z <= 1550.0) return 1;   // St. 2
    if (z >= 2150.0 && z <= 2750.0) return 2;  // St. 3
    return -1; // IFT or outside
}

inline EventEval evaluateEvent(
    int lt,
    const std::vector<double>& c2,
    const std::vector<double>& nd,
    const std::vector<double>& x0,
    const std::vector<double>& y0,
    const std::vector<double>& pz0,
    const std::vector<double>& x1,
    const std::vector<double>& y1,
    const std::vector<double>& og_x0,
    const std::vector<double>& og_y0,
    const std::vector<double>& og_z0,
    const std::vector<double>& og_px0,
    const std::vector<double>& og_py0,
    const std::vector<double>& og_pz0,
    const std::vector<double>& og_x1,
    const std::vector<double>& og_y1,
    const std::vector<double>& og_z1,
    const std::vector<double>& og_px1,
    const std::vector<double>& og_py1,
    const std::vector<double>& og_pz1,
    unsigned int nSeg_u,
    const std::vector<double>& seg_x,
    const std::vector<double>& seg_y,
    const std::vector<double>& seg_z,
    unsigned int nGBSeg_u,
    const std::vector<double>& gb_x,
    const std::vector<double>& gb_y,
    const std::vector<double>& gb_z) {

    EventEval ev;
    ev.lt = lt;
    int n3_seg = 0;
    for (double z : seg_z) if (stationOfZ(z) >= 0) n3_seg++;
    int n3_gb = 0;
    for (double z : gb_z) if (stationOfZ(z) >= 0) n3_gb++;
    ev.nSeg = n3_seg;
    ev.nGBSeg = n3_gb;
    int ntr = static_cast<int>(x0.size());

    auto trackPassC2R = [&](int j) -> bool {
        if (j < 0 || j >= ntr) return false;
        if (nd[j] <= 0.0 || (c2[j] / nd[j]) > 25.0) return false;
        if (std::hypot(x0[j], y0[j]) >= 95.0) return false;
        if (std::hypot(x1[j], y1[j]) >= 95.0) return false;
        return true;
    };

    ev.pass[MET_LT2] = (lt == 2);
    ev.pass[MET_GE2] = (lt >= 2);
    if (lt == 2 && ntr == 2 && trackPassC2R(0) && trackPassC2R(1)) {
        ev.pass[MET_PAWAN_C2R] = true;
    } else if (lt == 2) {
        ev.lt2_fail_c2r = true;
    }

    // 1-to-1 spatial matching (< 5 mm at St. 1 & St. 3) + c2R for eff_matched & momentum residual
    if (lt >= 2 && ntr >= 2) {
        double best_sum = 1e30;
        int best_a = -1, best_b = -1;
        for (int a = 0; a < ntr; ++a) {
            if (!trackPassC2R(a)) continue;
            double da0 = std::hypot(x0[a] - og_x0[0], y0[a] - og_y0[0]);
            double da1 = std::hypot(x1[a] - og_x1[0], y1[a] - og_y1[0]);
            if (da0 >= 5.0 || da1 >= 5.0) continue;
            for (int b = 0; b < ntr; ++b) {
                if (b == a || !trackPassC2R(b)) continue;
                double db0 = std::hypot(x0[b] - og_x0[1], y0[b] - og_y0[1]);
                double db1 = std::hypot(x1[b] - og_x1[1], y1[b] - og_y1[1]);
                if (db0 >= 5.0 || db1 >= 5.0) continue;
                double s = da0 + da1 + db0 + db1;
                if (s < best_sum) {
                    best_sum = s;
                    best_a = a;
                    best_b = b;
                }
            }
        }
        if (best_a >= 0 && best_b >= 0) {
            ev.pass[MET_MATCHED] = true;
            if (og_pz0[0] > 0.0) ev.matched_res[0] = (pz0[best_a] / og_pz0[0]) - 1.0;
            if (og_pz0[1] > 0.0) ev.matched_res[1] = (pz0[best_b] / og_pz0[1]) - 1.0;
        }
    }

    // Stage-of-loss classification
    if (ev.pass[MET_PAWAN_C2R]) {
        ev.stage_code = 0;
    } else if (lt == 2) {
        ev.stage_code = 1; // lt2_fail_c2R
    } else if (lt >= 3) {
        ev.stage_code = 2; // lt_ge3
    } else if (lt <= 0) {
        ev.stage_code = 3; // lt0
    } else {
        // lt == 1: identify which original track was kept and which was lost
        if (ntr >= 1) {
            double d0 = std::hypot(x0[0] - og_x0[0], y0[0] - og_y0[0]) +
                        std::hypot(x1[0] - og_x1[0], y1[0] - og_y1[0]);
            double d1 = std::hypot(x0[0] - og_x0[1], y0[0] - og_y0[1]) +
                        std::hypot(x1[0] - og_x1[1], y1[0] - og_y1[1]);
            if (std::min(d0, d1) < 15.0) {
                ev.lt1_kept_idx = (d0 <= d1) ? 0 : 1;
            }
        }
        if (ev.lt1_kept_idx < 0) {
            ev.stage_code = 10; // lt1_unmatched_kept
        } else {
            int lost = 1 - ev.lt1_kept_idx;
            double z0l = og_z0[lost], x0l = og_x0[lost], y0l = og_y0[lost];
            double z1l = og_z1[lost], x1l = og_x1[lost], y1l = og_y1[lost];
            double tx0 = (std::abs(og_pz0[lost]) > 1e-3) ? (og_px0[lost] / og_pz0[lost]) : 0.0;
            double ty0 = (std::abs(og_pz0[lost]) > 1e-3) ? (og_py0[lost] / og_pz0[lost]) : 0.0;
            double tx1 = (std::abs(og_pz1[lost]) > 1e-3) ? (og_px1[lost] / og_pz1[lost]) : 0.0;
            double ty1 = (std::abs(og_pz1[lost]) > 1e-3) ? (og_py1[lost] / og_pz1[lost]) : 0.0;

            bool seg_has[3] = {false, false, false};
            for (size_t k = 0; k < seg_z.size(); ++k) {
                int st = stationOfZ(seg_z[k]);
                if (st < 0 || seg_has[st]) continue;
                double xp = 0.0, yp = 0.0;
                hermitePredict(seg_z[k], z0l, x0l, y0l, tx0, ty0, z1l, x1l, y1l, tx1, ty1, xp, yp);
                if (std::hypot(seg_x[k] - xp, seg_y[k] - yp) < 5.0) {
                    seg_has[st] = true;
                }
            }
            int n_miss_seg = (!seg_has[0]) + (!seg_has[1]) + (!seg_has[2]);
            if (n_miss_seg >= 2) {
                ev.stage_code = 7; // lt1_miss_seg_multi
            } else if (!seg_has[0]) {
                ev.stage_code = 4; // lt1_miss_seg_st1
            } else if (!seg_has[1]) {
                ev.stage_code = 5; // lt1_miss_seg_st2
            } else if (!seg_has[2]) {
                ev.stage_code = 6; // lt1_miss_seg_st3
            } else {
                // All 3 stations had a SegmentFit segment; check GhostBustedTrackSegments
                bool gb_has[3] = {false, false, false};
                for (size_t k = 0; k < gb_z.size(); ++k) {
                    int st = stationOfZ(gb_z[k]);
                    if (st < 0 || gb_has[st]) continue;
                    double xp = 0.0, yp = 0.0;
                    hermitePredict(gb_z[k], z0l, x0l, y0l, tx0, ty0, z1l, x1l, y1l, tx1, ty1, xp, yp);
                    if (std::hypot(gb_x[k] - xp, gb_y[k] - yp) < 5.0) {
                        gb_has[st] = true;
                    }
                }
                if (gb_has[0] && gb_has[1] && gb_has[2]) {
                    ev.stage_code = 9; // lt1_all3gb_lost_in_ckf
                } else {
                    ev.stage_code = 8; // lt1_killed_by_gb
                }
            }
        }
    }
    return ev;
}

inline void addStageCount(StageLossCounts& s, int code) {
    s.n_den++;
    switch (code) {
        case 0: s.pass_c2R++; break;
        case 1: s.lt2_fail_c2R++; break;
        case 2: s.lt_ge3++; break;
        case 3: s.lt0++; break;
        case 4: s.lt1_total++; s.lt1_miss_seg_st1++; break;
        case 5: s.lt1_total++; s.lt1_miss_seg_st2++; break;
        case 6: s.lt1_total++; s.lt1_miss_seg_st3++; break;
        case 7: s.lt1_total++; s.lt1_miss_seg_multi++; break;
        case 8: s.lt1_total++; s.lt1_killed_by_gb++; break;
        case 9: s.lt1_total++; s.lt1_all3gb_lost_in_ckf++; break;
        default: s.lt1_total++; s.lt1_unmatched_kept++; break;
    }
}

void fillVersionStats(VersionDatasetStats& st,
                      const EventEval& ev,
                      const double seps[N_CURVES],
                      const bool in_reg[N_REGIMES],
                      double og_x0_0, double og_y0_0, double og_x1_0, double og_y1_0,
                      double og_x0_1, double og_y0_1, double og_x1_1, double og_y1_1,
                      double og_p0_0_gev, double og_p0_1_gev,
                      double og_pz0_0, double og_pz0_1,
                      const std::vector<double>& reco_pz0) {
    st.n_events++;
    st.n_og2++;
    int ltb = std::min(std::max(ev.lt, 0), 5);
    st.lt_counts[ltb]++;
    if (ev.lt2_fail_c2r) st.lt2_fail_c2R++;

    for (int r = 0; r < N_REGIMES; ++r) {
        if (!in_reg[r]) continue;
        st.reg_den[r]++;
        for (int m = 0; m < N_METRICS; ++m) {
            if (ev.pass[m]) st.reg_pass[r][m]++;
        }
        addStageCount(st.stage_reg[r], ev.stage_code);
    }

    for (int c = 0; c < N_CURVES; ++c) {
        st.curves[c].fill(seps[c], ev.pass);
    }
    for (int m = 0; m < N_MULT_CURVES; ++m) {
        int cid = MULT_CURVE_IDS[m];
        st.mult_curves[m].fill(seps[cid], ev.lt, ev.pass[MET_PAWAN_C2R]);
    }

    // Reconstructed Track_pz0
    for (double v : reco_pz0) {
        if (!std::isfinite(v) || v <= 1.0) { st.pz_reco_uf++; continue; }
        double lv = std::log10(v);
        if (lv < PZ_LOG_MIN) { st.pz_reco_uf++; continue; }
        if (lv >= PZ_LOG_MAX) { st.pz_reco_of++; continue; }
        int b = static_cast<int>((lv - PZ_LOG_MIN) / (PZ_LOG_MAX - PZ_LOG_MIN) * N_PZ_BINS);
        if (b < 0) b = 0;
        if (b >= N_PZ_BINS) b = N_PZ_BINS - 1;
        st.pz_reco[b]++;
        st.pz_reco_tot++;
    }
    // Original og_Track_pz0 & 40-bin og_p0
    double og_pzs[2] = {og_pz0_0, og_pz0_1};
    double og_p0s[2] = {og_p0_0_gev, og_p0_1_gev};
    for (int k = 0; k < 2; ++k) {
        double v = og_pzs[k];
        if (!std::isfinite(v) || v <= 1.0) { st.pz_og_uf++; }
        else {
            double lv = std::log10(v);
            if (lv < PZ_LOG_MIN) st.pz_og_uf++;
            else if (lv >= PZ_LOG_MAX) st.pz_og_of++;
            else {
                int b = static_cast<int>((lv - PZ_LOG_MIN) / (PZ_LOG_MAX - PZ_LOG_MIN) * N_PZ_BINS);
                if (b < 0) b = 0;
                if (b >= N_PZ_BINS) b = N_PZ_BINS - 1;
                st.pz_og[b]++;
                st.pz_og_tot++;
            }
        }
        if (og_p0s[k] >= 50.0 && og_p0s[k] <= 5000.0) {
            double lp = std::log10(og_p0s[k]);
            int b = static_cast<int>((lp - RW_P0_LOG_MIN) / (RW_P0_LOG_MAX - RW_P0_LOG_MIN) * N_RW_P0_BINS);
            if (b < 0) b = 0;
            if (b >= N_RW_P0_BINS) b = N_RW_P0_BINS - 1;
            st.og_p0_40[b]++;
            st.og_p0_40_tot++;
        }
    }
    // Matched relative residual
    for (int k = 0; k < 2; ++k) {
        double r = ev.matched_res[k];
        if (!std::isfinite(r)) continue;
        if (r < RES_MIN) { st.pz_res_uf++; continue; }
        if (r >= RES_MAX) { st.pz_res_of++; continue; }
        int b = static_cast<int>((r - RES_MIN) / (RES_MAX - RES_MIN) * N_RES_BINS);
        if (b < 0) b = 0;
        if (b >= N_RES_BINS) b = N_RES_BINS - 1;
        st.pz_res[b]++;
        st.pz_res_tot++;
        st.pz_res_sum += r;
        st.pz_res_sum2 += r * r;
    }

    // Plateau diagnostics
    double min_p0 = std::min(og_p0_0_gev, og_p0_1_gev);
    double max_R0 = std::max(std::hypot(og_x0_0, og_y0_0), std::hypot(og_x0_1, og_y0_1));
    if (in_reg[REG_PLATEAU]) {
        st.plat_nSeg.fillSingleMetric(std::min(ev.nSeg, 12) + 0.5, ev.pass[MET_PAWAN_C2R]);
        st.plat_nGBSeg.fillSingleMetric(std::min(ev.nGBSeg, 12) + 0.5, ev.pass[MET_PAWAN_C2R]);
        st.plat_min_p0.fillSingleMetric(min_p0, ev.pass[MET_PAWAN_C2R]);
        st.plat_max_R0.fillSingleMetric(max_R0, ev.pass[MET_PAWAN_C2R]);
    }

    // Large-|dX0| diagnostics
    double adx0 = seps[8];
    if (adx0 >= 0.0 && adx0 < 100.0) {
        int b = static_cast<int>(adx0 / 5.0);
        if (b >= 0 && b < N_DX_STRAT_BINS) addStageCount(st.abs_dX0_stage[b], ev.stage_code);
    }
    double max_abs_x0 = std::max(std::abs(og_x0_0), std::abs(og_x0_1));
    int sx = (max_abs_x0 < 33.0) ? 0 : ((max_abs_x0 < 66.0) ? 1 : 2);
    st.eff_dx0_by_x0[sx].fillSingleMetric(adx0, ev.pass[MET_PAWAN_C2R]);

    int sr = (max_R0 < 60.0) ? 0 : ((max_R0 < 80.0) ? 1 : 2);
    st.eff_dx0_by_R0[sr].fillSingleMetric(adx0, ev.pass[MET_PAWAN_C2R]);

    // 2D single-track loss maps: track 0 and track 1
    bool lost0 = (ev.lt <= 0) || (ev.lt == 1 && ev.lt1_kept_idx == 1);
    bool lost1 = (ev.lt <= 0) || (ev.lt == 1 && ev.lt1_kept_idx == 0);
    st.loss_map_st1.fill(og_x0_0, og_y0_0, lost0);
    st.loss_map_st1.fill(og_x0_1, og_y0_1, lost1);
    st.loss_map_st3.fill(og_x1_0, og_y1_0, lost0);
    st.loss_map_st3.fill(og_x1_1, og_y1_1, lost1);
}

void fillPairedStats(PairedDatasetStats& ps,
                     const EventEval& evA,
                     const EventEval& evB,
                     const double seps[N_CURVES],
                     const bool in_reg[N_REGIMES]) {
    int la = std::min(std::max(evA.lt, 0), 4);
    int lb = std::min(std::max(evB.lt, 0), 4);
    for (int r = 0; r < N_REGIMES; ++r) {
        if (!in_reg[r]) continue;
        ps.trans[r][la][lb]++;
        for (int m = 0; m < N_METRICS; ++m) {
            ps.reg_cells[r][m].add(evA.pass[m], evB.pass[m]);
        }
    }
    for (int c = 0; c < N_CURVES; ++c) {
        ps.curves[c].fill(seps[c],
                          evA.pass[MET_PAWAN_C2R], evB.pass[MET_PAWAN_C2R],
                          evA.pass[MET_LT2], evB.pass[MET_LT2]);
    }
}

struct ThreadAccum {
    std::array<std::array<VersionDatasetStats, N_RAW_DS>, N_MAX_VER> vstats;
    std::array<std::array<PairedDatasetStats, N_RAW_DS>, N_MAX_PAIRS> pstats;
    long long og_checks_passed = 0;
    long long og_checks_failed = 0;

    void merge(const ThreadAccum& o, int n_ver, int n_pairs) {
        og_checks_passed += o.og_checks_passed;
        og_checks_failed += o.og_checks_failed;
        for (int v = 0; v < n_ver; ++v)
            for (int d = 0; d < N_RAW_DS; ++d)
                vstats[v][d].merge(o.vstats[v][d]);
        for (int p = 0; p < n_pairs; ++p)
            for (int d = 0; d < N_RAW_DS; ++d)
                pstats[p][d].merge(o.pstats[p][d]);
    }
};

struct TreeBinding {
    std::unique_ptr<TFile> file;
    std::unique_ptr<TTreeReader> reader;
    std::unique_ptr<TTreeReaderValue<Int_t>> lt, og_lt;
    std::unique_ptr<TTreeReaderValue<std::vector<double>>> c2, nd, x0, y0, pz0, x1, y1;
    std::unique_ptr<TTreeReaderValue<std::vector<double>>> og_x0, og_y0, og_z0, og_px0, og_py0, og_pz0;
    std::unique_ptr<TTreeReaderValue<std::vector<double>>> og_x1, og_y1, og_z1, og_px1, og_py1, og_pz1;
    std::unique_ptr<TTreeReaderValue<UInt_t>> nSeg, nGBSeg;
    std::unique_ptr<TTreeReaderValue<std::vector<double>>> seg_x, seg_y, seg_z;
    std::unique_ptr<TTreeReaderValue<std::vector<double>>> gb_x, gb_y, gb_z;
    Long64_t n_entries = 0;

    void open(const std::string& path) {
        file.reset(TFile::Open(path.c_str(), "READ"));
        if (!file || file->IsZombie()) {
            throw std::runtime_error("Cannot open ROOT file: " + path);
        }
        TTree* t = dynamic_cast<TTree*>(file->Get("nt"));
        if (!t) throw std::runtime_error("Missing TTree 'nt' in: " + path);
        n_entries = t->GetEntries();
        reader = std::make_unique<TTreeReader>(t);
        lt = std::make_unique<TTreeReaderValue<Int_t>>(*reader, "longTracks");
        og_lt = std::make_unique<TTreeReaderValue<Int_t>>(*reader, "og_longTracks");
        c2 = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "Track_Chi2");
        nd = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "Track_nDoF");
        x0 = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "Track_x0");
        y0 = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "Track_y0");
        pz0 = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "Track_pz0");
        x1 = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "Track_x1");
        y1 = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "Track_y1");
        og_x0 = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "og_Track_x0");
        og_y0 = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "og_Track_y0");
        og_z0 = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "og_Track_z0");
        og_px0 = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "og_Track_px0");
        og_py0 = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "og_Track_py0");
        og_pz0 = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "og_Track_pz0");
        og_x1 = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "og_Track_x1");
        og_y1 = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "og_Track_y1");
        og_z1 = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "og_Track_z1");
        og_px1 = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "og_Track_px1");
        og_py1 = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "og_Track_py1");
        og_pz1 = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "og_Track_pz1");
        nSeg = std::make_unique<TTreeReaderValue<UInt_t>>(*reader, "TrackSegments");
        seg_x = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "TrackSegment_x");
        seg_y = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "TrackSegment_y");
        seg_z = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "TrackSegment_z");
        nGBSeg = std::make_unique<TTreeReaderValue<UInt_t>>(*reader, "GhostBustedTrackSegments");
        gb_x = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "GhostBustedTrackSegment_x");
        gb_y = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "GhostBustedTrackSegment_y");
        gb_z = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "GhostBustedTrackSegment_z");
    }
};

std::string findLatestFile(const std::string& latest_dir,
                           const std::map<std::string, std::string>& scanned_map,
                           const std::string& ds_name,
                           int pbin, int nbin) {
    std::string key = ds_name + "_PosBin" + std::to_string(pbin) + "_NegBin" + std::to_string(nbin);
    auto it = scanned_map.find(key);
    if (it != scanned_map.end()) return it->second;

    std::vector<std::string> candidates = {
        latest_dir + "/" + key + "/" + key + "_nt.root",
        latest_dir + "/" + key + "/nt_faser.root",
        latest_dir + "/" + ds_name + "/myOverlay_PosBin" + std::to_string(pbin) + "_NegBin" + std::to_string(nbin) + "-RDO.pool-ntuple.root",
        latest_dir + "/" + ds_name + "/" + key + "_nt.root",
        latest_dir + "/" + key + "_nt.root"
    };
    for (const auto& c : candidates) {
        if (fs::exists(c)) return c;
    }
    return "";
}

std::map<std::string, std::string> indexLatestDir(const std::string& latest_dir) {
    std::map<std::string, std::string> m;
    if (latest_dir.empty() || !fs::exists(latest_dir)) return m;
    std::regex re("(008749|010698|Pos_008749-Neg_010698|10004X).*PosBin([0-8])_NegBin([0-8])");
    for (const auto& e : fs::recursive_directory_iterator(latest_dir)) {
        if (!e.is_regular_file()) continue;
        std::string p = e.path().string();
        if (p.size() < 5 || p.substr(p.size() - 5) != ".root") continue;
        if (p.find("raw") != std::string::npos) continue;
        std::smatch sm;
        if (std::regex_search(p, sm, re)) {
            std::string key = sm[1].str() + "_PosBin" + sm[2].str() + "_NegBin" + sm[3].str();
            m[key] = p;
        }
    }
    return m;
}

// JSON helper writers
void writeEffScalarJson(std::ostream& os, long long pass, long long den) {
    double eff, lo, hi;
    clopperPearson1Sigma(pass, den, eff, lo, hi);
    os << "{\"pass\": " << pass << ", \"den\": " << den
       << ", \"eff\": " << fmtDoubleOrNull(eff)
       << ", \"err_lo\": " << fmtDoubleOrNull(lo)
       << ", \"err_hi\": " << fmtDoubleOrNull(hi) << "}";
}

void writePairedScalarJson(std::ostream& os, const PairedCell& c) {
    double d, e;
    mcnemarDelta(c.bp, c.g, c.l, c.bf, d, e);
    long long n = c.bp + c.g + c.l + c.bf;
    os << "{\"n_den\": " << n
       << ", \"n_both_pass\": " << c.bp
       << ", \"n_gain\": " << c.g
       << ", \"n_loss\": " << c.l
       << ", \"n_both_fail\": " << c.bf
       << ", \"delta\": " << fmtDoubleOrNull(d)
       << ", \"delta_err\": " << fmtDoubleOrNull(e) << "}";
}

inline std::string jsonEscape(const std::string& s) {
    std::string out;
    out.reserve(s.size() + 8);
    for (char c : s) {
        if (c == '\\') out += "\\\\";
        else if (c == '"') out += "\\\"";
        else out += c;
    }
    return out;
}

void writeEffCurveJson(std::ostream& os, const std::string& name, const std::string& xlabel,
                       const std::string& unit, const EffBinSet& bs, int metric_idx) {
    int nb = bs.nbins;
    os << "{\"name\": \"" << jsonEscape(name) << "\", \"xlabel\": \"" << jsonEscape(xlabel)
       << "\", \"unit\": \"" << jsonEscape(unit) << "\", ";
    std::vector<double> edges(nb + 1), centers(nb), hws(nb);
    if (bs.log_x) {
        double llo = std::log10(bs.xmin), lhi = std::log10(bs.xmax);
        for (int i = 0; i <= nb; ++i) edges[i] = std::pow(10.0, llo + (lhi - llo) * i / nb);
    } else {
        double bw = (bs.xmax - bs.xmin) / nb;
        for (int i = 0; i <= nb; ++i) edges[i] = bs.xmin + i * bw;
    }
    for (int i = 0; i < nb; ++i) {
        centers[i] = 0.5 * (edges[i] + edges[i + 1]);
        hws[i] = 0.5 * (edges[i + 1] - edges[i]);
    }
    auto writeArrD = [&](const char* k, const std::vector<double>& a) {
        os << "\"" << k << "\": [";
        for (size_t i = 0; i < a.size(); ++i) {
            if (i > 0) os << ", ";
            os << fmtDoubleOrNull(a[i]);
        }
        os << "], ";
    };
    auto writeArrL = [&](const char* k, const std::vector<long long>& a) {
        os << "\"" << k << "\": [";
        for (size_t i = 0; i < a.size(); ++i) {
            if (i > 0) os << ", ";
            os << a[i];
        }
        os << "], ";
    };
    writeArrD("edges", edges);
    writeArrD("centers", centers);
    writeArrD("half_widths", hws);
    os << "\"underflow_den\": " << bs.underflow_den
       << ", \"overflow_den\": " << bs.overflow_den
       << ", \"invalid_den\": " << bs.invalid_den << ", ";
    writeArrL("den", bs.den);
    writeArrL("pass", bs.pass[metric_idx]);
    std::vector<double> eff(nb), elo(nb), ehi(nb);
    for (int i = 0; i < nb; ++i) clopperPearson1Sigma(bs.pass[metric_idx][i], bs.den[i], eff[i], elo[i], ehi[i]);
    writeArrD("eff", eff);
    writeArrD("err_lo", elo);
    os << "\"err_hi\": [";
    for (int i = 0; i < nb; ++i) { if (i > 0) os << ", "; os << fmtDoubleOrNull(ehi[i]); }
    os << "]}";
}

void writeRatioCurveJson(std::ostream& os, const EffBinSet& mc, const EffBinSet& data, int metric_idx) {
    int nb = mc.nbins;
    double bw = (mc.xmax - mc.xmin) / nb;
    os << "{\"edges\": [";
    for (int i = 0; i <= nb; ++i) { if (i > 0) os << ", "; os << fmtDoubleOrNull(mc.xmin + i * bw); }
    os << "], \"centers\": [";
    for (int i = 0; i < nb; ++i) { if (i > 0) os << ", "; os << fmtDoubleOrNull(mc.xmin + (i + 0.5) * bw); }
    os << "], \"half_widths\": [";
    for (int i = 0; i < nb; ++i) { if (i > 0) os << ", "; os << fmtDoubleOrNull(0.5 * bw); }
    os << "], ";
    std::vector<double> rat(nb), rlo(nb), rhi(nb);
    for (int i = 0; i < nb; ++i) {
        double em, lom, him, ed, lod, hid;
        clopperPearson1Sigma(mc.pass[metric_idx][i], mc.den[i], em, lom, him);
        clopperPearson1Sigma(data.pass[metric_idx][i], data.den[i], ed, lod, hid);
        if (std::isfinite(em) && std::isfinite(ed) && em > 0.0 && ed > 0.0) {
            rat[i] = em / ed;
            rlo[i] = rat[i] * std::hypot(lom / em, hid / ed);
            rhi[i] = rat[i] * std::hypot(him / em, lod / ed);
        } else {
            rat[i] = rlo[i] = rhi[i] = std::numeric_limits<double>::quiet_NaN();
        }
    }
    auto writeArr = [&](const char* k, const std::vector<double>& a, bool comma) {
        os << "\"" << k << "\": [";
        for (int i = 0; i < nb; ++i) { if (i > 0) os << ", "; os << fmtDoubleOrNull(a[i]); }
        os << "]" << (comma ? ", " : "");
    };
    writeArr("ratio", rat, true);
    writeArr("err_lo", rlo, true);
    writeArr("err_hi", rhi, false);
    os << "}";
}

void writeWeightedCurveJson(std::ostream& os, const std::string& name, const WeightedBinSet& ws) {
    int nb = ws.nbins;
    double bw = (ws.xmax - ws.xmin) / nb;
    os << "{\"name\": \"" << name << "\", \"edges\": [";
    for (int i = 0; i <= nb; ++i) { if (i > 0) os << ", "; os << fmtDoubleOrNull(ws.xmin + i * bw); }
    os << "], \"centers\": [";
    for (int i = 0; i < nb; ++i) { if (i > 0) os << ", "; os << fmtDoubleOrNull(ws.xmin + (i + 0.5) * bw); }
    os << "], \"half_widths\": [";
    for (int i = 0; i < nb; ++i) { if (i > 0) os << ", "; os << fmtDoubleOrNull(0.5 * bw); }
    os << "], ";
    auto writeVec = [&](const char* k, const std::vector<double>& v, bool comma) {
        os << "\"" << k << "\": [";
        for (int i = 0; i < nb; ++i) { if (i > 0) os << ", "; os << fmtDoubleOrNull(v[i]); }
        os << "]" << (comma ? ", " : "");
    };
    writeVec("sumw_den", ws.sumw_den, true);
    writeVec("sumw2_den", ws.sumw2_den, true);
    writeVec("sumw_pass", ws.sumw_pass, true);
    writeVec("sumw2_pass", ws.sumw2_pass, true);
    writeVec("sumw2_fail", ws.sumw2_fail, true);
    std::vector<double> eff(nb), err(nb);
    for (int i = 0; i < nb; ++i) {
        if (ws.sumw_den[i] > 0.0) {
            eff[i] = ws.sumw_pass[i] / ws.sumw_den[i];
            double var = (ws.sumw2_pass[i] * (1.0 - eff[i]) * (1.0 - eff[i]) +
                          ws.sumw2_fail[i] * eff[i] * eff[i]) /
                         (ws.sumw_den[i] * ws.sumw_den[i]);
            err[i] = std::sqrt(std::max(0.0, var));
        } else {
            eff[i] = err[i] = std::numeric_limits<double>::quiet_NaN();
        }
    }
    writeVec("eff", eff, true);
    writeVec("err", err, false);
    os << "}";
}

void writeWeightedRatioCurveJson(std::ostream& os, const WeightedBinSet& wmc, const EffBinSet& data) {
    int nb = wmc.nbins;
    double bw = (wmc.xmax - wmc.xmin) / nb;
    os << "{\"edges\": [";
    for (int i = 0; i <= nb; ++i) { if (i > 0) os << ", "; os << fmtDoubleOrNull(wmc.xmin + i * bw); }
    os << "], \"centers\": [";
    for (int i = 0; i < nb; ++i) { if (i > 0) os << ", "; os << fmtDoubleOrNull(wmc.xmin + (i + 0.5) * bw); }
    os << "], \"half_widths\": [";
    for (int i = 0; i < nb; ++i) { if (i > 0) os << ", "; os << fmtDoubleOrNull(0.5 * bw); }
    os << "], ";
    std::vector<double> rat(nb), rlo(nb), rhi(nb);
    for (int i = 0; i < nb; ++i) {
        double ed, lod, hid;
        clopperPearson1Sigma(data.pass[MET_PAWAN_C2R][i], data.den[i], ed, lod, hid);
        if (wmc.sumw_den[i] > 0.0 && std::isfinite(ed) && ed > 0.0) {
            double em = wmc.sumw_pass[i] / wmc.sumw_den[i];
            double errm = std::sqrt(std::max(0.0,
                (wmc.sumw2_pass[i] * (1.0 - em) * (1.0 - em) + wmc.sumw2_fail[i] * em * em))) / wmc.sumw_den[i];
            if (em > 0.0) {
                rat[i] = em / ed;
                rlo[i] = rat[i] * std::hypot(errm / em, hid / ed);
                rhi[i] = rat[i] * std::hypot(errm / em, lod / ed);
            } else {
                rat[i] = rlo[i] = rhi[i] = std::numeric_limits<double>::quiet_NaN();
            }
        } else {
            rat[i] = rlo[i] = rhi[i] = std::numeric_limits<double>::quiet_NaN();
        }
    }
    auto writeArr = [&](const char* k, const std::vector<double>& a, bool comma) {
        os << "\"" << k << "\": [";
        for (int i = 0; i < nb; ++i) { if (i > 0) os << ", "; os << fmtDoubleOrNull(a[i]); }
        os << "]" << (comma ? ", " : "");
    };
    writeArr("ratio", rat, true);
    writeArr("err_lo", rlo, true);
    writeArr("err_hi", rhi, false);
    os << "}";
}

void writePairedCurveJson(std::ostream& os,
                          const std::string& pair_key, const std::string& ver_a, const std::string& ver_b,
                          const std::string& dataset, const std::string& metric,
                          const PairedBinSet& pb, int midx) {
    int nb = pb.nbins;
    double bw = (pb.xmax - pb.xmin) / nb;
    os << "{\"pair\": \"" << pair_key << "\", \"ver_a\": \"" << ver_a << "\", \"ver_b\": \"" << ver_b
       << "\", \"dataset\": \"" << dataset << "\", \"metric\": \"" << metric << "\", \"edges\": [";
    for (int i = 0; i <= nb; ++i) { if (i > 0) os << ", "; os << fmtDoubleOrNull(pb.xmin + i * bw); }
    os << "], \"centers\": [";
    for (int i = 0; i < nb; ++i) { if (i > 0) os << ", "; os << fmtDoubleOrNull(pb.xmin + (i + 0.5) * bw); }
    os << "], \"half_widths\": [";
    for (int i = 0; i < nb; ++i) { if (i > 0) os << ", "; os << fmtDoubleOrNull(0.5 * bw); }
    os << "], ";
    const auto& cs = pb.cells[midx];
    auto writeL = [&](const char* k, auto fn) {
        os << "\"" << k << "\": [";
        for (int i = 0; i < nb; ++i) { if (i > 0) os << ", "; os << fn(cs[i]); }
        os << "], ";
    };
    writeL("n_den", [](const PairedCell& c){ return c.bp + c.g + c.l + c.bf; });
    writeL("n_both_pass", [](const PairedCell& c){ return c.bp; });
    writeL("n_gain", [](const PairedCell& c){ return c.g; });
    writeL("n_loss", [](const PairedCell& c){ return c.l; });
    writeL("n_both_fail", [](const PairedCell& c){ return c.bf; });
    std::vector<double> d(nb), e(nb);
    for (int i = 0; i < nb; ++i) mcnemarDelta(cs[i].bp, cs[i].g, cs[i].l, cs[i].bf, d[i], e[i]);
    os << "\"delta\": [";
    for (int i = 0; i < nb; ++i) { if (i > 0) os << ", "; os << fmtDoubleOrNull(d[i]); }
    os << "], \"delta_err\": [";
    for (int i = 0; i < nb; ++i) { if (i > 0) os << ", "; os << fmtDoubleOrNull(e[i]); }
    os << "]}";
}

void writeStageCountsJson(std::ostream& os, const StageLossCounts& s) {
    os << "{\"n_den\": " << s.n_den
       << ", \"pass_c2R\": " << s.pass_c2R
       << ", \"lt2_fail_c2R\": " << s.lt2_fail_c2R
       << ", \"lt_ge3\": " << s.lt_ge3
       << ", \"lt0\": " << s.lt0
       << ", \"lt1_total\": " << s.lt1_total
       << ", \"lt1_miss_seg_st1\": " << s.lt1_miss_seg_st1
       << ", \"lt1_miss_seg_st2\": " << s.lt1_miss_seg_st2
       << ", \"lt1_miss_seg_st3\": " << s.lt1_miss_seg_st3
       << ", \"lt1_miss_seg_multi\": " << s.lt1_miss_seg_multi
       << ", \"lt1_killed_by_gb\": " << s.lt1_killed_by_gb
       << ", \"lt1_all3gb_lost_in_ckf\": " << s.lt1_all3gb_lost_in_ckf
       << ", \"lt1_unmatched_kept\": " << s.lt1_unmatched_kept << "}";
}

void writeMap2DJson(std::ostream& os, const Map2D& m) {
    os << "{\"x_edges\": [";
    for (int i = 0; i <= N_MAP_BINS; ++i) { if (i > 0) os << ", "; os << (-100.0 + 5.0 * i); }
    os << "], \"y_edges\": [";
    for (int i = 0; i <= N_MAP_BINS; ++i) { if (i > 0) os << ", "; os << (-100.0 + 5.0 * i); }
    os << "], \"den\": [";
    for (int iy = 0; iy < N_MAP_BINS; ++iy) {
        if (iy > 0) os << ", ";
        os << "[";
        for (int ix = 0; ix < N_MAP_BINS; ++ix) { if (ix > 0) os << ", "; os << m.den[iy][ix]; }
        os << "]";
    }
    os << "], \"lost\": [";
    for (int iy = 0; iy < N_MAP_BINS; ++iy) {
        if (iy > 0) os << ", ";
        os << "[";
        for (int ix = 0; ix < N_MAP_BINS; ++ix) { if (ix > 0) os << ", "; os << m.lost[iy][ix]; }
        os << "]";
    }
    os << "], \"p_loss\": [";
    for (int iy = 0; iy < N_MAP_BINS; ++iy) {
        if (iy > 0) os << ", ";
        os << "[";
        for (int ix = 0; ix < N_MAP_BINS; ++ix) {
            if (ix > 0) os << ", ";
            double v = (m.den[iy][ix] > 0) ? static_cast<double>(m.lost[iy][ix]) / m.den[iy][ix]
                                           : std::numeric_limits<double>::quiet_NaN();
            os << fmtDoubleOrNull(v);
        }
        os << "]";
    }
    os << "]}";
}

int main(int argc, char** argv) {
    auto t_start = std::chrono::steady_clock::now();
    std::string eos_base = "/eos/experiment/faser/data0/overlay";
    std::string mode = "baseline-only";
    std::string latest_dir;
    std::string aux_file;
    std::string out_dir = "overlay_analysis_out";
    int n_threads = 16;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--eos-base" && i + 1 < argc) eos_base = argv[++i];
        else if (arg == "--mode" && i + 1 < argc) mode = argv[++i];
        else if (arg == "--baseline-only") mode = "baseline-only";
        else if (arg == "--full") mode = "full";
        else if (arg == "--latest-dir" && i + 1 < argc) latest_dir = argv[++i];
        else if (arg == "--aux-file" && i + 1 < argc) aux_file = argv[++i];
        else if (arg == "--out-dir" && i + 1 < argc) out_dir = argv[++i];
        else if (arg == "--threads" && i + 1 < argc) n_threads = std::stoi(argv[++i]);
        else {
            std::cerr << "Unknown argument: " << arg << "\n";
            return 2;
        }
    }

    const bool is_full = (mode == "full");
    if (is_full && latest_dir.empty()) {
        std::cerr << "ERROR: --mode full requires --latest-dir <dir>\n";
        return 2;
    }
    const int n_ver = is_full ? 3 : 2;
    const int n_pairs = is_full ? 3 : 1;

    ROOT::EnableThreadSafety();
    omp_set_num_threads(n_threads);

    const std::array<std::string, N_RAW_DS> raw_ds_names = {
        "008749", "010698", "Pos_008749-Neg_010698", "10004X"
    };

    auto scanned_latest = is_full ? indexLatestDir(latest_dir) : std::map<std::string, std::string>{};

    // Build exact 322 matched file list
    std::vector<FileTriple> triples;
    std::array<int, N_RAW_DS> ds_file_counts{};
    for (int d = 0; d < N_RAW_DS; ++d) {
        const auto& dn = raw_ds_names[d];
        for (int pb = 0; pb < 9; ++pb) {
            for (int nb = 0; nb < 9; ++nb) {
                std::string fname = "myOverlay_PosBin" + std::to_string(pb) + "_NegBin" + std::to_string(nb) + "-RDO.pool-ntuple.root";
                std::string p_bug = eos_base + "/" + dn + "/OverlaidNTPs_bug/" + fname;
                std::string p_jul = eos_base + "/" + dn + "/OverlaidNTPs/" + fname;
                if (!fs::exists(p_bug) || !fs::exists(p_jul)) continue;
                FileTriple ft;
                ft.ds_idx = d;
                ft.ds_name = dn;
                ft.pos_bin = pb;
                ft.neg_bin = nb;
                ft.path_may30 = p_bug;
                ft.path_jul8 = p_jul;
                if (is_full) {
                    ft.path_latest = findLatestFile(latest_dir, scanned_latest, dn, pb, nb);
                    if (ft.path_latest.empty()) {
                        std::cerr << "ERROR: Missing Latest GeoRotation file for " << dn
                                  << " PosBin" << pb << "_NegBin" << nb << " under " << latest_dir << "\n";
                        return 1;
                    }
                }
                triples.push_back(std::move(ft));
                ds_file_counts[d]++;
            }
        }
    }

    std::cout << "[overlay_analyze] Mode=" << mode << " threads=" << n_threads
              << " matched_triples=" << triples.size()
              << " (008749=" << ds_file_counts[0]
              << ", 010698=" << ds_file_counts[1]
              << ", Cross=" << ds_file_counts[2]
              << ", 10004X=" << ds_file_counts[3] << ")\n";

    // Pass 1: Paired lockstep scan over all 322 matched file triples
    std::vector<ThreadAccum> locals(n_threads);
    std::string fatal_error;
    std::atomic<bool> has_fatal{false};

    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t fi = 0; fi < triples.size(); ++fi) {
        int tid = omp_get_thread_num();
        if (has_fatal.load(std::memory_order_relaxed)) continue;

        try {
            const auto& ft = triples[fi];
            TreeBinding tb[N_MAX_VER];
            tb[VER_MAY30].open(ft.path_may30);
            tb[VER_JUL8].open(ft.path_jul8);
            if (is_full) tb[VER_LATEST].open(ft.path_latest);

            Long64_t N = tb[VER_JUL8].n_entries;
            if (tb[VER_MAY30].n_entries != N || (is_full && tb[VER_LATEST].n_entries != N)) {
                throw std::runtime_error("Entry count mismatch on " + ft.ds_name +
                                         " PosBin" + std::to_string(ft.pos_bin) + "_NegBin" + std::to_string(ft.neg_bin));
            }

            auto& acc = locals[tid];
            for (int v = 0; v < n_ver; ++v) acc.vstats[v][ft.ds_idx].n_files++;

            for (Long64_t e = 0; e < N; ++e) {
                bool ok = tb[VER_MAY30].reader->Next() && tb[VER_JUL8].reader->Next();
                if (is_full) ok = ok && tb[VER_LATEST].reader->Next();
                if (!ok) throw std::runtime_error("Premature EOF at entry " + std::to_string(e));

                // Hard per-event og_Track identity assertion across all versions
                const auto& ox0_ref = *tb[VER_JUL8].og_x0->Get();
                const auto& oy0_ref = *tb[VER_JUL8].og_y0->Get();
                const auto& oz0_ref = *tb[VER_JUL8].og_z0->Get();
                const auto& opx0_ref = *tb[VER_JUL8].og_px0->Get();
                const auto& opy0_ref = *tb[VER_JUL8].og_py0->Get();
                const auto& opz0_ref = *tb[VER_JUL8].og_pz0->Get();
                const auto& ox1_ref = *tb[VER_JUL8].og_x1->Get();
                const auto& oy1_ref = *tb[VER_JUL8].og_y1->Get();
                const auto& oz1_ref = *tb[VER_JUL8].og_z1->Get();
                const auto& opx1_ref = *tb[VER_JUL8].og_px1->Get();
                const auto& opy1_ref = *tb[VER_JUL8].og_py1->Get();
                const auto& opz1_ref = *tb[VER_JUL8].og_pz1->Get();

                bool id_ok = (*tb[VER_JUL8].og_lt->Get() == 2 && ox0_ref.size() == 2 && ox1_ref.size() == 2);
                for (int v = 0; v < n_ver && id_ok; ++v) {
                    if (v == VER_JUL8) continue;
                    const auto& ox0_v = *tb[v].og_x0->Get();
                    const auto& oy0_v = *tb[v].og_y0->Get();
                    const auto& opz0_v = *tb[v].og_pz0->Get();
                    if (*tb[v].og_lt->Get() != 2 || ox0_v.size() != 2 ||
                        std::abs(ox0_v[0] - ox0_ref[0]) > 1e-4 ||
                        std::abs(oy0_v[0] - oy0_ref[0]) > 1e-4 ||
                        std::abs(opz0_v[0] - opz0_ref[0]) > 1e-2 ||
                        std::abs(ox0_v[1] - ox0_ref[1]) > 1e-4 ||
                        std::abs(oy0_v[1] - oy0_ref[1]) > 1e-4 ||
                        std::abs(opz0_v[1] - opz0_ref[1]) > 1e-2) {
                        id_ok = false;
                    }
                }
                if (!id_ok) {
                    acc.og_checks_failed++;
                    throw std::runtime_error("og_Track identity mismatch in " + ft.ds_name +
                                             " PosBin" + std::to_string(ft.pos_bin) +
                                             "_NegBin" + std::to_string(ft.neg_bin) +
                                             " entry " + std::to_string(e));
                }
                acc.og_checks_passed++;

                double dx0 = ox0_ref[0] - ox0_ref[1];
                double dy0 = oy0_ref[0] - oy0_ref[1];
                double dr0 = std::hypot(dx0, dy0);
                double dx1 = ox1_ref[0] - ox1_ref[1];
                double dy1 = oy1_ref[0] - oy1_ref[1];
                double dr1 = std::hypot(dx1, dy1);
                double adx0 = std::abs(dx0), ady0 = std::abs(dy0);
                double adx1 = std::abs(dx1), ady1 = std::abs(dy1);

                double seps[N_CURVES] = {
                    dr1, dr1, dr0, dr0,
                    ady0, ady0, ady1, ady1,
                    adx0, dx0, adx1, dx1,
                    dy1, dy0
                };
                bool in_reg[N_REGIMES] = {
                    true,
                    (dr1 < 5.0),
                    (dr1 >= 30.0 && dr1 <= 70.0),
                    (adx0 >= 70.0 && adx0 <= 100.0)
                };

                double p0_0_gev = std::sqrt(opx0_ref[0]*opx0_ref[0] + opy0_ref[0]*opy0_ref[0] + opz0_ref[0]*opz0_ref[0]) * 1e-3;
                double p0_1_gev = std::sqrt(opx0_ref[1]*opx0_ref[1] + opy0_ref[1]*opy0_ref[1] + opz0_ref[1]*opz0_ref[1]) * 1e-3;

                EventEval evs[N_MAX_VER];
                for (int v = 0; v < n_ver; ++v) {
                    evs[v] = evaluateEvent(
                        *tb[v].lt->Get(),
                        *tb[v].c2->Get(), *tb[v].nd->Get(),
                        *tb[v].x0->Get(), *tb[v].y0->Get(), *tb[v].pz0->Get(),
                        *tb[v].x1->Get(), *tb[v].y1->Get(),
                        ox0_ref, oy0_ref, oz0_ref, opx0_ref, opy0_ref, opz0_ref,
                        ox1_ref, oy1_ref, oz1_ref, opx1_ref, opy1_ref, opz1_ref,
                        *tb[v].nSeg->Get(), *tb[v].seg_x->Get(), *tb[v].seg_y->Get(), *tb[v].seg_z->Get(),
                        *tb[v].nGBSeg->Get(), *tb[v].gb_x->Get(), *tb[v].gb_y->Get(), *tb[v].gb_z->Get());

                    fillVersionStats(acc.vstats[v][ft.ds_idx], evs[v], seps, in_reg,
                                     ox0_ref[0], oy0_ref[0], ox1_ref[0], oy1_ref[0],
                                     ox0_ref[1], oy0_ref[1], ox1_ref[1], oy1_ref[1],
                                     p0_0_gev, p0_1_gev, opz0_ref[0], opz0_ref[1],
                                     *tb[v].pz0->Get());
                }

                fillPairedStats(acc.pstats[PAIR_JUL8_MAY30][ft.ds_idx], evs[VER_MAY30], evs[VER_JUL8], seps, in_reg);
                if (is_full) {
                    fillPairedStats(acc.pstats[PAIR_LATEST_JUL8][ft.ds_idx], evs[VER_JUL8], evs[VER_LATEST], seps, in_reg);
                    fillPairedStats(acc.pstats[PAIR_LATEST_MAY30][ft.ds_idx], evs[VER_MAY30], evs[VER_LATEST], seps, in_reg);
                }
            }
        } catch (const std::exception& ex) {
            #pragma omp critical
            {
                if (fatal_error.empty()) fatal_error = ex.what();
                has_fatal.store(true, std::memory_order_relaxed);
            }
        }
    }

    if (!fatal_error.empty()) {
        std::cerr << "FATAL ERROR during Pass 1: " << fatal_error << "\n";
        return 1;
    }

    ThreadAccum merged;
    for (int t = 0; t < n_threads; ++t) merged.merge(locals[t], n_ver, n_pairs);

    // Build combined dataset maps (6 named datasets)
    const std::vector<std::string> ds_order = {
        "All3Data", "10004X", "008749", "010698", "Pos_008749-Neg_010698", "008749+010698"
    };
    std::vector<std::string> ver_names = is_full
        ? std::vector<std::string>{"baseline_may30", "baseline_jul8", "latest_georot"}
        : std::vector<std::string>{"baseline_may30", "baseline_jul8"};
    std::vector<std::string> pair_names = is_full
        ? std::vector<std::string>{"baseline_jul8_minus_baseline_may30", "latest_georot_minus_baseline_jul8", "latest_georot_minus_baseline_may30"}
        : std::vector<std::string>{"baseline_jul8_minus_baseline_may30"};
    std::vector<std::pair<std::string, std::string>> pair_ab = {
        {"baseline_may30", "baseline_jul8"},
        {"baseline_jul8", "latest_georot"},
        {"baseline_may30", "latest_georot"}
    };

    std::map<std::string, std::map<std::string, VersionDatasetStats>> all_vstats;
    for (int v = 0; v < n_ver; ++v) {
        const auto& vn = ver_names[v];
        all_vstats[vn]["008749"] = merged.vstats[v][DS_008749];
        all_vstats[vn]["010698"] = merged.vstats[v][DS_010698];
        all_vstats[vn]["Pos_008749-Neg_010698"] = merged.vstats[v][DS_CROSS];
        all_vstats[vn]["10004X"] = merged.vstats[v][DS_10004X];
        VersionDatasetStats d8710 = merged.vstats[v][DS_008749];
        d8710.merge(merged.vstats[v][DS_010698]);
        all_vstats[vn]["008749+010698"] = d8710;
        VersionDatasetStats all3 = d8710;
        all3.merge(merged.vstats[v][DS_CROSS]);
        all_vstats[vn]["All3Data"] = all3;
    }

    std::map<std::string, std::map<std::string, PairedDatasetStats>> all_pstats;
    for (int p = 0; p < n_pairs; ++p) {
        const auto& pn = pair_names[p];
        all_pstats[pn]["008749"] = merged.pstats[p][DS_008749];
        all_pstats[pn]["010698"] = merged.pstats[p][DS_010698];
        all_pstats[pn]["Pos_008749-Neg_010698"] = merged.pstats[p][DS_CROSS];
        all_pstats[pn]["10004X"] = merged.pstats[p][DS_10004X];
        PairedDatasetStats p8710 = merged.pstats[p][DS_008749];
        p8710.merge(merged.pstats[p][DS_010698]);
        all_pstats[pn]["008749+010698"] = p8710;
        PairedDatasetStats pall3 = p8710;
        pall3.merge(merged.pstats[p][DS_CROSS]);
        all_pstats[pn]["All3Data"] = pall3;
    }

    // Pass 2: MC -> Data og_p0 momentum reweighting (40 log10 bins [50, 5000] GeV)
    std::array<double, N_RW_P0_BINS> data_p0_frac{}, mc_p0_frac{}, rw_weights{};
    const auto& data_ref = all_vstats["baseline_jul8"]["All3Data"];
    const auto& mc_ref = all_vstats["baseline_jul8"]["10004X"];
    for (int b = 0; b < N_RW_P0_BINS; ++b) {
        data_p0_frac[b] = (data_ref.og_p0_40_tot > 0) ? static_cast<double>(data_ref.og_p0_40[b]) / data_ref.og_p0_40_tot : 0.0;
        mc_p0_frac[b] = (mc_ref.og_p0_40_tot > 0) ? static_cast<double>(mc_ref.og_p0_40[b]) / mc_ref.og_p0_40_tot : 0.0;
        rw_weights[b] = (mc_p0_frac[b] > 0.0) ? (data_p0_frac[b] / mc_p0_frac[b]) : 0.0;
    }

    std::vector<const FileTriple*> mc_triples;
    for (const auto& ft : triples) {
        if (ft.ds_idx == DS_10004X) mc_triples.push_back(&ft);
    }

    std::vector<std::array<ReweightStats, N_MAX_VER>> rw_locals(n_threads);
    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t fi = 0; fi < mc_triples.size(); ++fi) {
        int tid = omp_get_thread_num();
        const auto& ft = *mc_triples[fi];
        std::string paths[N_MAX_VER] = {ft.path_may30, ft.path_jul8, ft.path_latest};
        for (int v = 0; v < n_ver; ++v) {
            std::unique_ptr<TFile> f(TFile::Open(paths[v].c_str(), "READ"));
            TTree* t = dynamic_cast<TTree*>(f->Get("nt"));
            TTreeReader rd(t);
            TTreeReaderValue<Int_t> lt(rd, "longTracks");
            TTreeReaderValue<std::vector<double>> c2(rd, "Track_Chi2"), nd(rd, "Track_nDoF");
            TTreeReaderValue<std::vector<double>> x0(rd, "Track_x0"), y0(rd, "Track_y0");
            TTreeReaderValue<std::vector<double>> x1(rd, "Track_x1"), y1(rd, "Track_y1");
            TTreeReaderValue<std::vector<double>> ox0(rd, "og_Track_x0"), oy0(rd, "og_Track_y0");
            TTreeReaderValue<std::vector<double>> opx0(rd, "og_Track_px0"), opy0(rd, "og_Track_py0"), opz0(rd, "og_Track_pz0");
            TTreeReaderValue<std::vector<double>> ox1(rd, "og_Track_x1"), oy1(rd, "og_Track_y1");

            auto& rw = rw_locals[tid][v];
            while (rd.Next()) {
                double p0 = std::sqrt((*opx0)[0]*(*opx0)[0] + (*opy0)[0]*(*opy0)[0] + (*opz0)[0]*(*opz0)[0]) * 1e-3;
                double p1 = std::sqrt((*opx0)[1]*(*opx0)[1] + (*opy0)[1]*(*opy0)[1] + (*opz0)[1]*(*opz0)[1]) * 1e-3;
                auto binP = [](double p) {
                    double lp = std::log10(std::clamp(p, 50.0, 5000.0));
                    int b = static_cast<int>((lp - RW_P0_LOG_MIN) / (RW_P0_LOG_MAX - RW_P0_LOG_MIN) * N_RW_P0_BINS);
                    return std::clamp(b, 0, N_RW_P0_BINS - 1);
                };
                int b0 = binP(p0), b1 = binP(p1);
                double w = rw_weights[b0] * rw_weights[b1];

                bool pass_c2r = false;
                if (*lt == 2 && x0->size() == 2) {
                    bool ok0 = ((*nd)[0] > 0.0 && (*c2)[0]/(*nd)[0] <= 25.0 &&
                                std::hypot((*x0)[0], (*y0)[0]) < 95.0 && std::hypot((*x1)[0], (*y1)[0]) < 95.0);
                    bool ok1 = ((*nd)[1] > 0.0 && (*c2)[1]/(*nd)[1] <= 25.0 &&
                                std::hypot((*x0)[1], (*y0)[1]) < 95.0 && std::hypot((*x1)[1], (*y1)[1]) < 95.0);
                    pass_c2r = (ok0 && ok1);
                }

                double dx0 = (*ox0)[0] - (*ox0)[1], dy0 = (*oy0)[0] - (*oy0)[1];
                double dx1 = (*ox1)[0] - (*ox1)[1], dy1 = (*oy1)[0] - (*oy1)[1];
                double dr0 = std::hypot(dx0, dy0), dr1 = std::hypot(dx1, dy1);
                double vals[N_RW_CURVES] = {dr1, dr1, dr0, std::abs(dx0), std::abs(dy0)};

                rw.n_events++;
                rw.sumw += w;
                rw.sumw2 += w * w;
                rw.closure_sumw[b0] += w;
                rw.closure_sumw[b1] += w;
                rw.closure_tot += 2.0 * w;
                for (int c = 0; c < N_RW_CURVES; ++c) rw.curves[c].fill(vals[c], w, pass_c2r);
            }
        }
    }

    std::array<ReweightStats, N_MAX_VER> rw_merged;
    for (int t = 0; t < n_threads; ++t)
        for (int v = 0; v < n_ver; ++v) rw_merged[v].merge(rw_locals[t][v]);

    // Optional auxiliary shard (010698 PosBin2_NegBin3)
    if (is_full && aux_file.empty()) {
        aux_file = findLatestFile(latest_dir, scanned_latest, "010698", 2, 3);
    }
    bool has_aux = (!aux_file.empty() && fs::exists(aux_file));
    VersionDatasetStats aux_stats;
    if (has_aux) {
        TreeBinding tb;
        tb.open(aux_file);
        aux_stats.n_files = 1;
        while (tb.reader->Next()) {
            const auto& ox0 = *tb.og_x0->Get();
            const auto& oy0 = *tb.og_y0->Get();
            const auto& oz0 = *tb.og_z0->Get();
            const auto& opx0 = *tb.og_px0->Get();
            const auto& opy0 = *tb.og_py0->Get();
            const auto& opz0 = *tb.og_pz0->Get();
            const auto& ox1 = *tb.og_x1->Get();
            const auto& oy1 = *tb.og_y1->Get();
            const auto& oz1 = *tb.og_z1->Get();
            const auto& opx1 = *tb.og_px1->Get();
            const auto& opy1 = *tb.og_py1->Get();
            const auto& opz1 = *tb.og_pz1->Get();
            if (*tb.og_lt->Get() != 2 || ox0.size() < 2) continue;
            double dx0 = ox0[0]-ox0[1], dy0 = oy0[0]-oy0[1], dr0 = std::hypot(dx0, dy0);
            double dx1 = ox1[0]-ox1[1], dy1 = oy1[0]-oy1[1], dr1 = std::hypot(dx1, dy1);
            double seps[N_CURVES] = {dr1, dr1, dr0, dr0, std::abs(dy0), std::abs(dy0), std::abs(dy1), std::abs(dy1),
                                     std::abs(dx0), dx0, std::abs(dx1), dx1, dy1, dy0};
            bool in_reg[N_REGIMES] = {true, dr1 < 5.0, dr1 >= 30.0 && dr1 <= 70.0, std::abs(dx0) >= 70.0 && std::abs(dx0) <= 100.0};
            double p0_0 = std::sqrt(opx0[0]*opx0[0]+opy0[0]*opy0[0]+opz0[0]*opz0[0])*1e-3;
            double p0_1 = std::sqrt(opx0[1]*opx0[1]+opy0[1]*opy0[1]+opz0[1]*opz0[1])*1e-3;
            auto ev = evaluateEvent(*tb.lt->Get(), *tb.c2->Get(), *tb.nd->Get(),
                                    *tb.x0->Get(), *tb.y0->Get(), *tb.pz0->Get(), *tb.x1->Get(), *tb.y1->Get(),
                                    ox0, oy0, oz0, opx0, opy0, opz0, ox1, oy1, oz1, opx1, opy1, opz1,
                                    *tb.nSeg->Get(), *tb.seg_x->Get(), *tb.seg_y->Get(), *tb.seg_z->Get(),
                                    *tb.nGBSeg->Get(), *tb.gb_x->Get(), *tb.gb_y->Get(), *tb.gb_z->Get());
            fillVersionStats(aux_stats, ev, seps, in_reg,
                             ox0[0], oy0[0], ox1[0], oy1[0], ox0[1], oy0[1], ox1[1], oy1[1],
                             p0_0, p0_1, opz0[0], opz0[1], *tb.pz0->Get());
        }
    }

    auto t_end = std::chrono::steady_clock::now();
    double elapsed_s = std::chrono::duration<double>(t_end - t_start).count();

    // Create output directories
    fs::create_directories(out_dir + "/bin_counts");

    // Write bin_counts/*.csv
    for (int v = 0; v < n_ver; ++v) {
        const auto& vn = ver_names[v];
        for (const std::string dn : {"All3Data", "10004X"}) {
            const auto& vs = all_vstats[vn][dn];
            for (int c = 0; c < N_CURVES; ++c) {
                const auto& sp = CURVE_SPECS[c];
                const auto& bs = vs.curves[c];
                double bw = (sp.xmax - sp.xmin) / sp.nbins;
                std::ofstream f(out_dir + "/bin_counts/" + vn + "_" + dn + "_" + sp.key + ".csv");
                f << "bin_lo,bin_hi,bin_center,den,pass_c2R,eff_c2R,err_lo_c2R,err_hi_c2R,"
                     "pass_lt2,eff_lt2,pass_ge2,eff_ge2,pass_matched,eff_matched\n";
                for (int i = 0; i < sp.nbins; ++i) {
                    double lo = sp.xmin + i * bw, hi = lo + bw, cen = 0.5 * (lo + hi);
                    double e0, l0, h0, e1, l1, h1, e2, l2, h2, e3, l3, h3;
                    clopperPearson1Sigma(bs.pass[0][i], bs.den[i], e0, l0, h0);
                    clopperPearson1Sigma(bs.pass[1][i], bs.den[i], e1, l1, h1);
                    clopperPearson1Sigma(bs.pass[2][i], bs.den[i], e2, l2, h2);
                    clopperPearson1Sigma(bs.pass[3][i], bs.den[i], e3, l3, h3);
                    f << std::setprecision(8) << lo << "," << hi << "," << cen << "," << bs.den[i] << ","
                      << bs.pass[0][i] << "," << (std::isfinite(e0)?e0:-1.0) << "," << l0 << "," << h0 << ","
                      << bs.pass[1][i] << "," << (std::isfinite(e1)?e1:-1.0) << ","
                      << bs.pass[2][i] << "," << (std::isfinite(e2)?e2:-1.0) << ","
                      << bs.pass[3][i] << "," << (std::isfinite(e3)?e3:-1.0) << "\n";
                }
            }
            // Momentum CSV
            std::ofstream fpz(out_dir + "/bin_counts/momentum_" + vn + "_" + dn + "_track_pz0.csv");
            fpz << "bin_lo_mev,bin_hi_mev,bin_center_mev,counts,norm,norm_err\n";
            for (int i = 0; i < N_PZ_BINS; ++i) {
                double lo = std::pow(10.0, PZ_LOG_MIN + (PZ_LOG_MAX - PZ_LOG_MIN) * i / N_PZ_BINS);
                double hi = std::pow(10.0, PZ_LOG_MIN + (PZ_LOG_MAX - PZ_LOG_MIN) * (i + 1) / N_PZ_BINS);
                double cen = std::sqrt(lo * hi);
                double norm = (vs.pz_reco_tot > 0) ? static_cast<double>(vs.pz_reco[i]) / vs.pz_reco_tot : 0.0;
                double err = (vs.pz_reco_tot > 0) ? std::sqrt(static_cast<double>(vs.pz_reco[i])) / vs.pz_reco_tot : 0.0;
                fpz << std::setprecision(8) << lo << "," << hi << "," << cen << ","
                    << vs.pz_reco[i] << "," << norm << "," << err << "\n";
            }
        }
        // Reweighted MC CSVs
        for (int rc = 0; rc < N_RW_CURVES; ++rc) {
            int cid = RW_CURVE_IDS[rc];
            const auto& sp = CURVE_SPECS[cid];
            const auto& ws = rw_merged[v].curves[rc];
            double bw = (sp.xmax - sp.xmin) / sp.nbins;
            std::ofstream frw(out_dir + "/bin_counts/reweighted_" + vn + "_" + sp.key + ".csv");
            frw << "bin_lo,bin_hi,bin_center,sumw_den,sumw2_den,sumw_pass,sumw2_pass,sumw2_fail,eff,err\n";
            for (int i = 0; i < sp.nbins; ++i) {
                double lo = sp.xmin + i * bw, hi = lo + bw, cen = 0.5 * (lo + hi);
                double eff = (ws.sumw_den[i] > 0.0) ? ws.sumw_pass[i] / ws.sumw_den[i] : -1.0;
                double err = (ws.sumw_den[i] > 0.0)
                                 ? std::sqrt(std::max(0.0, ws.sumw2_pass[i]*(1.0-eff)*(1.0-eff) + ws.sumw2_fail[i]*eff*eff)) / ws.sumw_den[i]
                                 : 0.0;
                frw << std::setprecision(8) << lo << "," << hi << "," << cen << ","
                    << ws.sumw_den[i] << "," << ws.sumw2_den[i] << ","
                    << ws.sumw_pass[i] << "," << ws.sumw2_pass[i] << "," << ws.sumw2_fail[i] << ","
                    << eff << "," << err << "\n";
            }
        }
    }
    for (int p = 0; p < n_pairs; ++p) {
        const auto& pn = pair_names[p];
        for (const std::string dn : {"All3Data", "10004X"}) {
            const auto& ps = all_pstats[pn][dn];
            for (int c = 0; c < N_CURVES; ++c) {
                const auto& sp = CURVE_SPECS[c];
                const auto& pb = ps.curves[c];
                double bw = (sp.xmax - sp.xmin) / sp.nbins;
                std::ofstream fp(out_dir + "/bin_counts/paired_" + pn + "_" + dn + "_" + sp.key + ".csv");
                fp << "bin_lo,bin_hi,bin_center,n_den,n_both_pass,n_gain,n_loss,n_both_fail,delta_c2R,delta_err_c2R,delta_lt2,delta_err_lt2\n";
                for (int i = 0; i < sp.nbins; ++i) {
                    double lo = sp.xmin + i * bw, hi = lo + bw, cen = 0.5 * (lo + hi);
                    const auto& c0 = pb.cells[0][i];
                    const auto& c1 = pb.cells[1][i];
                    double d0, e0, d1, e1;
                    mcnemarDelta(c0.bp, c0.g, c0.l, c0.bf, d0, e0);
                    mcnemarDelta(c1.bp, c1.g, c1.l, c1.bf, d1, e1);
                    fp << std::setprecision(8) << lo << "," << hi << "," << cen << ","
                       << (c0.bp+c0.g+c0.l+c0.bf) << "," << c0.bp << "," << c0.g << "," << c0.l << "," << c0.bf << ","
                       << (std::isfinite(d0)?d0:0.0) << "," << (std::isfinite(e0)?e0:0.0) << ","
                       << (std::isfinite(d1)?d1:0.0) << "," << (std::isfinite(e1)?e1:0.0) << "\n";
                }
            }
        }
    }

    // Write overlay_hists.root
    {
        std::unique_ptr<TFile> hf(TFile::Open((out_dir + "/overlay_hists.root").c_str(), "RECREATE"));
        for (int v = 0; v < n_ver; ++v) {
            const auto& vn = ver_names[v];
            for (const std::string dn : {"All3Data", "10004X"}) {
                const auto& vs = all_vstats[vn][dn];
                for (int c = 0; c < N_CURVES; ++c) {
                    const auto& sp = CURVE_SPECS[c];
                    std::string hname = vn + "_" + dn + "_" + sp.key + "_eff_c2R";
                    TH1D* h = new TH1D(hname.c_str(), sp.xlabel, sp.nbins, sp.xmin, sp.xmax);
                    for (int i = 0; i < sp.nbins; ++i) {
                        double e, lo, hi;
                        clopperPearson1Sigma(vs.curves[c].pass[MET_PAWAN_C2R][i], vs.curves[c].den[i], e, lo, hi);
                        if (std::isfinite(e)) {
                            h->SetBinContent(i + 1, e);
                            h->SetBinError(i + 1, 0.5 * (lo + hi));
                        }
                    }
                    h->Write();
                }
                std::string m1name = vn + "_" + dn + "_loss_map_st1";
                TH2D* h2 = new TH2D(m1name.c_str(), "Single-track loss rate St.1;og_x0 [mm];og_y0 [mm]",
                                    N_MAP_BINS, -100.0, 100.0, N_MAP_BINS, -100.0, 100.0);
                for (int iy = 0; iy < N_MAP_BINS; ++iy)
                    for (int ix = 0; ix < N_MAP_BINS; ++ix)
                        if (vs.loss_map_st1.den[iy][ix] > 0)
                            h2->SetBinContent(ix + 1, iy + 1, static_cast<double>(vs.loss_map_st1.lost[iy][ix]) / vs.loss_map_st1.den[iy][ix]);
                h2->Write();
            }
        }
        hf->Close();
    }

    // Write overlay_summary.json matching overlay_analysis_schema.md verbatim
    std::ofstream js(out_dir + "/overlay_summary.json");
    js << "{\n";
    js << "  \"schema_version\": \"1.0\",\n";
    js << "  \"mode\": \"" << (is_full ? "full" : "baseline_only") << "\",\n";
    js << "  \"primary_version\": \"" << (is_full ? "latest_georot" : "baseline_jul8") << "\",\n";
    js << "  \"reference_version\": \"" << (is_full ? "baseline_jul8" : "baseline_may30") << "\",\n";
    js << "  \"versions\": [";
    for (int v = 0; v < n_ver; ++v) { if (v > 0) js << ", "; js << "\"" << ver_names[v] << "\""; }
    js << "],\n";
    js << "  \"datasets\": [";
    for (size_t i = 0; i < ds_order.size(); ++i) { if (i > 0) js << ", "; js << "\"" << ds_order[i] << "\""; }
    js << "],\n";
    js << "  \"provenance\": {\n"
       << "    \"frame\": \"FASER\",\n"
       << "    \"cond_tag\": \"OFLCOND-FASER-05\",\n"
       << "    \"geom_tag\": \"TI12MC04 / FASERNU-04\",\n"
       << "    \"matched_bin_pairs\": {\"008749\": " << ds_file_counts[0]
       << ", \"010698\": " << ds_file_counts[1]
       << ", \"Pos_008749-Neg_010698\": " << ds_file_counts[2]
       << ", \"10004X\": " << ds_file_counts[3]
       << ", \"total\": " << triples.size() << "},\n"
       << "    \"og_identity_checks_passed\": " << merged.og_checks_passed << ",\n"
       << "    \"og_identity_checks_failed\": " << merged.og_checks_failed << ",\n"
       << "    \"runtime_seconds\": " << std::fixed << std::setprecision(3) << elapsed_s << ",\n"
       << "    \"n_threads\": " << n_threads << "\n"
       << "  },\n";

    // global_summary
    const char* reg_keys[N_REGIMES] = {"global", "small_dR1_lt5", "plateau_dR1_30_70", "large_abs_dX0_70_100"};
    const char* met_keys[N_METRICS] = {"eff_pawan_c2R", "eff_lt2", "eff_ge2", "eff_matched"};
    js << "  \"global_summary\": {\n";
    for (int v = 0; v < n_ver; ++v) {
        const auto& vn = ver_names[v];
        js << "    \"" << vn << "\": {\n";
        for (size_t di = 0; di < ds_order.size(); ++di) {
            const auto& dn = ds_order[di];
            const auto& s = all_vstats[vn][dn];
            js << "      \"" << dn << "\": {\n"
               << "        \"n_files\": " << s.n_files << ", \"n_events\": " << s.n_events << ", \"n_og2\": " << s.n_og2 << ",\n"
               << "        \"lt_counts\": {\"0\": " << s.lt_counts[0] << ", \"1\": " << s.lt_counts[1]
               << ", \"2\": " << s.lt_counts[2] << ", \"3\": " << s.lt_counts[3]
               << ", \"4\": " << s.lt_counts[4] << ", \"5plus\": " << s.lt_counts[5] << "},\n"
               << "        \"lt2_fail_c2R\": " << s.lt2_fail_c2R << ",\n";
            for (int m = 0; m < N_METRICS; ++m) {
                js << "        \"" << met_keys[m] << "\": ";
                writeEffScalarJson(js, s.reg_pass[REG_GLOBAL][m], s.reg_den[REG_GLOBAL]);
                js << ",\n";
            }
            js << "        \"regimes\": {\n";
            for (int r = 1; r < N_REGIMES; ++r) {
                js << "          \"" << reg_keys[r] << "\": {";
                for (int m = 0; m < N_METRICS; ++m) {
                    if (m > 0) js << ", ";
                    js << "\"" << met_keys[m] << "\": ";
                    writeEffScalarJson(js, s.reg_pass[r][m], s.reg_den[r]);
                }
                js << "}" << (r + 1 < N_REGIMES ? ",\n" : "\n");
            }
            js << "        }\n      }" << (di + 1 < ds_order.size() ? ",\n" : "\n");
        }
        js << "    }" << (v + 1 < n_ver ? ",\n" : "\n");
    }
    js << "  },\n";

    // paired_global_deltas
    js << "  \"paired_global_deltas\": {\n";
    for (int p = 0; p < n_pairs; ++p) {
        const auto& pn = pair_names[p];
        js << "    \"" << pn << "\": {\n";
        for (size_t di = 0; di < ds_order.size(); ++di) {
            const auto& dn = ds_order[di];
            const auto& ps = all_pstats[pn][dn];
            js << "      \"" << dn << "\": {\n";
            for (int r = 0; r < N_REGIMES; ++r) {
                js << "        \"" << reg_keys[r] << "\": {";
                for (int m = 0; m < N_METRICS; ++m) {
                    if (m > 0) js << ", ";
                    js << "\"" << met_keys[m] << "\": ";
                    writePairedScalarJson(js, ps.reg_cells[r][m]);
                }
                js << "}" << (r + 1 < N_REGIMES ? ",\n" : "\n");
            }
            js << "      }" << (di + 1 < ds_order.size() ? ",\n" : "\n");
        }
        js << "    }" << (p + 1 < n_pairs ? ",\n" : "\n");
    }
    js << "  },\n";

    // transition_matrices
    js << "  \"transition_matrices\": {\n";
    for (int p = 0; p < n_pairs; ++p) {
        const auto& pn = pair_names[p];
        js << "    \"" << pn << "\": {\n";
        for (size_t di = 0; di < ds_order.size(); ++di) {
            const auto& dn = ds_order[di];
            const auto& ps = all_pstats[pn][dn];
            js << "      \"" << dn << "\": {\"labels\": [\"0\", \"1\", \"2\", \"3\", \"4+\"], ";
            for (int r = 0; r < N_REGIMES; ++r) {
                js << "\"" << reg_keys[r] << "\": [";
                for (int a = 0; a < 5; ++a) {
                    if (a > 0) js << ", ";
                    js << "[";
                    for (int b = 0; b < 5; ++b) { if (b > 0) js << ", "; js << ps.trans[r][a][b]; }
                    js << "]";
                }
                js << "]" << (r + 1 < N_REGIMES ? ", " : "");
            }
            js << "}" << (di + 1 < ds_order.size() ? ",\n" : "\n");
        }
        js << "    }" << (p + 1 < n_pairs ? ",\n" : "\n");
    }
    js << "  },\n";

    // curves
    const std::vector<int> sub_curve_ids = {0, 1, 2, 8, 4}; // subset datasets
    js << "  \"curves\": {\n";
    for (int v = 0; v < n_ver; ++v) {
        const auto& vn = ver_names[v];
        js << "    \"" << vn << "\": {\n";
        for (size_t di = 0; di < ds_order.size(); ++di) {
            const auto& dn = ds_order[di];
            const auto& vs = all_vstats[vn][dn];
            bool full_ds = (dn == "All3Data" || dn == "10004X");
            std::vector<int> cids;
            if (full_ds) { for (int c = 0; c < N_CURVES; ++c) cids.push_back(c); }
            else cids = sub_curve_ids;
            js << "      \"" << dn << "\": {\n";
            for (size_t ci = 0; ci < cids.size(); ++ci) {
                int c = cids[ci];
                const auto& sp = CURVE_SPECS[c];
                js << "        \"" << sp.key << "\": {\n";
                for (int m = 0; m < N_METRICS; ++m) {
                    js << "          \"" << met_keys[m] << "\": ";
                    writeEffCurveJson(js, sp.key, sp.xlabel, "mm", vs.curves[c], m);
                    js << (m + 1 < N_METRICS ? ",\n" : "\n");
                }
                js << "        }" << (ci + 1 < cids.size() ? ",\n" : "\n");
            }
            js << "      }" << (di + 1 < ds_order.size() ? ",\n" : "\n");
        }
        js << "    }" << (v + 1 < n_ver ? ",\n" : "\n");
    }
    js << "  },\n";

    // mc_data_ratios
    js << "  \"mc_data_ratios\": {\n";
    for (int v = 0; v < n_ver; ++v) {
        const auto& vn = ver_names[v];
        const auto& vmc = all_vstats[vn]["10004X"];
        const auto& vda = all_vstats[vn]["All3Data"];
        js << "    \"" << vn << "\": {\n";
        for (int c = 0; c < N_CURVES; ++c) {
            const auto& sp = CURVE_SPECS[c];
            js << "      \"" << sp.key << "\": {\"eff_pawan_c2R\": ";
            writeRatioCurveJson(js, vmc.curves[c], vda.curves[c], MET_PAWAN_C2R);
            js << ", \"eff_lt2\": ";
            writeRatioCurveJson(js, vmc.curves[c], vda.curves[c], MET_LT2);
            js << "}" << (c + 1 < N_CURVES ? ",\n" : "\n");
        }
        js << "    }" << (v + 1 < n_ver ? ",\n" : "\n");
    }
    js << "  },\n";

    // paired_curves
    js << "  \"paired_curves\": {\n";
    for (int p = 0; p < n_pairs; ++p) {
        const auto& pn = pair_names[p];
        js << "    \"" << pn << "\": {\n";
        std::vector<std::string> pds = {"All3Data", "10004X"};
        for (size_t di = 0; di < pds.size(); ++di) {
            const auto& dn = pds[di];
            const auto& ps = all_pstats[pn][dn];
            js << "      \"" << dn << "\": {\n";
            for (int c = 0; c < N_CURVES; ++c) {
                const auto& sp = CURVE_SPECS[c];
                js << "        \"" << sp.key << "\": {\"eff_pawan_c2R\": ";
                writePairedCurveJson(js, pn, pair_ab[p].first, pair_ab[p].second, dn, "eff_pawan_c2R", ps.curves[c], 0);
                js << ", \"eff_lt2\": ";
                writePairedCurveJson(js, pn, pair_ab[p].first, pair_ab[p].second, dn, "eff_lt2", ps.curves[c], 1);
                js << "}" << (c + 1 < N_CURVES ? ",\n" : "\n");
            }
            js << "      }" << (di + 1 < pds.size() ? ",\n" : "\n");
        }
        js << "    }" << (p + 1 < n_pairs ? ",\n" : "\n");
    }
    js << "  },\n";

    // multiplicity_curves
    js << "  \"multiplicity_curves\": {\n";
    for (int v = 0; v < n_ver; ++v) {
        const auto& vn = ver_names[v];
        js << "    \"" << vn << "\": {\n";
        std::vector<std::string> mds = {"All3Data", "10004X"};
        for (size_t di = 0; di < mds.size(); ++di) {
            const auto& dn = mds[di];
            const auto& vs = all_vstats[vn][dn];
            js << "      \"" << dn << "\": {\n";
            for (int mi = 0; mi < N_MULT_CURVES; ++mi) {
                int cid = MULT_CURVE_IDS[mi];
                const auto& sp = CURVE_SPECS[cid];
                const auto& mb = vs.mult_curves[mi];
                double bw = (sp.xmax - sp.xmin) / sp.nbins;
                js << "        \"" << sp.key << "\": {\"edges\": [";
                for (int i = 0; i <= sp.nbins; ++i) { if (i > 0) js << ", "; js << fmtDoubleOrNull(sp.xmin + i * bw); }
                js << "], \"centers\": [";
                for (int i = 0; i < sp.nbins; ++i) { if (i > 0) js << ", "; js << fmtDoubleOrNull(sp.xmin + (i + 0.5) * bw); }
                js << "], \"half_widths\": [";
                for (int i = 0; i < sp.nbins; ++i) { if (i > 0) js << ", "; js << fmtDoubleOrNull(0.5 * bw); }
                js << "], \"den\": [";
                for (int i = 0; i < sp.nbins; ++i) { if (i > 0) js << ", "; js << mb.den[i]; }
                auto writeProb = [&](const char* k, const std::vector<long long>& num, bool comma) {
                    js << "], \"" << k << "\": [";
                    for (int i = 0; i < sp.nbins; ++i) {
                        if (i > 0) js << ", ";
                        double pval = (mb.den[i] > 0) ? static_cast<double>(num[i]) / mb.den[i]
                                                      : std::numeric_limits<double>::quiet_NaN();
                        js << fmtDoubleOrNull(pval);
                    }
                    if (!comma) js << "]";
                };
                writeProb("p_lt0", mb.c_lt0, true);
                writeProb("p_lt1", mb.c_lt1, true);
                writeProb("p_lt2_pass_c2R", mb.c_lt2_pass, true);
                writeProb("p_lt2_fail_c2R", mb.c_lt2_fail, true);
                writeProb("p_lt_ge3", mb.c_ge3, false);
                js << "}" << (mi + 1 < N_MULT_CURVES ? ",\n" : "\n");
            }
            js << "      }" << (di + 1 < mds.size() ? ",\n" : "\n");
        }
        js << "    }" << (v + 1 < n_ver ? ",\n" : "\n");
    }
    js << "  },\n";

    // momentum
    std::vector<double> pz_edges(N_PZ_BINS + 1), pz_centers(N_PZ_BINS);
    for (int i = 0; i <= N_PZ_BINS; ++i) pz_edges[i] = std::pow(10.0, PZ_LOG_MIN + (PZ_LOG_MAX - PZ_LOG_MIN) * i / N_PZ_BINS);
    for (int i = 0; i < N_PZ_BINS; ++i) pz_centers[i] = std::sqrt(pz_edges[i] * pz_edges[i + 1]);

    js << "  \"momentum\": {\n";
    for (int v = 0; v < n_ver; ++v) {
        const auto& vn = ver_names[v];
        js << "    \"" << vn << "\": {\n";
        std::vector<std::string> mds = {"All3Data", "10004X"};
        for (size_t di = 0; di < mds.size(); ++di) {
            const auto& dn = mds[di];
            const auto& vs = all_vstats[vn][dn];
            auto writePzHist = [&](const std::array<long long, N_PZ_BINS>& cnts, long long tot, long long uf, long long of) {
                js << "{\"edges\": [";
                for (int i = 0; i <= N_PZ_BINS; ++i) { if (i > 0) js << ", "; js << fmtDoubleOrNull(pz_edges[i]); }
                js << "], \"centers\": [";
                for (int i = 0; i < N_PZ_BINS; ++i) { if (i > 0) js << ", "; js << fmtDoubleOrNull(pz_centers[i]); }
                js << "], \"counts\": [";
                for (int i = 0; i < N_PZ_BINS; ++i) { if (i > 0) js << ", "; js << cnts[i]; }
                js << "], \"norm\": [";
                for (int i = 0; i < N_PZ_BINS; ++i) {
                    if (i > 0) js << ", ";
                    js << fmtDoubleOrNull(tot > 0 ? static_cast<double>(cnts[i]) / tot : 0.0);
                }
                js << "], \"norm_err\": [";
                for (int i = 0; i < N_PZ_BINS; ++i) {
                    if (i > 0) js << ", ";
                    js << fmtDoubleOrNull(tot > 0 ? std::sqrt(static_cast<double>(cnts[i])) / tot : 0.0);
                }
                js << "], \"total_tracks\": " << tot << ", \"underflow\": " << uf << ", \"overflow\": " << of << "}";
            };
            // Residual stats (mean, std, median, iqr68 from histogram)
            double rmean = (vs.pz_res_tot > 0) ? vs.pz_res_sum / vs.pz_res_tot : 0.0;
            double rstd = (vs.pz_res_tot > 0) ? std::sqrt(std::max(0.0, vs.pz_res_sum2 / vs.pz_res_tot - rmean * rmean)) : 0.0;
            auto histQuantile = [&](double q) {
                if (vs.pz_res_tot <= 0) return 0.0;
                double target = q * vs.pz_res_tot, cum = 0.0;
                double bw = (RES_MAX - RES_MIN) / N_RES_BINS;
                for (int i = 0; i < N_RES_BINS; ++i) {
                    double prev = cum;
                    cum += vs.pz_res[i];
                    if (cum >= target && vs.pz_res[i] > 0) {
                        double frac = (target - prev) / vs.pz_res[i];
                        return RES_MIN + (i + frac) * bw;
                    }
                }
                return RES_MAX;
            };
            double rmed = histQuantile(0.5);
            double riqr68 = 0.5 * (histQuantile(0.8413447) - histQuantile(0.1586553));

            js << "      \"" << dn << "\": {\n"
               << "        \"track_pz0\": "; writePzHist(vs.pz_reco, vs.pz_reco_tot, vs.pz_reco_uf, vs.pz_reco_of); js << ",\n"
               << "        \"og_track_pz0\": "; writePzHist(vs.pz_og, vs.pz_og_tot, vs.pz_og_uf, vs.pz_og_of); js << ",\n"
               << "        \"pz0_rel_res\": {\"edges\": [";
            double rbw = (RES_MAX - RES_MIN) / N_RES_BINS;
            for (int i = 0; i <= N_RES_BINS; ++i) { if (i > 0) js << ", "; js << fmtDoubleOrNull(RES_MIN + i * rbw); }
            js << "], \"centers\": [";
            for (int i = 0; i < N_RES_BINS; ++i) { if (i > 0) js << ", "; js << fmtDoubleOrNull(RES_MIN + (i + 0.5) * rbw); }
            js << "], \"counts\": [";
            for (int i = 0; i < N_RES_BINS; ++i) { if (i > 0) js << ", "; js << vs.pz_res[i]; }
            js << "], \"norm\": [";
            for (int i = 0; i < N_RES_BINS; ++i) {
                if (i > 0) js << ", ";
                js << fmtDoubleOrNull(vs.pz_res_tot > 0 ? static_cast<double>(vs.pz_res[i]) / vs.pz_res_tot : 0.0);
            }
            js << "], \"total_matched\": " << vs.pz_res_tot
               << ", \"mean\": " << fmtDoubleOrNull(rmean)
               << ", \"std\": " << fmtDoubleOrNull(rstd)
               << ", \"median\": " << fmtDoubleOrNull(rmed)
               << ", \"iqr68\": " << fmtDoubleOrNull(riqr68)
               << ", \"underflow\": " << vs.pz_res_uf
               << ", \"overflow\": " << vs.pz_res_of << "}\n"
               << "      }" << (di + 1 < mds.size() ? ",\n" : "\n");
        }
        js << "    },\n";
    }
    js << "    \"mc_data_pz0_ratio\": {\n";
    for (int v = 0; v < n_ver; ++v) {
        const auto& vn = ver_names[v];
        const auto& vmc = all_vstats[vn]["10004X"];
        const auto& vda = all_vstats[vn]["All3Data"];
        js << "      \"" << vn << "\": {\"edges\": [";
        for (int i = 0; i <= N_PZ_BINS; ++i) { if (i > 0) js << ", "; js << fmtDoubleOrNull(pz_edges[i]); }
        js << "], \"centers\": [";
        for (int i = 0; i < N_PZ_BINS; ++i) { if (i > 0) js << ", "; js << fmtDoubleOrNull(pz_centers[i]); }
        js << "], \"ratio\": [";
        std::vector<double> rerr(N_PZ_BINS);
        for (int i = 0; i < N_PZ_BINS; ++i) {
            if (i > 0) js << ", ";
            double fm = (vmc.pz_reco_tot > 0) ? static_cast<double>(vmc.pz_reco[i]) / vmc.pz_reco_tot : 0.0;
            double fd = (vda.pz_reco_tot > 0) ? static_cast<double>(vda.pz_reco[i]) / vda.pz_reco_tot : 0.0;
            if (fm > 0.0 && fd > 0.0 && vmc.pz_reco[i] > 0 && vda.pz_reco[i] > 0) {
                double rat = fm / fd;
                rerr[i] = rat * std::sqrt(1.0 / vmc.pz_reco[i] + 1.0 / vda.pz_reco[i]);
                js << fmtDoubleOrNull(rat);
            } else {
                rerr[i] = std::numeric_limits<double>::quiet_NaN();
                js << "null";
            }
        }
        js << "], \"err\": [";
        for (int i = 0; i < N_PZ_BINS; ++i) { if (i > 0) js << ", "; js << fmtDoubleOrNull(rerr[i]); }
        js << "]}" << (v + 1 < n_ver ? ",\n" : "\n");
    }
    js << "    }\n  },\n";

    // stage_of_loss
    js << "  \"stage_of_loss\": {\n";
    for (int v = 0; v < n_ver; ++v) {
        const auto& vn = ver_names[v];
        js << "    \"" << vn << "\": {\n";
        std::vector<std::string> sds = {"All3Data", "10004X"};
        for (size_t di = 0; di < sds.size(); ++di) {
            const auto& dn = sds[di];
            const auto& vs = all_vstats[vn][dn];
            js << "      \"" << dn << "\": {\n        \"regimes\": {\n";
            for (int r = 0; r < N_REGIMES; ++r) {
                js << "          \"" << reg_keys[r] << "\": ";
                writeStageCountsJson(js, vs.stage_reg[r]);
                js << (r + 1 < N_REGIMES ? ",\n" : "\n");
            }
            js << "        },\n";
            js << "        \"plateau_vs_nSeg\": ";
            writeEffCurveJson(js, "plateau_vs_nSeg", "Total TrackSegments", "count", vs.plat_nSeg, 0); js << ",\n";
            js << "        \"plateau_vs_nGBSeg\": ";
            writeEffCurveJson(js, "plateau_vs_nGBSeg", "GhostBusted TrackSegments", "count", vs.plat_nGBSeg, 0); js << ",\n";
            js << "        \"plateau_vs_min_og_p0\": ";
            writeEffCurveJson(js, "plateau_vs_min_og_p0", "min(og_p0) [GeV]", "GeV", vs.plat_min_p0, 0); js << ",\n";
            js << "        \"plateau_vs_max_og_R0\": ";
            writeEffCurveJson(js, "plateau_vs_max_og_R0", "max(og_R0) [mm]", "mm", vs.plat_max_R0, 0); js << ",\n";
            js << "        \"abs_dX0_stage_fractions\": {\"edges\": [";
            for (int i = 0; i <= N_DX_STRAT_BINS; ++i) { if (i > 0) js << ", "; js << (5.0 * i); }
            js << "], \"centers\": [";
            for (int i = 0; i < N_DX_STRAT_BINS; ++i) { if (i > 0) js << ", "; js << (5.0 * i + 2.5); }
            js << "], \"den\": [";
            for (int i = 0; i < N_DX_STRAT_BINS; ++i) { if (i > 0) js << ", "; js << vs.abs_dX0_stage[i].n_den; }
            auto writeSFrac = [&](const char* k, auto fn, bool comma) {
                js << "], \"" << k << "\": [";
                for (int i = 0; i < N_DX_STRAT_BINS; ++i) {
                    if (i > 0) js << ", ";
                    long long den = vs.abs_dX0_stage[i].n_den;
                    double f = (den > 0) ? static_cast<double>(fn(vs.abs_dX0_stage[i])) / den
                                         : std::numeric_limits<double>::quiet_NaN();
                    js << fmtDoubleOrNull(f);
                }
                if (!comma) js << "]";
            };
            writeSFrac("frac_pass_c2R", [](const StageLossCounts& s){ return s.pass_c2R; }, true);
            writeSFrac("frac_lt1_miss_seg", [](const StageLossCounts& s){ return s.lt1_miss_seg_st1 + s.lt1_miss_seg_st2 + s.lt1_miss_seg_st3 + s.lt1_miss_seg_multi; }, true);
            writeSFrac("frac_lt1_killed_gb", [](const StageLossCounts& s){ return s.lt1_killed_by_gb; }, true);
            writeSFrac("frac_lt1_lost_ckf", [](const StageLossCounts& s){ return s.lt1_all3gb_lost_in_ckf + s.lt1_unmatched_kept; }, true);
            writeSFrac("frac_lt0", [](const StageLossCounts& s){ return s.lt0; }, true);
            writeSFrac("frac_lt_ge3", [](const StageLossCounts& s){ return s.lt_ge3; }, true);
            writeSFrac("frac_lt2_fail_c2R", [](const StageLossCounts& s){ return s.lt2_fail_c2R; }, false);
            js << "}\n      }" << (di + 1 < sds.size() ? ",\n" : "\n");
        }
        js << "    }" << (v + 1 < n_ver ? ",\n" : "\n");
    }
    js << "  },\n";

    // large_dX_diagnostics
    js << "  \"large_dX_diagnostics\": {\n";
    for (int v = 0; v < n_ver; ++v) {
        const auto& vn = ver_names[v];
        js << "    \"" << vn << "\": {\n";
        std::vector<std::string> lds = {"All3Data", "10004X"};
        for (size_t di = 0; di < lds.size(); ++di) {
            const auto& dn = lds[di];
            const auto& vs = all_vstats[vn][dn];
            js << "      \"" << dn << "\": {\n"
               << "        \"eff_abs_dX0_by_max_abs_og_x0\": {\n"
               << "          \"0_33\": "; writeEffCurveJson(js, "abs_dX0_x0_0_33", "|\\Delta X_0| [mm]", "mm", vs.eff_dx0_by_x0[0], 0); js << ",\n"
               << "          \"33_66\": "; writeEffCurveJson(js, "abs_dX0_x0_33_66", "|\\Delta X_0| [mm]", "mm", vs.eff_dx0_by_x0[1], 0); js << ",\n"
               << "          \"66_95\": "; writeEffCurveJson(js, "abs_dX0_x0_66_95", "|\\Delta X_0| [mm]", "mm", vs.eff_dx0_by_x0[2], 0); js << "\n"
               << "        },\n"
               << "        \"eff_abs_dX0_by_max_og_R0\": {\n"
               << "          \"0_60\": "; writeEffCurveJson(js, "abs_dX0_R0_0_60", "|\\Delta X_0| [mm]", "mm", vs.eff_dx0_by_R0[0], 0); js << ",\n"
               << "          \"60_80\": "; writeEffCurveJson(js, "abs_dX0_R0_60_80", "|\\Delta X_0| [mm]", "mm", vs.eff_dx0_by_R0[1], 0); js << ",\n"
               << "          \"80_95\": "; writeEffCurveJson(js, "abs_dX0_R0_80_95", "|\\Delta X_0| [mm]", "mm", vs.eff_dx0_by_R0[2], 0); js << "\n"
               << "        },\n"
               << "        \"loss_map_st1\": "; writeMap2DJson(js, vs.loss_map_st1); js << ",\n"
               << "        \"loss_map_st3\": "; writeMap2DJson(js, vs.loss_map_st3); js << "\n"
               << "      }" << (di + 1 < lds.size() ? ",\n" : "\n");
        }
        js << "    }" << (v + 1 < n_ver ? ",\n" : "\n");
    }
    js << "  },\n";

    // momentum_reweighting
    js << "  \"momentum_reweighting\": {\n"
       << "    \"og_p0_bins_gev\": {\"edges\": [";
    std::vector<double> rpe(N_RW_P0_BINS + 1), rpc(N_RW_P0_BINS);
    for (int i = 0; i <= N_RW_P0_BINS; ++i) rpe[i] = std::pow(10.0, RW_P0_LOG_MIN + (RW_P0_LOG_MAX - RW_P0_LOG_MIN) * i / N_RW_P0_BINS);
    for (int i = 0; i < N_RW_P0_BINS; ++i) rpc[i] = std::sqrt(rpe[i] * rpe[i + 1]);
    for (int i = 0; i <= N_RW_P0_BINS; ++i) { if (i > 0) js << ", "; js << fmtDoubleOrNull(rpe[i]); }
    js << "], \"centers\": [";
    for (int i = 0; i < N_RW_P0_BINS; ++i) { if (i > 0) js << ", "; js << fmtDoubleOrNull(rpc[i]); }
    js << "], \"data_frac\": [";
    for (int i = 0; i < N_RW_P0_BINS; ++i) { if (i > 0) js << ", "; js << fmtDoubleOrNull(data_p0_frac[i]); }
    js << "], \"mc_frac\": [";
    for (int i = 0; i < N_RW_P0_BINS; ++i) { if (i > 0) js << ", "; js << fmtDoubleOrNull(mc_p0_frac[i]); }
    js << "], \"weights\": [";
    for (int i = 0; i < N_RW_P0_BINS; ++i) { if (i > 0) js << ", "; js << fmtDoubleOrNull(rw_weights[i]); }
    js << "]},\n";
    for (int v = 0; v < n_ver; ++v) {
        const auto& vn = ver_names[v];
        const auto& rw = rw_merged[v];
        const auto& vda = all_vstats[vn]["All3Data"];
        double neff = (rw.sumw2 > 0.0) ? (rw.sumw * rw.sumw) / rw.sumw2 : 0.0;
        js << "    \"" << vn << "\": {\n"
           << "      \"n_events\": " << rw.n_events << ", \"sumw\": " << fmtDoubleOrNull(rw.sumw)
           << ", \"sumw2\": " << fmtDoubleOrNull(rw.sumw2) << ", \"n_eff\": " << fmtDoubleOrNull(neff) << ",\n"
           << "      \"closure_mc_reweighted_frac\": [";
        for (int i = 0; i < N_RW_P0_BINS; ++i) {
            if (i > 0) js << ", ";
            js << fmtDoubleOrNull(rw.closure_tot > 0.0 ? rw.closure_sumw[i] / rw.closure_tot : 0.0);
        }
        js << "],\n      \"curves\": {\n";
        for (int rc = 0; rc < N_RW_CURVES; ++rc) {
            int cid = RW_CURVE_IDS[rc];
            const auto& sp = CURVE_SPECS[cid];
            js << "        \"" << sp.key << "\": ";
            writeWeightedCurveJson(js, std::string(sp.key) + "_mc_reweighted", rw.curves[rc]);
            js << (rc + 1 < N_RW_CURVES ? ",\n" : "\n");
        }
        js << "      },\n      \"mc_data_ratios_reweighted\": {\n";
        for (int rc = 0; rc < N_RW_CURVES; ++rc) {
            int cid = RW_CURVE_IDS[rc];
            const auto& sp = CURVE_SPECS[cid];
            js << "        \"" << sp.key << "\": ";
            writeWeightedRatioCurveJson(js, rw.curves[rc], vda.curves[cid]);
            js << (rc + 1 < N_RW_CURVES ? ",\n" : "\n");
        }
        js << "      }\n    }" << (v + 1 < n_ver ? ",\n" : "\n");
    }
    js << "  }";

    if (has_aux) {
        js << ",\n  \"auxiliary\": {\n    \"010698_PosBin2_NegBin3\": {\n"
           << "      \"file\": \"" << aux_file << "\", \"n_events\": " << aux_stats.n_events
           << ", \"eff_pawan_c2R\": ";
        writeEffScalarJson(js, aux_stats.reg_pass[REG_GLOBAL][MET_PAWAN_C2R], aux_stats.reg_den[REG_GLOBAL]);
        js << ", \"eff_lt2\": ";
        writeEffScalarJson(js, aux_stats.reg_pass[REG_GLOBAL][MET_LT2], aux_stats.reg_den[REG_GLOBAL]);
        js << "\n    }\n  }\n";
    } else {
        js << "\n";
    }
    js << "}\n";
    js.close();

    // Concise stdout summary
    std::cout << "[overlay_analyze] DONE in " << std::fixed << std::setprecision(2) << elapsed_s
              << " s | og_checks_passed=" << merged.og_checks_passed
              << " failed=" << merged.og_checks_failed << "\n";
    for (int v = 0; v < n_ver; ++v) {
        const auto& vn = ver_names[v];
        for (const std::string dn : {"All3Data", "10004X"}) {
            const auto& s = all_vstats[vn][dn];
            double ec, lo, hi, ep, lpp, hpp;
            clopperPearson1Sigma(s.reg_pass[REG_GLOBAL][MET_PAWAN_C2R], s.reg_den[REG_GLOBAL], ec, lo, hi);
            clopperPearson1Sigma(s.reg_pass[REG_PLATEAU][MET_PAWAN_C2R], s.reg_den[REG_PLATEAU], ep, lpp, hpp);
            const auto& sr = s.stage_reg[REG_PLATEAU];
            double nd = static_cast<double>(sr.n_den);
            std::cout << "  " << std::setw(15) << vn << " | " << std::setw(8) << dn
                      << " | global_c2R=" << std::setprecision(4) << ec
                      << " | plateau_c2R=" << ep
                      << " | plateau_loss: miss_seg=" << std::setprecision(4)
                      << (sr.lt1_miss_seg_st1 + sr.lt1_miss_seg_st2 + sr.lt1_miss_seg_st3 + sr.lt1_miss_seg_multi) / nd
                      << " (st1=" << sr.lt1_miss_seg_st1 / nd << ",st2=" << sr.lt1_miss_seg_st2 / nd
                      << ",st3=" << sr.lt1_miss_seg_st3 / nd << ",multi=" << sr.lt1_miss_seg_multi / nd
                      << ") gb_kill=" << sr.lt1_killed_by_gb / nd
                      << " ckf_lost=" << sr.lt1_all3gb_lost_in_ckf / nd
                      << " unm=" << sr.lt1_unmatched_kept / nd
                      << " lt0=" << sr.lt0 / nd
                      << " ge3=" << sr.lt_ge3 / nd
                      << " fail_c2R=" << sr.lt2_fail_c2R / nd << "\n";
        }
    }
    return 0;
}
