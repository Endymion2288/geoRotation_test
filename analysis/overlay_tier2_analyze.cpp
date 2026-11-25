// Compiled C++17 16-thread OpenMP 4-way Tier-2 Attribution Analyzer for Stage 08.
// Compares 4 arms on the exact 36 matched quadruples (18 008749 + 18 10004X, ~0.95M events x 4):
//   0: baseline_jul8   (Unrotated EOS OverlaidNTPs)
//   1: t2b_7fd40df8    (Pre-fix GeoRotation 7fd40df8)
//   2: latest_218813c2 (Latest GeoRotation fix_v1 218813c2, covPhi=0.1089)
//   3: t2a_covphi      (Latest GeoRotation fix_v1 218813c2 + covPhi=8.41e-6)

#include <Math/QuantFuncMathCore.h>
#include <TFile.h>
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
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

static constexpr double CL_1SIGMA = 0.6826894921370859;
static constexpr double ALPHA_1SIGMA = 1.0 - CL_1SIGMA;

enum ArmIdx { ARM_BASE = 0, ARM_T2B = 1, ARM_LATEST = 2, ARM_T2A = 3, N_ARMS = 4 };
enum PairIdx {
    PAIR_T2B_BASE = 0,
    PAIR_LATEST_BASE = 1,
    PAIR_T2A_BASE = 2,
    PAIR_T2A_LATEST = 3,
    N_PAIRS = 4
};
enum MetricIdx { MET_C2R = 0, MET_LT2 = 1, MET_GE2 = 2, MET_MATCHED = 3, N_METRICS = 4 };
enum SplitIdx {
    SP_8749_DIAG = 0,
    SP_8749_OFF = 1,
    SP_MC_DIAG = 2,
    SP_MC_OFF = 3,
    N_RAW_SPLITS = 4
};
enum RegimeIdx {
    REG_GLOBAL = 0,
    REG_CRATER = 1,
    REG_SMALL_Y0 = 2,
    REG_LARGE_Y0 = 3,
    REG_SMALL_DR1 = 4,
    REG_PLATEAU = 5,
    REG_LARGE_DX0 = 6,
    N_REGIMES = 7
};

struct CurveSpec {
    const char* key;
    const char* xlabel;
    double xmin;
    double xmax;
    int nbins;
};

static constexpr int N_CURVES = 5;
static constexpr CurveSpec CURVE_SPECS[N_CURVES] = {
    {"abs_dY0_zoom",    "Track Separation (|\\Delta Y_0|) [mm]",    0.0,  20.0,  40},
    {"dR1_zoom",        "Track Separation (\\Delta R_1) [mm]",      0.0,  20.0,  40},
    {"dR1_wide",        "Track Separation (\\Delta R_1) [mm]",      0.0, 100.0, 200},
    {"abs_dX0_wide",    "Track Separation (|\\Delta X_0|) [mm]",    0.0, 100.0, 100},
    {"dY1_zoom_signed", "Track Separation (\\Delta Y_1) [mm]",    -25.0,  25.0,  50},
};

static constexpr int N_RES_BINS = 100;
static constexpr double RES_MIN = -0.5;
static constexpr double RES_MAX = 0.5;

inline void clopperPearson1Sigma(long long k, long long n, double& eff, double& lo, double& hi) {
    if (n <= 0) {
        eff = std::numeric_limits<double>::quiet_NaN();
        lo = 0.0; hi = 0.0;
        return;
    }
    eff = static_cast<double>(k) / static_cast<double>(n);
    double l = (k == 0) ? 0.0 : ROOT::Math::beta_quantile(0.5 * ALPHA_1SIGMA, static_cast<double>(k), static_cast<double>(n - k + 1));
    double h = (k == n) ? 1.0 : ROOT::Math::beta_quantile(1.0 - 0.5 * ALPHA_1SIGMA, static_cast<double>(k + 1), static_cast<double>(n - k));
    lo = std::max(0.0, eff - l);
    hi = std::max(0.0, h - eff);
}

inline void mcnemarDelta(long long bp, long long g, long long l, long long bf, double& d, double& e) {
    long long n = bp + g + l + bf;
    if (n <= 0) {
        d = std::numeric_limits<double>::quiet_NaN();
        e = std::numeric_limits<double>::quiet_NaN();
        return;
    }
    double dn = static_cast<double>(n);
    double diff = static_cast<double>(g - l);
    d = diff / dn;
    e = std::sqrt(std::max(0.0, static_cast<double>(g + l) - (diff * diff) / dn)) / dn;
}

inline std::string fmtD(double v) {
    if (!std::isfinite(v)) return "null";
    std::ostringstream oss;
    oss << std::setprecision(10) << v;
    return oss.str();
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

struct PairedCell {
    long long bp = 0, g = 0, l = 0, bf = 0;
    void add(bool a, bool b) {
        if (a && b) bp++;
        else if (!a && b) g++;
        else if (a && !b) l++;
        else bf++;
    }
    void merge(const PairedCell& o) { bp += o.bp; g += o.g; l += o.l; bf += o.bf; }
};

struct StageLossCounts {
    long long n_den = 0, pass_c2R = 0, lt2_fail_c2R = 0, lt_ge3 = 0, lt0 = 0, lt1_total = 0;
    long long lt1_miss_seg_st1 = 0, lt1_miss_seg_st2 = 0, lt1_miss_seg_st3 = 0, lt1_miss_seg_multi = 0;
    long long lt1_killed_by_gb = 0, lt1_all3gb_lost_in_ckf = 0, lt1_unmatched_kept = 0;
    void merge(const StageLossCounts& o) {
        n_den += o.n_den; pass_c2R += o.pass_c2R; lt2_fail_c2R += o.lt2_fail_c2R;
        lt_ge3 += o.lt_ge3; lt0 += o.lt0; lt1_total += o.lt1_total;
        lt1_miss_seg_st1 += o.lt1_miss_seg_st1; lt1_miss_seg_st2 += o.lt1_miss_seg_st2;
        lt1_miss_seg_st3 += o.lt1_miss_seg_st3; lt1_miss_seg_multi += o.lt1_miss_seg_multi;
        lt1_killed_by_gb += o.lt1_killed_by_gb; lt1_all3gb_lost_in_ckf += o.lt1_all3gb_lost_in_ckf;
        lt1_unmatched_kept += o.lt1_unmatched_kept;
    }
};

struct EffBins2M {
    int nbins = 0;
    double xmin = 0.0, xmax = 1.0;
    long long underflow = 0, overflow = 0, invalid = 0;
    std::vector<long long> den, pass_c2r, pass_ge2;
    EffBins2M() = default;
    EffBins2M(int nb, double lo, double hi)
        : nbins(nb), xmin(lo), xmax(hi), den(nb, 0), pass_c2r(nb, 0), pass_ge2(nb, 0) {}
    void fill(double x, bool p_c2r, bool p_ge2) {
        if (!std::isfinite(x)) { invalid++; return; }
        if (x < xmin) { underflow++; return; }
        if (x >= xmax) { overflow++; return; }
        int b = std::clamp(static_cast<int>((x - xmin) / (xmax - xmin) * nbins), 0, nbins - 1);
        den[b]++;
        if (p_c2r) pass_c2r[b]++;
        if (p_ge2) pass_ge2[b]++;
    }
    void merge(const EffBins2M& o) {
        underflow += o.underflow; overflow += o.overflow; invalid += o.invalid;
        for (int i = 0; i < nbins; ++i) {
            den[i] += o.den[i]; pass_c2r[i] += o.pass_c2r[i]; pass_ge2[i] += o.pass_ge2[i];
        }
    }
};

struct PairedBins2M {
    int nbins = 0;
    double xmin = 0.0, xmax = 1.0;
    std::vector<PairedCell> c_c2r, c_ge2;
    PairedBins2M() = default;
    PairedBins2M(int nb, double lo, double hi) : nbins(nb), xmin(lo), xmax(hi), c_c2r(nb), c_ge2(nb) {}
    void fill(double x, bool a0, bool b0, bool a2, bool b2) {
        if (!std::isfinite(x) || x < xmin || x >= xmax) return;
        int b = std::clamp(static_cast<int>((x - xmin) / (xmax - xmin) * nbins), 0, nbins - 1);
        c_c2r[b].add(a0, b0);
        c_ge2[b].add(a2, b2);
    }
    void merge(const PairedBins2M& o) {
        for (int i = 0; i < nbins; ++i) { c_c2r[i].merge(o.c_c2r[i]); c_ge2[i].merge(o.c_ge2[i]); }
    }
};

struct ArmSplitStats {
    int n_files = 0;
    long long n_events = 0;
    std::array<long long, 6> lt_counts{};
    long long lt2_fail_c2R = 0;
    std::array<long long, N_REGIMES> reg_den{};
    std::array<std::array<long long, N_METRICS>, N_REGIMES> reg_pass{};
    std::array<StageLossCounts, N_REGIMES> stage_reg{};
    std::array<EffBins2M, N_CURVES> curves;
    std::array<long long, N_RES_BINS> res_counts{};
    long long res_tot = 0, res_uf = 0, res_of = 0;
    double res_sum = 0.0, res_sum2 = 0.0;

    ArmSplitStats() {
        for (int c = 0; c < N_CURVES; ++c)
            curves[c] = EffBins2M(CURVE_SPECS[c].nbins, CURVE_SPECS[c].xmin, CURVE_SPECS[c].xmax);
    }

    void merge(const ArmSplitStats& o) {
        n_files += o.n_files; n_events += o.n_events; lt2_fail_c2R += o.lt2_fail_c2R;
        for (int i = 0; i < 6; ++i) lt_counts[i] += o.lt_counts[i];
        for (int r = 0; r < N_REGIMES; ++r) {
            reg_den[r] += o.reg_den[r];
            for (int m = 0; m < N_METRICS; ++m) reg_pass[r][m] += o.reg_pass[r][m];
            stage_reg[r].merge(o.stage_reg[r]);
        }
        for (int c = 0; c < N_CURVES; ++c) curves[c].merge(o.curves[c]);
        for (int i = 0; i < N_RES_BINS; ++i) res_counts[i] += o.res_counts[i];
        res_tot += o.res_tot; res_uf += o.res_uf; res_of += o.res_of;
        res_sum += o.res_sum; res_sum2 += o.res_sum2;
    }
};

struct PairSplitStats {
    std::array<std::array<PairedCell, N_METRICS>, N_REGIMES> reg_cells{};
    std::array<PairedBins2M, N_CURVES> curves;
    PairSplitStats() {
        for (int c = 0; c < N_CURVES; ++c)
            curves[c] = PairedBins2M(CURVE_SPECS[c].nbins, CURVE_SPECS[c].xmin, CURVE_SPECS[c].xmax);
    }
    void merge(const PairSplitStats& o) {
        for (int r = 0; r < N_REGIMES; ++r)
            for (int m = 0; m < N_METRICS; ++m) reg_cells[r][m].merge(o.reg_cells[r][m]);
        for (int c = 0; c < N_CURVES; ++c) curves[c].merge(o.curves[c]);
    }
};

struct ThreadAccum {
    std::array<std::array<ArmSplitStats, N_RAW_SPLITS>, N_ARMS> arms;
    std::array<std::array<PairSplitStats, N_RAW_SPLITS>, N_PAIRS> pairs;
    long long og_passed = 0, og_failed = 0;
    void merge(const ThreadAccum& o) {
        og_passed += o.og_passed; og_failed += o.og_failed;
        for (int a = 0; a < N_ARMS; ++a)
            for (int s = 0; s < N_RAW_SPLITS; ++s) arms[a][s].merge(o.arms[a][s]);
        for (int p = 0; p < N_PAIRS; ++p)
            for (int s = 0; s < N_RAW_SPLITS; ++s) pairs[p][s].merge(o.pairs[p][s]);
    }
};

inline int stationOfZ(double z) {
    if (z >= -100.0 && z <= 350.0) return 0;
    if (z >= 950.0 && z <= 1550.0) return 1;
    if (z >= 2150.0 && z <= 2750.0) return 2;
    return -1;
}

inline void hermitePredict(double sz,
                           double z0, double x0, double y0, double tx0, double ty0,
                           double z1, double x1, double y1, double tx1, double ty1,
                           double& xp, double& yp) {
    double L = z1 - z0;
    if (std::abs(L) < 1e-6) { xp = x0; yp = y0; return; }
    double t = (sz - z0) / L, t2 = t * t, t3 = t2 * t;
    xp = (2*t3 - 3*t2 + 1)*x0 + (t3 - 2*t2 + t)*L*tx0 + (-2*t3 + 3*t2)*x1 + (t3 - t2)*L*tx1;
    yp = (2*t3 - 3*t2 + 1)*y0 + (t3 - 2*t2 + t)*L*ty0 + (-2*t3 + 3*t2)*y1 + (t3 - t2)*L*ty1;
}

struct EvalEv {
    int lt = 0;
    bool pass[N_METRICS] = {false, false, false, false};
    bool lt2_fail = false;
    int stage_code = 0;
    double res[2] = {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
};

inline EvalEv evaluateArm(
    int lt,
    const std::vector<double>& c2, const std::vector<double>& nd,
    const std::vector<double>& x0, const std::vector<double>& y0, const std::vector<double>& pz0,
    const std::vector<double>& x1, const std::vector<double>& y1,
    const std::vector<double>& ox0, const std::vector<double>& oy0, const std::vector<double>& oz0,
    const std::vector<double>& opx0, const std::vector<double>& opy0, const std::vector<double>& opz0,
    const std::vector<double>& ox1, const std::vector<double>& oy1, const std::vector<double>& oz1,
    const std::vector<double>& opx1, const std::vector<double>& opy1, const std::vector<double>& opz1,
    const std::vector<double>& sx, const std::vector<double>& sy, const std::vector<double>& sz,
    const std::vector<double>& gx, const std::vector<double>& gy, const std::vector<double>& gz) {

    EvalEv ev;
    ev.lt = lt;
    int ntr = static_cast<int>(x0.size());
    auto okC2R = [&](int j) {
        return (j >= 0 && j < ntr && nd[j] > 0.0 && (c2[j] / nd[j]) <= 25.0 &&
                std::hypot(x0[j], y0[j]) < 95.0 && std::hypot(x1[j], y1[j]) < 95.0);
    };
    ev.pass[MET_LT2] = (lt == 2);
    ev.pass[MET_GE2] = (lt >= 2);
    if (lt == 2 && ntr == 2 && okC2R(0) && okC2R(1)) ev.pass[MET_C2R] = true;
    else if (lt == 2) ev.lt2_fail = true;

    if (lt >= 2 && ntr >= 2) {
        double best = 1e30; int ba = -1, bb = -1;
        for (int a = 0; a < ntr; ++a) {
            if (!okC2R(a)) continue;
            double da0 = std::hypot(x0[a]-ox0[0], y0[a]-oy0[0]), da1 = std::hypot(x1[a]-ox1[0], y1[a]-oy1[0]);
            if (da0 >= 5.0 || da1 >= 5.0) continue;
            for (int b = 0; b < ntr; ++b) {
                if (b == a || !okC2R(b)) continue;
                double db0 = std::hypot(x0[b]-ox0[1], y0[b]-oy0[1]), db1 = std::hypot(x1[b]-ox1[1], y1[b]-oy1[1]);
                if (db0 >= 5.0 || db1 >= 5.0) continue;
                double s = da0 + da1 + db0 + db1;
                if (s < best) { best = s; ba = a; bb = b; }
            }
        }
        if (ba >= 0 && bb >= 0) {
            ev.pass[MET_MATCHED] = true;
            if (opz0[0] > 0.0) ev.res[0] = (pz0[ba] / opz0[0]) - 1.0;
            if (opz0[1] > 0.0) ev.res[1] = (pz0[bb] / opz0[1]) - 1.0;
        }
    }

    if (ev.pass[MET_C2R]) ev.stage_code = 0;
    else if (lt == 2) ev.stage_code = 1;
    else if (lt >= 3) ev.stage_code = 2;
    else if (lt <= 0) ev.stage_code = 3;
    else {
        int kept = -1;
        if (ntr >= 1) {
            double d0 = std::hypot(x0[0]-ox0[0], y0[0]-oy0[0]) + std::hypot(x1[0]-ox1[0], y1[0]-oy1[0]);
            double d1 = std::hypot(x0[0]-ox0[1], y0[0]-oy0[1]) + std::hypot(x1[0]-ox1[1], y1[0]-oy1[1]);
            if (std::min(d0, d1) < 15.0) kept = (d0 <= d1) ? 0 : 1;
        }
        if (kept < 0) ev.stage_code = 10;
        else {
            int lost = 1 - kept;
            double tx0 = std::abs(opz0[lost]) > 1e-3 ? opx0[lost]/opz0[lost] : 0.0;
            double ty0 = std::abs(opz0[lost]) > 1e-3 ? opy0[lost]/opz0[lost] : 0.0;
            double tx1 = std::abs(opz1[lost]) > 1e-3 ? opx1[lost]/opz1[lost] : 0.0;
            double ty1 = std::abs(opz1[lost]) > 1e-3 ? opy1[lost]/opz1[lost] : 0.0;
            bool sh[3] = {false, false, false};
            for (size_t k = 0; k < sz.size(); ++k) {
                int st = stationOfZ(sz[k]); if (st < 0 || sh[st]) continue;
                double xp = 0, yp = 0;
                hermitePredict(sz[k], oz0[lost], ox0[lost], oy0[lost], tx0, ty0,
                               oz1[lost], ox1[lost], oy1[lost], tx1, ty1, xp, yp);
                if (std::hypot(sx[k]-xp, sy[k]-yp) < 5.0) sh[st] = true;
            }
            int nm = (!sh[0]) + (!sh[1]) + (!sh[2]);
            if (nm >= 2) ev.stage_code = 7;
            else if (!sh[0]) ev.stage_code = 4;
            else if (!sh[1]) ev.stage_code = 5;
            else if (!sh[2]) ev.stage_code = 6;
            else {
                bool gh[3] = {false, false, false};
                for (size_t k = 0; k < gz.size(); ++k) {
                    int st = stationOfZ(gz[k]); if (st < 0 || gh[st]) continue;
                    double xp = 0, yp = 0;
                    hermitePredict(gz[k], oz0[lost], ox0[lost], oy0[lost], tx0, ty0,
                                   oz1[lost], ox1[lost], oy1[lost], tx1, ty1, xp, yp);
                    if (std::hypot(gx[k]-xp, gy[k]-yp) < 5.0) gh[st] = true;
                }
                ev.stage_code = (gh[0] && gh[1] && gh[2]) ? 9 : 8;
            }
        }
    }
    return ev;
}

inline void addStage(StageLossCounts& s, int code) {
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

struct TreeBind {
    std::unique_ptr<TFile> file;
    std::unique_ptr<TTreeReader> reader;
    std::unique_ptr<TTreeReaderValue<Int_t>> lt, og_lt;
    std::unique_ptr<TTreeReaderValue<std::vector<double>>> c2, nd, x0, y0, pz0, x1, y1;
    std::unique_ptr<TTreeReaderValue<std::vector<double>>> ox0, oy0, oz0, opx0, opy0, opz0;
    std::unique_ptr<TTreeReaderValue<std::vector<double>>> ox1, oy1, oz1, opx1, opy1, opz1;
    std::unique_ptr<TTreeReaderValue<std::vector<double>>> sx, sy, sz, gx, gy, gz;
    Long64_t n_entries = 0;

    void open(const std::string& path) {
        file.reset(TFile::Open(path.c_str(), "READ"));
        if (!file || file->IsZombie()) throw std::runtime_error("Cannot open " + path);
        TTree* t = dynamic_cast<TTree*>(file->Get("nt"));
        if (!t) throw std::runtime_error("No 'nt' tree in " + path);
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
        ox0 = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "og_Track_x0");
        oy0 = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "og_Track_y0");
        oz0 = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "og_Track_z0");
        opx0 = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "og_Track_px0");
        opy0 = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "og_Track_py0");
        opz0 = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "og_Track_pz0");
        ox1 = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "og_Track_x1");
        oy1 = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "og_Track_y1");
        oz1 = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "og_Track_z1");
        opx1 = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "og_Track_px1");
        opy1 = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "og_Track_py1");
        opz1 = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "og_Track_pz1");
        sx = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "TrackSegment_x");
        sy = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "TrackSegment_y");
        sz = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "TrackSegment_z");
        gx = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "GhostBustedTrackSegment_x");
        gy = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "GhostBustedTrackSegment_y");
        gz = std::make_unique<TTreeReaderValue<std::vector<double>>>(*reader, "GhostBustedTrackSegment_z");
    }
};

std::map<std::string, std::string> indexShardDir(const std::string& dir) {
    std::map<std::string, std::string> m;
    if (dir.empty() || !fs::exists(dir)) return m;
    std::regex re("/(008749|10004X)/.*PosBin([0-8])_NegBin([0-8])");
    for (const auto& e : fs::recursive_directory_iterator(dir)) {
        if (!e.is_regular_file()) continue;
        std::string p = e.path().string();
        if (p.size() < 5 || p.substr(p.size() - 5) != ".root") continue;
        if (p.find("raw") != std::string::npos || p.find("Pos_008749-Neg_010698") != std::string::npos) continue;
        std::smatch sm;
        if (std::regex_search(p, sm, re)) {
            m[sm[1].str() + "_PosBin" + sm[2].str() + "_NegBin" + sm[3].str()] = p;
        }
    }
    return m;
}

struct Quad {
    int split_idx = 0;
    std::string ds_name;
    int pb = 0, nb = 0;
    std::string paths[N_ARMS];
};

void writeEffScalar(std::ostream& os, long long pass, long long den) {
    double e, l, h; clopperPearson1Sigma(pass, den, e, l, h);
    os << "{\"pass\": " << pass << ", \"den\": " << den
       << ", \"eff\": " << fmtD(e) << ", \"err_lo\": " << fmtD(l) << ", \"err_hi\": " << fmtD(h) << "}";
}

void writePairScalar(std::ostream& os, const PairedCell& c) {
    double d, e; mcnemarDelta(c.bp, c.g, c.l, c.bf, d, e);
    os << "{\"n_den\": " << (c.bp+c.g+c.l+c.bf)
       << ", \"n_both_pass\": " << c.bp << ", \"n_gain\": " << c.g
       << ", \"n_loss\": " << c.l << ", \"n_both_fail\": " << c.bf
       << ", \"delta\": " << fmtD(d) << ", \"delta_err\": " << fmtD(e) << "}";
}

void writeStageJson(std::ostream& os, const StageLossCounts& s) {
    os << "{\"n_den\": " << s.n_den << ", \"pass_c2R\": " << s.pass_c2R
       << ", \"lt2_fail_c2R\": " << s.lt2_fail_c2R << ", \"lt_ge3\": " << s.lt_ge3
       << ", \"lt0\": " << s.lt0 << ", \"lt1_total\": " << s.lt1_total
       << ", \"lt1_miss_seg_st1\": " << s.lt1_miss_seg_st1
       << ", \"lt1_miss_seg_st2\": " << s.lt1_miss_seg_st2
       << ", \"lt1_miss_seg_st3\": " << s.lt1_miss_seg_st3
       << ", \"lt1_miss_seg_multi\": " << s.lt1_miss_seg_multi
       << ", \"lt1_killed_by_gb\": " << s.lt1_killed_by_gb
       << ", \"lt1_all3gb_lost_in_ckf\": " << s.lt1_all3gb_lost_in_ckf
       << ", \"lt1_unmatched_kept\": " << s.lt1_unmatched_kept << "}";
}

void writeEffCurve2M(std::ostream& os, const CurveSpec& sp, const EffBins2M& bs, bool use_ge2) {
    int nb = bs.nbins;
    double bw = (bs.xmax - bs.xmin) / nb;
    const auto& pv = use_ge2 ? bs.pass_ge2 : bs.pass_c2r;
    os << "{\"name\": \"" << sp.key << "\", \"xlabel\": \"" << jsonEscape(sp.xlabel)
       << "\", \"unit\": \"mm\", \"edges\": [";
    for (int i = 0; i <= nb; ++i) { if (i) os << ", "; os << fmtD(bs.xmin + i * bw); }
    os << "], \"centers\": [";
    for (int i = 0; i < nb; ++i) { if (i) os << ", "; os << fmtD(bs.xmin + (i + 0.5) * bw); }
    os << "], \"half_widths\": [";
    for (int i = 0; i < nb; ++i) { if (i) os << ", "; os << fmtD(0.5 * bw); }
    os << "], \"underflow_den\": " << bs.underflow << ", \"overflow_den\": " << bs.overflow
       << ", \"invalid_den\": " << bs.invalid << ", \"den\": [";
    for (int i = 0; i < nb; ++i) { if (i) os << ", "; os << bs.den[i]; }
    os << "], \"pass\": [";
    for (int i = 0; i < nb; ++i) { if (i) os << ", "; os << pv[i]; }
    std::vector<double> eff(nb), lo(nb), hi(nb);
    for (int i = 0; i < nb; ++i) clopperPearson1Sigma(pv[i], bs.den[i], eff[i], lo[i], hi[i]);
    os << "], \"eff\": [";
    for (int i = 0; i < nb; ++i) { if (i) os << ", "; os << fmtD(eff[i]); }
    os << "], \"err_lo\": [";
    for (int i = 0; i < nb; ++i) { if (i) os << ", "; os << fmtD(lo[i]); }
    os << "], \"err_hi\": [";
    for (int i = 0; i < nb; ++i) { if (i) os << ", "; os << fmtD(hi[i]); }
    os << "]}";
}

void writePairCurve2M(std::ostream& os, const std::string& pkey, const std::string& va, const std::string& vb,
                      const std::string& ds, const std::string& met, const PairedBins2M& pb, bool use_ge2) {
    int nb = pb.nbins;
    double bw = (pb.xmax - pb.xmin) / nb;
    const auto& cs = use_ge2 ? pb.c_ge2 : pb.c_c2r;
    os << "{\"pair\": \"" << pkey << "\", \"ver_a\": \"" << va << "\", \"ver_b\": \"" << vb
       << "\", \"dataset\": \"" << ds << "\", \"metric\": \"" << met << "\", \"edges\": [";
    for (int i = 0; i <= nb; ++i) { if (i) os << ", "; os << fmtD(pb.xmin + i * bw); }
    os << "], \"centers\": [";
    for (int i = 0; i < nb; ++i) { if (i) os << ", "; os << fmtD(pb.xmin + (i + 0.5) * bw); }
    os << "], \"half_widths\": [";
    for (int i = 0; i < nb; ++i) { if (i) os << ", "; os << fmtD(0.5 * bw); }
    os << "], \"n_den\": [";
    for (int i = 0; i < nb; ++i) { if (i) os << ", "; os << (cs[i].bp+cs[i].g+cs[i].l+cs[i].bf); }
    os << "], \"n_both_pass\": [";
    for (int i = 0; i < nb; ++i) { if (i) os << ", "; os << cs[i].bp; }
    os << "], \"n_gain\": [";
    for (int i = 0; i < nb; ++i) { if (i) os << ", "; os << cs[i].g; }
    os << "], \"n_loss\": [";
    for (int i = 0; i < nb; ++i) { if (i) os << ", "; os << cs[i].l; }
    os << "], \"n_both_fail\": [";
    for (int i = 0; i < nb; ++i) { if (i) os << ", "; os << cs[i].bf; }
    std::vector<double> d(nb), e(nb);
    for (int i = 0; i < nb; ++i) mcnemarDelta(cs[i].bp, cs[i].g, cs[i].l, cs[i].bf, d[i], e[i]);
    os << "], \"delta\": [";
    for (int i = 0; i < nb; ++i) { if (i) os << ", "; os << fmtD(d[i]); }
    os << "], \"delta_err\": [";
    for (int i = 0; i < nb; ++i) { if (i) os << ", "; os << fmtD(e[i]); }
    os << "]}";
}

void writeRatioCurve(std::ostream& os, const EffBins2M& mc, const EffBins2M& da) {
    int nb = mc.nbins;
    double bw = (mc.xmax - mc.xmin) / nb;
    os << "{\"edges\": [";
    for (int i = 0; i <= nb; ++i) { if (i) os << ", "; os << fmtD(mc.xmin + i * bw); }
    os << "], \"centers\": [";
    for (int i = 0; i < nb; ++i) { if (i) os << ", "; os << fmtD(mc.xmin + (i + 0.5) * bw); }
    os << "], \"half_widths\": [";
    for (int i = 0; i < nb; ++i) { if (i) os << ", "; os << fmtD(0.5 * bw); }
    os << "], ";
    std::vector<double> rat(nb), rlo(nb), rhi(nb);
    for (int i = 0; i < nb; ++i) {
        double em, lm, hm, ed, ld, hd;
        clopperPearson1Sigma(mc.pass_c2r[i], mc.den[i], em, lm, hm);
        clopperPearson1Sigma(da.pass_c2r[i], da.den[i], ed, ld, hd);
        if (std::isfinite(em) && std::isfinite(ed) && em > 0.0 && ed > 0.0) {
            rat[i] = em / ed;
            rlo[i] = rat[i] * std::hypot(lm / em, hd / ed);
            rhi[i] = rat[i] * std::hypot(hm / em, ld / ed);
        } else rat[i] = rlo[i] = rhi[i] = std::numeric_limits<double>::quiet_NaN();
    }
    os << "\"ratio\": ["; for (int i = 0; i < nb; ++i) { if (i) os << ", "; os << fmtD(rat[i]); }
    os << "], \"err_lo\": ["; for (int i = 0; i < nb; ++i) { if (i) os << ", "; os << fmtD(rlo[i]); }
    os << "], \"err_hi\": ["; for (int i = 0; i < nb; ++i) { if (i) os << ", "; os << fmtD(rhi[i]); }
    os << "]}";
}

int main(int argc, char** argv) {
    auto t0 = std::chrono::steady_clock::now();
    std::string eos_base = "/eos/experiment/faser/data0/overlay";
    std::string latest_dir = "/eos/home-x/xcheng/FASER/geoRotation_test_results/stage08_overlay_2track/run_20261125_v1/shards";
    std::string t2a_dir = "/eos/home-x/xcheng/FASER/geoRotation_test_results/stage08_overlay_2track/run_20261125_v1_t2a_covphi/shards";
    std::string t2b_dir = "/eos/home-x/xcheng/FASER/geoRotation_test_results/stage08_overlay_2track/run_20261125_v1_t2b_7fd40df8/shards";
    std::string out_dir = "/eos/home-x/xcheng/FASER/geoRotation_test_results/stage08_overlay_2track/run_20261125_v1/tier2_attribution";
    int n_threads = 16;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--eos-base" && i + 1 < argc) eos_base = argv[++i];
        else if (arg == "--latest-dir" && i + 1 < argc) latest_dir = argv[++i];
        else if (arg == "--t2a-dir" && i + 1 < argc) t2a_dir = argv[++i];
        else if (arg == "--t2b-dir" && i + 1 < argc) t2b_dir = argv[++i];
        else if (arg == "--out-dir" && i + 1 < argc) out_dir = argv[++i];
        else if (arg == "--threads" && i + 1 < argc) n_threads = std::stoi(argv[++i]);
        else { std::cerr << "Unknown argument: " << arg << "\n"; return 2; }
    }

    ROOT::EnableThreadSafety();
    omp_set_num_threads(n_threads);

    auto map_latest = indexShardDir(latest_dir);
    auto map_t2a = indexShardDir(t2a_dir);
    auto map_t2b = indexShardDir(t2b_dir);

    std::vector<Quad> quads;
    int n_8749 = 0, n_mc = 0;
    for (const std::string ds : {"008749", "10004X"}) {
        for (int pb = 0; pb < 9; ++pb) {
            for (int nb = 0; nb < 9; ++nb) {
                std::string key = ds + "_PosBin" + std::to_string(pb) + "_NegBin" + std::to_string(nb);
                auto ia = map_t2a.find(key);
                auto ib = map_t2b.find(key);
                auto il = map_latest.find(key);
                if (ia == map_t2a.end() || ib == map_t2b.end() || il == map_latest.end()) continue;
                std::string p_base = eos_base + "/" + ds + "/OverlaidNTPs/myOverlay_PosBin" +
                                     std::to_string(pb) + "_NegBin" + std::to_string(nb) + "-RDO.pool-ntuple.root";
                if (!fs::exists(p_base)) continue;
                Quad q;
                q.ds_name = ds; q.pb = pb; q.nb = nb;
                q.split_idx = (ds == "008749") ? (pb == nb ? SP_8749_DIAG : SP_8749_OFF)
                                               : (pb == nb ? SP_MC_DIAG : SP_MC_OFF);
                q.paths[ARM_BASE] = p_base;
                q.paths[ARM_T2B] = ib->second;
                q.paths[ARM_LATEST] = il->second;
                q.paths[ARM_T2A] = ia->second;
                quads.push_back(std::move(q));
                if (ds == "008749") n_8749++; else n_mc++;
            }
        }
    }

    std::cout << "[overlay_tier2_analyze] Found matched quadruples: total=" << quads.size()
              << " (008749=" << n_8749 << ", 10004X=" << n_mc << ")\n";
    if (quads.empty()) {
        std::cerr << "ERROR: Zero matched quadruples found!\n";
        return 1;
    }

    std::vector<ThreadAccum> locals(n_threads);
    std::string fatal_err;
    std::atomic<bool> has_fatal{false};

    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t qi = 0; qi < quads.size(); ++qi) {
        if (has_fatal.load(std::memory_order_relaxed)) continue;
        int tid = omp_get_thread_num();
        try {
            const auto& q = quads[qi];
            TreeBind tb[N_ARMS];
            for (int a = 0; a < N_ARMS; ++a) tb[a].open(q.paths[a]);
            Long64_t N = tb[ARM_BASE].n_entries;
            for (int a = 1; a < N_ARMS; ++a) {
                if (tb[a].n_entries != N)
                    throw std::runtime_error("Entry mismatch on " + q.ds_name + " PosBin" + std::to_string(q.pb) + "_NegBin" + std::to_string(q.nb));
            }
            auto& acc = locals[tid];
            for (int a = 0; a < N_ARMS; ++a) acc.arms[a][q.split_idx].n_files++;

            for (Long64_t e = 0; e < N; ++e) {
                for (int a = 0; a < N_ARMS; ++a)
                    if (!tb[a].reader->Next()) throw std::runtime_error("Premature EOF at entry " + std::to_string(e));

                const auto& ox0 = *tb[ARM_BASE].ox0->Get();
                const auto& oy0 = *tb[ARM_BASE].oy0->Get();
                const auto& oz0 = *tb[ARM_BASE].oz0->Get();
                const auto& opx0 = *tb[ARM_BASE].opx0->Get();
                const auto& opy0 = *tb[ARM_BASE].opy0->Get();
                const auto& opz0 = *tb[ARM_BASE].opz0->Get();
                const auto& ox1 = *tb[ARM_BASE].ox1->Get();
                const auto& oy1 = *tb[ARM_BASE].oy1->Get();
                const auto& oz1 = *tb[ARM_BASE].oz1->Get();
                const auto& opx1 = *tb[ARM_BASE].opx1->Get();
                const auto& opy1 = *tb[ARM_BASE].opy1->Get();
                const auto& opz1 = *tb[ARM_BASE].opz1->Get();

                bool id_ok = (*tb[ARM_BASE].og_lt->Get() == 2 && ox0.size() == 2 && ox1.size() == 2);
                for (int a = 1; a < N_ARMS && id_ok; ++a) {
                    const auto& vx = *tb[a].ox0->Get();
                    const auto& vy = *tb[a].oy0->Get();
                    const auto& vpz = *tb[a].opz0->Get();
                    if (*tb[a].og_lt->Get() != 2 || vx.size() != 2 ||
                        std::abs(vx[0]-ox0[0]) > 1e-4 || std::abs(vy[0]-oy0[0]) > 1e-4 || std::abs(vpz[0]-opz0[0]) > 1e-2 ||
                        std::abs(vx[1]-ox0[1]) > 1e-4 || std::abs(vy[1]-oy0[1]) > 1e-4 || std::abs(vpz[1]-opz0[1]) > 1e-2) {
                        id_ok = false;
                    }
                }
                if (!id_ok) {
                    acc.og_failed++;
                    throw std::runtime_error("og_Track mismatch at " + q.ds_name + " entry " + std::to_string(e));
                }
                acc.og_passed++;

                double dx0 = ox0[0] - ox0[1], dy0 = oy0[0] - oy0[1];
                double dx1 = ox1[0] - ox1[1], dy1 = oy1[0] - oy1[1];
                double adx0 = std::abs(dx0), ady0 = std::abs(dy0);
                double dr1 = std::hypot(dx1, dy1);
                double cvals[N_CURVES] = {ady0, dr1, dr1, adx0, dy1};
                bool in_reg[N_REGIMES] = {
                    true,
                    (ady0 >= 12.5 && ady0 <= 18.5),
                    (ady0 < 2.0),
                    (ady0 >= 19.5),
                    (dr1 < 5.0),
                    (dr1 >= 30.0 && dr1 <= 70.0),
                    (adx0 >= 70.0 && adx0 <= 100.0)
                };

                EvalEv evs[N_ARMS];
                for (int a = 0; a < N_ARMS; ++a) {
                    evs[a] = evaluateArm(
                        *tb[a].lt->Get(), *tb[a].c2->Get(), *tb[a].nd->Get(),
                        *tb[a].x0->Get(), *tb[a].y0->Get(), *tb[a].pz0->Get(),
                        *tb[a].x1->Get(), *tb[a].y1->Get(),
                        ox0, oy0, oz0, opx0, opy0, opz0, ox1, oy1, oz1, opx1, opy1, opz1,
                        *tb[a].sx->Get(), *tb[a].sy->Get(), *tb[a].sz->Get(),
                        *tb[a].gx->Get(), *tb[a].gy->Get(), *tb[a].gz->Get());

                    auto& st = acc.arms[a][q.split_idx];
                    st.n_events++;
                    st.lt_counts[std::clamp(evs[a].lt, 0, 5)]++;
                    if (evs[a].lt2_fail) st.lt2_fail_c2R++;
                    for (int r = 0; r < N_REGIMES; ++r) {
                        if (!in_reg[r]) continue;
                        st.reg_den[r]++;
                        for (int m = 0; m < N_METRICS; ++m) if (evs[a].pass[m]) st.reg_pass[r][m]++;
                        addStage(st.stage_reg[r], evs[a].stage_code);
                    }
                    for (int c = 0; c < N_CURVES; ++c) st.curves[c].fill(cvals[c], evs[a].pass[MET_C2R], evs[a].pass[MET_GE2]);
                    for (int k = 0; k < 2; ++k) {
                        double rv = evs[a].res[k];
                        if (!std::isfinite(rv)) continue;
                        if (rv < RES_MIN) st.res_uf++;
                        else if (rv >= RES_MAX) st.res_of++;
                        else {
                            int b = std::clamp(static_cast<int>((rv - RES_MIN) / (RES_MAX - RES_MIN) * N_RES_BINS), 0, N_RES_BINS - 1);
                            st.res_counts[b]++; st.res_tot++; st.res_sum += rv; st.res_sum2 += rv * rv;
                        }
                    }
                }

                const int pa[N_PAIRS] = {ARM_BASE, ARM_BASE, ARM_BASE, ARM_LATEST};
                const int pb_arm[N_PAIRS] = {ARM_T2B, ARM_LATEST, ARM_T2A, ARM_T2A};
                for (int p = 0; p < N_PAIRS; ++p) {
                    auto& ps = acc.pairs[p][q.split_idx];
                    const auto& eA = evs[pa[p]];
                    const auto& eB = evs[pb_arm[p]];
                    for (int r = 0; r < N_REGIMES; ++r) {
                        if (!in_reg[r]) continue;
                        for (int m = 0; m < N_METRICS; ++m) ps.reg_cells[r][m].add(eA.pass[m], eB.pass[m]);
                    }
                    for (int c = 0; c < N_CURVES; ++c)
                        ps.curves[c].fill(cvals[c], eA.pass[MET_C2R], eB.pass[MET_C2R], eA.pass[MET_GE2], eB.pass[MET_GE2]);
                }
            }
        } catch (const std::exception& ex) {
            #pragma omp critical
            { if (fatal_err.empty()) fatal_err = ex.what(); has_fatal.store(true, std::memory_order_relaxed); }
        }
    }

    if (!fatal_err.empty()) { std::cerr << "FATAL: " << fatal_err << "\n"; return 1; }

    ThreadAccum merged;
    for (int t = 0; t < n_threads; ++t) merged.merge(locals[t]);

    const std::array<std::string, N_ARMS> arm_names = {"baseline_jul8", "t2b_7fd40df8", "latest_218813c2", "t2a_covphi"};
    const std::array<std::string, N_PAIRS> pair_names = {
        "t2b_7fd40df8_minus_baseline_jul8", "latest_218813c2_minus_baseline_jul8",
        "t2a_covphi_minus_baseline_jul8", "t2a_covphi_minus_latest_218813c2"
    };
    const std::array<std::pair<std::string, std::string>, N_PAIRS> pair_ab = {{
        {"baseline_jul8", "t2b_7fd40df8"}, {"baseline_jul8", "latest_218813c2"},
        {"baseline_jul8", "t2a_covphi"}, {"latest_218813c2", "t2a_covphi"}
    }};
    const std::vector<std::string> ds_names = {
        "008749", "008749_diag", "008749_offdiag", "10004X", "10004X_diag", "10004X_offdiag"
    };
    const std::array<std::string, N_REGIMES> reg_names = {
        "global", "crater_abs_dY0_12p5_18p5", "small_abs_dY0_lt2", "large_abs_dY0_ge19p5",
        "small_dR1_lt5", "plateau_dR1_30_70", "large_abs_dX0_70_100"
    };
    const std::array<std::string, N_METRICS> met_names = {"eff_pawan_c2R", "eff_lt2", "eff_ge2", "eff_matched"};

    std::map<std::string, std::map<std::string, ArmSplitStats>> amap;
    for (int a = 0; a < N_ARMS; ++a) {
        const auto& an = arm_names[a];
        amap[an]["008749_diag"] = merged.arms[a][SP_8749_DIAG];
        amap[an]["008749_offdiag"] = merged.arms[a][SP_8749_OFF];
        ArmSplitStats s87 = merged.arms[a][SP_8749_DIAG]; s87.merge(merged.arms[a][SP_8749_OFF]);
        amap[an]["008749"] = s87;
        amap[an]["10004X_diag"] = merged.arms[a][SP_MC_DIAG];
        amap[an]["10004X_offdiag"] = merged.arms[a][SP_MC_OFF];
        ArmSplitStats smc = merged.arms[a][SP_MC_DIAG]; smc.merge(merged.arms[a][SP_MC_OFF]);
        amap[an]["10004X"] = smc;
    }

    std::map<std::string, std::map<std::string, PairSplitStats>> pmap;
    for (int p = 0; p < N_PAIRS; ++p) {
        const auto& pn = pair_names[p];
        pmap[pn]["008749_diag"] = merged.pairs[p][SP_8749_DIAG];
        pmap[pn]["008749_offdiag"] = merged.pairs[p][SP_8749_OFF];
        PairSplitStats p87 = merged.pairs[p][SP_8749_DIAG]; p87.merge(merged.pairs[p][SP_8749_OFF]);
        pmap[pn]["008749"] = p87;
        pmap[pn]["10004X_diag"] = merged.pairs[p][SP_MC_DIAG];
        pmap[pn]["10004X_offdiag"] = merged.pairs[p][SP_MC_OFF];
        PairSplitStats pmc = merged.pairs[p][SP_MC_DIAG]; pmc.merge(merged.pairs[p][SP_MC_OFF]);
        pmap[pn]["10004X"] = pmc;
    }

    double elapsed_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    fs::create_directories(out_dir + "/bin_counts");

    // Write CSVs
    for (int a = 0; a < N_ARMS; ++a) {
        for (const std::string dn : {"008749", "10004X"}) {
            const auto& st = amap[arm_names[a]][dn];
            for (int c = 0; c < N_CURVES; ++c) {
                const auto& sp = CURVE_SPECS[c];
                const auto& bs = st.curves[c];
                double bw = (sp.xmax - sp.xmin) / sp.nbins;
                std::ofstream f(out_dir + "/bin_counts/" + arm_names[a] + "_" + dn + "_" + sp.key + ".csv");
                f << "bin_lo,bin_hi,bin_center,den,pass_c2R,eff_c2R,err_lo_c2R,err_hi_c2R,pass_ge2,eff_ge2\n";
                for (int i = 0; i < sp.nbins; ++i) {
                    double lo = sp.xmin + i * bw, hi = lo + bw, cen = 0.5 * (lo + hi);
                    double e0, l0, h0, e2, l2, h2;
                    clopperPearson1Sigma(bs.pass_c2r[i], bs.den[i], e0, l0, h0);
                    clopperPearson1Sigma(bs.pass_ge2[i], bs.den[i], e2, l2, h2);
                    f << std::setprecision(8) << lo << "," << hi << "," << cen << "," << bs.den[i] << ","
                      << bs.pass_c2r[i] << "," << (std::isfinite(e0)?e0:-1.0) << "," << l0 << "," << h0 << ","
                      << bs.pass_ge2[i] << "," << (std::isfinite(e2)?e2:-1.0) << "\n";
                }
            }
        }
    }
    for (int p = 0; p < N_PAIRS; ++p) {
        for (const std::string dn : {"008749", "10004X"}) {
            const auto& ps = pmap[pair_names[p]][dn];
            for (int c = 0; c < N_CURVES; ++c) {
                const auto& sp = CURVE_SPECS[c];
                const auto& pb = ps.curves[c];
                double bw = (sp.xmax - sp.xmin) / sp.nbins;
                std::ofstream f(out_dir + "/bin_counts/paired_" + pair_names[p] + "_" + dn + "_" + sp.key + ".csv");
                f << "bin_lo,bin_hi,bin_center,n_den,n_both_pass,n_gain,n_loss,n_both_fail,delta_c2R,delta_err_c2R,delta_ge2,delta_err_ge2\n";
                for (int i = 0; i < sp.nbins; ++i) {
                    double lo = sp.xmin + i * bw, hi = lo + bw, cen = 0.5 * (lo + hi);
                    const auto& c0 = pb.c_c2r[i]; const auto& c2 = pb.c_ge2[i];
                    double d0, e0, d2, e2;
                    mcnemarDelta(c0.bp, c0.g, c0.l, c0.bf, d0, e0);
                    mcnemarDelta(c2.bp, c2.g, c2.l, c2.bf, d2, e2);
                    f << std::setprecision(8) << lo << "," << hi << "," << cen << ","
                      << (c0.bp+c0.g+c0.l+c0.bf) << "," << c0.bp << "," << c0.g << "," << c0.l << "," << c0.bf << ","
                      << (std::isfinite(d0)?d0:0.0) << "," << (std::isfinite(e0)?e0:0.0) << ","
                      << (std::isfinite(d2)?d2:0.0) << "," << (std::isfinite(e2)?e2:0.0) << "\n";
                }
            }
        }
    }

    // Write tier2_attribution_summary.json
    std::ofstream js(out_dir + "/tier2_attribution_summary.json");
    js << "{\n"
       << "  \"schema_version\": \"1.0\",\n"
       << "  \"arms\": [\"baseline_jul8\", \"t2b_7fd40df8\", \"latest_218813c2\", \"t2a_covphi\"],\n"
       << "  \"arm_labels\": {\n"
       << "    \"baseline_jul8\": \"Baseline-Jul8 (Unrotated)\",\n"
       << "    \"t2b_7fd40df8\": \"Pre-fix GeoRot (7fd40df8)\",\n"
       << "    \"latest_218813c2\": \"Latest GeoRot (218813c2, covPhi=0.1089)\",\n"
       << "    \"t2a_covphi\": \"Latest GeoRot + covPhi=8.41e-6\"\n"
       << "  },\n"
       << "  \"pairs\": [\"t2b_7fd40df8_minus_baseline_jul8\", \"latest_218813c2_minus_baseline_jul8\", \"t2a_covphi_minus_baseline_jul8\", \"t2a_covphi_minus_latest_218813c2\"],\n"
       << "  \"datasets\": [\"008749\", \"008749_diag\", \"008749_offdiag\", \"10004X\", \"10004X_diag\", \"10004X_offdiag\"],\n"
       << "  \"regimes\": [\"global\", \"crater_abs_dY0_12p5_18p5\", \"small_abs_dY0_lt2\", \"large_abs_dY0_ge19p5\", \"small_dR1_lt5\", \"plateau_dR1_30_70\", \"large_abs_dX0_70_100\"],\n"
       << "  \"provenance\": {\n"
       << "    \"frame\": \"FASER\", \"cond_tag\": \"OFLCOND-FASER-05\",\n"
       << "    \"matched_quadruples\": {\"008749\": " << n_8749 << ", \"10004X\": " << n_mc << ", \"total\": " << quads.size() << "},\n"
       << "    \"og_identity_checks_passed\": " << merged.og_passed << ", \"og_identity_checks_failed\": " << merged.og_failed << ",\n"
       << "    \"runtime_seconds\": " << std::fixed << std::setprecision(3) << elapsed_s << ", \"n_threads\": " << n_threads << "\n"
       << "  },\n";

    // global_summary
    js << "  \"global_summary\": {\n";
    for (int a = 0; a < N_ARMS; ++a) {
        const auto& an = arm_names[a];
        js << "    \"" << an << "\": {\n";
        for (size_t di = 0; di < ds_names.size(); ++di) {
            const auto& dn = ds_names[di];
            const auto& s = amap[an][dn];
            js << "      \"" << dn << "\": {\"n_files\": " << s.n_files << ", \"n_events\": " << s.n_events
               << ", \"lt_counts\": {\"0\": " << s.lt_counts[0] << ", \"1\": " << s.lt_counts[1]
               << ", \"2\": " << s.lt_counts[2] << ", \"3\": " << s.lt_counts[3]
               << ", \"4\": " << s.lt_counts[4] << ", \"5plus\": " << s.lt_counts[5]
               << "}, \"lt2_fail_c2R\": " << s.lt2_fail_c2R << ", \"regimes\": {\n";
            for (int r = 0; r < N_REGIMES; ++r) {
                js << "        \"" << reg_names[r] << "\": {";
                for (int m = 0; m < N_METRICS; ++m) {
                    if (m) js << ", ";
                    js << "\"" << met_names[m] << "\": ";
                    writeEffScalar(js, s.reg_pass[r][m], s.reg_den[r]);
                }
                js << "}" << (r + 1 < N_REGIMES ? ",\n" : "\n");
            }
            js << "      }}" << (di + 1 < ds_names.size() ? ",\n" : "\n");
        }
        js << "    }" << (a + 1 < N_ARMS ? ",\n" : "\n");
    }
    js << "  },\n";

    // paired_deltas
    js << "  \"paired_deltas\": {\n";
    for (int p = 0; p < N_PAIRS; ++p) {
        const auto& pn = pair_names[p];
        js << "    \"" << pn << "\": {\n";
        for (size_t di = 0; di < ds_names.size(); ++di) {
            const auto& dn = ds_names[di];
            const auto& ps = pmap[pn][dn];
            js << "      \"" << dn << "\": {\n";
            for (int r = 0; r < N_REGIMES; ++r) {
                js << "        \"" << reg_names[r] << "\": {";
                for (int m = 0; m < N_METRICS; ++m) {
                    if (m) js << ", ";
                    js << "\"" << met_names[m] << "\": ";
                    writePairScalar(js, ps.reg_cells[r][m]);
                }
                js << "}" << (r + 1 < N_REGIMES ? ",\n" : "\n");
            }
            js << "      }" << (di + 1 < ds_names.size() ? ",\n" : "\n");
        }
        js << "    }" << (p + 1 < N_PAIRS ? ",\n" : "\n");
    }
    js << "  },\n";

    // stage_of_loss
    js << "  \"stage_of_loss\": {\n";
    for (int a = 0; a < N_ARMS; ++a) {
        const auto& an = arm_names[a];
        js << "    \"" << an << "\": {\n";
        for (size_t di = 0; di < ds_names.size(); ++di) {
            const auto& dn = ds_names[di];
            const auto& s = amap[an][dn];
            js << "      \"" << dn << "\": {\n";
            for (int r = 0; r < N_REGIMES; ++r) {
                js << "        \"" << reg_names[r] << "\": ";
                writeStageJson(js, s.stage_reg[r]);
                js << (r + 1 < N_REGIMES ? ",\n" : "\n");
            }
            js << "      }" << (di + 1 < ds_names.size() ? ",\n" : "\n");
        }
        js << "    }" << (a + 1 < N_ARMS ? ",\n" : "\n");
    }
    js << "  },\n";

    // curves
    js << "  \"curves\": {\n";
    for (int a = 0; a < N_ARMS; ++a) {
        const auto& an = arm_names[a];
        js << "    \"" << an << "\": {\n";
        std::vector<std::string> cds = {"008749", "10004X"};
        for (size_t di = 0; di < cds.size(); ++di) {
            const auto& dn = cds[di];
            const auto& s = amap[an][dn];
            js << "      \"" << dn << "\": {\n";
            for (int c = 0; c < N_CURVES; ++c) {
                js << "        \"" << CURVE_SPECS[c].key << "\": {\"eff_pawan_c2R\": ";
                writeEffCurve2M(js, CURVE_SPECS[c], s.curves[c], false);
                js << ", \"eff_ge2\": ";
                writeEffCurve2M(js, CURVE_SPECS[c], s.curves[c], true);
                js << "}" << (c + 1 < N_CURVES ? ",\n" : "\n");
            }
            js << "      }" << (di + 1 < cds.size() ? ",\n" : "\n");
        }
        js << "    }" << (a + 1 < N_ARMS ? ",\n" : "\n");
    }
    js << "  },\n";

    // paired_curves
    js << "  \"paired_curves\": {\n";
    for (int p = 0; p < N_PAIRS; ++p) {
        const auto& pn = pair_names[p];
        js << "    \"" << pn << "\": {\n";
        std::vector<std::string> cds = {"008749", "10004X"};
        for (size_t di = 0; di < cds.size(); ++di) {
            const auto& dn = cds[di];
            const auto& ps = pmap[pn][dn];
            js << "      \"" << dn << "\": {\n";
            for (int c = 0; c < N_CURVES; ++c) {
                js << "        \"" << CURVE_SPECS[c].key << "\": {\"eff_pawan_c2R\": ";
                writePairCurve2M(js, pn, pair_ab[p].first, pair_ab[p].second, dn, "eff_pawan_c2R", ps.curves[c], false);
                js << ", \"eff_ge2\": ";
                writePairCurve2M(js, pn, pair_ab[p].first, pair_ab[p].second, dn, "eff_ge2", ps.curves[c], true);
                js << "}" << (c + 1 < N_CURVES ? ",\n" : "\n");
            }
            js << "      }" << (di + 1 < cds.size() ? ",\n" : "\n");
        }
        js << "    }" << (p + 1 < N_PAIRS ? ",\n" : "\n");
    }
    js << "  },\n";

    // mc_data_ratios
    js << "  \"mc_data_ratios\": {\n";
    for (int a = 0; a < N_ARMS; ++a) {
        const auto& an = arm_names[a];
        js << "    \"" << an << "\": {\n";
        for (int c = 0; c < N_CURVES; ++c) {
            js << "      \"" << CURVE_SPECS[c].key << "\": {\"eff_pawan_c2R\": ";
            writeRatioCurve(js, amap[an]["10004X"].curves[c], amap[an]["008749"].curves[c]);
            js << "}" << (c + 1 < N_CURVES ? ",\n" : "\n");
        }
        js << "    }" << (a + 1 < N_ARMS ? ",\n" : "\n");
    }
    js << "  },\n";

    // momentum_residual
    js << "  \"momentum_residual\": {\n";
    double rbw = (RES_MAX - RES_MIN) / N_RES_BINS;
    for (int a = 0; a < N_ARMS; ++a) {
        const auto& an = arm_names[a];
        js << "    \"" << an << "\": {\n";
        std::vector<std::string> cds = {"008749", "10004X"};
        for (size_t di = 0; di < cds.size(); ++di) {
            const auto& dn = cds[di];
            const auto& s = amap[an][dn];
            double rmean = s.res_tot > 0 ? s.res_sum / s.res_tot : 0.0;
            double rstd = s.res_tot > 0 ? std::sqrt(std::max(0.0, s.res_sum2 / s.res_tot - rmean * rmean)) : 0.0;
            auto hq = [&](double q) {
                if (s.res_tot <= 0) return 0.0;
                double tgt = q * s.res_tot, cum = 0.0;
                for (int i = 0; i < N_RES_BINS; ++i) {
                    double prev = cum; cum += s.res_counts[i];
                    if (cum >= tgt && s.res_counts[i] > 0) return RES_MIN + (i + (tgt - prev) / s.res_counts[i]) * rbw;
                }
                return RES_MAX;
            };
            double rmed = hq(0.5), riqr = 0.5 * (hq(0.8413447) - hq(0.1586553));
            js << "      \"" << dn << "\": {\"edges\": [";
            for (int i = 0; i <= N_RES_BINS; ++i) { if (i) js << ", "; js << fmtD(RES_MIN + i * rbw); }
            js << "], \"centers\": [";
            for (int i = 0; i < N_RES_BINS; ++i) { if (i) js << ", "; js << fmtD(RES_MIN + (i + 0.5) * rbw); }
            js << "], \"counts\": [";
            for (int i = 0; i < N_RES_BINS; ++i) { if (i) js << ", "; js << s.res_counts[i]; }
            js << "], \"norm\": [";
            for (int i = 0; i < N_RES_BINS; ++i) { if (i) js << ", "; js << fmtD(s.res_tot > 0 ? static_cast<double>(s.res_counts[i]) / s.res_tot : 0.0); }
            js << "], \"total_matched\": " << s.res_tot << ", \"mean\": " << fmtD(rmean)
               << ", \"std\": " << fmtD(rstd) << ", \"median\": " << fmtD(rmed) << ", \"iqr68\": " << fmtD(riqr)
               << ", \"underflow\": " << s.res_uf << ", \"overflow\": " << s.res_of << "}"
               << (di + 1 < cds.size() ? ",\n" : "\n");
        }
        js << "    }" << (a + 1 < N_ARMS ? ",\n" : "\n");
    }
    js << "  }\n}\n";
    js.close();

    std::cout << "[overlay_tier2_analyze] DONE in " << std::fixed << std::setprecision(2) << elapsed_s
              << " s | og_passed=" << merged.og_passed << " failed=" << merged.og_failed << "\n";
    for (const std::string dn : {"008749", "10004X"}) {
        std::cout << "  === " << dn << " ===\n";
        for (int a = 0; a < N_ARMS; ++a) {
            const auto& s = amap[arm_names[a]][dn];
            auto e = [&](int r){ return s.reg_den[r] > 0 ? static_cast<double>(s.reg_pass[r][MET_C2R])/s.reg_den[r] : 0.0; };
            std::cout << "    " << std::setw(15) << arm_names[a]
                      << " | glob=" << std::setprecision(4) << e(REG_GLOBAL)
                      << " | crater=" << e(REG_CRATER)
                      << " | dR1<5=" << e(REG_SMALL_DR1)
                      << " | plat=" << e(REG_PLATEAU)
                      << " | |dX0|>70=" << e(REG_LARGE_DX0)
                      << " | lt0=" << s.lt_counts[0] << "\n";
        }
    }
    return 0;
}
