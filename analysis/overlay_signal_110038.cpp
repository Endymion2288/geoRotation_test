// Compiled C++17 / ROOT RDataFrame extractor for A' signal MC (DSID 110038, m = 100 MeV, eps = 1e-5).
// Extracts BOTH truth-level unbiased A' daughter (e+, e-) separation distributions at Station 1 & Station 3
// (primary right-axis overlay curves) AND reconstructed 2-track distributions & efficiency vs truth fiducial.

#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TROOT.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

struct Hist1D {
    std::string name;
    double xmin = 0.0;
    double xmax = 100.0;
    int nbins = 100;
    std::vector<long long> counts;
    long long underflow = 0;
    long long overflow = 0;
    long long total_in_range = 0;

    Hist1D() = default;
    Hist1D(std::string n, double lo, double hi, int nb)
        : name(std::move(n)), xmin(lo), xmax(hi), nbins(nb), counts(nb, 0) {}

    void fill(double x) {
        if (!std::isfinite(x)) return;
        if (x < xmin) { underflow++; return; }
        if (x >= xmax) { overflow++; return; }
        int b = static_cast<int>((x - xmin) / (xmax - xmin) * nbins);
        if (b < 0) b = 0;
        if (b >= nbins) b = nbins - 1;
        counts[b]++;
        total_in_range++;
    }

    void merge(const Hist1D& o) {
        for (int i = 0; i < nbins; ++i) counts[i] += o.counts[i];
        underflow += o.underflow;
        overflow += o.overflow;
        total_in_range += o.total_in_range;
    }
};

struct SignalResult {
    std::string input_dir;
    int n_files = 0;
    long long n_events_total = 0;
    long long n_events_truth_fiducial = 0;
    long long n_events_reco_2trk_os_r95 = 0;
    long long n_events_reco_2trk_os_c2r = 0;
    long long n_events_fid_and_reco_c2r = 0;
    long long n_truth_dy1_pos_minus_neg_gt0 = 0;
    long long n_truth_dy1_neg_minus_pos_gt0 = 0;
    std::vector<std::string> curve_order;
    std::map<std::string, Hist1D> truth_hists;
    std::map<std::string, Hist1D> reco_hists;
    Hist1D truth_dy1_neg_minus_pos{"dY1_signed_neg_minus_pos", -100.0, 100.0, 200};
    Hist1D truth_dy0_neg_minus_pos{"dY0_signed_neg_minus_pos", -100.0, 100.0, 200};
    Hist1D truth_dr1_pass_c2r{"truth_dR1_zoom_pass_c2r", 0.0, 20.0, 40};
};

std::map<std::string, Hist1D> makeSignalHists(std::vector<std::string>& order) {
    order = {
        "dR1_wide", "dR1_zoom", "dR0_wide", "dR0_zoom",
        "abs_dY0_wide", "abs_dY0_zoom", "abs_dY1_wide", "abs_dY1_zoom",
        "abs_dX0_wide", "dX0_signed", "abs_dX1_wide", "dX1_signed",
        "dY1_signed", "dY0_signed"
    };
    std::map<std::string, Hist1D> m;
    m["dR1_wide"]     = Hist1D("dR1_wide", 0.0, 100.0, 200);
    m["dR1_zoom"]     = Hist1D("dR1_zoom", 0.0, 20.0, 40);
    m["dR0_wide"]     = Hist1D("dR0_wide", 0.0, 100.0, 200);
    m["dR0_zoom"]     = Hist1D("dR0_zoom", 0.0, 20.0, 40);
    m["abs_dY0_wide"] = Hist1D("abs_dY0_wide", 0.0, 100.0, 100);
    m["abs_dY0_zoom"] = Hist1D("abs_dY0_zoom", 0.0, 20.0, 40);
    m["abs_dY1_wide"] = Hist1D("abs_dY1_wide", 0.0, 100.0, 100);
    m["abs_dY1_zoom"] = Hist1D("abs_dY1_zoom", 0.0, 20.0, 40);
    m["abs_dX0_wide"] = Hist1D("abs_dX0_wide", 0.0, 100.0, 100);
    m["dX0_signed"]   = Hist1D("dX0_signed", -100.0, 100.0, 200);
    m["abs_dX1_wide"] = Hist1D("abs_dX1_wide", 0.0, 100.0, 100);
    m["dX1_signed"]   = Hist1D("dX1_signed", -100.0, 100.0, 200);
    m["dY1_signed"]   = Hist1D("dY1_signed", -100.0, 100.0, 200);
    m["dY0_signed"]   = Hist1D("dY0_signed", -100.0, 100.0, 200);
    return m;
}

std::vector<std::string> collectRootFiles(const std::string& dir) {
    std::vector<std::string> files;
    if (!fs::exists(dir)) return files;
    for (const auto& entry : fs::recursive_directory_iterator(dir)) {
        if (!entry.is_regular_file()) continue;
        std::string p = entry.path().string();
        if (p.size() >= 5 && p.substr(p.size() - 5) == ".root") files.push_back(p);
    }
    std::sort(files.begin(), files.end());
    return files;
}

inline void fillSet(std::map<std::string, Hist1D>& h, double dx0, double dy0, double dx1, double dy1) {
    double dr0 = std::hypot(dx0, dy0);
    double dr1 = std::hypot(dx1, dy1);
    h["dR1_wide"].fill(dr1);
    h["dR1_zoom"].fill(dr1);
    h["dR0_wide"].fill(dr0);
    h["dR0_zoom"].fill(dr0);
    h["abs_dY0_wide"].fill(std::abs(dy0));
    h["abs_dY0_zoom"].fill(std::abs(dy0));
    h["abs_dY1_wide"].fill(std::abs(dy1));
    h["abs_dY1_zoom"].fill(std::abs(dy1));
    h["abs_dX0_wide"].fill(std::abs(dx0));
    h["dX0_signed"].fill(dx0);
    h["abs_dX1_wide"].fill(std::abs(dx1));
    h["dX1_signed"].fill(dx1);
    h["dY1_signed"].fill(dy1);
    h["dY0_signed"].fill(dy0);
}

SignalResult processSignalDir(const std::string& dir) {
    SignalResult res;
    res.input_dir = dir;
    res.truth_hists = makeSignalHists(res.curve_order);
    res.reco_hists = makeSignalHists(res.curve_order);

    auto files = collectRootFiles(dir);
    res.n_files = static_cast<int>(files.size());
    if (files.empty()) throw std::runtime_error("No ROOT files found in signal directory: " + dir);

    bool has_truth_st = false;
    {
        std::unique_ptr<TFile> f(TFile::Open(files[0].c_str(), "READ"));
        if (!f || f->IsZombie()) throw std::runtime_error("Failed to open " + files[0]);
        TTree* t = dynamic_cast<TTree*>(f->Get("nt"));
        if (!t) throw std::runtime_error("Missing TTree 'nt' in " + files[0]);
        has_truth_st = (t->GetBranch("truth_pdg") != nullptr &&
                        t->GetBranch("truth_st1_x") != nullptr &&
                        t->GetBranch("truth_st3_x") != nullptr);
    }

    ROOT::RDataFrame df("nt", files);
    const unsigned int nslots = std::max(1u, df.GetNSlots());
    std::vector<SignalResult> slot_res(nslots);
    for (unsigned int s = 0; s < nslots; ++s) {
        slot_res[s].truth_hists = makeSignalHists(slot_res[s].curve_order);
        slot_res[s].reco_hists = makeSignalHists(slot_res[s].curve_order);
    }

    auto eval_reco = [&](SignalResult& sr,
                         int longTracks,
                         const std::vector<double>& x0, const std::vector<double>& y0,
                         const std::vector<double>& x1, const std::vector<double>& y1,
                         const std::vector<int>& charge,
                         const std::vector<double>& chi2, const std::vector<double>& ndof) -> bool {
        if (longTracks != 2) return false;
        if (x0.size() < 2 || y0.size() < 2 || x1.size() < 2 || y1.size() < 2 || charge.size() < 2) return false;
        if (charge[0] * charge[1] >= 0) return false;
        if (std::hypot(x0[0], y0[0]) >= 95.0 || std::hypot(x1[0], y1[0]) >= 95.0 ||
            std::hypot(x0[1], y0[1]) >= 95.0 || std::hypot(x1[1], y1[1]) >= 95.0) return false;
        sr.n_events_reco_2trk_os_r95++;
        bool ok_c2 = (chi2.size() >= 2 && ndof.size() >= 2 &&
                      ndof[0] > 0.0 && ndof[1] > 0.0 &&
                      (chi2[0] / ndof[0]) <= 25.0 && (chi2[1] / ndof[1]) <= 25.0);
        if (ok_c2) sr.n_events_reco_2trk_os_c2r++;
        int ipos = (charge[0] > 0) ? 0 : 1;
        int ineg = 1 - ipos;
        fillSet(sr.reco_hists, x0[ipos] - x0[ineg], y0[ipos] - y0[ineg],
                x1[ipos] - x1[ineg], y1[ipos] - y1[ineg]);
        return ok_c2;
    };

    if (has_truth_st) {
        df.ForeachSlot(
            [&](unsigned int slot,
                int longTracks,
                const std::vector<double>& x0, const std::vector<double>& y0,
                const std::vector<double>& x1, const std::vector<double>& y1,
                const std::vector<int>& charge,
                const std::vector<double>& chi2, const std::vector<double>& ndof,
                const std::vector<int>& tpdg,
                const std::vector<double>& ts1x, const std::vector<double>& ts1y,
                const std::vector<double>& ts3x, const std::vector<double>& ts3y) {
                auto& sr = slot_res[slot];
                sr.n_events_total++;
                bool reco_c2r = eval_reco(sr, longTracks, x0, y0, x1, y1, charge, chi2, ndof);

                int ipos = -1, ineg = -1;
                size_t nt = std::min({tpdg.size(), ts1x.size(), ts1y.size(), ts3x.size(), ts3y.size()});
                for (size_t i = 0; i < nt; ++i) {
                    if (!std::isfinite(ts1x[i]) || !std::isfinite(ts1y[i]) ||
                        !std::isfinite(ts3x[i]) || !std::isfinite(ts3y[i])) continue;
                    if (tpdg[i] == -11 && ipos < 0) ipos = static_cast<int>(i);
                    else if (tpdg[i] == 11 && ineg < 0) ineg = static_cast<int>(i);
                }
                if (ipos >= 0 && ineg >= 0) {
                    double r0p = std::hypot(ts1x[ipos], ts1y[ipos]);
                    double r1p = std::hypot(ts3x[ipos], ts3y[ipos]);
                    double r0m = std::hypot(ts1x[ineg], ts1y[ineg]);
                    double r1m = std::hypot(ts3x[ineg], ts3y[ineg]);
                    if (std::max({r0p, r1p, r0m, r1m}) < 95.0) {
                        sr.n_events_truth_fiducial++;
                        double dx0 = ts1x[ipos] - ts1x[ineg];
                        double dy0 = ts1y[ipos] - ts1y[ineg];
                        double dx1 = ts3x[ipos] - ts3x[ineg];
                        double dy1 = ts3y[ipos] - ts3y[ineg];
                        double dr1 = std::hypot(dx1, dy1);
                        if (dy1 > 0.0) sr.n_truth_dy1_pos_minus_neg_gt0++;
                        if (-dy1 > 0.0) sr.n_truth_dy1_neg_minus_pos_gt0++;
                        fillSet(sr.truth_hists, dx0, dy0, dx1, dy1);
                        sr.truth_dy1_neg_minus_pos.fill(-dy1);
                        sr.truth_dy0_neg_minus_pos.fill(-dy0);
                        if (reco_c2r) {
                            sr.n_events_fid_and_reco_c2r++;
                            sr.truth_dr1_pass_c2r.fill(dr1);
                        }
                    }
                }
            },
            {"longTracks", "Track_x0", "Track_y0", "Track_x1", "Track_y1",
             "Track_charge", "Track_Chi2", "Track_nDoF",
             "truth_pdg", "truth_st1_x", "truth_st1_y", "truth_st3_x", "truth_st3_y"});
    } else {
        df.ForeachSlot(
            [&](unsigned int slot,
                int longTracks,
                const std::vector<double>& x0, const std::vector<double>& y0,
                const std::vector<double>& x1, const std::vector<double>& y1,
                const std::vector<int>& charge,
                const std::vector<double>& chi2, const std::vector<double>& ndof) {
                auto& sr = slot_res[slot];
                sr.n_events_total++;
                eval_reco(sr, longTracks, x0, y0, x1, y1, charge, chi2, ndof);
            },
            {"longTracks", "Track_x0", "Track_y0", "Track_x1", "Track_y1",
             "Track_charge", "Track_Chi2", "Track_nDoF"});
    }

    for (unsigned int s = 0; s < nslots; ++s) {
        res.n_events_total += slot_res[s].n_events_total;
        res.n_events_truth_fiducial += slot_res[s].n_events_truth_fiducial;
        res.n_events_reco_2trk_os_r95 += slot_res[s].n_events_reco_2trk_os_r95;
        res.n_events_reco_2trk_os_c2r += slot_res[s].n_events_reco_2trk_os_c2r;
        res.n_events_fid_and_reco_c2r += slot_res[s].n_events_fid_and_reco_c2r;
        res.n_truth_dy1_pos_minus_neg_gt0 += slot_res[s].n_truth_dy1_pos_minus_neg_gt0;
        res.n_truth_dy1_neg_minus_pos_gt0 += slot_res[s].n_truth_dy1_neg_minus_pos_gt0;
        for (const auto& k : res.curve_order) {
            res.truth_hists[k].merge(slot_res[s].truth_hists[k]);
            res.reco_hists[k].merge(slot_res[s].reco_hists[k]);
        }
        res.truth_dy1_neg_minus_pos.merge(slot_res[s].truth_dy1_neg_minus_pos);
        res.truth_dy0_neg_minus_pos.merge(slot_res[s].truth_dy0_neg_minus_pos);
        res.truth_dr1_pass_c2r.merge(slot_res[s].truth_dr1_pass_c2r);
    }
    return res;
}

void writeHistJson(std::ostream& os, const Hist1D& h, const std::string& indent) {
    double bw = (h.xmax - h.xmin) / h.nbins;
    double tot = static_cast<double>(h.total_in_range);
    os << "{\n";
    os << indent << "  \"edges\": [";
    for (int i = 0; i <= h.nbins; ++i) { if (i > 0) os << ", "; os << std::setprecision(8) << (h.xmin + i * bw); }
    os << "],\n" << indent << "  \"centers\": [";
    for (int i = 0; i < h.nbins; ++i) { if (i > 0) os << ", "; os << std::setprecision(8) << (h.xmin + (i + 0.5) * bw); }
    os << "],\n" << indent << "  \"half_widths\": [";
    for (int i = 0; i < h.nbins; ++i) { if (i > 0) os << ", "; os << std::setprecision(8) << (0.5 * bw); }
    os << "],\n" << indent << "  \"counts\": [";
    for (int i = 0; i < h.nbins; ++i) { if (i > 0) os << ", "; os << h.counts[i]; }
    os << "],\n" << indent << "  \"norm\": [";
    for (int i = 0; i < h.nbins; ++i) {
        if (i > 0) os << ", ";
        os << std::setprecision(8) << ((tot > 0.0) ? static_cast<double>(h.counts[i]) / tot : 0.0);
    }
    os << "],\n" << indent << "  \"norm_err\": [";
    for (int i = 0; i < h.nbins; ++i) {
        if (i > 0) os << ", ";
        os << std::setprecision(8) << ((tot > 0.0) ? std::sqrt(static_cast<double>(h.counts[i])) / tot : 0.0);
    }
    os << "],\n" << indent << "  \"total_in_range\": " << h.total_in_range
       << ",\n" << indent << "  \"underflow\": " << h.underflow
       << ",\n" << indent << "  \"overflow\": " << h.overflow << "\n" << indent << "}";
}

void writeSignalCsv(const std::string& csv_dir, const std::string& src_name, const SignalResult& sr) {
    if (csv_dir.empty()) return;
    fs::create_directories(csv_dir);
    const auto& primary = (sr.n_events_truth_fiducial > 0) ? sr.truth_hists : sr.reco_hists;
    for (const auto& key : sr.curve_order) {
        const auto& h = primary.at(key);
        double bw = (h.xmax - h.xmin) / h.nbins;
        double tot = static_cast<double>(h.total_in_range);
        std::ofstream f(csv_dir + "/signal_110038_" + src_name + "_" + key + ".csv");
        f << "bin_lo,bin_hi,bin_center,counts,norm,norm_err\n";
        for (int i = 0; i < h.nbins; ++i) {
            double lo = h.xmin + i * bw, hi = lo + bw, c = 0.5 * (lo + hi);
            double norm = (tot > 0.0) ? static_cast<double>(h.counts[i]) / tot : 0.0;
            double err = (tot > 0.0) ? std::sqrt(static_cast<double>(h.counts[i])) / tot : 0.0;
            f << std::setprecision(8) << lo << "," << hi << "," << c << ","
              << h.counts[i] << "," << norm << "," << err << "\n";
        }
    }
}

void writeSourceJson(std::ostream& os, const SignalResult& sr) {
    double den_tot = static_cast<double>(sr.n_events_total);
    double eff_tot = (den_tot > 0.0) ? static_cast<double>(sr.n_events_reco_2trk_os_r95) / den_tot : 0.0;
    double err_tot = (den_tot > 0.0) ? std::sqrt(eff_tot * (1.0 - eff_tot) / den_tot) : 0.0;
    double den_fid = static_cast<double>(sr.n_events_truth_fiducial);
    double eff_fid = (den_fid > 0.0) ? static_cast<double>(sr.n_events_fid_and_reco_c2r) / den_fid : 0.0;
    double err_fid = (den_fid > 0.0) ? std::sqrt(eff_fid * (1.0 - eff_fid) / den_fid) : 0.0;

    const auto& primary = (sr.n_events_truth_fiducial > 0) ? sr.truth_hists : sr.reco_hists;
    os << "    {\n"
       << "      \"input_dir\": \"" << sr.input_dir << "\",\n"
       << "      \"n_files\": " << sr.n_files << ",\n"
       << "      \"n_events_total\": " << sr.n_events_total << ",\n"
       << "      \"n_events_truth_fiducial\": " << sr.n_events_truth_fiducial << ",\n"
       << "      \"n_events_reco_2trk_os_r95\": " << sr.n_events_reco_2trk_os_r95 << ",\n"
       << "      \"n_events_reco_2trk_os_c2r\": " << sr.n_events_reco_2trk_os_c2r << ",\n"
       << "      \"n_events_fid_and_reco_c2r\": " << sr.n_events_fid_and_reco_c2r << ",\n"
       << "      \"two_track_eff_total\": " << std::setprecision(8) << eff_tot << ",\n"
       << "      \"two_track_eff_err\": " << std::setprecision(8) << err_tot << ",\n"
       << "      \"two_track_eff_vs_truth_fid_c2r\": " << std::setprecision(8) << eff_fid << ",\n"
       << "      \"two_track_eff_vs_truth_fid_err\": " << std::setprecision(8) << err_fid << ",\n"
       << "      \"primary_curve_level\": \"" << (sr.n_events_truth_fiducial > 0 ? "truth_st1_st3_fiducial_R95" : "reco_2trk_os_R95") << "\",\n"
       << "      \"sign_convention\": \"pos - neg (e+ minus e-); in FASER +X dipole B-field e- bends +Y and e+ bends -Y so dY(pos-neg)<0 and dY(neg-pos)>0\",\n"
       << "      \"n_truth_dy1_pos_minus_neg_gt0\": " << sr.n_truth_dy1_pos_minus_neg_gt0 << ",\n"
       << "      \"n_truth_dy1_neg_minus_pos_gt0\": " << sr.n_truth_dy1_neg_minus_pos_gt0 << ",\n"
       << "      \"curves\": {\n";
    for (size_t i = 0; i < sr.curve_order.size(); ++i) {
        const auto& k = sr.curve_order[i];
        os << "        \"" << k << "\": ";
        writeHistJson(os, primary.at(k), "        ");
        os << ",\n";
    }
    os << "        \"dY1_signed_neg_minus_pos\": ";
    writeHistJson(os, sr.truth_dy1_neg_minus_pos, "        ");
    os << ",\n        \"dY0_signed_neg_minus_pos\": ";
    writeHistJson(os, sr.truth_dy0_neg_minus_pos, "        ");
    os << "\n      },\n      \"reco_curves\": {\n";
    for (size_t i = 0; i < sr.curve_order.size(); ++i) {
        const auto& k = sr.curve_order[i];
        os << "        \"" << k << "\": ";
        writeHistJson(os, sr.reco_hists.at(k), "        ");
        os << (i + 1 < sr.curve_order.size() ? ",\n" : "\n");
    }
    os << "      },\n      \"truth_dR1_zoom_pass_c2r\": ";
    writeHistJson(os, sr.truth_dr1_pass_c2r, "      ");
    os << "\n    }";
}

int main(int argc, char** argv) {
    std::string master_dev_dir = "/eos/experiment/faser/data0/sim/mc24/foresee/110038/phy/master_dev";
    std::string latest_dir;
    std::string out_json = "signal_110038.json";
    std::string out_csv_dir;
    int threads = 16;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--master-dev" && i + 1 < argc) master_dev_dir = argv[++i];
        else if (arg == "--latest-georot" && i + 1 < argc) latest_dir = argv[++i];
        else if (arg == "--out-json" && i + 1 < argc) out_json = argv[++i];
        else if (arg == "--out-csv-dir" && i + 1 < argc) out_csv_dir = argv[++i];
        else if (arg == "--threads" && i + 1 < argc) threads = std::stoi(argv[++i]);
        else { std::cerr << "Unknown argument: " << arg << "\n"; return 2; }
    }

    ROOT::EnableImplicitMT(threads);
    auto master_res = processSignalDir(master_dev_dir);
    writeSignalCsv(out_csv_dir, "master_dev", master_res);

    bool has_latest = !latest_dir.empty();
    SignalResult latest_res;
    if (has_latest) {
        latest_res = processSignalDir(latest_dir);
        writeSignalCsv(out_csv_dir, "latest_georot", latest_res);
    }

    fs::path out_p(out_json);
    if (out_p.has_parent_path()) fs::create_directories(out_p.parent_path());
    std::ofstream os(out_json);
    os << "{\n"
       << "  \"schema_version\": \"1.0\",\n"
       << "  \"sample_label\": \"A' MC (m = 100 MeV, eps = 1e-5)\",\n"
       << "  \"dsid\": 110038,\n"
       << "  \"sources\": {\n"
       << "    \"master_dev\":\n";
    writeSourceJson(os, master_res);
    os << ",\n    \"latest_georot\": ";
    if (has_latest) { os << "\n"; writeSourceJson(os, latest_res); os << "\n"; }
    else os << "null\n";
    os << "  }\n}\n";
    os.close();

    const auto& dr0 = master_res.truth_hists.at("dR0_wide");
    const auto& dr1 = master_res.truth_hists.at("dR1_wide");
    double p0 = (dr0.total_in_range > 0) ? static_cast<double>(dr0.counts[0]) / dr0.total_in_range : 0.0;
    double p1 = (dr1.total_in_range > 0) ? static_cast<double>(dr1.counts[0]) / dr1.total_in_range : 0.0;
    std::cout << "[overlay_signal_110038] master_dev: files=" << master_res.n_files
              << " total=" << master_res.n_events_total
              << " truth_fid=" << master_res.n_events_truth_fiducial
              << " reco_2trk_c2r=" << master_res.n_events_reco_2trk_os_c2r
              << " eff_vs_fid=" << std::fixed << std::setprecision(4)
              << (master_res.n_events_truth_fiducial > 0 ? static_cast<double>(master_res.n_events_fid_and_reco_c2r) / master_res.n_events_truth_fiducial : 0.0)
              << " truth_dR0[0,0.5]=" << p0 << " truth_dR1[0,0.5]=" << p1
              << " dy1(neg-pos)>0=" << master_res.n_truth_dy1_neg_minus_pos_gt0 << "\n";
    return 0;
}
