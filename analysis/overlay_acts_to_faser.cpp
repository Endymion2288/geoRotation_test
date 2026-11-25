// Fast C++17 frame normalizer & verifier for TwoTrackOverlayNTuplerAlg output.
//
// Applies FaserActs::actsToFaser (x_F = -z_A, y_F = y_A, z_F = x_A) exclusively
// to the 12 reconstructed Track_{x,y,z,px,py,pz}{0,1} branches when --rotate 1
// is set, copying all remaining 36 branches (og_Track_*, TrackSegment_*,
// GhostBustedTrackSegment_*, Track_Chi2, Track_nDoF, longTracks, etc.)
// bit-for-bit. Attaches self-describing TNamed provenance metadata and enforces
// hard coordinate-frame & completeness assertions.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TNamed.h"
#include "TTree.h"

namespace {

struct Config {
  std::string input_path;
  std::string output_path;
  std::string rdo_path;
  std::string stats_json_path;
  int64_t expected_entries = -1;
  int rotate = 1;
  std::string meta_frame = "FASER";
  std::string meta_calypso_commit;
  std::string meta_cond_tag;
  std::string meta_material_map;
  std::string meta_source_rdo;
  std::string meta_dataset;
  std::string meta_pos_bin;
  std::string meta_neg_bin;
};

std::string json_escape(const std::string& s) {
  std::ostringstream oss;
  for (char c : s) {
    switch (c) {
      case '"': oss << "\\\""; break;
      case '\\': oss << "\\\\"; break;
      case '\n': oss << "\\n"; break;
      case '\r': oss << "\\r"; break;
      case '\t': oss << "\\t"; break;
      default: oss << c; break;
    }
  }
  return oss.str();
}

double compute_median(std::vector<double>& vals) {
  if (vals.empty()) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const size_t n = vals.size();
  const size_t mid = n / 2;
  std::nth_element(vals.begin(), vals.begin() + mid, vals.end());
  double m = vals[mid];
  if (n % 2 == 0) {
    const double lower = *std::max_element(vals.begin(), vals.begin() + mid);
    m = 0.5 * (lower + m);
  }
  return m;
}

std::string format_double_json(double v) {
  if (!std::isfinite(v)) {
    return "null";
  }
  std::ostringstream oss;
  oss << std::setprecision(10) << v;
  return oss.str();
}

bool parse_args(int argc, char** argv, Config& cfg) {
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    auto need_val = [&](const std::string& name) -> std::string {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << name << "\n";
        std::exit(2);
      }
      return argv[++i];
    };
    if (arg == "--input") {
      cfg.input_path = need_val(arg);
    } else if (arg == "--output") {
      cfg.output_path = need_val(arg);
    } else if (arg == "--rdo") {
      cfg.rdo_path = need_val(arg);
    } else if (arg == "--expected-entries") {
      cfg.expected_entries = std::stoll(need_val(arg));
    } else if (arg == "--rotate") {
      cfg.rotate = std::stoi(need_val(arg));
    } else if (arg == "--stats-json") {
      cfg.stats_json_path = need_val(arg);
    } else if (arg == "--meta-frame") {
      cfg.meta_frame = need_val(arg);
    } else if (arg == "--meta-calypso-commit") {
      cfg.meta_calypso_commit = need_val(arg);
    } else if (arg == "--meta-cond-tag") {
      cfg.meta_cond_tag = need_val(arg);
    } else if (arg == "--meta-material-map") {
      cfg.meta_material_map = need_val(arg);
    } else if (arg == "--meta-source-rdo") {
      cfg.meta_source_rdo = need_val(arg);
    } else if (arg == "--meta-dataset") {
      cfg.meta_dataset = need_val(arg);
    } else if (arg == "--meta-pos-bin") {
      cfg.meta_pos_bin = need_val(arg);
    } else if (arg == "--meta-neg-bin") {
      cfg.meta_neg_bin = need_val(arg);
    } else {
      std::cerr << "Unknown argument: " << arg << "\n";
      return false;
    }
  }
  if (cfg.input_path.empty() || cfg.output_path.empty()) {
    std::cerr << "Usage: overlay_acts_to_faser --input <raw.root> --output <faser.root> "
                 "[--rdo <rdo.root>] [--expected-entries N] [--rotate 0|1] "
                 "[--stats-json <stats.json>] [--meta-* ...]\n";
    return false;
  }
  if (cfg.rotate != 0 && cfg.rotate != 1) {
    std::cerr << "--rotate must be 0 or 1, got " << cfg.rotate << "\n";
    return false;
  }
  return true;
}

}  // namespace

int main(int argc, char** argv) {
  Config cfg;
  if (!parse_args(argc, argv, cfg)) {
    return 2;
  }

  int64_t expected_entries = cfg.expected_entries;
  if (expected_entries < 0 && !cfg.rdo_path.empty()) {
    std::unique_ptr<TFile> rdo_file(TFile::Open(cfg.rdo_path.c_str(), "READ"));
    if (!rdo_file || rdo_file->IsZombie()) {
      std::cerr << "ERROR: Failed to open input RDO file: " << cfg.rdo_path << "\n";
      return 3;
    }
    TTree* col_tree = dynamic_cast<TTree*>(rdo_file->Get("CollectionTree"));
    if (!col_tree) {
      std::cerr << "ERROR: Missing CollectionTree in RDO file: " << cfg.rdo_path << "\n";
      return 3;
    }
    expected_entries = static_cast<int64_t>(col_tree->GetEntries());
  }

  std::unique_ptr<TFile> in_file(TFile::Open(cfg.input_path.c_str(), "READ"));
  if (!in_file || in_file->IsZombie()) {
    std::cerr << "ERROR: Failed to open input raw ntuple: " << cfg.input_path << "\n";
    return 4;
  }
  TTree* in_tree = dynamic_cast<TTree*>(in_file->Get("nt"));
  if (!in_tree) {
    std::cerr << "ERROR: Missing TTree 'nt' in " << cfg.input_path << "\n";
    return 4;
  }

  const int64_t n_entries = static_cast<int64_t>(in_tree->GetEntries());

  int long_tracks = 0;
  int og_long_tracks = 0;
  std::vector<double>* x0 = nullptr;
  std::vector<double>* y0 = nullptr;
  std::vector<double>* z0 = nullptr;
  std::vector<double>* px0 = nullptr;
  std::vector<double>* py0 = nullptr;
  std::vector<double>* pz0 = nullptr;
  std::vector<double>* x1 = nullptr;
  std::vector<double>* y1 = nullptr;
  std::vector<double>* z1 = nullptr;
  std::vector<double>* px1 = nullptr;
  std::vector<double>* py1 = nullptr;
  std::vector<double>* pz1 = nullptr;

  auto bind_branch = [&](const char* name, auto* ptr) {
    if (!in_tree->GetBranch(name)) {
      std::cerr << "ERROR: Missing required branch '" << name << "' in input nt tree\n";
      std::exit(5);
    }
    in_tree->SetBranchAddress(name, ptr);
  };

  bind_branch("longTracks", &long_tracks);
  bind_branch("og_longTracks", &og_long_tracks);
  bind_branch("Track_x0", &x0);
  bind_branch("Track_y0", &y0);
  bind_branch("Track_z0", &z0);
  bind_branch("Track_px0", &px0);
  bind_branch("Track_py0", &py0);
  bind_branch("Track_pz0", &pz0);
  bind_branch("Track_x1", &x1);
  bind_branch("Track_y1", &y1);
  bind_branch("Track_z1", &z1);
  bind_branch("Track_px1", &px1);
  bind_branch("Track_py1", &py1);
  bind_branch("Track_pz1", &pz1);

  std::unique_ptr<TFile> out_file(TFile::Open(cfg.output_path.c_str(), "RECREATE"));
  if (!out_file || out_file->IsZombie()) {
    std::cerr << "ERROR: Failed to create output ROOT file: " << cfg.output_path << "\n";
    return 6;
  }
  out_file->SetCompressionLevel(1);
  TTree* out_tree = in_tree->CloneTree(0);

  int64_t lt_counts[5] = {0, 0, 0, 0, 0};  // 0, 1, 2, 3, >=4
  int64_t og_not_two = 0;
  int64_t total_tracks = 0;
  int64_t non_finite_tracks = 0;
  int64_t abs_z1_lt_500 = 0;
  int64_t in_faser_bounds = 0;
  int64_t out_z0_bounds = 0;
  int64_t out_z1_bounds = 0;
  int64_t non_pos_pz0 = 0;
  int64_t non_pos_pz1 = 0;

  double min_z0 = std::numeric_limits<double>::infinity();
  double max_z0 = -std::numeric_limits<double>::infinity();
  double min_z1 = std::numeric_limits<double>::infinity();
  double max_z1 = -std::numeric_limits<double>::infinity();
  double min_pz0 = std::numeric_limits<double>::infinity();
  double max_pz0 = -std::numeric_limits<double>::infinity();
  double min_pz1 = std::numeric_limits<double>::infinity();
  double max_pz1 = -std::numeric_limits<double>::infinity();

  std::vector<double> all_z0;
  std::vector<double> all_z1;
  all_z0.reserve(static_cast<size_t>(std::max<int64_t>(0, n_entries) * 2));
  all_z1.reserve(static_cast<size_t>(std::max<int64_t>(0, n_entries) * 2));

  for (int64_t i = 0; i < n_entries; ++i) {
    in_tree->GetEntry(i);

    if (og_long_tracks != 2) {
      ++og_not_two;
    }
    const int lt_bin = (long_tracks <= 0) ? 0 : (long_tracks >= 4 ? 4 : long_tracks);
    ++lt_counts[lt_bin];

    const size_t ntrk = x0->size();
    if (y0->size() != ntrk || z0->size() != ntrk || px0->size() != ntrk ||
        py0->size() != ntrk || pz0->size() != ntrk || x1->size() != ntrk ||
        y1->size() != ntrk || z1->size() != ntrk || px1->size() != ntrk ||
        py1->size() != ntrk || pz1->size() != ntrk) {
      std::cerr << "ERROR: Track vector size mismatch at entry " << i << "\n";
      return 7;
    }

    for (size_t k = 0; k < ntrk; ++k) {
      if (cfg.rotate == 1) {
        // Exact FaserActs::actsToFaser: (x_F, y_F, z_F) = (-z_A, y_A, x_A)
        const double ax0 = (*x0)[k];
        const double az0 = (*z0)[k];
        (*x0)[k] = -az0;
        (*z0)[k] = ax0;

        const double apx0 = (*px0)[k];
        const double apz0 = (*pz0)[k];
        (*px0)[k] = -apz0;
        (*pz0)[k] = apx0;

        const double ax1 = (*x1)[k];
        const double az1 = (*z1)[k];
        (*x1)[k] = -az1;
        (*z1)[k] = ax1;

        const double apx1 = (*px1)[k];
        const double apz1 = (*pz1)[k];
        (*px1)[k] = -apz1;
        (*pz1)[k] = apx1;
      }

      const double fx0 = (*x0)[k];
      const double fy0 = (*y0)[k];
      const double fz0 = (*z0)[k];
      const double fpx0 = (*px0)[k];
      const double fpy0 = (*py0)[k];
      const double fpz0 = (*pz0)[k];
      const double fx1 = (*x1)[k];
      const double fy1 = (*y1)[k];
      const double fz1 = (*z1)[k];
      const double fpx1 = (*px1)[k];
      const double fpy1 = (*py1)[k];
      const double fpz1 = (*pz1)[k];

      ++total_tracks;
      const bool finite = std::isfinite(fx0) && std::isfinite(fy0) && std::isfinite(fz0) &&
                          std::isfinite(fpx0) && std::isfinite(fpy0) && std::isfinite(fpz0) &&
                          std::isfinite(fx1) && std::isfinite(fy1) && std::isfinite(fz1) &&
                          std::isfinite(fpx1) && std::isfinite(fpy1) && std::isfinite(fpz1);
      if (!finite) {
        ++non_finite_tracks;
        continue;
      }

      all_z0.push_back(fz0);
      all_z1.push_back(fz1);
      min_z0 = std::min(min_z0, fz0);
      max_z0 = std::max(max_z0, fz0);
      min_z1 = std::min(min_z1, fz1);
      max_z1 = std::max(max_z1, fz1);
      min_pz0 = std::min(min_pz0, fpz0);
      max_pz0 = std::max(max_pz0, fpz0);
      min_pz1 = std::min(min_pz1, fpz1);
      max_pz1 = std::max(max_pz1, fpz1);

      if (std::abs(fz1) < 500.0) {
        ++abs_z1_lt_500;
      }
      const bool ok_z0 = (fz0 >= -100.0 && fz0 <= 200.0);
      const bool ok_z1 = (fz1 >= 2200.0 && fz1 <= 2650.0);
      const bool ok_pz0 = (fpz0 > 0.0);
      const bool ok_pz1 = (fpz1 > 0.0);
      if (!ok_z0) ++out_z0_bounds;
      if (!ok_z1) ++out_z1_bounds;
      if (!ok_pz0) ++non_pos_pz0;
      if (!ok_pz1) ++non_pos_pz1;
      if (ok_z0 && ok_z1 && ok_pz0 && ok_pz1) {
        ++in_faser_bounds;
      }
    }

    out_tree->Fill();
  }

  out_file->cd();
  out_tree->Write("", TObject::kOverwrite);

  auto write_meta = [&](const char* key, const std::string& val) {
    TNamed named(key, val.c_str());
    named.Write("", TObject::kOverwrite);
  };
  write_meta("meta_frame", cfg.meta_frame);
  write_meta("meta_calypso_commit", cfg.meta_calypso_commit);
  write_meta("meta_cond_tag", cfg.meta_cond_tag);
  write_meta("meta_material_map", cfg.meta_material_map);
  write_meta("meta_source_rdo", cfg.meta_source_rdo);
  write_meta("meta_dataset", cfg.meta_dataset);
  write_meta("meta_pos_bin", cfg.meta_pos_bin);
  write_meta("meta_neg_bin", cfg.meta_neg_bin);
  out_file->Close();
  in_file->Close();

  const double med_z0 = compute_median(all_z0);
  const double med_z1 = compute_median(all_z1);
  const double in_bounds_frac =
      (total_tracks > 0) ? (static_cast<double>(in_faser_bounds) / static_cast<double>(total_tracks)) : 1.0;

  std::vector<std::string> errors;
  if (expected_entries >= 0 && n_entries != expected_entries) {
    std::ostringstream oss;
    oss << "Entry count mismatch: nt_entries=" << n_entries << " != expected=" << expected_entries;
    errors.push_back(oss.str());
  }
  if (og_not_two != 0) {
    std::ostringstream oss;
    oss << "og_longTracks != 2 on " << og_not_two << " / " << n_entries << " events";
    errors.push_back(oss.str());
  }
  if (non_finite_tracks != 0) {
    std::ostringstream oss;
    oss << "Found " << non_finite_tracks << " tracks with non-finite coordinates/momenta";
    errors.push_back(oss.str());
  }
  if (total_tracks > 0) {
    if (!(med_z0 >= -20.0 && med_z0 <= 80.0)) {
      std::ostringstream oss;
      oss << "Median Track_z0=" << med_z0 << " outside [-20, 80] mm";
      errors.push_back(oss.str());
    }
    if (!(med_z1 >= 2350.0 && med_z1 <= 2500.0)) {
      std::ostringstream oss;
      oss << "Median Track_z1=" << med_z1 << " outside [2350, 2500] mm";
      errors.push_back(oss.str());
    }
    if (abs_z1_lt_500 != 0) {
      std::ostringstream oss;
      oss << "Found " << abs_z1_lt_500 << " tracks with |Track_z1| < 500 mm (unrotated ACTS frame?)";
      errors.push_back(oss.str());
    }
    if (in_bounds_frac < 0.999) {
      std::ostringstream oss;
      oss << "FASER bounds fraction " << in_bounds_frac << " < 0.999 ("
          << in_faser_bounds << "/" << total_tracks << ")";
      errors.push_back(oss.str());
    }
  }

  const bool passed = errors.empty();

  std::ostringstream js;
  js << "{\n"
     << "  \"passed\": " << (passed ? "true" : "false") << ",\n"
     << "  \"rotate\": " << cfg.rotate << ",\n"
     << "  \"meta_frame\": \"" << json_escape(cfg.meta_frame) << "\",\n"
     << "  \"nt_entries\": " << n_entries << ",\n"
     << "  \"expected_entries\": " << expected_entries << ",\n"
     << "  \"og_longTracks_not_2\": " << og_not_two << ",\n"
     << "  \"longTracks_counts\": {\n"
     << "    \"0\": " << lt_counts[0] << ",\n"
     << "    \"1\": " << lt_counts[1] << ",\n"
     << "    \"2\": " << lt_counts[2] << ",\n"
     << "    \"3\": " << lt_counts[3] << ",\n"
     << "    \"4plus\": " << lt_counts[4] << "\n"
     << "  },\n"
     << "  \"total_tracks\": " << total_tracks << ",\n"
     << "  \"non_finite_tracks\": " << non_finite_tracks << ",\n"
     << "  \"abs_z1_lt_500\": " << abs_z1_lt_500 << ",\n"
     << "  \"in_faser_bounds\": " << in_faser_bounds << ",\n"
     << "  \"in_faser_bounds_fraction\": " << format_double_json(in_bounds_frac) << ",\n"
     << "  \"outliers\": {\n"
     << "    \"z0_outside_m100_200\": " << out_z0_bounds << ",\n"
     << "    \"z1_outside_2200_2650\": " << out_z1_bounds << ",\n"
     << "    \"pz0_le_0\": " << non_pos_pz0 << ",\n"
     << "    \"pz1_le_0\": " << non_pos_pz1 << "\n"
     << "  },\n"
     << "  \"track_z0_mm\": {\"min\": " << format_double_json(min_z0)
     << ", \"median\": " << format_double_json(med_z0)
     << ", \"max\": " << format_double_json(max_z0) << "},\n"
     << "  \"track_z1_mm\": {\"min\": " << format_double_json(min_z1)
     << ", \"median\": " << format_double_json(med_z1)
     << ", \"max\": " << format_double_json(max_z1) << "},\n"
     << "  \"track_pz0_mev\": {\"min\": " << format_double_json(min_pz0)
     << ", \"max\": " << format_double_json(max_pz0) << "},\n"
     << "  \"track_pz1_mev\": {\"min\": " << format_double_json(min_pz1)
     << ", \"max\": " << format_double_json(max_pz1) << "},\n"
     << "  \"errors\": [";
  for (size_t i = 0; i < errors.size(); ++i) {
    if (i > 0) js << ", ";
    js << "\"" << json_escape(errors[i]) << "\"";
  }
  js << "]\n}\n";

  if (!cfg.stats_json_path.empty()) {
    std::ofstream ofs(cfg.stats_json_path);
    ofs << js.str();
  }
  std::cout << js.str();

  if (!passed) {
    for (const auto& e : errors) {
      std::cerr << "ASSERTION FAILED: " << e << "\n";
    }
    return 10;
  }
  return 0;
}
