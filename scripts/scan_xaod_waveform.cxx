// Fast standalone xAOD waveform scanner for 2022-2026 data and MC.
// Reads VetoNuWaveformHits, VetoWaveformHits, and PreshowerWaveformHits (v1/v2/v3)
// directly via ROOT TTree + setStore and emits a single JSON summary line.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <TBranch.h>
#include <TFile.h>
#include <TTree.h>

#include "xAODCore/AuxContainerBase.h"
#include "xAODFaserWaveform/versions/WaveformHitAuxContainer_v1.h"
#include "xAODFaserWaveform/versions/WaveformHitAuxContainer_v2.h"
#include "xAODFaserWaveform/versions/WaveformHitAuxContainer_v3.h"
#include "xAODFaserWaveform/versions/WaveformHitContainer_v1.h"
#include "xAODFaserWaveform/versions/WaveformHitContainer_v2.h"
#include "xAODFaserWaveform/versions/WaveformHitContainer_v3.h"

struct HitLite {
  int ch;
  unsigned int st;
  double lt;
  double raw_q;
  bool has_wave;
};

struct WaveHandle {
  TBranch* b_vec = nullptr;
  TBranch* b_aux = nullptr;
  xAOD::WaveformHitContainer_v1* vec1 = nullptr;
  xAOD::WaveformHitContainer_v2* vec2 = nullptr;
  xAOD::WaveformHitContainer_v3* vec3 = nullptr;
  xAOD::WaveformHitAuxContainer_v1* aux1 = nullptr;
  xAOD::WaveformHitAuxContainer_v2* aux2 = nullptr;
  xAOD::WaveformHitAuxContainer_v3* aux3 = nullptr;
  int ver_v = 0;
  int ver_a = 0;
  std::string cls;

  bool bind(TTree* t, const char* name) {
    b_vec = t->GetBranch(name);
    std::string aux_name = std::string(name) + "Aux.";
    b_aux = t->GetBranch(aux_name.c_str());
    if (!b_vec || !b_aux) return false;
    std::string vcls = b_vec->GetClassName();
    cls = b_aux->GetClassName();
    t->SetBranchStatus(name, 1);
    t->SetBranchStatus(aux_name.c_str(), 1);
    if (vcls.find("_v3") != std::string::npos) {
      ver_v = 3;
      t->SetBranchAddress(name, &vec3);
    } else if (vcls.find("_v2") != std::string::npos) {
      ver_v = 2;
      t->SetBranchAddress(name, &vec2);
    } else {
      ver_v = 1;
      t->SetBranchAddress(name, &vec1);
    }
    if (cls.find("_v3") != std::string::npos) {
      ver_a = 3;
      t->SetBranchAddress(aux_name.c_str(), &aux3);
    } else if (cls.find("_v2") != std::string::npos) {
      ver_a = 2;
      t->SetBranchAddress(aux_name.c_str(), &aux2);
    } else {
      ver_a = 1;
      t->SetBranchAddress(aux_name.c_str(), &aux1);
    }
    return true;
  }

  void get_hits(long long ie, std::vector<HitLite>& out) {
    out.clear();
    if (!b_vec || !b_aux) return;
    b_vec->GetEntry(ie);
    b_aux->GetEntry(ie);
    const SG::IConstAuxStore* store = nullptr;
    if (ver_a == 3) store = aux3;
    else if (ver_a == 2) store = aux2;
    else store = aux1;
    if (!store) return;
    auto fill = [&](auto* vec) {
      if (!vec) return;
      vec->setStore(store);
      const auto& cvec = *vec;
      for (size_t k = 0; k < cvec.size(); ++k) {
        const auto* h = cvec.at(k);
        out.push_back({(int)h->channel(), (unsigned int)h->hit_status(),
                       (double)h->localtime(), (double)h->raw_integral() / 50.0,
                       !h->wave_vector().empty()});
      }
    };
    if (ver_v == 3) fill(vec3);
    else if (ver_v == 2) fill(vec2);
    else fill(vec1);
  }
};

static bool has_val(const std::vector<int>& v, int x) {
  return std::find(v.begin(), v.end(), x) != v.end();
}

struct PeakStat {
  long long n = 0;
  bool valid = false;
  double peak = 0.0;
  double median = 0.0;
  double core_std = 0.0;
};

static PeakStat compute_peak(std::vector<double>& v, double bin_w = 0.2, double win = 2.0) {
  PeakStat s;
  v.erase(std::remove_if(v.begin(), v.end(), [](double x) { return !std::isfinite(x); }), v.end());
  s.n = (long long)v.size();
  if (s.n < 20) return s;
  std::sort(v.begin(), v.end());
  s.median = v[s.n / 2];
  double lo = s.median - 25.0, hi = s.median + 25.0;
  int nb = std::max(10, (int)std::llround((hi - lo) / bin_w));
  std::vector<double> cnt(nb, 0.0);
  for (double x : v) {
    if (x >= lo && x < hi) {
      int b = (int)((x - lo) / bin_w);
      if (b >= 0 && b < nb) cnt[b] += 1.0;
    }
  }
  int imax = (int)(std::max_element(cnt.begin(), cnt.end()) - cnt.begin());
  if (cnt[imax] < 5.0) return s;
  double mode = lo + (imax + 0.5) * bin_w;
  // Weighted quadratic fit to log(counts) in [imax-4, imax+4]
  int i0 = std::max(0, imax - 4);
  int i1 = std::min(nb - 1, imax + 4);
  // 3x3 normal equations for y = a*x^2 + b*x + c with x = center - mode
  double S[5] = {0, 0, 0, 0, 0};
  double Y[3] = {0, 0, 0};
  int np = 0;
  for (int i = i0; i <= i1; ++i) {
    if (cnt[i] <= std::max(3.0, 0.25 * cnt[imax])) continue;
    double x = (lo + (i + 0.5) * bin_w) - mode;
    double y = std::log(cnt[i]);
    double w = cnt[i];
    double x2 = x * x;
    S[0] += w; S[1] += w * x; S[2] += w * x2; S[3] += w * x2 * x; S[4] += w * x2 * x2;
    Y[0] += w * y; Y[1] += w * x * y; Y[2] += w * x2 * y;
    np++;
  }
  double peak = mode;
  if (np >= 3) {
    // Solve symmetric 3x3: [S4 S3 S2; S3 S2 S1; S2 S1 S0] [a; b; c] = [Y2; Y1; Y0]
    double A[3][4] = {{S[4], S[3], S[2], Y[2]}, {S[3], S[2], S[1], Y[1]}, {S[2], S[1], S[0], Y[0]}};
    for (int p = 0; p < 3; ++p) {
      int piv = p;
      for (int r = p + 1; r < 3; ++r) if (std::fabs(A[r][p]) > std::fabs(A[piv][p])) piv = r;
      if (std::fabs(A[piv][p]) < 1e-12) { piv = -1; break; }
      for (int c = p; c < 4; ++c) std::swap(A[p][c], A[piv][c]);
      for (int r = 0; r < 3; ++r) {
        if (r == p) continue;
        double f = A[r][p] / A[p][p];
        for (int c = p; c < 4; ++c) A[r][c] -= f * A[p][c];
      }
    }
    double a = (std::fabs(A[0][0]) > 1e-12) ? A[0][3] / A[0][0] : 0.0;
    double b = (std::fabs(A[1][1]) > 1e-12) ? A[1][3] / A[1][1] : 0.0;
    if (a < -1e-6) {
      double dx = -b / (2.0 * a);
      if (std::fabs(dx) <= 1.5) peak = mode + dx;
    }
  }
  double sum = 0.0, sum2 = 0.0;
  long long nc = 0;
  for (double x : v) {
    if (std::fabs(x - peak) <= win) {
      sum += x; sum2 += x * x; nc++;
    }
  }
  s.valid = true;
  s.peak = peak;
  if (nc >= 10) {
    double m = sum / nc;
    s.core_std = std::sqrt(std::max(0.0, sum2 / nc - m * m));
  }
  return s;
}

static void print_stat(const char* name, const PeakStat& s, bool comma = true) {
  std::cout << "\"" << name << "\":{\"n\":" << s.n << ",\"peak\":";
  if (s.valid) std::cout << s.peak << ",\"median\":" << s.median << ",\"core_std\":" << s.core_std;
  else std::cout << "null,\"median\":null,\"core_std\":null";
  std::cout << "}" << (comma ? "," : "");
}

static void print_ivec(const char* name, std::vector<int> v, bool comma = true) {
  std::sort(v.begin(), v.end());
  std::cout << "\"" << name << "\":[";
  for (size_t i = 0; i < v.size(); ++i) std::cout << (i ? "," : "") << v[i];
  std::cout << "]" << (comma ? "," : "");
}

int main(int argc, char** argv) {
  if (argc < 2) return 2;
  const char* path = argv[1];
  long long max_events = (argc >= 3) ? std::atoll(argv[2]) : 40000;

  std::unique_ptr<TFile> f(TFile::Open(path, "READ"));
  if (!f || f->IsZombie()) {
    std::cout << "{\"error\":\"open_failed\"}\n";
    return 0;
  }
  TTree* t = dynamic_cast<TTree*>(f->Get("CollectionTree"));
  if (!t) {
    std::cout << "{\"error\":\"no_CollectionTree\"}\n";
    return 0;
  }
  t->SetBranchStatus("*", 0);
  WaveHandle h_vn, h_v2, h_ps;
  if (!h_vn.bind(t, "VetoNuWaveformHits")) {
    std::cout << "{\"error\":\"no_VetoNuWaveformHits\"}\n";
    return 0;
  }
  bool has_v2 = h_v2.bind(t, "VetoWaveformHits");
  bool has_ps = h_ps.bind(t, "PreshowerWaveformHits");

  long long n = t->GetEntries();
  if (max_events > 0 && n > max_events) n = max_events;

  long long n_vn_hits = 0, n_vn_nonempty = 0, n_vn_sec = 0, n_vn_ge3 = 0, n_vn_swap = 0, n_vn_512 = 0;
  long long n_v2_0 = 0, n_v2_512 = 0, n_pst_0 = 0, n_pst_512 = 0, n_psb_0 = 0, n_psb_512 = 0;
  std::vector<int> vn_ch, v2_ch, ps_ch;
  std::vector<double> dt_v2, dt_v2_0, dt_v2_512c, dt_v2_512r;
  std::vector<double> dt_ps, dt_pst, dt_psb, dt_psb_0, dt_psb_512c, dt_psb_512r;

  std::vector<HitLite> vn, v2, ps;
  for (long long ie = 0; ie < n; ++ie) {
    h_vn.get_hits(ie, vn);
    if (vn.empty()) continue;
    if (vn.size() >= 3) n_vn_ge3++;
    int i4 = -1, i5 = -1;
    double vn_sum = 0.0;
    int vn_cnt = 0;
    bool vn_hq[2] = {false, false};
    for (size_t i = 0; i < vn.size(); ++i) {
      const auto& h = vn[i];
      n_vn_hits++;
      if (!has_val(vn_ch, h.ch)) vn_ch.push_back(h.ch);
      if (h.ch == 4 && i4 < 0) i4 = (int)i;
      if (h.ch == 5 && i5 < 0) i5 = (int)i;
      if (h.has_wave) n_vn_nonempty++;
      if (h.st & 2U) n_vn_sec++;
      if (h.st & 512U) n_vn_512++;
      if ((h.st & ~512U) == 0 && (h.ch == 4 || h.ch == 5) && h.raw_q > 40.0) {
        vn_hq[h.ch == 4 ? 0 : 1] = true;
        vn_sum += h.lt + ((h.st & 512U) ? 16.0 : 0.0);
        vn_cnt++;
      }
    }
    if (i4 >= 0 && i5 >= 0 && i5 < i4) n_vn_swap++;
    if (vn_cnt < 2 || !vn_hq[0] || !vn_hq[1]) continue;
    double t_vn = vn_sum / vn_cnt;

    if (has_v2) {
      h_v2.get_hits(ie, v2);
      double sa = 0, s0 = 0, s5c = 0, s5r = 0;
      int na = 0, n0 = 0, n5 = 0;
      for (const auto& h : v2) {
        if (!has_val(v2_ch, h.ch)) v2_ch.push_back(h.ch);
        if (h.ch != 6 && h.ch != 7 && h.ch != 22 && h.ch != 23) continue;
        if ((h.st & ~512U) != 0) continue;
        double tc = h.lt + ((h.st & 512U) ? 16.0 : 0.0);
        sa += tc; na++;
        if (h.st == 0) { s0 += h.lt; n0++; n_v2_0++; }
        else if (h.st == 512U) { s5c += tc; s5r += h.lt; n5++; n_v2_512++; }
      }
      if (na > 0) dt_v2.push_back(sa / na - t_vn);
      if (n0 > 0) dt_v2_0.push_back(s0 / n0 - t_vn);
      if (n5 > 0) { dt_v2_512c.push_back(s5c / n5 - t_vn); dt_v2_512r.push_back(s5r / n5 - t_vn); }
    }

    if (has_ps) {
      h_ps.get_hits(ie, ps);
      double sa = 0, stp = 0, sbm = 0, sb0 = 0, sb5c = 0, sb5r = 0;
      int na = 0, ntp = 0, nbm = 0, nb0 = 0, nb5 = 0;
      for (const auto& h : ps) {
        if (!has_val(ps_ch, h.ch)) ps_ch.push_back(h.ch);
        if ((h.st & ~512U) != 0) continue;
        double tc = h.lt + ((h.st & 512U) ? 16.0 : 0.0);
        sa += tc; na++;
        if (h.ch == 12 || h.ch == 13) {
          stp += tc; ntp++;
          if (h.st == 0) n_pst_0++; else if (h.st == 512U) n_pst_512++;
        } else if (h.ch == 24 || h.ch == 25) {
          sbm += tc; nbm++;
          if (h.st == 0) { sb0 += h.lt; nb0++; n_psb_0++; }
          else if (h.st == 512U) { sb5c += tc; sb5r += h.lt; nb5++; n_psb_512++; }
        }
      }
      if (na > 0) dt_ps.push_back(sa / na - t_vn);
      if (ntp > 0) dt_pst.push_back(stp / ntp - t_vn);
      if (nbm > 0) dt_psb.push_back(sbm / nbm - t_vn);
      if (nb0 > 0) dt_psb_0.push_back(sb0 / nb0 - t_vn);
      if (nb5 > 0) { dt_psb_512c.push_back(sb5c / nb5 - t_vn); dt_psb_512r.push_back(sb5r / nb5 - t_vn); }
    }
  }

  std::cout << "{"
            << "\"aux_class\":\"" << h_vn.cls << "\","
            << "\"n_events\":" << n << ","
            << "\"n_vn_hits\":" << n_vn_hits << ","
            << "\"n_vn_nonempty_wave\":" << n_vn_nonempty << ","
            << "\"n_vn_secondary\":" << n_vn_sec << ","
            << "\"n_vn_ge3_hits\":" << n_vn_ge3 << ","
            << "\"n_vn_ch5_before_ch4\":" << n_vn_swap << ","
            << "\"n_vn_512\":" << n_vn_512 << ","
            << "\"n_v2_st0\":" << n_v2_0 << ","
            << "\"n_v2_st512\":" << n_v2_512 << ","
            << "\"n_ps_top_st0\":" << n_pst_0 << ","
            << "\"n_ps_top_st512\":" << n_pst_512 << ","
            << "\"n_ps_bot_st0\":" << n_psb_0 << ","
            << "\"n_ps_bot_st512\":" << n_psb_512 << ",";
  print_ivec("vn_channels", vn_ch);
  print_ivec("v2_channels", v2_ch);
  print_ivec("ps_channels", ps_ch);
  print_stat("veto2_minus_vetonu", compute_peak(dt_v2));
  print_stat("veto2_st0_minus_vetonu", compute_peak(dt_v2_0));
  print_stat("veto2_st512_corr_minus_vetonu", compute_peak(dt_v2_512c));
  print_stat("veto2_st512_raw_minus_vetonu", compute_peak(dt_v2_512r));
  print_stat("ps_all_minus_vetonu", compute_peak(dt_ps));
  print_stat("ps_top_minus_vetonu", compute_peak(dt_pst));
  print_stat("ps_bot_minus_vetonu", compute_peak(dt_psb));
  print_stat("ps_bot_st0_minus_vetonu", compute_peak(dt_psb_0));
  print_stat("ps_bot_st512_corr_minus_vetonu", compute_peak(dt_psb_512c));
  print_stat("ps_bot_st512_raw_minus_vetonu", compute_peak(dt_psb_512r), false);
  std::cout << "}\n";
  return 0;
}
