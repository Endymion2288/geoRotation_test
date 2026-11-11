"""Small independent final-track observer; coordinates are written in the ACTS frame."""
import json
import math

import AthenaPython.PyAthena as PyAthena


def vector(v):
    return [float(v.x()), float(v.y()), float(v.z())]


class CompactEventWriter(PyAthena.Alg):
    def __init__(self, name, output):
        super().__init__(name=name)
        self.output = output

    def initialize(self):
        self.sg = PyAthena.py_svc("StoreGateSvc")
        self.stream = open(self.output, "w", buffering=1024 * 1024)
        return PyAthena.StatusCode.Success

    def execute(self):
        ei = self.sg.retrieve("xAOD::EventInfo", "EventInfo")
        truth = self.sg.retrieve("xAOD::TruthParticleContainer", "TruthParticles")
        hits = self.sg.retrieve("FaserSiHitCollection", "SCT_Hits")
        tracks = self.sg.retrieve("TrackCollection", "CKFTrackCollectionWithoutIFT")
        particles = []
        hit_counts = {}
        for hit in hits:
            barcode = int(hit.trackNumber())
            station = int(hit.getStation())
            counts = hit_counts.setdefault(barcode, [0] * 4)
            if 0 <= station < 4:
                counts[station] += 1
        for p in truth:
            if abs(p.pdgId()) not in (11, 13):
                continue
            barcode = int(p.barcode())
            # Keep generated particles, their G4 continuation vertices, and hit
            # secondaries; omit unrelated shower particles with no tracker hit.
            if barcode % 1000000 >= 200000 and barcode not in hit_counts:
                continue
            particles.append(dict(barcode=int(p.barcode()), pdg=int(p.pdgId()), status=int(p.status()),
                                  signal=0 < barcode < 200000 and int(p.status()) == 1,
                                  charge=float(p.charge()), p_MeV=math.sqrt(p.px()**2 + p.py()**2 + p.pz()**2),
                                  momentum_F_MeV=[float(p.px()), float(p.py()), float(p.pz())],
                                  vertex_F_mm=vector(p.prodVtx()) if p.hasProdVtx() else None,
                                  hit_counts=hit_counts.get(int(p.barcode()), [0] * 4)))
        rows = []
        for index, track in enumerate(tracks):
            parameters = list(track.trackParameters())
            if not parameters:
                rows.append(dict(index=index, valid=False))
                continue
            p = parameters[0]
            fq = track.fitQuality()
            rows.append(dict(index=index, valid=True, position_A_mm=vector(p.position()),
                             momentum_A_MeV=vector(p.momentum()), charge=float(p.charge()),
                             chi2=float(fq.chiSquared()), ndof=float(fq.numberDoF()),
                             measurements=int(track.measurementsOnTrack().size()),
                             rmax_mm=max(math.hypot(q.position().y(), q.position().z()) for q in parameters)))
        record = dict(run=int(ei.runNumber()), event=int(ei.eventNumber()), lumi=int(ei.lumiBlock()), truth=particles,
                      n_si_hits=int(hits.size()), tracks=rows)
        self.stream.write(json.dumps(record, allow_nan=False, separators=(",", ":")) + "\n")
        return PyAthena.StatusCode.Success

    def finalize(self):
        self.stream.close()
        return PyAthena.StatusCode.Success
