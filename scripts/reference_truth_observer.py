"""Read native Calypso production truth and output order without altering tracking."""
import json
import math
import AthenaPython.PyAthena as PyAthena


class ReferenceTruthObserver(PyAthena.Alg):
    def __init__(self, name, skip=0, rotated=True):
        super().__init__(name=name)
        self.skip = skip
        self.rotated = rotated

    def initialize(self):
        self.sg = PyAthena.py_svc('StoreGateSvc')
        self.stream = open('production_truth.jsonl', 'w', buffering=1024*1024)
        self.ordinal = 0
        return PyAthena.StatusCode.Success

    def execute(self):
        ei = self.sg.retrieve('xAOD::EventInfo', 'EventInfo')
        truth = self.sg.retrieve('xAOD::TruthParticleContainer', 'TruthParticles')
        hits = self.sg.retrieve('FaserSiHitCollection', 'SCT_Hits')
        tracks = self.sg.retrieve('TrackCollection', 'CKFTrackCollectionWithoutIFT')
        hit_codes = {int(hit.trackNumber()) for hit in hits}
        particles = {}
        for particle in truth:
            barcode = int(particle.barcode())
            particles[barcode] = dict(p_MeV=math.sqrt(particle.px()**2+particle.py()**2+particle.pz()**2),
                charge=float(particle.charge()), pdg=int(particle.pdgId()), status=int(particle.status()),
                has_sihit=barcode in hit_codes)
        final = []
        for track in tracks:
            params = list(track.trackParameters())
            quality = track.fitQuality()
            def radius(parameter):
                pos = parameter.position()
                return math.hypot(pos.y(), pos.z()) if self.rotated else math.hypot(pos.x(), pos.y())
            final.append(dict(chi2=float(quality.chiSquared()), ndof=float(quality.numberDoF()),
                measurements=int(track.measurementsOnTrack().size()),
                rmax_mm=max((radius(q) for q in params), default=None)))
        row = dict(source_entry=self.skip+self.ordinal, run=int(ei.runNumber()), event=int(ei.eventNumber()),
            ntracks=len(final), truth=particles, final=final)
        self.stream.write(json.dumps(row, allow_nan=False, separators=(',', ':'))+'\n')
        self.ordinal += 1
        return PyAthena.StatusCode.Success

    def finalize(self):
        self.stream.close()
        return PyAthena.StatusCode.Success
