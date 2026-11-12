"""Supply empty waveform products for a tracking-only Ntuple branch check."""
import AthenaPython.PyAthena as PyAthena
import ROOT


class EmptyWaveforms(PyAthena.Alg):
    def initialize(self):
        self.sg = PyAthena.py_svc('StoreGateSvc')
        return PyAthena.StatusCode.Success

    def execute(self):
        for detector in ['VetoNu', 'Veto', 'Trigger', 'Preshower', 'Muon',
                         'VetoGm', 'Calo', 'Calo2', 'CaloNu']:
            name = detector + 'WaveformHits'
            hits = ROOT.xAOD.WaveformHitContainer()
            aux = ROOT.xAOD.WaveformHitAuxContainer()
            hits.setStore(aux)
            if not self.sg.record(aux, name + 'Aux.') or not self.sg.record(hits, name):
                return PyAthena.StatusCode.Failure
            ROOT.SetOwnership(aux, False)
            ROOT.SetOwnership(hits, False)
        return PyAthena.StatusCode.Success
