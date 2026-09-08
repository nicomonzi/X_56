"""Read-only source checks; generated decks live in disposable temporary folders."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from model import ROOT, PID_NAMES, build, load_config
from analyse import release_check


class ModelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory(prefix='bff_dust_test_')
        cls.path = Path(cls.tmp.name)/'case'
        cls.config = build(cls.path, load_config(ROOT/'case.json'))
        cls.audit = json.loads((cls.path/'model_audit.json').read_text())
        cls.main = (cls.path/'main.mbd').read_text()

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_modes(self):
        self.assertEqual(self.config['mode_ids'], list(range(7,32)))
        self.assertAlmostEqual(self.config['mode_frequencies_hz'][-1],36.970058,places=5)
        self.assertGreaterEqual(self.audit['samples_per_highest_period'],12)

    def test_only_yaw_controllers_removed(self):
        self.assertEqual(set(self.audit['invariants']['controllers']),set(PID_NAMES))
        for name in ('YAW_PID','R_PID','YAW_DRIVE','R_DRIVE','NASTRAN_DLM_ROM'):
            self.assertNotIn(name,self.main)
        self.assertIn('force: SAS_MODAL_DAMPER_FORCE',self.main)
        self.assertTrue(self.audit['invariants']['yaw_locked'])
        self.assertTrue(self.audit['invariants']['no_mbdyn_aero'])
        self.assertEqual(self.audit['invariants']['hold_surfaces'],10)

    def test_filters_stable_unit_dc(self):
        for coefficients in self.audit['filters'].values():
            a,b = coefficients['a'],coefficients['b']
            self.assertTrue(np.all(np.abs(np.roots(a))<1))
            self.assertAlmostEqual(sum(b)/sum(a),1.,places=10)
        self.assertAlmostEqual(self.audit['filters']['ACT']['a'][1],-9/11)

    def test_hinge_probes_unique_and_small_change(self):
        points=np.loadtxt(self.path/'coupling_nodes.in')
        self.assertEqual(len(np.unique(points,axis=0)),59)
        offsets=np.asarray(self.audit['coupling_probe_adjustments_in'])
        np.testing.assert_array_equal(offsets[:39],0)
        self.assertLessEqual(np.linalg.norm(offsets,axis=1).max(),.010000001)

    def test_velocity_units_and_single_point(self):
        self.assertIn('u_inf = (/2611.811023622047,0.,0./)',(self.path/'dust.in').read_text())
        self.assertIn('max-time value="15.5"',(self.path/'precice.xml').read_text())

    def test_release_gate(self):
        time=np.arange(9.5,10.5,.002)
        state={'time':time,'v':np.zeros((len(time),3)),'omega':np.zeros((len(time),3))}
        self.assertTrue(release_check(state,self.config)['passed'])
        state['v'][:,2]=10
        self.assertFalse(release_check(state,self.config)['passed'])

    def test_bounded_wake(self):
        dust=(self.path/'dust.in').read_text()
        self.assertIn('n_wake_panels = 2\n',dust)
        self.assertIn('particles_box_max = (/581.0,1500.,1500./)',dust)
        self.assertIn('n_box = (/4,6,6/)',dust)
        self.assertIn('n_wake_particles = 50000\n',dust)
        self.assertEqual(self.config['wake_length_chords']*self.config['wake_reference_chord_in'],360)


if __name__=='__main__':
    unittest.main()
