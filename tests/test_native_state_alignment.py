import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'analysis'))
from slope_native_stages import align_transport


class NativeStateAlignmentTest(unittest.TestCase):
    def test_transport_event_order_preserves_stage_payload(self):
        state=dict(entry=np.array([0,1,0,1]),track_index=np.array([0,0,0,0]),
                   stage=np.array([3,3,11,11]),barcode=np.array([9,8,9,8]),
                   layer=np.array([0,0,-1,-1]))
        order=np.array([0,2,1,3])
        transport={b:state[a][order] for a,b in [('entry','source_entry'),
            ('track_index','seed'),('stage','stage'),('barcode','barcode'),('layer','layer')]}
        transport['payload']=np.array([30,110,31,111])
        actual=align_transport(state,transport)
        np.testing.assert_array_equal(actual['payload'],[30,31,110,111])

    def test_same_row_count_does_not_hide_missing_identity(self):
        state=dict(entry=np.array([0,1]),track_index=np.array([0,0]),
                   stage=np.array([3,3]),barcode=np.array([9,8]),layer=np.array([0,0]))
        transport=dict(source_entry=np.array([0,2]),seed=np.array([0,0]),
                       stage=state['stage'],barcode=state['barcode'],layer=state['layer'])
        with self.assertRaisesRegex(AssertionError,'coverage differs'):
            align_transport(state,transport)


if __name__=='__main__':unittest.main()
