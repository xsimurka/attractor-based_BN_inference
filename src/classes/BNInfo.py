from typing import Tuple, Dict, List


State = Tuple[bool]


class BNInfo:

    def __init__(self):
        self.wt_sinks: List[State] = []
        self.ko_sinks: Dict[int, List[State]] = {}
        self.oe_sinks: Dict[int, List[State]] = {}