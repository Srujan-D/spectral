class LongitudinalState:
    def __init__(self) -> None:
        self.s = 0.0
        self.ds = 0.0
        self.dds = 0.0
        self.t = 0.0
    
    def __init__(self, s, ds, dds, t) -> None:
        self.s = s
        self.ds = ds
        self.dds = dds
        self.t = t