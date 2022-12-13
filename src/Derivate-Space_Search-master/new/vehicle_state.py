class VehicleState:
    def __init__(self) -> None:
        self.s = 0.0
        self.ds = 0.0
        self.dds = 0.0

        self.l = 0.0
        self.dl = 0.0
        self.ddl = 0.0

        self.t = 0.0
    
    def __init__(self, s, ds, dds, l, dl, ddl, t) -> None:
        self.s = s
        self.ds = ds
        self.dds = dds

        self.l = l
        self.dl = dl
        self.ddl = ddl

        self.t = t

