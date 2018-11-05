class TreeConfig:
    def __init__(self, hbmin, r, multiplier, stop):
        self.stop = stop
        self.multiplier = multiplier
        self.cutoff_search_range = r
        self.hbmin = hbmin
