class TreeConfig:
    def __init__(self, hbmin, mincomp, r, multiplier, stop, re_consensus):
        self.re_consensus = re_consensus
        self.stop = stop
        self.multiplier = multiplier
        self.cutoff_search_range = r
        self.mincomp = mincomp
        self.hbmin = hbmin
