class Node:
    def __init__(self, label, coords, bcs):
        self.label=label
        self.coordinates = coords
        self.bcs = set(bcs)
