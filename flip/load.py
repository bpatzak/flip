class Load:
    def __init__(self):
        self.nodalLabels = []
        self.elementLabels=[]
    def actsOnElement (self, elem):
        return (elem in self.elementLabels)
    # returns list of location arrays
    def getLocationArrays(self, domain):
        pass
    # returns list of value arrays
    def getValues(self, domain):
        pass
