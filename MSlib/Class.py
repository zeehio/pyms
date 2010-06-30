import sys
sys.path.append("/x/PyMS/")

class msrecord(object):

    def __init__(self):
        self.name = ''
        self.regno = ''
        self.mi = {}

    def set(self, name, regno, mi):
        self.name = name
        self.regno = regno
        self.mi = mi

class mslib(object):

    def __init__(self):
        self.msrecord_list = []

    def addrecord(self, msrecord):
        self.msrecord_list.append(msrecord)

