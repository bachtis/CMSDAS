import ROOT
ROOT.gROOT.ProcessLine(".x tdrstyle.C")
from CMSDAS.FourLeptons.sources import *
from CMSDAS.FourLeptons.FourLeptonAnalyzer  import *


analyzer = FourLeptonAnalyzer()

analyzer.declareHistos()
analyzer.optimize=False
analyzer.calculateFakeRates(fakeRate)

for sample in sources:
    analyzer.processSample(sample,5000)
analyzer.renormalizeFakes()
analyzer.exportData()
