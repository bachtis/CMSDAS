from CMSDAS.FourLeptons.Analyzer import *

class FourLeptonAnalyzer(Analyzer):
    def __init__(self):
        super(FourLeptonAnalyzer,self).__init__()






    #####CHANGE THIS METHOD TO CHANGE MUON ID######    
    def muonID(self,muon,vertex):
        if muon.pt()<5 or abs(muon.eta())>2.4:
            return False
        if muon.innerTrack().dxy( vertex.position() )>0.02:
            return False
        if muon.innerTrack().dz( vertex.position() )>0.2:
            return False
        if not (muon.isPFMuon() and \
               ( muon.isGlobalMuon() or muon.isTrackerMuon() )):
            return False
        # muon ISO variable
        if (muon.chargedHadronIso()+max(0.0,muon.photonIso()+muon.neutralHadronIso()-0.5*muon.puChargedHadronIso()))/muon.pt()>0.6:
            return False
        # muon SIP variable
        if (muon.dB(2)/muon.edB(2))>8:
            return False

        return True
    

    
    #####CHANGE THIS METHOD TO CHANGE ELECTRON ID######    
    def electronID(self,electron,vertex):
        if electron.pt()<7 or abs(electron.eta())>2.5:
            return False


        mvaRegions = [{'ptMin':0,'ptMax':10, 'etaMin':0.0, 'etaMax':0.8,'mva':0.47},\
                      {'ptMin':0,'ptMax':10, 'etaMin':0.8 ,'etaMax':1.479,'mva':0.004},\
                      {'ptMin':0,'ptMax':10, 'etaMin':1.479, 'etaMax':3.0,'mva':0.295},\
                      {'ptMin':10,'ptMax':99999999, 'etaMin':0.0, 'etaMax':0.8,'mva':-0.34},\
                      {'ptMin':10,'ptMax':99999999, 'etaMin':0.8, 'etaMax':1.479,'mva':-0.65},\
                      {'ptMin':10,'ptMax':99999999, 'etaMin':1.479, 'etaMax':3.0,'mva':0.6}]
        ID=False 
        for element in mvaRegions:
            if electron.pt()>= element['ptMin'] and \
               electron.pt()< element['ptMax'] and \
               abs(electron.superCluster().eta())>=element['etaMin'] and \
               abs(electron.superCluster().eta())<element['etaMax'] and \
               electron.electronID("mvaNonTrigV0")> element['mva']: 
                ID=True
        if not ID:
            return False

        if electron.gsfTrack().trackerExpectedHitsInner().numberOfHits()>1:
            return False

        # electron ISO variable
        if (electron.chargedHadronIso()+max(0.0,electron.photonIso()+electron.neutralHadronIso()-0.5*electron.puChargedHadronIso()))/electron.pt()>0.6:
            return False
        # electron SIP variable
        if (electron.dB(2)/electron.edB(2))>8:
            return False


        return True


    def analyze(self,box):

        #####START FROM A bOX CONTAINING SELECTED MUONS AND ELECTRONS and MAKE
        #FOUR LEPTON CANDIDATES
        
        
        #Now check if there are at least four leptons:
        box.leptons=set(box.selectedMuons+box.selectedElectrons)
        if len(box.leptons)<4:
            return False

        #Now create Z candidates and apply cuts:
        box.zcandidates=[]

        for l1,l2 in itertools.combinations(box.leptons,2):
            #ask that lepton ID passed for the Z1 ONLY
            #for data driven estimation
            #they need to have same flavour and OS
            if abs(l1.pdgId()) != abs(l2.pdgId()):
                continue
            if l1.charge() +l2.charge() !=0:
                continue
            #now create a di lepton object and check mass
            z=DiObject(l1,l2)
            if not (z.mass()>12 and z.mass()<120):
                continue
            box.zcandidates.append(DiObject(l1,l2))
        # OK if there are more than one Z candidates
        #pick the one with the best mass.
        if len(box.zcandidates)==0:
            return False



        sortedZs=sorted(box.zcandidates,key=lambda x: abs(x.mass()-91.118))
        box.Z1 = sortedZs[0]
        
        #now remove the used leptons from the list and make Z2 pairs
        box.leptons.remove(box.Z1.l1)        
        box.leptons.remove(box.Z1.l2)


        #now the same thing with the second
        box.zcandidates2=[]

        for l1,l2 in itertools.combinations(box.leptons,2):
            #they need to have same flavour and OS
            if abs(l1.pdgId()) != abs(l2.pdgId()):
                continue
            if l1.charge() +l2.charge() !=0:
                continue
            #now create a di lepton object and check mass
            z=DiObject(l1,l2)
            if not (z.mass()>4 and z.mass()<120):
                continue
            box.zcandidates2.append(DiObject(l1,l2))
        # OK if there are more than one Z candidates
       #pick the one with the highest lepton pt sum
        if len(box.zcandidates2)==0:
            return False


        
        sortedZ2s=sorted(box.zcandidates2,key=lambda x: x.l1.pt()+x.l2.pt(),reverse=True)
        box.Z2 = sortedZ2s[0]

        #kill the candidate if a OS pair has mll<4 GeV
        for l1,l2 in itertools.combinations([box.Z1.l1,box.Z1.l2,box.Z2.l1,box.Z2.l2],2):
            ll =DiObject(l1,l2)
            if (l1.charge()+l2.charge()) ==0:
                if ll.mass()<4 :
                    return False


        
        #create the ZZ
        box.ZZ = DiObject(box.Z1,box.Z2)
        
        return True



    def declareHistos(self):
        super(FourLeptonAnalyzer,self).declareHistos()
        ###ADD YOUR HISTOGRAMS AFTER THIS LINE AS AbOVE#####
        self.declareHisto('mass',30,70,150,"m_{4l} [GeV]")
        self.declareHisto('massFull',100,70,570,"m_{4l} [GeV]")
        self.declareHisto('massZ1',20,12,120,"m_{Z1} [GeV]")
        self.declareHisto('massZ2',20,4,74,"m_{Z2} [GeV]")


    def fillHistos(self,box,sample,weight = 1):
        super(FourLeptonAnalyzer,self).fillHistos(box,sample,weight)

        self.fillHisto('mass',sample,box.ZZ.mass(),weight)        
        self.fillHisto('massFull',sample,box.ZZ.mass(),weight)        
        self.fillHisto('massZ1',sample,box.ZZ.l1.mass(),weight)        
        self.fillHisto('massZ2',sample,box.ZZ.l2.mass(),weight)        



        
        
