from DataFormats.FWLite import Events, Handle,Lumis
import itertools
import ROOT
import os
import json
import math

def getFullPath(path):
    return os.environ['CMSSW_BASE']+"/src/CMSDAS/FourLeptons/"+path



class EventBox(object):
    def __init__(self):

        self.Z1_l1=None
        self.Z1_l2=None
        self.Z2_l1=None
        self.Z2_l2=None



class DiObject(object):
    def __init__(self,l1,l2):
        self.l1=l1
        self.l2=l2
        self.p4 = ROOT.TLorentzVector(self.l1.px()+self.l2.px(),\
                                           self.l1.py()+self.l2.py(),\
                                           self.l1.pz()+self.l2.pz(),\
                                           self.l1.energy()+self.l2.energy())
    def mass(self):
        return self.p4.M()
    
    def pdgId(self):
        return 23
            
    def p4(self):
        return self.p4
        
    def pt(self):
        return self.p4.Pt() 

    def px(self):
        return self.p4.Px() 

    def py(self):
        return self.p4.Py() 

    def pz(self):
        return self.p4.Pz() 

    def energy(self):
        return self.p4.Energy() 
    
        
class Analyzer (object):
    def __init__(self):
        self.vertexHandle  = Handle ('std::vector<reco::Vertex>')
        self.muonHandle    = Handle('std::vector<pat::Muon>')
        self.electronHandle    = Handle('std::vector<pat::Electron>')
        self.histograms={}
        self.samples = ['fakes','ZZ','higgs','data','ZZ3P1F','fakes3P1F']
        for sample in self.samples:
            self.histograms[sample]={}

        puFileMC = getFullPath('data/puProfile_Summer12_53X.root')
        puFileData = getFullPath('data/puProfile_Data12.root')
        self.fPUMC = ROOT.TFile(puFileMC)
        self.puHistoMC = self.fPUMC.Get('pileup')
        self.fPUData = ROOT.TFile(puFileData)
        self.puHistoData = self.fPUData.Get('pileup')
        self.puHistoMC.Scale(1./self.puHistoMC.Integral())
        self.puHistoData.Scale(1./self.puHistoData.Integral())
        self.data=[]
        
        self.muonFRD = ROOT.TH2F("muonFRD","muonFR",3,5,50,2,0,2.4) 
        self.electronFRD = ROOT.TH2F("electronFRD","electronFR",3,7,50,2,0,2.5) 
        self.muonFR = ROOT.TH2F("muonFR","muonFR",3,5,50,2,0,2.4) 
        self.electronFR = ROOT.TH2F("electronFR","electronFR",3,7,50,2,0,2.5) 
        self.optimize=True
        

    def muonID(self,muon,vertex):
        return True

    def electronID(self,electron,vertex):
        return True

    def muonPreID(self,muon,vertex):
        if muon.pt()<5 or abs(muon.eta())>2.4:
            return False
        if not ( muon.isGlobalMuon() or muon.isTrackerMuon() ):
            return False
        if (muon.chargedHadronIso()+max(0.0,muon.photonIso()+muon.neutralHadronIso()-0.5*muon.puChargedHadronIso()))/muon.pt()>1:
            return False
        return True

    def electronPreID(self,electron,vertex):
        if electron.pt()<7 or abs(electron.eta())>2.5:
            return False
        if electron.electronID("mvaNonTrigV0")<-0.65:
            return False
        if (electron.chargedHadronIso()+max(0.0,electron.photonIso()+electron.neutralHadronIso()-0.5*electron.puChargedHadronIso()))/electron.pt()>1.0:
            return False

        return True



    def leptonID(self,lepton,vertex):
        if abs(lepton.pdgId())==11:
            return self.electronID(lepton,vertex)
        elif abs(lepton.pdgId())==13:
            return self.muonID(lepton,vertex)

    def readCollections(self,event,box,isFake = False):
        event.getByLabel('offlinePrimaryVertices',self.vertexHandle)
        event.getByLabel('muons',self.muonHandle)
        event.getByLabel('electrons',self.electronHandle)

        box.muons = self.muonHandle.product()
        box.electrons = self.electronHandle.product()
        box.vertex = self.vertexHandle.product()[0]  

        #first select muons and electrons
        box.selectedMuons = []
        for mu in box.muons:
            if mu.pt()<5 or abs(mu.eta())>2.4:
                continue
            if (self.muonID(mu,box.vertex) or (isFake and self.muonPreID(mu,box.vertex))):
                box.selectedMuons.append(mu)
                
        box.selectedElectrons = []
        for ele in box.electrons:
            if ele.pt()<5 or abs(ele.eta())>2.5:
                continue
            if self.electronID(ele,box.vertex) or (isFake and self.electronPreID(ele,box.vertex)):
                box.selectedElectrons.append(ele)


    def analyze(self,box):
        return True



    def declareHistos(self):
        self.declareHisto('signalR110150',30,110,150)
        self.declareHisto('signalR120130',10,120,130)

    def fillHistos(self,box,sample,weight = 1):
        self.fillHisto('signalR110150',sample,box.ZZ.mass(),weight)        
        self.fillHisto('signalR120130',sample,box.ZZ.mass(),weight)        


    def fakeRate(self,lepton):
        histo=None
        if abs(lepton.pdgId())==11:
            histo=self.electronFR
        else:    
            histo=self.muonFR

        binx = histo.GetXaxis().FindBin(lepton.pt())
        if binx>3:
            binx=3
        biny = histo.GetYaxis().FindBin(abs(lepton.eta()))
        if biny>3:
            biny=3
        return histo.GetBinContent(histo.GetBin(binx,biny))

    
    def calculateFakeRates(self,sample,maxEv = -1):
        events=Events(sample.files)
        for event in events:
            box = EventBox()
            self.readCollections(event,box,True)
            box.leptons=set(box.selectedMuons+box.selectedElectrons)
            if len(box.leptons)!=3:
                continue
            
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
                if not (z.mass()>70 and z.mass()<120):
                    continue
                box.zcandidates.append(DiObject(l1,l2))
            # OK if there are more than one Z candidates
            #pick the one with the best mass.
            if len(box.zcandidates)==0:
                continue


            sortedZs=sorted(box.zcandidates,key=lambda x: abs(x.mass()-91.118))
            box.Z1 = sortedZs[0]
            
            box.leptons.remove(box.Z1.l1)        
            box.leptons.remove(box.Z1.l2)

            if len(box.leptons)!=1:
                continue

            lepton=list(box.leptons)[0]
            
            if abs(lepton.pdgId())==11:
                self.electronFRD.Fill(lepton.pt(),abs(lepton.eta()))
                if self.electronID(lepton,box.vertex):
                    self.electronFR.Fill(lepton.pt(),abs(lepton.eta())) 
            if abs(lepton.pdgId())==13:
                self.muonFRD.Fill(lepton.pt(),abs(lepton.eta()))
                if self.muonID(lepton,box.vertex):
                    self.muonFR.Fill(lepton.pt(),abs(lepton.eta())) 
                   
        self.muonFR.Divide(self.muonFRD)
        self.electronFR.Divide(self.electronFRD)
            

    def processSample(self,sample,maxEv = -1):
        
        print 'Processing Files'
        print sample.files

        events=Events(sample.files)
        counterHandle = Handle("edm::MergeableCounter")

        #loop and find the #events before skim
        if sample.isMC:
            lumis = Lumis(sample.files)
            print 'Counting the events for normalization' 
            totalEv=0.0

            for lumi in lumis:
                lumi.getByLabel('prePathCounter',counterHandle)
                totalEv+=counterHandle.product().value
            print 'Continuing on the event Loop' 


        #Now setup PU reweighting
        puHandle = Handle("std::vector<PileupSummaryInfo>")

        #Allow to process less events to run faster
        weightFactor=1
        if sample.isMC:
            if maxEv>0 and maxEv<events.size():
                totalEv=totalEv*(float(maxEv))/(float(events.size()))
        N=0        
        for event in events:
            if sample.isMC and maxEv>0 and N>maxEv:
                break
            N=N+1

            weight=1
            if sample.isMC:
                weight=sample.crossSection/totalEv
                #now read the PU info
                inter=0
                event.getByLabel('addPileupInfo',puHandle)
                for puInfo in puHandle.product():
                    if puInfo.getBunchCrossing()==0:
                        inter = puInfo.getTrueNumInteractions()
                data=self.puHistoData.GetBinContent(self.puHistoData.FindBin(inter))        
                mc=self.puHistoMC.GetBinContent(self.puHistoMC.FindBin(inter))
                weight=weight*data/mc
                
            box = EventBox()
            self.readCollections(event,box)
            if (not sample.isMC): 
                boxNoID = EventBox()
                self.readCollections(event,boxNoID,True)
                #Now select the candidate on that
                if self.analyze(boxNoID):
                    if self.leptonID(boxNoID.ZZ.l1.l1,boxNoID.vertex) and self.leptonID(boxNoID.ZZ.l1.l2,boxNoID.vertex):
                        if (not self.leptonID(boxNoID.ZZ.l2.l1,boxNoID.vertex)) and (not self.leptonID(boxNoID.ZZ.l2.l2,boxNoID.vertex)):
                            fr1=self.fakeRate(boxNoID.ZZ.l2.l1)
                            fr2=self.fakeRate(boxNoID.ZZ.l2.l2)
                            self.fillHistos(boxNoID,'fakes',fr1*fr2/((1-fr1)*(1-fr2)))
            if not self.analyze(box):
                continue

            if not sample.isMC and box.ZZ.mass()>110 and box.ZZ.mass()<150 and self.optimize:
                continue
            self.fillHistos(box,sample.type,weight)
            if not sample.isMC:
                if box.ZZ.mass()>110 and box.ZZ.mass()<150 :
                    self.data.append(box.ZZ.mass())
#                    print 'EVENT SELECTED:', event.eventAuxiliary().id().run(),event.eventAuxiliary().id().luminosityBlock(),event.eventAuxiliary().id().event()

        #Now apply the weights to the histopgrams for the fakes

    def renormalizeFakes(self):
            for name,histo in self.histograms['fakes'].iteritems():
                histo.Smooth()
                histo.Smooth()


    def exportData(self):
        f = open('data.json','w')
        info={'data':self.data , 'higgs':self.histograms['higgs']['signalR110150'].Integral(), 'ZZ':self.histograms['ZZ']['signalR110150'].Integral(), 'ZX':self.histograms['fakes']['signalR110150'].Integral()}
        json.dump(info,f)
        f.close()
    
    def declareHisto(self,name,bins,min,max,xlabel=''):
        for sample in self.samples:
            self.histograms[sample][name]=ROOT.TH1F(sample+'_'+name,name,bins,min,max)
            self.histograms[sample][name].GetXaxis().SetTitle(xlabel)
            self.histograms[sample][name].GetYaxis().SetTitle('events')
            


    def fillHisto(self,name,sample,value,weight = 1):
        self.histograms[sample][name].Fill(value,weight)




    def convertToPoisson(self,h):
        graph = ROOT.TGraphAsymmErrors()
        q = (1-0.6827)/2.

        for i in range(1,h.GetNbinsX()+1):
            x=h.GetXaxis().GetBinCenter(i)
            xLow =h.GetXaxis().GetBinLowEdge(i) 
            xHigh =h.GetXaxis().GetBinUpEdge(i) 
            y=h.GetBinContent(i)
            yLow=0
            yHigh=0
            if y !=0.0:
                yLow = y-ROOT.Math.chisquared_quantile_c(1-q,2*y)/2.
                yHigh = ROOT.Math.chisquared_quantile_c(q,2*(y+1))/2.-y
                graph.SetPoint(i-1,x,y)
                graph.SetPointEYlow(i-1,yLow)
                graph.SetPointEYhigh(i-1,yHigh)
                #            graph.SetPointEXlow(i-1,x-xLow)
                #            graph.SetPointEXhigh(i-1,xHigh-x)
                graph.SetPointEXlow(i-1,0.0)
                graph.SetPointEXhigh(i-1,0.0)


        graph.SetMarkerStyle(20)
        graph.SetLineWidth(2)
        graph.SetMarkerSize(1.)
        graph.SetMarkerColor(ROOT.kBlack)
    

        return graph    


    def makeStack(self,histogram):
        sandbox = []
        canvas = ROOT.TCanvas("stack"+histogram)
        canvas.cd()
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptTitle(0)
        canvas.Range(-68.75,-7.5,856.25,42.5)
        canvas.SetFillColor(0)
        canvas.SetBorderMode(0)
        canvas.SetBorderSize(2)
        canvas.SetTickx(1)
        canvas.SetTicky(1)
        canvas.SetLeftMargin(0.15)
        canvas.SetRightMargin(0.05)
        canvas.SetTopMargin(0.05)
        canvas.SetBottomMargin(0.15)
        canvas.SetFrameFillStyle(0)
        canvas.SetFrameBorderMode(0)
        canvas.SetFrameFillStyle(0)
        canvas.SetFrameBorderMode(0)


        canvas.cd()

        stack = ROOT.THStack("stackH"+histogram,"")
        self.histograms['fakes'][histogram].SetFillStyle(1001)
        self.histograms['fakes'][histogram].SetFillColor(ROOT.kGreen-5)
        stack.Add(self.histograms['fakes'][histogram])
        self.histograms['ZZ'][histogram].SetFillStyle(1001)
        self.histograms['ZZ'][histogram].SetFillColor(ROOT.kAzure-9)
        stack.Add(self.histograms['ZZ'][histogram])
        self.histograms['higgs'][histogram].SetFillStyle(0)
        self.histograms['higgs'][histogram].SetFillColor(ROOT.kWhite)
        self.histograms['higgs'][histogram].SetLineStyle(1)
        self.histograms['higgs'][histogram].SetLineWidth(3)
        self.histograms['higgs'][histogram].SetLineColor(ROOT.kOrange+10)
        stack.Add(self.histograms['higgs'][histogram])

        self.histograms['data'][histogram].SetMarkerStyle(20)

        
        datamax = max(ROOT.Math.chisquared_quantile_c((1-0.6827)/2.,2*(self.histograms['data'][histogram].GetMaximum()+1))/2.,stack.GetMaximum())*1.2

        frame = canvas.DrawFrame(self.histograms['higgs'][histogram].GetXaxis().GetXmin(),0.0,self.histograms['higgs'][histogram].GetXaxis().GetXmax(),datamax)
        frame.GetXaxis().SetLabelFont(42)
        frame.GetXaxis().SetLabelOffset(0.007)
        frame.GetXaxis().SetLabelSize(0.045)
        frame.GetXaxis().SetTitleSize(0.05)
        frame.GetXaxis().SetTitleOffset(1.15)
        frame.GetXaxis().SetTitleFont(42)
        frame.GetYaxis().SetLabelFont(42)
        frame.GetYaxis().SetLabelOffset(0.007)
        frame.GetYaxis().SetLabelSize(0.045)
        frame.GetYaxis().SetTitleSize(0.05)
        frame.GetYaxis().SetTitleOffset(1.4)
        frame.GetYaxis().SetTitleFont(42)
        frame.GetZaxis().SetLabelFont(42)
        frame.GetZaxis().SetLabelOffset(0.007)
        frame.GetZaxis().SetLabelSize(0.045)
        frame.GetZaxis().SetTitleSize(0.05)
        frame.GetZaxis().SetTitleFont(42)


##         if len(units)>0:
##             frame.GetXaxis().SetTitle(titlex + " (" +units+")")
##             frame.GetYaxis().SetTitle("Events / "+str((maxi-mini)/bins)+ " "+units)
##         else:    
##             frame.GetXaxis().SetTitle(titlex)
##             frame.GetYaxis().SetTitle("Events")

        frame.Draw()
        stack.Draw("A,HIST,SAME")
        self.histograms['data'][histogram].Sumw2()
        dataG = None
        if self.histograms['data'][histogram].Integral()>0:
            dataG = self.convertToPoisson(self.histograms['data'][histogram])
            dataG.Draw("Psame")
            

        
        legend = ROOT.TLegend(0.62,0.6,0.92,0.90,"","brNDC")
	legend.SetBorderSize(0)
	legend.SetLineColor(1)
	legend.SetLineStyle(1)
	legend.SetLineWidth(1)
	legend.SetFillColor(0)
	legend.SetFillStyle(0)
	legend.SetTextFont(42)
        legend.AddEntry(self.histograms['fakes'][histogram],"Z+X","f")
        legend.AddEntry(self.histograms['ZZ'][histogram],"Z#gamma^{*}","f")
        legend.AddEntry(self.histograms['higgs'][histogram],"m_{H} = 126 GeV","f")
        legend.AddEntry(self.histograms['data'][histogram],"Data","p")
        legend.Draw()

        pt =ROOT.TPaveText(0.1577181,0.9562937,0.9580537,0.9947552,"brNDC")
	pt.SetBorderSize(0)
	pt.SetTextAlign(12)
	pt.SetFillStyle(0)
	pt.SetTextFont(42)
	pt.SetTextSize(0.03)
	text = pt.AddText(0.01,0.3,"CMS")
	text = pt.AddText(0.25,0.3,"                             #sqrt{s} = 8 TeV, L = 19.7 fb^{-1}")
	pt.Draw()   

        plot={'canvas':canvas,'stack':stack,'legend':legend,'dataG':dataG,'latex1':pt}



        print 'STATISTICS at the region 110-150'
        print '--------------------------------'
        S = self.histograms['higgs']['signalR110150'].Integral()
        ZZ = self.histograms['ZZ']['signalR110150'].Integral()
        ZX = self.histograms['fakes']['signalR110150'].Integral()
        DATA= self.histograms['data']['signalR110150'].Integral()

        print 'Expected Signal = ',S
        print 'Expected ZZ = ',ZZ
        print 'Expected reducible = ',ZX
        print 'Total background = ',ZX+ZZ
        print 'Total expected = ',ZX+ZZ+S
        print 'Observed Data = ',DATA

        print 'expected s/sqrt(b)',S/math.sqrt(ZX+ZZ)
        print '                                ' 
        print 'STATISTICS at the region 120-130'
        print '--------------------------------'
        S = self.histograms['higgs']['signalR120130'].Integral()
        ZZ = self.histograms['ZZ']['signalR120130'].Integral()
        ZX = self.histograms['fakes']['signalR120130'].Integral()
        DATA= self.histograms['data']['signalR120130'].Integral()

        print 'Expected Signal = ',S
        print 'Expected ZZ = ',ZZ
        print 'Expected reducible = ',ZX
        print 'Total background = ',ZX+ZZ
        print 'Total expected = ',ZX+ZZ+S
        print 'Observed Data = ',DATA
        print 'expected s/sqrt(b)',S/math.sqrt(ZX+ZZ)
        
        
        canvas.RedrawAxis()
        canvas.Update()
        


        return plot




    def makeComparison(self,histogram):
        sandbox = []
        canvas = ROOT.TCanvas("stack"+histogram)
        canvas.cd()
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptTitle(0)
        canvas.Range(-68.75,-7.5,856.25,42.5)
        canvas.SetFillColor(0)
        canvas.SetBorderMode(0)
        canvas.SetBorderSize(2)
        canvas.SetTickx(1)
        canvas.SetTicky(1)
        canvas.SetLeftMargin(0.15)
        canvas.SetRightMargin(0.05)
        canvas.SetTopMargin(0.05)
        canvas.SetBottomMargin(0.15)
        canvas.SetFrameFillStyle(0)
        canvas.SetFrameBorderMode(0)
        canvas.SetFrameFillStyle(0)
        canvas.SetFrameBorderMode(0)


        canvas.cd()

        stack = ROOT.THStack("stackH"+histogram,"")


        self.histograms['fakes'][histogram].SetFillStyle(3017)
        self.histograms['fakes'][histogram].SetFillColor(ROOT.kGreen-5)
        fakesClone = self.histograms['fakes'][histogram].Clone()
        norm = fakesClone.Integral()
        if norm>0:
            fakesClone.Scale(1./norm)
        stack.Add(fakesClone)
        


        self.histograms['ZZ'][histogram].SetFillStyle(3018)
        self.histograms['ZZ'][histogram].SetFillColor(ROOT.kBlue)
        zzClone = self.histograms['ZZ'][histogram].Clone()
        norm = zzClone.Integral()
        if norm>0:
            
            zzClone.Scale(1./norm)
        stack.Add(zzClone)

        self.histograms['higgs'][histogram].SetFillStyle(0)
        self.histograms['higgs'][histogram].SetFillColor(ROOT.kWhite)
        self.histograms['higgs'][histogram].SetLineStyle(1)
        self.histograms['higgs'][histogram].SetLineWidth(3)
        self.histograms['higgs'][histogram].SetLineColor(ROOT.kOrange+10)

        higgsClone = self.histograms['higgs'][histogram].Clone()
        norm=higgsClone.Integral()
        if norm>0:
            
            higgsClone.Scale(1./norm)
        stack.Add(higgsClone)

        stack.Draw("HIST,NOSTACK")
        legend = ROOT.TLegend(0.62,0.6,0.92,0.90,"","brNDC")
	legend.SetBorderSize(0)
	legend.SetLineColor(1)
	legend.SetLineStyle(1)
	legend.SetLineWidth(1)
	legend.SetFillColor(0)
	legend.SetFillStyle(0)
	legend.SetTextFont(42)
        legend.AddEntry(self.histograms['fakes'][histogram],"Z+X","f")
        legend.AddEntry(self.histograms['ZZ'][histogram],"Z#gamma^{*}","f")
        legend.AddEntry(self.histograms['higgs'][histogram],"m_{H} = 126 GeV","f")
        legend.Draw()

        pt =ROOT.TPaveText(0.1577181,0.9562937,0.9580537,0.9947552,"brNDC")
	pt.SetBorderSize(0)
	pt.SetTextAlign(12)
	pt.SetFillStyle(0)
	pt.SetTextFont(42)
	pt.SetTextSize(0.03)
	text = pt.AddText(0.01,0.3,"CMS")
	text = pt.AddText(0.25,0.3,"                             #sqrt{s} = 8 TeV, L = 19.7 fb^{-1}")
	pt.Draw()   

        plot={'canvas':canvas,'stack':stack,'legend':legend,'latex1':pt,'clones':[zzClone,fakesClone,higgsClone]}


        canvas.RedrawAxis()
        canvas.Update()
        
        return plot




        
        
