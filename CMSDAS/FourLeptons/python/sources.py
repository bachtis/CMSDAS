from CMSDAS.FourLeptons.Sample import Sample
lumi=19.8



ZZTo2e2mu = Sample('ZZ',True,['root://eoscms//eos/cms//store/cmst3/group/das2014/FourLeptons/ZZTo2e2mu.root'],lumi*1000*0.1767)

ZZTo4e = Sample('ZZ',True,['root://eoscms//eos/cms//store/cmst3/group/das2014/FourLeptons/ZZTo4e.root'],lumi*1000*0.07691)
ZZTo4mu = Sample('ZZ',True,['root://eoscms//eos/cms//store/cmst3/group/das2014/FourLeptons/ZZTo4mu.root'],lumi*1000*0.07691)
GGZZTo2L2L = Sample('ZZ',True,['root://eoscms//eos/cms//store/cmst3/group/das2014/FourLeptons/GGZZ2L2L.root'],lumi*1000*0.01203)
GGZZTo4L = Sample('ZZ',True,['root://eoscms//eos/cms//store/cmst3/group/das2014/FourLeptons/GGZZ4L.root'],lumi*1000*0.0048)
higgs = Sample('higgs',True,['root://eoscms//eos/cms//store/cmst3/group/das2014/FourLeptons/higgs.root'],lumi*1000*21.749*3.02e-4)



fakeRate = Sample('fakeRate',False,['root://eoscms//eos/cms//store/cmst3/group/das2014/FourLeptons/fakeRate.root'],1)



data_files=[
'root://eoscms//eos/cms//store/cmst3/group/das2014/FourLeptons/DoubleElectron.root',
'root://eoscms//eos/cms//store/cmst3/group/das2014/FourLeptons/DoubleMu.root',
'root://eoscms//eos/cms//store/cmst3/group/das2014/FourLeptons/MuEG.root'
]

data = Sample('data',False,data_files,1)


sources = [ZZTo2e2mu,ZZTo4e,ZZTo4mu,GGZZTo2L2L,GGZZTo4L,higgs,data]







