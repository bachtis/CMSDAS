import ROOT
import json

ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")


f=ROOT.TFile("workspace.root","RECREATE")

w = ROOT.RooWorkspace('w')

w.factory('MH[126,110,150]')
w.factory('m4l[126,110,150]')
w.factory('scale[0.003,0,1]')
w.factory('resolution[0.2,0,1]')
w.factory('ggH_norm[1]')



w.factory("expr::m0('MH-675.028021767+(28.4058912104*MH)+(-0.471030234007*MH*MH)+(0.00385523878002*MH*MH*MH)+(-1.56033149e-05*MH*MH*MH*MH)+(2.50152559617e-08*MH*MH*MH*MH*MH)*(1+scale)',MH,scale)")

w.factory("expr::sigma('(5893.99328641+(-223.357259613*MH)++(3.3783014563*MH*MH)+(-0.0254887585829*MH*MH*MH)+(9.59351485306e-05*MH*MH*MH*MH)+(-1.44104967891e-07*MH*MH*MH*MH*MH))*(1+resolution)',MH,resolution)")

w.factory("expr::alphaR('(2405.13381227+(-89.695355535*MH)++(1.334717089*MH*MH)+(-0.00990329417255*MH*MH*MH)+(3.6642501963e-05*MH*MH*MH*MH)+(-5.40906340625e-08*MH*MH*MH*MH*MH))',MH)")

w.factory("expr::nL('(84.2437375249+(-0.0619644765407*MH)++(-0.0387295344451*MH*MH)+(0.000546309838433*MH*MH*MH)+(-2.84401727179e-06*MH*MH*MH*MH)+(5.20107205018e-09*MH*MH*MH*MH*MH))',MH)")

w.factory("expr::nR('(20.0042185053+(-0.000158974754724*MH)++(2.39215331567e-06*MH*MH)+(-1.79657858317e-08*MH*MH*MH)+(6.73435367775e-11*MH*MH*MH*MH)+(-1.00791862259e-13*MH*MH*MH*MH*MH))',MH)")

w.factory("expr::alphaL('(3521.36524827+(-134.348137936*MH)+(2.04466932799*MH*MH)+(-0.0155122272155*MH*MH*MH)+(5.86665610984e-05*MH*MH*MH*MH)+(-8.84845540231e-08*MH*MH*MH*MH*MH))',MH)")



pdf    = ROOT.RooDoubleCB('ggH','signal',w.var('m4l'),w.function('m0'),w.function('sigma'),w.function('alphaL'),w.function('nL'),w.function('alphaR'),w.function('nR'))
getattr(w,'importClassCode')(ROOT.RooDoubleCB.Class(),1)
getattr(w,'import')(pdf)


w.factory('RooUniform::ZZ(m4l)')
w.factory('RooUniform::ZX(m4l)')

f.cd()

#now import data
ff=open('data.json')
info=json.load(ff)
events = info['data']

data = ROOT.RooDataSet("data_obs","data",ROOT.RooArgSet(w.var('m4l')))
for event in events:
    w.var('m4l').setVal(event)
    data.add(ROOT.RooArgSet(w.var('m4l')))


getattr(w,'import')(data)
w.Write()
ff.close()
f.Close()


d = open('datacard.txt','w')

d.write("""
imax 1 number of bins
jmax 2 number of processes minus 1
kmax * number of nuisance parameters
--------------------------------------------------------- 
shapes *                hzz4l  workspace.root w:$PROCESS
---------------------------------------------------------
bin          hzz4l  
""")
d.write('observation '+str(len(events)))
d.write("""
---------------------------------------------------------
bin                   hzz4l  hzz4l  hzz4l  
process               ggH    ZZ      ZX    
process               0       1       2
""")
d.write('rate    '+str(info['higgs'])+"  "+str(info['ZZ'])+"   "+str(info['ZX']))
d.write("""
---------------------------------------------------------
BRhiggs_hzz4l  lnN    1.02   -        - 
leptonEff      lnN    1.04   1.04     -
QCDscale_ggH   lnN    1.08    -       - 
lumi_8TeV      lnN    1.044  1.044    -
pdf_gg         lnN    1.071   -       -
flatApprox     lnN    -       1.15    -
flatApprox     lnN    -       -       1.15
fakeRate       lnN    -       -       1.90
scale        param    0.0 0.003
resolution   param    0.0 0.2
--------------------------------------------------------
""")

d.close()

       







