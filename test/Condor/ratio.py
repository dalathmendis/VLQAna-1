#! /usr/bin/env python                                                          
from array import array
from ROOT import TH1D,TFile,TCanvas, THStack,TF1, TH1

from ROOT import TH1D,TH2D,TFile,TMath,TCanvas,THStack,TLegend,TPave,TLine,TLatex, TF1, TGraph, TMultiGraph
from ROOT import gROOT,gStyle,gPad,gStyle


#path ='/uscms_data/d3/dmendis/80x/CMSSW_8_0_15/src/Analysis/VLQAna/test/Condor/Histo5Mu/'
path = '/uscms_data/d3/dmendis/80x/CMSSW_8_0_15/src/Analysis/VLQAna/test/Condor/Macro/electroncntnob/'
path2 = '/uscms_data/d3/dmendis/80x/CMSSW_8_0_15/src/Analysis/VLQAna/test/Condor/Macro/electroncntnob/'

c1= TCanvas()
f1 = TFile(path+'nob_ht.root')
f2 = TFile(path2+'nob_ht.root')
h1= TH1D()
h2=TH1D()
h3=TH1D()
h4=TH1D()
#h1 = f1.Get('ana/pre/ht_pre')
h1= f1.Get('eeT0Z0H0b2__DATA')
h2 = f2.Get('eeT0Z0H0b2__DY')

c1.Divide(1,2)
c1.cd(1)
#h1.Draw()
#h2.Draw()

hs = THStack()
hs.Add(h1)
hs.Add(h2)

hs.Draw('BAR')


c1.cd(2)
h3=h1.Clone('h3')
h3.Divide(h2)
h3.GetYaxis().SetTitle('Ratio Data/DY')
h3.GetYaxis().SetTitleFont(43)
h3.GetYaxis().SetTitleOffset(0.5)
h3.GetYaxis().SetTitleSize(20)
h3.SetMarkerStyle(20)

#for i in range(1,100):
  #  nbin = h3.GetBinContent(i)
  #  nbinX= h3.GetNbinsX()
  #  bincen = h3.GetBinCenter(i)
   # binlab=h3.GetBinLabel(i)
#f3 = TF1("myfit"," [0]+ ([1]*x) +([2]*x^2)+([3]*x^3)", 200,1500)
#f3 = TF1("myfit"," [0]", 1500,3000) 
f3 = TF1("myfit"," [0]+ [1]*x", 201,1200) 

h3.Fit("myfit","r")
h3.Draw()

raw_input("hod_on")
'''
    if bincen>2000:
        #print nbin,",",bincen #nbin, i
        #f3 = TF1("myfit"," [0]+ ([1]*x) +([2]*x^2)+([3]*x^3)", 200,2500)
        f4 = TF1("myfit1"," [0]", 2000,3000)
#f3 = TF1("myfit"," [0]+ ([1]*x) +([2]*x^2)", 200,2500)
#f3 = TF1("myfit"," [0]+ ([1]*x)", 201,1200) 

        h3.Fit("myfit1","r") 
        h3.Draw('')

raw_input("hod_on")
'''
