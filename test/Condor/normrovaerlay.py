#! /usr/bin/env python                                                                                                                                                       
from array import array
from ROOT import TH1D,TFile,TCanvas, THStack,TF1, TH1,TLegend,kRed,kBlue,TPad,gPad,TLine,kBlack


path1 = '/uscms_data/d3/dmendis/80x/CMSSW_8_0_15/src/Analysis/VLQAna/test/Condor/Macro/mutest1/CR_Zmumu/'
path2 = '/uscms_data/d3/dmendis/80x/CMSSW_8_0_15/src/Analysis/VLQAna/test/Condor/Macro/mutest2/CR_Zmumu/'

f1 =TFile(path1+'nob_ht.root')
f2= TFile(path2+'1b_ht.root')
#c1 = TCanvas("c1")
c1 = TCanvas('c1', 'c1', 800, 600)
c1.Divide(2,1)
pad = c1.cd(1)
#pad.SetPad(0, 0.3, 1, 1)
#pad.SetTopMargin(0.1)
#pad.SetBottomMargin(0.005)
gPad.SetLogy()
h1= TH1D ()
h2 = TH1D()
h3=TH1D() 
h1= f1.Get('eeT0Z0H0b2__DY')
h3 = f2.Get('eeT0Z0H0b2__DY')

nbins1= h1.GetNbinsX()
nbins2 = h3.GetNbinsX()
print "nbins",nbins1
print "nbins2", nbins2
print " ********** 0 btag befor normalization*********************"
for i in range (1,nbins1+1):
    print h1.GetBinContent(i)

#for i in range(1,nbins+1):
 #   if i>2:
  #      x=h3.GetBinContent(i)
   #     h1.SetBinContent(i,x)
print " ****************** 1 BTAG before normalixation*************"

for i in range (1,nbins2+1):
    print h3.GetBinContent(i)

h2=h1.Clone('h2')
for i in range(1,nbins1+1):                                                                                                                                      
    if i==6:
        h2.SetBinContent(i,0)
    elif i>6:         
        x=h3.GetBinContent(i)                                                                                                                                   
        h2.SetBinContent(i,x)   


print " ****************** 1 BTAG ater settin bin vales*************"
for i in range (1,nbins2+1):
    print h2.GetBinContent(i)
h1.Scale(100/h1.Integral())
h2.Scale(100/h2.Integral())

print " ****************** o BTAG after normalixation*************"
for i in range (1,nbins2+1):
    print h1.GetBinContent(i)

print " ****************** 1 BTAG after  normalixation*************"
for i in range (1,nbins2+1):
    print h2.GetBinContent(i)




h1.GetYaxis().SetTitle("nob ht / 1b ht ");

#h1.SetOption("bar")
#h2.SetOption("bar")

h2.SetLineColor(kBlue);
h2.SetMarkerStyle(22)
h2.SetMarkerColorAlpha(kBlue, 0.35);
h2.Draw()

h1.SetLineColor(kRed);
h1.SetMarkerStyle(20)
h1.SetMarkerColorAlpha(kRed, 0.35);
h1.Draw("SAME")

h2.SetTitle("")
#print h1.GetOption();
#print h2.GetOption();
leg =TLegend(0.1,0.7,0.48,0.9)
#); // option "C" allows to center the header
leg.AddEntry(h1,"no btag","l")
leg.AddEntry(h2.SetMarkerStyle(20),"1 btag","l")
leg.Draw()



c1.cd(2)
h4= TH1D()
h4=h1.Clone('h4')
h4.Divide(h2)

h4.Draw()

line = TLine(0, 1, 3000, 1)
line.SetLineColor(kBlack)
line.Draw()

h4.SetTitle("ratio 0 btag/ >=1 btag")
h4.GetYaxis().SetTitle("ratio 0 btag/ >=1 btag")
raw_input("hold_on")
