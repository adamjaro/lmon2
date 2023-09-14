
import ROOT as rt
from ROOT import TLatex

class parameter_descriptor(object):
   #inserts list of parameters and values aligned to the equal sign
   #_____________________________________________________________________________
   def __init__(self, frame, x, y, sep, nsep=0.01):
      self.frame = frame
      #position and separations relative to the frame coordinates
      gxmin = self.frame.GetXaxis().GetXmin()
      gxmax = self.frame.GetXaxis().GetXmax()
      if frame.InheritsFrom(rt.TGraph.Class()) == True:
        gymin = self.frame.GetYaxis().GetXmin()
        gymax = self.frame.GetYaxis().GetXmax()
      else:
        gymin = self.frame.GetMinimum()
        gymax = self.frame.GetMaximum()
      self.xpos = gxmin + (gxmax-gxmin)*x # horizontal position of first equal sign
      self.ypos = gymin + (gymax-gymin)*y # vertical position of first equal sign
      self.itsep = (gymax-gymin)*sep # vertical separation of items
      self.eqsep = (gxmax-gxmin)*nsep # space between name of the variable and equal sign
      self.tsiz = 0.035 #text size
      self.prec = 3 # precision
      self.fmt = "f" # notation, fixed point on exponent
      self.lnam = TLatex() # put names
      self.lnam.SetTextFont(42)
      self.lnam.SetTextSize(self.tsiz)
      self.lnam.SetTextAlign(31)
      self.lval = TLatex() # put values
      self.lval.SetTextFont(42)
      self.lval.SetTextSize(self.tsiz)
      self.item_list = []

   #_____________________________________________________________________________
   def itemRes(self, nam, r1, ipar, col=rt.kBlack):
      #from TFitResult
      self.itemD(nam, r1.Parameter(ipar), r1.ParError(ipar), col)

   #_____________________________________________________________________________
   def itemR(self, nam, x, col=rt.kBlack):
      #from RooRealVar
      self.itemD(nam, x.getVal(), x.getError(), col)

   #_____________________________________________________________________________
   def itemD(self, nam, val, err=-1., col=rt.kBlack):
      #sval = "{0:.{1:d}f}".format(val, self.prec)
      sval = "{0:.{1:d}{2:s}}".format(val, self.prec, self.fmt)
      if err > -1.: sval = self.put_err(sval, err)
      self.item(nam, sval, col)

   #_____________________________________________________________________________
   def item(self, nam, val="", col=rt.kBlack):
      if val != "":
         val = "= " + val
      self.item_list.append([col, nam, val])

   #_____________________________________________________________________________
   def draw(self):
      for icnt in range(len(self.item_list)):
         putY = self.ypos - self.itsep*icnt # linear y axis
         self.lnam.SetTextColor(self.item_list[icnt][0])
         self.lnam.DrawLatex(self.xpos-self.eqsep, putY, self.item_list[icnt][1])
         if self.item_list[icnt][2] != "":
            self.lval.DrawLatex(self.xpos, putY, self.item_list[icnt][2])

   #_____________________________________________________________________________
   def put_err(self, val, err):
      val += " #pm {0:.{1:d}{2:s}}".format(err, self.prec, self.fmt)
      return val

   #_____________________________________________________________________________
   def set_text_size(self, siz):
      self.tsiz = siz
      self.lnam.SetTextSize(self.tsiz)
      self.lval.SetTextSize(self.tsiz)

























