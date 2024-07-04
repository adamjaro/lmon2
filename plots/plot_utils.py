
from ctypes import c_double

import ROOT as rt
from ROOT import TMath, TH1D, TCanvas, TLegend, TLine, TIter, TH1, TH2D, TH2, TF2, TGraph
from ROOT import RooHist, TLatex, gROOT, TIter, TGraphErrors, TGaxis, TF1, TFrame, TH3D
from ROOT.Fit import FitResult
from ROOT import std, vector

#_____________________________________________________________________________
def prepare_TH1D(name, binsiz, xmin, xmax):

  #fixed bin size binsiz

  nbins, xmax = get_nbins(binsiz, xmin, xmax)

  return prepare_TH1D_n(name, nbins, xmin, xmax)

#_____________________________________________________________________________
def prepare_TH1D_mid(name, bin1, bin2, xmin, xmax, xmid):

  #bin size bin1 below xmid and bin2 above xmid

  return prepare_TH1D_vec(name, get_bins_vec_2pt(bin1, bin2, xmin, xmax, xmid))

#_____________________________________________________________________________
def prepare_TH1D_3pt(name, bin1, bin2, bin3, xmin, xmax, xmid1, xmid2):

  #binning with bin1 below xmid1, bin2 between xmid1 and xmid2 and bin3 above xmid2

  return prepare_TH1D_vec(name, get_bins_vec_3pt(bin1, bin2, bin3, xmin, xmax, xmid1, xmid2))

#_____________________________________________________________________________
def get_nbins(binsiz, xmin, xmax):

  nbins = int(TMath.Ceil( (xmax-xmin)/binsiz )) #round-up value
  xmax = xmin + float(binsiz*nbins) # move max up to pass the bins

  return nbins, xmax

#_____________________________________________________________________________
def get_bins_vec_2pt(bin1, bin2, xmin, xmax, xmid):

    #evaluate binning with bin1 width below xmid and bin2 above xmid
    bins = vector(rt.double)()
    bins.push_back(xmin)
    while True:
        if bins[bins.size()-1] < xmid:
            increment = bin1
        else:
            increment = bin2
        bins.push_back( bins[bins.size()-1] + increment )
        if bins[bins.size()-1] > xmax: break

    return bins

#_____________________________________________________________________________
def get_bins_vec_3pt(bin1, bin2, bin3, xmin, xmax, xmid1, xmid2):

    #evaluate binning with bin1 below xmid1, bin2 between xmid1 and xmid2
    #and bin3 above xmid2

    bins = vector(rt.double)()
    bins.push_back(xmin)

    while True:
        xpos = bins[bins.size()-1]

        increment = bin1
        if xpos > xmid1: increment = bin2
        if xpos > xmid2: increment = bin3

        bins.push_back( xpos + increment )
        if bins[bins.size()-1] > xmax: break

    return bins

#_____________________________________________________________________________
def prepare_TH1D_n(name, nbins, xmin, xmax):

  hx = TH1D(name, name, nbins, xmin, xmax)
  set_H1D(hx)
  hx.SetTitle("");

  return hx

#_____________________________________________________________________________
def prepare_TH1D_vec(name, vec):

    #th1d from vector<double>

    hx = TH1D(name, name, vec.size()-1, vec.data())
    set_H1D(hx)
    hx.SetTitle("");

    return hx

#_____________________________________________________________________________
def set_H1D(hx):

  hx.SetOption("E1");
  hx.SetMarkerStyle(rt.kFullCircle);
  #hx.SetMarkerColor(rt.kBlack)
  #hx.SetLineColor(rt.kBlack);
  set_H1D_col(hx, rt.kBlack)
  hx.SetLineWidth(2);
  hx.SetYTitle("Counts");
  siz = 0.035;
  hx.SetTitleSize(siz)
  hx.SetLabelSize(siz)
  hx.SetTitleSize(siz, "Y")
  hx.SetLabelSize(siz, "Y")
  hx.SetTitle("")

#_____________________________________________________________________________
def set_H1_text_size(hx, siz):

    hx.SetTitleSize(siz)
    hx.SetLabelSize(siz)
    hx.SetTitleSize(siz, "Y")
    hx.SetLabelSize(siz, "Y")

#_____________________________________________________________________________
def set_H1D_col(hx, col):

    hx.SetMarkerColor(col)
    hx.SetLineColor(col);

#_____________________________________________________________________________
def set_F1(fx, col=rt.kRed):

    fx.SetLineWidth(3)
    fx.SetNpx(1000)
    fx.SetTitle("")

    fx.SetLineColor(col)

#set_F1

#_____________________________________________________________________________
def set_axis(axis):

    axis.SetTextFont(42)
    axis.SetLabelFont(42)
    siz = 0.035
    axis.SetTitleSize(siz)
    axis.SetLabelSize(siz)

#_____________________________________________________________________________
def set_graph(tx, col=rt.kBlack, style=rt.kFullCircle):

    tx.SetMarkerStyle(style)
    tx.SetMarkerColor(col)
    tx.SetLineColor(col)
    tx.SetLineWidth(2)

#_____________________________________________________________________________
def set_line_graph(tx, col=rt.kBlack, style=rt.kSolid, width=2):

    tx.SetLineColor(col)
    tx.SetLineStyle(style)
    tx.SetLineWidth(width)

#_____________________________________________________________________________
def h1_to_graph(hx, col=None, style=rt.kSolid, width=2):

    tx = TGraphErrors(hx.GetNbinsX())
    for ibin in range(1,hx.GetNbinsX()+1):
        #print ibin, hx.GetBinContent(ibin)
        tx.SetPoint(ibin-1, hx.GetBinCenter(ibin), hx.GetBinContent(ibin))

    if col is not None:
        set_line_graph(tx, col, style, width)

    return tx

#_____________________________________________________________________________
def h1_to_graph_nz(hx, delt=0.001):

    #skip bins with zero content

    points = []
    for ibin in range(1,hx.GetNbinsX()+1):
        if hx.GetBinContent(ibin) < delt: continue

        points.append( (hx.GetBinCenter(ibin), hx.GetBinContent(ibin)) )

    tx = TGraphErrors(len(points))
    i = 0
    for p in points:
        tx.SetPoint(i, p[0], p[1])
        i += 1

    return tx

#_____________________________________________________________________________
def h1_to_arrays(hx):

    ifirst = -1
    ilast = -1
    for ibin in range(1, hx.GetNbinsX()+1):
        if( hx.GetBinContent(ibin) < 1e-12 ): continue

        if ifirst < 0: ifirst = ibin
        ilast = ibin

    xp = []
    yp = []
    #for ibin in range(1, hx.GetNbinsX()+1):
        #if( hx.GetBinContent(ibin) < 1e-12 ): continue
    for ibin in range(ifirst, ilast+1):

        x0 = hx.GetBinLowEdge(ibin)
        x1 = x0 + hx.GetBinWidth(ibin)

        xp.append( x0 )
        xp.append( x1 )

        yp.append( hx.GetBinContent(ibin) )
        yp.append( hx.GetBinContent(ibin) )

    xp = [xp[0]] + xp
    yp = [0] + yp

    xp.append(xp[-1])
    yp.append(0)

    return xp, yp

#_____________________________________________________________________________
def graph_to_arrays(gx):

    xp = []
    yp = []

    xv = c_double(0)
    yv = c_double(0)

    for i in range(gx.GetN()):

        gx.GetPoint(i, xv, yv)

        xp.append( xv.value - gx.GetErrorXlow(i) )
        xp.append( xv.value + gx.GetErrorXhigh(i) )

        yp.append( yv.value )
        yp.append( yv.value )

    return xp, yp

#_____________________________________________________________________________
def h1_to_arrays(hx):

    xp = []
    yp = []

    for ibin in range(1,hx.GetNbinsX()+1):

        xp.append( hx.GetBinLowEdge(ibin) )
        xp.append( hx.GetBinLowEdge(ibin)+hx.GetBinWidth(ibin) )

        yp.append( hx.GetBinContent(ibin) )
        yp.append( hx.GetBinContent(ibin) )

    return xp, yp

#_____________________________________________________________________________
def h1_to_arrays_centers(hx):

    xp = []
    yp = []

    for ibin in range(1,hx.GetNbinsX()+1):

        xp.append( hx.GetBinCenter(ibin) )
        #xp.append( hx.GetBinLowEdge(ibin)+hx.GetBinWidth(ibin) )

        yp.append( hx.GetBinContent(ibin) )
        #yp.append( hx.GetBinContent(ibin) )

    return xp, yp

#_____________________________________________________________________________
def prepare_TH2D(name, xbin, xmin, xmax, ybin, ymin, ymax):

  #bins along x and y
  nbinsX, xmax = get_nbins(xbin, xmin, xmax)
  nbinsY, ymax = get_nbins(ybin, ymin, ymax)

  return prepare_TH2D_n(name, nbinsX, xmin, xmax, nbinsY, ymin, ymax)

#_____________________________________________________________________________
def prepare_TH2D_n(name, nbinsX, xmin, xmax, nbinsY, ymin, ymax):

  hx = TH2D(name, name, nbinsX, xmin, xmax, nbinsY, ymin, ymax)
  hx.SetOption("COLZ");
  siz = 0.035
  hx.SetTitleSize(siz)
  hx.SetLabelSize(siz)
  hx.SetTitleSize(siz, "Y")
  hx.SetLabelSize(siz, "Y")
  hx.SetTitleSize(siz, "Z")
  hx.SetLabelSize(siz, "Z")

  hx.SetTitle("")

  return hx

#_____________________________________________________________________________
def prepare_TH3D(name, xbin, xmin, xmax, ybin, ymin, ymax, zbin, zmin, zmax):

    #bins along x, y and z
    nx, xmax = get_nbins(xbin, xmin, xmax)
    ny, ymax = get_nbins(ybin, ymin, ymax)
    nz, zmax = get_nbins(zbin, zmin, zmax)

    hx = TH3D(name, name, nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax)

    return hx

#_____________________________________________________________________________
def put_yx_tit(hx, ytit, xtit, yofs=1.5, xofs=1.1):

    hx.SetYTitle(ytit)
    hx.SetXTitle(xtit)

    hx.SetTitleOffset(yofs, "Y")
    hx.SetTitleOffset(xofs, "X")

#_____________________________________________________________________________
def put_frame_yx_tit(frame, ytit, xtit, yofs=1.7, xofs=1.2):

    frame.SetTitle("")

    frame.SetYTitle(ytit)
    frame.SetXTitle(xtit)

    frame.GetYaxis().SetTitleOffset(yofs);
    frame.GetXaxis().SetTitleOffset(xofs);


#_____________________________________________________________________________
def norm_to_data(hMC, hDat, col=rt.kBlue, lo=0., hi=-1.):

    #normalize MC hMC to data hDat, suppress errors to draw as line and set color,
    #optional range for data to be used to normalize the MC

    hMC.Sumw2()
    #norm in full or restricted range
    if hi < lo:
        hMC.Scale(hDat.Integral("width")/hMC.Integral("width"))
    else:
        b_lo = hDat.FindBin(lo)
        b_hi = hDat.FindBin(hi)
        hMC.Scale(hDat.Integral(b_lo, b_hi, "width")/hMC.Integral("width"))
    
    for ibin in range(hMC.GetNbinsX()+1):
        hMC.SetBinError(ibin, 0)

    hMC.SetLineColor(col)
    hMC.SetMarkerColor(col)

#_____________________________________________________________________________
def norm_to_num(hMC, num, col=rt.kBlue):

    #normalize hMC to a given number of entries num
    #suppress errors to draw as line and set color,

    hMC.Sumw2()
    hMC.Scale(num/hMC.Integral())

    for ibin in range(hMC.GetNbinsX()+1):
        hMC.SetBinError(ibin, 0)

    hMC.SetLineColor(col)
    hMC.SetMarkerColor(col)

#_____________________________________________________________________________
def norm_to_integral(hMC, ival, col=rt.kBlue, keep_err=False):

    #normalize hMC to integral given by ival
    #suppress errors to draw as line and set color,

    hMC.Sumw2()
    hMC.Scale(ival/hMC.Integral("width"))

    if not keep_err:
        for ibin in range(hMC.GetNbinsX()+1):
            hMC.SetBinError(ibin, 0)

    hMC.SetLineColor(col)
    hMC.SetMarkerColor(col)

#_____________________________________________________________________________
def norm_to_den_w(hx, den):

    #divide bins by denominator den and bin width

    for ibin in range(hx.GetNbinsX()+1):
        hx.SetBinContent(ibin, hx.GetBinContent(ibin)/(den*hx.GetBinWidth(ibin)))
        hx.SetBinError(ibin, hx.GetBinError(ibin)/(den*hx.GetBinWidth(ibin)))

    set_H1D(hx)

#_____________________________________________________________________________
def norm_to_integral_list(hlist, ival):

    for i in hlist:
        norm_to_den_w(i, 1)
        norm_to_integral(i, ival)

#_____________________________________________________________________________
def fill_h1_tf(hx, func, col=rt.kBlue):

    #fill h1 histogram from function

    for ibin in xrange(1,hx.GetNbinsX()+1):
        edge = hx.GetBinLowEdge(ibin)
        w = hx.GetBinWidth(ibin)
        hx.SetBinContent(ibin, func.Integral(edge, edge+w))
        hx.SetBinError(ibin, 0.)

    hx.SetLineColor(col)
    hx.SetMarkerColor(col)

#_____________________________________________________________________________
def line_h1(hx, col=rt.kBlue, wdt=-1):

    #set H1 to show as a line of bin content

    for ibin in range(1,hx.GetNbinsX()+1):
        hx.SetBinError(ibin, 0.)

    hx.SetLineColor(col)
    hx.SetMarkerColor(col)

    if wdt > 0: hx.SetLineWidth(wdt)

#_____________________________________________________________________________
def box_canvas(dx=768, dy=768):

    can = TCanvas("c3", "Analysis", dx, dy)
    rt.gStyle.SetOptStat("")
    rt.gStyle.SetPalette(1)
    rt.gStyle.SetLineWidth(2)
    rt.TGaxis.SetMaxDigits(4)

    can.cd(1)

    return can

#_____________________________________________________________________________
def set_margin_lbtr(gPad, lm, bm, tm, rm):

    gPad.SetLeftMargin(lm)
    gPad.SetBottomMargin(bm)
    gPad.SetTopMargin(tm)
    gPad.SetRightMargin(rm)

#_____________________________________________________________________________
def prepare_leg(xl, yl, dxl, dyl, tsiz=0.045):

  leg = TLegend(xl, yl, xl+dxl, yl+dyl)
  leg.SetFillStyle(0)
  leg.SetBorderSize(0)
  leg.SetTextSize(tsiz)

  return leg

#_____________________________________________________________________________
def add_leg_y_pt(leg, ymin, ymax, ptmax):

    leg.AddEntry(None, "#bf{%2.1f < #it{y} < %2.1f}" % (ymin, ymax), "")
    leg.AddEntry(None, "#bf{#it{p}_{T} < "+"{0:.2f}".format(ptmax)+" GeV}", "")

#_____________________________________________________________________________
def add_leg_pt_mass(leg, ptmax, mmin, mmax):

    leg.AddEntry(None, "#bf{#it{p}_{T} < "+"{0:.2f}".format(ptmax)+" GeV/c}", "")
    add_leg_mass(leg, mmin, mmax)

#_____________________________________________________________________________
def add_leg_mass(leg, mmin, mmax):

    mmin_fmt = "{0:.1f}".format(mmin)
    mmax_fmt = "{0:.1f}".format(mmax)
    leg.AddEntry(None, "#bf{"+mmin_fmt+" < #it{m}_{e^{+}e^{-}} < "+mmax_fmt+" GeV/c^{2}}", "")

#_____________________________________________________________________________
def make_uo_leg(hx, xl, yl, dxl, dyl, tsiz=0.03):

    leg = prepare_leg(xl, yl, dxl, dyl, tsiz)

    uolin = "Underflow: {0:.0f}, overflow: {1:.0f}".format(hx.GetBinContent(0), hx.GetBinContent(hx.GetNbinsX()+1))
    leg.AddEntry(None, uolin, "")

    return leg

#_____________________________________________________________________________
def col_lin(col, w=4, st=rt.kSolid):

  #create line of a given color
  lin = TLine()
  lin.SetLineColor(col)
  lin.SetLineWidth(w)
  lin.SetLineStyle(st)

  return lin

#_____________________________________________________________________________
def cut_line(cut_val, yl, hx, logy=False):

    #vertical line representing a cut value

    if logy == False:
        tsel_pos = yl*hx.GetMaximum()
    else:
        if hx.GetMinimum() > 0.:
            tsel_pos = TMath.Log10(hx.GetMinimum()) + yl*(TMath.Log10(hx.GetMaximum())-TMath.Log10(hx.GetMinimum()))
        else:
            tsel_pos = yl*TMath.Log10(hx.GetMaximum())
        tsel_pos = TMath.Power(10, tsel_pos)

    lin = TLine(cut_val, 0., cut_val, tsel_pos)
    lin.SetLineColor(rt.kViolet)
    lin.SetLineStyle(rt.kDashed)
    lin.SetLineWidth(4)

    return lin

#_____________________________________________________________________________
def log_fit_result(r1, lmg=6):

  result = ""

  ss = std.stringstream()
  r1.printMultiline(ss, 0, rt.kTRUE)
  stream = ss.str()

  #put parameter numbers
  line_stream = stream.split("\n")
  pnum = -2
  for i in range(len(line_stream)):
    if pnum >= 0 and line_stream[i] != "":
      numlin = "  " + str(pnum) + line_stream[i][2+len(str(pnum)):]
      line_stream[i] = numlin
      pnum += 1
    if pnum == -1:
      pnum += 1
    if line_stream[i].find("InitialValue") > -1:
      pnum += 1
    result += line_stream[i] + "\n"

  cor = r1.correlationMatrix()
  result += print_matrix(cor)

  return insert_left_margin(result, lmg)

#_____________________________________________________________________________
def log_tfit_result(r1, lmg=6):

    result = "Minimizer status: " + str(r1.Status()) + ", "
    result += "cov matrix status: " + str(r1.CovMatrixStatus()) + "\n"

    ss = std.stringstream()
    #call to base class FitResult::Print
    FitResult.Print(r1, ss, True)
    result += ss.str()

    return insert_left_margin(result, lmg)

#_____________________________________________________________________________
def insert_left_margin(res, lmg):

    #put left margin
    rline = res.split("\n")
    result = ""
    for line in rline:
        result += " ".ljust(lmg) + line + "\n"

    return result

#_____________________________________________________________________________
def log_fit_parameters(r1, lmg=6, prec=3):

    #direct access to fit parameters
    result = ""
    fmt = "{0:."+str(prec)+"f}"
    arglist = r1.floatParsFinal()
    idx = 0
    arg = arglist.at(idx)
    while arg != None:
        result += "".ljust(lmg) + arg.GetName() + " = "
        result += fmt.format(arg.getVal()) + " +/- "
        result += fmt.format(arg.getError())
        result += "\n"
        #move to next
        idx += 1
        arg = arglist.at(idx)

    return result

#_____________________________________________________________________________
def table_fit_parameters(r1):

    #LaTex table with fit parameters
    result = ""
    arglist = r1.floatParsFinal()
    idx = 0
    arg = arglist.at(idx)
    while arg != None:
        result += "$" + arg.GetName() + "$ & "
        result += "{0:.3f}".format(arg.getVal()) + " $\pm$ "
        result += "{0:.3f}".format(arg.getError()) + " \\\\"
        result += "\n"
        #move to next
        idx += 1
        arg = arglist.at(idx)

    return result

#_____________________________________________________________________________
def log_results(out, msg, lmg=6):

    for line in msg.split("\n"):
        out.write(" ".ljust(lmg) + line + "\n")

#_____________________________________________________________________________
def print_matrix(mat):

  result = "Correlation matrix:\n"

  line = "".rjust(4) + "| "
  for col in range(mat.GetNcols()):
    line += "{0:5d}  |".format(col)
  result += line + "\n"

  delim = "".ljust(5,"-")
  for col in range(mat.GetNcols()):
    delim += "".ljust(8,"-")
  delim += "-"
  result += delim + "\n"

  for row in range(mat.GetNrows()):
    line = "{0:3d} |".format(row)
    for col in range(mat.GetNcols()):
      line += "{0:8.3f}".format(mat(row,col))
    result += line + "\n"

  result += delim

  return result

#_____________________________________________________________________________
def make_log_string(*names):

    #prepare string to log parameter names and values

    strlog=""
    for iset in xrange(len(names)):
        for nam, val in names[iset]:
            strlog += nam + " " + str(val) + " "
        if iset < len(names)-1:
            strlog += "\n"

    return strlog

#_____________________________________________________________________________
def print_pad(pad):
  
  next = TIter(pad.GetListOfPrimitives())

  print("#####################")
  obj = next()
  while obj != None:
    print(obj.GetName(), obj.ClassName())
    obj = next()
  print("#####################")

#_____________________________________________________________________________
def invert_col(pad, bgcol=rt.kBlack):

   #set foreground and background color
   #fgcol = rt.kGreen
   fgcol = rt.kOrange-3

   pad.SetFillColor(bgcol)
   pad.SetFrameLineColor(fgcol)

   next = TIter(pad.GetListOfPrimitives())
   obj = next()
   while obj != None:
      #H1
      if obj.InheritsFrom(TH1.Class()) == True:
         if obj.GetLineColor() == rt.kBlack:
            obj.SetLineColor(fgcol)
            obj.SetFillColor(bgcol)
         if obj.GetMarkerColor() == rt.kBlack: obj.SetMarkerColor(fgcol)
         obj.SetAxisColor(fgcol, "X")
         obj.SetAxisColor(fgcol, "Y")
         obj.SetLabelColor(fgcol, "X")
         obj.SetLabelColor(fgcol, "Y")
         obj.GetXaxis().SetTitleColor(fgcol)
         obj.GetYaxis().SetTitleColor(fgcol)
      #Legend
      if obj.InheritsFrom(TLegend.Class()) == True:
         obj.SetTextColor(fgcol)
         #obj.SetFillStyle(1000)
         #obj.SetFillColor(fgcol)
         #obj.SetTextColor(bgcol)
         #ln = TIter(obj.GetListOfPrimitives())
         #lo = ln.Next()
         #while lo != None:
           #if lo.GetObject() == None:
             #lo = ln.Next()
             #continue
           #if lo.GetObject().InheritsFrom(TH1.Class()) == True:
             #hx = lo.GetObject()
             #hx.SetFillColor(bgcol)
             #if hx.GetMarkerColor() == rt.kBlack:
               #hx.SetMarkerColor(fgcol)
               #hx.SetLineColor(fgcol)
               #pass
           #lo = ln.Next()
      #RooHist
      if obj.InheritsFrom(RooHist.Class()) == True:
         if obj.GetMarkerColor() == rt.kBlack:
            obj.SetLineColor(fgcol)
            obj.SetMarkerColor(fgcol)
      #H2
      if obj.InheritsFrom(TH2.Class()) == True:
         obj.SetAxisColor(fgcol, "Z")
         obj.SetLabelColor(fgcol, "Z")
         obj.GetZaxis().SetTitleColor(fgcol)
         #obj.SetLineColor(fgcol)
         #obj.SetMarkerColor(fgcol)
      #TLatex
      if obj.InheritsFrom(TLatex.Class()) == True:
         if obj.GetTextColor() == rt.kBlack:
            obj.SetTextColor( fgcol )
      #F2
      if obj.InheritsFrom(TF2.Class()) == True:
        axes = [obj.GetXaxis(), obj.GetYaxis(), obj.GetZaxis()]
        for i in range(len(axes)):
            axes[i].SetAxisColor(fgcol)
            axes[i].SetLabelColor(fgcol)
            axes[i].SetTitleColor(fgcol)
      #F1
      if obj.InheritsFrom(TF1.Class()) == True:
        axes = [obj.GetXaxis(), obj.GetYaxis()]
        for i in range(len(axes)):
            axes[i].SetAxisColor(fgcol)
            axes[i].SetLabelColor(fgcol)
            axes[i].SetTitleColor(fgcol)
      #TGraph
      if obj.InheritsFrom(TGraph.Class()) == True:
            if obj.GetFillColor() == rt.kWhite:
                obj.SetFillColor(bgcol)
            ax = obj.GetXaxis()
            ay = obj.GetYaxis()
            ax.SetAxisColor(fgcol)
            ay.SetAxisColor(fgcol)
            ax.SetLabelColor(fgcol)
            ay.SetLabelColor(fgcol)
            ax.SetTitleColor(fgcol)
            ay.SetTitleColor(fgcol)
            if obj.GetLineColor() == rt.kBlack:
                obj.SetLineColor(fgcol)
                obj.SetMarkerColor(fgcol)
      #TGaxis
      if obj.InheritsFrom(TGaxis.Class()) == True:
            obj.SetLineColor(fgcol)
            obj.SetLabelColor(fgcol)
            obj.SetTitleColor(fgcol)

      #TFrame
      if obj.InheritsFrom( TFrame.Class() ) == True:
            if obj.GetLineColor() == rt.kBlack:
                obj.SetLineColor(fgcol)
                obj.SetFillColor(bgcol)

      #TLine
      if obj.InheritsFrom( TLine.Class() ) == True:
            if obj.GetLineColor() == rt.kBlack:
                obj.SetLineColor(fgcol)

      #move to next item
      obj = next()

#_____________________________________________________________________________
def invert_col_can(can):

   #set foreground and background color for the entire TCanvas

    bgcol = rt.kBlack
    can.SetFillColor(bgcol)

    i = 1
    while True:
        pad = can.GetPad(i)
        if type(pad) is rt.TVirtualPad: break

        invert_col(pad, bgcol)
        i += 1

#_____________________________________________________________________________
def frame_pow10_labels(frame, minl, maxl, xy="x", step=1, offset=-1):

    #frame axis labels in powers of 10

    if xy == "x":
        ax = frame.GetXaxis()

    if xy == "y":
        ax = frame.GetYaxis()

    labels = range(minl, maxl+1, step)
    for i in range(len(labels)):
        ax.ChangeLabel(i+1, -1, -1, -1, -1, -1, "10^{"+str(labels[i])+"}")

    if offset > 0:
        ax.SetLabelOffset(offset)

#frame_pow10_labels

#_____________________________________________________________________________
def frame_pow10_labels_float(frame, minl, maxl, xy="x", step=1, offset=-1):

    #frame axis labels in powers of 10 with floating point steps

    if xy == "x":
        ax = frame.GetXaxis()

    if xy == "y":
        ax = frame.GetYaxis()

    x = minl
    i = 0
    while x < maxl+step:

        ax.ChangeLabel(i+1, -1, -1, -1, -1, -1, "10^{"+str(x)+"}")

        x += step
        i += 1

    #labels = range(minl, maxl+1, step)
    #for i in range(len(labels)):
        #ax.ChangeLabel(i+1, -1, -1, -1, -1, -1, "10^{"+str(labels[i])+"}")

    if offset > 0:
        ax.SetLabelOffset(offset)

#frame_pow10_labels_float















