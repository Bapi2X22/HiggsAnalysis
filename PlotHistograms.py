# import ROOT
# import numpy as np

# def plot_histograms(arrays, labels, colors=None, mode="single", title="",
#                     xlabel="", ylabel="", bins=100, xrange=(0, 100),
#                     outname="hist", outdir="default_dir", outfile="Higgs_plots.root"):

#     assert mode in ["single", "overlay", "ratio"], "Mode must be 'single', 'overlay', or 'ratio'"
#     assert len(arrays) == len(labels), "arrays and labels must be same length"

#     # Default color palette (cycled if arrays > len(colors))
#     default_colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2,
#                       ROOT.kMagenta+1, ROOT.kOrange+7, ROOT.kCyan+1]
#     if colors is None:
#         colors = [default_colors[i % len(default_colors)] for i in range(len(arrays))]
#     else:
#         assert len(colors) == len(arrays), "If provided, colors must match arrays count"

#     arrays = [np.asarray(a) for a in arrays]
#     nbins = bins
#     xmin, xmax = xrange
#     bin_edges = np.linspace(xmin, xmax, nbins + 1)

#     hists = []
#     for i, arr in enumerate(arrays):
#         hist = ROOT.TH1F(f"h_{i}", "", nbins, xmin, xmax)
#         for val in arr:
#             hist.Fill(val)
#         hist.SetLineColor(colors[i])
#         hist.SetLineWidth(2)
#         hist.SetStats(0)
#         hists.append(hist)

#     canvas = ROOT.TCanvas("canvas", title, 800, 800)

#     if mode == "ratio":
#         pad1 = ROOT.TPad("pad1", "Top pad", 0, 0.3, 1, 1.0)
#         pad1.SetBottomMargin(0.03)
#         pad1.Draw()
#         pad1.cd()
#     else:
#         canvas.cd()

#     hists[0].GetXaxis().SetTitle(xlabel if mode != "ratio" else "")
#     hists[0].GetYaxis().SetTitle(ylabel)
#     hists[0].Draw("HIST")
#     for h in hists[1:]:
#         h.Draw("HIST SAME")

#     legend = ROOT.TLegend(0.6, 0.75, 0.88, 0.88)
#     for h, l in zip(hists, labels):
#         legend.AddEntry(h, l, "l")
#     legend.Draw()

#     if mode == "ratio":
#         canvas.cd()
#         pad2 = ROOT.TPad("pad2", "Bottom pad", 0, 0.05, 1, 0.3)
#         pad2.SetTopMargin(0.03)
#         pad2.SetBottomMargin(0.3)
#         pad2.Draw()
#         pad2.cd()

#         h_ref = hists[0]
#         for i, h in enumerate(hists[1:], start=1):
#             h_ratio = h.Clone(f"h_ratio_{i}")
#             h_ratio.Divide(h_ref)
#             h_ratio.SetStats(0)
#             h_ratio.SetLineColor(colors[i])
#             h_ratio.GetYaxis().SetTitle(f"{labels[i]}/{labels[0]}")
#             h_ratio.GetYaxis().SetNdivisions(505)
#             h_ratio.GetYaxis().SetTitleSize(0.08)
#             h_ratio.GetYaxis().SetTitleOffset(0.4)
#             h_ratio.GetYaxis().SetLabelSize(0.08)
#             h_ratio.GetXaxis().SetTitle(xlabel)
#             h_ratio.GetXaxis().SetTitleSize(0.1)
#             h_ratio.GetXaxis().SetLabelSize(0.08)
#             drawopt = "HIST" if i == 1 else "HIST SAME"
#             h_ratio.Draw(drawopt)

#     canvas.Update()
#     canvas.SaveAs(f"{outname}.pdf")

#     # Save to ROOT file
#     f = ROOT.TFile(outfile, "UPDATE")
#     if not f.GetDirectory(outdir):
#         f.mkdir(outdir)
#     f.cd(outdir)
#     dir = ROOT.gDirectory
#     if dir.Get(outname):
#         dir.Delete(f"{outname};*")
#     canvas.Write(outname)
#     f.Close()

#     print(f"Saved canvas '{outname}' in directory '{outdir}' of {outfile}")


# import ROOT
# import numpy as np

# def plot_histograms(arrays, labels, colors=None, mode="single", title="",
#                     xlabel="", ylabel="", bins=100, xrange=(0, 100),
#                     outname="hist", outdir="default_dir", outfile="Higgs_plots.root"):
#     assert mode in ["single", "overlay", "ratio"], "Mode must be 'single', 'overlay', or 'ratio'"
#     assert len(arrays) == len(labels), "arrays and labels must be same length"

#     # Default color palette
#     default_colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2,
#                       ROOT.kMagenta+1, ROOT.kOrange+7, ROOT.kCyan+1]
#     if colors is None:
#         colors = [default_colors[i % len(default_colors)] for i in range(len(arrays))]
#     else:
#         assert len(colors) == len(arrays), "If provided, colors must match arrays count"

#     arrays = [np.asarray(a) for a in arrays]
#     nbins = bins
#     xmin, xmax = xrange

#     # Create histograms
#     hists = []
#     for i, arr in enumerate(arrays):
#         hist = ROOT.TH1F(f"h_{i}", "", nbins, xmin, xmax)
#         for val in arr:
#             hist.Fill(val)
#         hist.SetLineColor(colors[i])
#         hist.SetLineWidth(2)
#         hist.SetStats(1)  # Enable stat box
#         hists.append(hist)

#     canvas = ROOT.TCanvas("canvas", title, 800, 800)

#     if mode == "ratio":
#         pad1 = ROOT.TPad("pad1", "Top pad", 0, 0.3, 1, 1.0)
#         pad1.SetBottomMargin(0.03)
#         pad1.Draw()
#         pad1.cd()
#     else:
#         canvas.cd()

#     hists[0].GetXaxis().SetTitle(xlabel if mode != "ratio" else "")
#     hists[0].GetYaxis().SetTitle(ylabel)
#     hists[0].Draw("HIST")

#     for i, h in enumerate(hists[1:], start=1):
#         h.Draw("HIST SAME")

#     # Position stat boxes
#     canvas.Update()
#     # for i, h in enumerate(hists):
#     #     stat = h.FindObject("stats")
#     #     if stat:
#     #         stat.SetX1NDC(0.7)
#     #         stat.SetX2NDC(0.9)
#     #         y1 = 0.85 - i * 0.12
#     #         y2 = y1 + 0.1
#     #         stat.SetY1NDC(y1)
#     #         stat.SetY2NDC(y2)
#     #         stat.SetTextColor(h.GetLineColor())
#     #         stat.Draw()

#     # Draw legend
#     legend = ROOT.TLegend(0.6, 0.75, 0.88, 0.88)
#     for h, l in zip(hists, labels):
#         legend.AddEntry(h, l, "l")
#     legend.Draw()

#     # Ratio plot
#     if mode == "ratio":
#         canvas.cd()
#         pad2 = ROOT.TPad("pad2", "Bottom pad", 0, 0.05, 1, 0.3)
#         pad2.SetTopMargin(0.03)
#         pad2.SetBottomMargin(0.3)
#         pad2.Draw()
#         pad2.cd()

#         h_ref = hists[0]
#         for i, h in enumerate(hists[1:], start=1):
#             h_ratio = h.Clone(f"h_ratio_{i}")
#             h_ratio.Divide(h_ref)
#             h_ratio.SetStats(0)
#             h_ratio.SetLineColor(colors[i])
#             h_ratio.GetYaxis().SetTitle(f"{labels[i]}/{labels[0]}")
#             h_ratio.GetYaxis().SetNdivisions(505)
#             h_ratio.GetYaxis().SetTitleSize(0.08)
#             h_ratio.GetYaxis().SetTitleOffset(0.4)
#             h_ratio.GetYaxis().SetLabelSize(0.08)
#             h_ratio.GetXaxis().SetTitle(xlabel)
#             h_ratio.GetXaxis().SetTitleSize(0.1)
#             h_ratio.GetXaxis().SetLabelSize(0.08)
#             drawopt = "HIST" if i == 1 else "HIST SAME"
#             h_ratio.Draw(drawopt)

#     canvas.Update()
#     canvas.SaveAs(f"{outname}.pdf")

#     # Save to ROOT file
#     f = ROOT.TFile(outfile, "UPDATE")
#     if not f.GetDirectory(outdir):
#         f.mkdir(outdir)
#     f.cd(outdir)
#     dir = ROOT.gDirectory
#     if dir.Get(outname):
#         dir.Delete(f"{outname};*")
#     canvas.Write(outname)
#     f.Close()

#     print(f"Saved canvas '{outname}' in directory '{outdir}' of {outfile}")

import ROOT
import numpy as np

def plot_histograms(arrays, labels, colors=None, mode="single", title="",
                    xlabel="", ylabel="", bins=100, xrange=(0, 100),
                    outname="hist", outdir="default_dir", outfile="Higgs_plots.root"):
    assert mode in ["single", "overlay", "ratio"], "Mode must be 'single', 'overlay', or 'ratio'"
    assert len(arrays) == len(labels), "arrays and labels must be same length"

    # Default color palette
    default_colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2,
                      ROOT.kMagenta+1, ROOT.kOrange+7, ROOT.kCyan+1]
    if colors is None:
        colors = [default_colors[i % len(default_colors)] for i in range(len(arrays))]
    else:
        assert len(colors) == len(arrays), "If provided, colors must match arrays count"

    arrays = [np.asarray(a) for a in arrays]
    nbins = bins
    xmin, xmax = xrange

    ROOT.gStyle.SetOptStat(1110)  # Show stat box with entries, mean, RMS
    ROOT.gStyle.SetStatX(0.9)
    ROOT.gStyle.SetStatY(0.9)

    # Create histograms
    hists = []
    for i, arr in enumerate(arrays):
        hist = ROOT.TH1F(f"h_{i}", "", nbins, xmin, xmax)
        for val in arr:
            hist.Fill(val)
        hist.SetLineColor(colors[i])
        hist.SetLineWidth(2)
        hist.SetStats(1)  # Request stat box
        hists.append(hist)

    canvas = ROOT.TCanvas("canvas", title, 800, 800)

    if mode == "ratio":
        pad1 = ROOT.TPad("pad1", "Top pad", 0, 0.3, 1, 1.0)
        pad1.SetBottomMargin(0.03)
        pad1.Draw()
        pad1.cd()
    else:
        canvas.cd()

    hists[0].GetXaxis().SetTitle(xlabel if mode != "ratio" else "")
    hists[0].GetYaxis().SetTitle(ylabel)
    hists[0].Draw("HIST")

    # Draw additional histograms
    for h in hists[1:]:
        h.Draw("HIST SAME")

    # Redraw all histograms with stat boxes on top (trick)
    # for h in hists:
    #     h.Draw("SAME")
    hists[0].Draw("HIST")
    for h in hists[1:]:
        h.Draw("HIST SAME")

    # Position and style stat boxes
    canvas.Update()
    for i, h in enumerate(hists):
        stat = h.FindObject("stats")
        if stat:
            stat.SetX1NDC(0.7)
            stat.SetX2NDC(0.9)
            y1 = 0.85 - i * 0.15
            y2 = y1 + 0.12
            stat.SetY1NDC(y1)
            stat.SetY2NDC(y2)
            stat.SetTextColor(h.GetLineColor())
            stat.Draw()

    # Draw legend
    legend = ROOT.TLegend(0.6, 0.75, 0.88, 0.88)
    for h, l in zip(hists, labels):
        legend.AddEntry(h, l, "l")
    legend.Draw()

    # Ratio plot
    if mode == "ratio":
        canvas.cd()
        pad2 = ROOT.TPad("pad2", "Bottom pad", 0, 0.05, 1, 0.3)
        pad2.SetTopMargin(0.03)
        pad2.SetBottomMargin(0.3)
        pad2.Draw()
        pad2.cd()

        h_ref = hists[0]
        for i, h in enumerate(hists[1:], start=1):
            h_ratio = h.Clone(f"h_ratio_{i}")
            h_ratio.Divide(h_ref)
            h_ratio.SetStats(0)
            h_ratio.SetLineColor(colors[i])
            h_ratio.GetYaxis().SetTitle(f"{labels[i]}/{labels[0]}")
            h_ratio.GetYaxis().SetNdivisions(505)
            h_ratio.GetYaxis().SetTitleSize(0.08)
            h_ratio.GetYaxis().SetTitleOffset(0.4)
            h_ratio.GetYaxis().SetLabelSize(0.08)
            h_ratio.GetXaxis().SetTitle(xlabel)
            h_ratio.GetXaxis().SetTitleSize(0.1)
            h_ratio.GetXaxis().SetLabelSize(0.08)
            drawopt = "HIST" if i == 1 else "HIST SAME"
            h_ratio.Draw(drawopt)

    canvas.Update()
    canvas.SaveAs(f"{outname}.pdf")

    # Save to ROOT file
    f = ROOT.TFile(outfile, "UPDATE")
    if not f.GetDirectory(outdir):
        f.mkdir(outdir)
    f.cd(outdir)
    dir = ROOT.gDirectory
    if dir.Get(outname):
        dir.Delete(f"{outname};*")
    canvas.Write(outname)
    f.Close()

    print(f"Saved canvas '{outname}' in directory '{outdir}' of {outfile}")

