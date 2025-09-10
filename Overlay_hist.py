import ROOT
import sys
import argparse
import re
import os

def make_ratio(h1, h2):
    ratio = h1.Clone(f"{h1.GetName()}_ratio")
    ratio.Divide(h2)
    return ratio

def overlay_branches_from_dir(infile, outfile, directory, branches, ratio=False, rebin=None, xrange = None, title = None, save_name = None, normalize = True):
    dir_obj = infile.Get(directory)
    if not dir_obj:
        print(f"Directory '{directory}' not found in file.")
        sys.exit(1)

    # --- Collect histograms & find max y ---
    hists = []
    max_y = 0
    for branch in branches:
        hist = dir_obj.Get(branch)
        if not hist:
            print(f"Branch/histogram '{branch}' not found in directory '{directory}'.")
            sys.exit(1)
        if rebin:
            hist = hist.Rebin(rebin)
        if normalize:
            integral = hist.Integral()
            if integral > 0:
                hist.Scale(1.0 / integral)
        local_max = max(hist.GetBinContent(i) for i in range(1, hist.GetNbinsX()+1))
        if local_max > max_y:
            max_y = local_max
        hists.append(hist)

    max_y *= 1.2
    if max_y == 0:
        max_y = 1

    # --- Create canvas ---
    if ratio and len(hists) >= 2:
        c = ROOT.TCanvas(f"c_{directory}_branches", f"Overlay branches in {directory}", 800, 800)
        pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1)
        pad2 = ROOT.TPad("pad2", "pad2", 0, 0, 1, 0.3)
        pad1.SetBottomMargin(0.02)
        pad2.SetTopMargin(0.05)
        pad2.SetBottomMargin(0.3)
        pad1.SetRightMargin(0.03)
        pad2.SetRightMargin(0.03)
        pad2.SetLeftMargin(0.12)
        pad1.SetLeftMargin(0.12)
        pad1.Draw()
        pad2.Draw()
        pad1.cd()
    else:
        c = ROOT.TCanvas(f"c_{directory}_branches", f"Overlay branches in {directory}", 800, 600)
        pad1 = c

    # --- Colors & Legend ---
    colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2, ROOT.kMagenta, ROOT.kCyan+1, ROOT.kOrange+7]
    leg = ROOT.TLegend(0.7, 0.7, 0.88, 0.88)

    # --- Draw histograms ---
    hists[0].SetMaximum(max_y)
    if xrange:
        hists[0].GetXaxis().SetRangeUser(xrange[0], xrange[1])
    for i, hist in enumerate(hists):
        if xrange:
            hist.GetXaxis().SetRangeUser(xrange[0], xrange[1])
        hist.SetLineColor(colors[i % len(colors)])
        hist.SetMarkerColor(colors[i % len(colors)])
        hist.SetMarkerStyle(20 + i)
        hist.SetLineWidth(1)   # default is usually 2
        hist.SetMarkerSize(0.5)
        if title and i == 0:  # only set once
            hist.SetTitle(title)
        if i == 0:
            hist.Draw("HIST E")
        else:
            hist.Draw("HIST E SAME")
        leg.AddEntry(hist, hist.GetName(), "lep")

    leg.Draw()

    # --- Draw ratio plots ---
    if ratio and len(hists) >= 2:
        # Hide x-axis labels on upper pad
        for hist in hists:
            hist.GetXaxis().SetLabelSize(0)
            hist.GetXaxis().SetTitleSize(0)

        # leg_ratio = ROOT.TLegend(0.75, 0.15, 0.9, 0.35)  # adjust position as needed
        # leg_ratio.SetBorderSize(0)
        # leg_ratio.SetFillStyle(0)
        # leg_ratio.SetTextSize(0.08)

        ratio_colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2, ROOT.kMagenta, ROOT.kCyan+1, ROOT.kOrange+7]
        ratio_hists = []

        for i, hist in enumerate(hists[1:], start=1):
            rhist = hist.Clone(f"{hist.GetName()}_ratio")
            if xrange:
                rhist.GetXaxis().SetRangeUser(xrange[0], xrange[1])
            rhist.Divide(hists[0])
            rhist.SetLineColor(ratio_colors[(i-1) % len(ratio_colors)])
            rhist.SetMarkerColor(ratio_colors[(i-1) % len(ratio_colors)])
            rhist.SetMarkerStyle(20 + (i-1))
            ratio_hists.append(rhist)  # store without drawing
            # leg_ratio.AddEntry(rhist, hist.GetName(), "lep")

        # pad2.cd()
        # ROOT.gStyle.SetOptStat(0)

        # Draw empty frame first
        pad2.cd()
        ROOT.gStyle.SetOptStat(0)
        # Frame for ratio pad
        frame = hists[0].Clone("frame")
        frame.Reset()
        frame.SetTitle("")
        frame.GetYaxis().SetTitle("Ratio")
        frame.GetYaxis().SetRangeUser(0, 4)
        frame.GetYaxis().SetNdivisions(505)
        frame.GetYaxis().SetTitleSize(20)
        frame.GetYaxis().SetTitleFont(43)
        frame.GetYaxis().SetTitleOffset(1.55)
        frame.GetYaxis().SetLabelFont(43)
        frame.GetYaxis().SetLabelSize(20)
        frame.GetXaxis().SetTitleSize(20)
        frame.GetXaxis().SetTitleFont(43)
        frame.GetXaxis().SetTitleOffset(1.0)
        frame.GetXaxis().SetLabelFont(43)
        frame.GetXaxis().SetLabelSize(20)
        frame.Draw("AXIS")

        # Now draw stored ratios
        for rhist in ratio_hists:
            rhist.SetStats(0)
            rhist.Draw("E SAME")
        # leg_ratio.Draw()
        pad2.Update()

    # --- Save ---
    c.Update()
    outfile.cd()

    # add "_norm" suffix if normalization applied
    suffix = "_norm" if normalize else ""
    name_to_use = (save_name if save_name else c.GetName()) + suffix

    c.Write(name_to_use, ROOT.TObject.kOverwrite)
    print(f"Canvas '{name_to_use}' saved to output file '{outfile.GetName()}'.")

    outdir = "/eos/user/b/bbapi/www/ZHToAA/overlay_hist/"  # directory to store pdf/png
    pdf_path = os.path.join(outdir, f"{name_to_use}.pdf")
    png_path = os.path.join(outdir, f"{name_to_use}.png")

    c.SaveAs(pdf_path)
    c.SaveAs(png_path)

    print(f"Extra copies saved: {pdf_path}, {png_path}")


def overlay_same_branch_multiple_dirs(infile, outfile, directories, branch, ratio=False, rebin=None, xrange = None, title = None, save_name = None, normalize = True):
    # --- Collect histograms ---
    hists = []
    max_y = 0
    for d in directories:
        dir_obj = infile.Get(d)
        if not dir_obj:
            print(f"Directory '{d}' not found in file.")
            sys.exit(1)
        hist = dir_obj.Get(branch)
        if not hist:
            print(f"Branch '{branch}' not found in directory '{d}'.")
            sys.exit(1)

        hist.SetName(f"{branch}_{d}")  # rename so they are unique in memory
        if rebin:
            hist = hist.Rebin(rebin)
        if normalize:
            integral = hist.Integral()
            if integral > 0:
                hist.Scale(1.0 / integral)

        # Track max Y for axis scaling
        local_max = max(hist.GetBinContent(i) for i in range(1, hist.GetNbinsX() + 1))
        if local_max > max_y:
            max_y = local_max

        hists.append(hist)

    max_y *= 1.2
    if max_y == 0:
        max_y = 1

    # --- Create canvas & pads ---
    if ratio and len(hists) >= 2:
        c = ROOT.TCanvas(f"c_{branch}_dirs", f"Overlay {branch} from multiple directories", 800, 800)
        pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1)
        pad2 = ROOT.TPad("pad2", "pad2", 0, 0, 1, 0.3)
        pad1.SetBottomMargin(0.02)
        pad1.SetRightMargin(0.03)
        pad2.SetRightMargin(0.03)
        pad2.SetLeftMargin(0.12)
        pad1.SetLeftMargin(0.12)
        pad2.SetTopMargin(0.05)
        pad2.SetBottomMargin(0.3)
        pad1.Draw()
        pad2.Draw()
        pad1.cd()
    else:
        c = ROOT.TCanvas(f"c_{branch}_dirs", f"Overlay {branch} from multiple directories", 800, 600)
        pad1 = c

    # --- Colors & legend ---
    colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2, ROOT.kMagenta, ROOT.kCyan+1, ROOT.kOrange+7]
    leg = ROOT.TLegend(0.7, 0.7, 0.88, 0.88)

    # --- Draw histograms ---
    hists[0].SetMaximum(max_y)
    if xrange:
        hists[0].GetXaxis().SetRangeUser(xrange[0], xrange[1])
    for i, hist in enumerate(hists):
        # if xrange:
        #     hist.GetXaxis().SetRangeUser(xrange[0], xrange[1])
        hist.SetLineColor(colors[i % len(colors)])
        hist.SetMarkerColor(colors[i % len(colors)])
        hist.SetMarkerStyle(20 + i)
        hist.SetLineWidth(1)   # default is usually 2
        hist.SetMarkerSize(0.5)
        if title and i == 0:  # only set once
            hist.SetTitle(title)
        if i == 0:
            hist.Draw("HIST E")
        else:
            hist.Draw("HIST E SAME")
        # leg.AddEntry(hist, hist.GetName(), "lep")

        match = re.match(r"M(\d+).*?(20UL\d+)", directories[i])
        if match:
            mass_val = match.group(1)      # e.g., '25'
            year_ul  = match.group(2)      # e.g., '20UL18'
            # TLatex style: M_a = 25 GeV (20UL18)
            legend_name = f"#it{{M_{{a}}}} = {mass_val} GeV ({year_ul})"
        else:
            legend_name = directories[i]
        
        leg.AddEntry(hist, legend_name, "lep")
        leg.Draw()

    # --- Ratio plot ---
    if ratio and len(hists) >= 2:
        # Hide x labels in upper pad
        for hist in hists:
            hist.GetXaxis().SetLabelSize(0)
            hist.GetXaxis().SetTitleSize(0)

        pad2.cd()
        ROOT.gStyle.SetOptStat(0)

        # Draw empty frame for ratio axes
        # leg_ratio = ROOT.TLegend(0.75, 0.15, 0.9, 0.35)  # adjust position as needed
        # leg_ratio.SetBorderSize(0)
        # leg_ratio.SetFillStyle(0)
        # leg_ratio.SetTextSize(0.08)

        ratio_colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2, ROOT.kMagenta, ROOT.kCyan+1, ROOT.kOrange+7]
        ratio_hists = []

        for i, hist in enumerate(hists[1:], start=1):
            rhist = hist.Clone(f"{hist.GetName()}_ratio")
            if xrange:
                rhist.GetXaxis().SetRangeUser(xrange[0], xrange[1])
            rhist.Divide(hists[0])
            rhist.SetLineColor(ratio_colors[(i-1) % len(ratio_colors)])
            rhist.SetMarkerColor(ratio_colors[(i-1) % len(ratio_colors)])
            rhist.SetMarkerStyle(20 + (i-1))
            ratio_hists.append(rhist)  # store without drawing
            # leg_ratio.AddEntry(rhist, hist.GetName(), "lep")

        # pad2.cd()
        # ROOT.gStyle.SetOptStat(0)

        # # Draw empty frame first
        # frame = hists[0].Clone("frame")
        # frame.Reset()
        # frame.SetTitle("")
        # frame.GetYaxis().SetTitle("Ratio")
        # frame.GetYaxis().SetRangeUser(0, 6)
        # frame.Draw("AXIS")

        frame = hists[0].Clone("frame")
        frame.Reset()
        frame.SetTitle("")
        frame.GetYaxis().SetTitle("Ratio")
        frame.GetXaxis().SetTitle("pT [GeV]")
        frame.GetYaxis().SetRangeUser(0, 4)
        frame.GetYaxis().SetNdivisions(505)
        frame.GetYaxis().SetTitleSize(20)
        frame.GetYaxis().SetTitleFont(43)
        frame.GetYaxis().SetTitleOffset(1.55)
        frame.GetYaxis().SetLabelFont(43)
        frame.GetYaxis().SetLabelSize(20)
        frame.GetXaxis().SetTitleSize(20)
        frame.GetXaxis().SetTitleFont(43)
        frame.GetXaxis().SetTitleOffset(1.0)
        frame.GetXaxis().SetLabelFont(43)
        frame.GetXaxis().SetLabelSize(20)
        frame.Draw("AXIS")
        # frame.OptStat(0)

        # Now draw stored ratios
        for rhist in ratio_hists:
            rhist.SetStats(0)
            rhist.Draw("E SAME")
        # leg_ratio.Draw()

        pad2.Update()

    # --- Save ---
    c.Update()
    outfile.cd()

    # add "_norm" suffix if normalization applied
    suffix = "_norm" if normalize else ""
    name_to_use = (save_name if save_name else c.GetName()) + suffix

    c.Write(name_to_use, ROOT.TObject.kOverwrite)
    print(f"Canvas '{name_to_use}' saved to output file '{outfile.GetName()}'.")

    outdir = "/eos/user/b/bbapi/www/ZHToAA/overlay_hist/"  # directory to store pdf/png
    pdf_path = os.path.join(outdir, f"{name_to_use}.pdf")
    png_path = os.path.join(outdir, f"{name_to_use}.png")

    c.SaveAs(pdf_path)
    c.SaveAs(png_path)

    print(f"Extra copies saved: {pdf_path}, {png_path}")



def main():
    parser = argparse.ArgumentParser(description="Overlay ROOT histograms with optional ratio plot and save to separate output ROOT file.")
    parser.add_argument("rootfile", help="Input ROOT file name")
    parser.add_argument("mode", choices=["branches", "dirs"], help="Mode: overlay branches in one dir or same branch in multiple dirs")
    parser.add_argument("args", nargs="+", help="Arguments: For branches mode: directory branch1 branch2 ...; for dirs mode: dir1 dir2 ... branch")
    parser.add_argument("--ratio", action="store_true", help="Enable ratio plot (first two histograms)")
    parser.add_argument("--outfile", default="overlay_output_ZH_skipped.root", help="Output ROOT file to save canvases (default: overlay_output.root)")
    parser.add_argument("--rebin", type=int, default=None, help="Rebin histograms by this factor before overlaying")
    parser.add_argument("--xrange", type=float, nargs=2, metavar=("XMIN", "XMAX"),help="Set x-axis range (e.g., --xrange 0 200)")
    parser.add_argument("--title", type=str, default=None,help="Custom plot title (leave empty for no title)")
    parser.add_argument("--save_name", type=str, default=None,help="Custom plot save name")


    args = parser.parse_args()

    infile = ROOT.TFile.Open(args.rootfile)
    if not infile or infile.IsZombie():
        print(f"Cannot open input file {args.rootfile}")
        sys.exit(1)

    outfile = ROOT.TFile.Open(args.outfile, "UPDATE")
    if not outfile or outfile.IsZombie():
        print(f"Cannot open output file {args.outfile} in UPDATE mode")
        sys.exit(1)

    if args.mode == "branches":
        if len(args.args) < 2:
            print("For 'branches' mode, provide directory and at least one branch")
            sys.exit(1)
        directory = args.args[0]
        branches = args.args[1:]
        overlay_branches_from_dir(infile, outfile, directory, branches, args.ratio, rebin=args.rebin, xrange=args.xrange, title=args.title, save_name=args.save_name)

    elif args.mode == "dirs":
        if len(args.args) < 2:
            print("For 'dirs' mode, provide at least two directories and one branch name at the end")
            sys.exit(1)
        directories = args.args[:-1]
        branch = args.args[-1]
        overlay_same_branch_multiple_dirs(infile, outfile, directories, branch, args.ratio, rebin=args.rebin, xrange=args.xrange, title=args.title, save_name=args.save_name)

    outfile.Close()
    infile.Close()

if __name__ == "__main__":
    main()

