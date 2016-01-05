#!/usr/bin/env python

import matplotlib.pyplot as plt
from  matplotlib.ticker import ScalarFormatter
import matplotlib.gridspec as gridspec
import sys
import numpy as np

if len(sys.argv) < 2:
    print "Usage: %s file1.bed.mid [file2.bed.mid [...]]" % sys.argv[0]
    #./scripts/beaf_figure5.py data/*.bed.mid
    quit()

orderedNames = [
    "kib",  "ss",
    "mer",  "Otd",
    "aPKC", "Sens",
    "lgl",  "Pph13",
    "hpo",  "TJ",
    "sav",  "melt",
    "mats", "rh5",
    "wts",  "rh6",
    "yki",  "",
    "sd",     
]
peaks = {}
peaks_n = {}

for gene in orderedNames:
    peaks[gene] = []
    peaks_n[gene] = []

for g_name in sys.argv[1:]:
    print g_name
    g = open(g_name)
    f   = open("~/beaf32-hippo-chip/flyBase_genes_dm3.bed")
    goi = open("~/beaf32-hippo-chip/hippo_genes.txt")
    
    f_name = ""

    data = []
    titles = { "kib":"Hippo Pathay", "ss":"Non Hippo Pathway"}

    
    name2sym = {}

    for line1 in f:
        goi.seek(0)
        for gName in goi:
            geneName = gName.split()[0]
            geneSym = gName.split()[1]
            name2sym[geneName] = geneSym
            if geneSym == line1.split()[3].split("-")[0]:
                if line1.startswith("chr"):
                    fields1 = line1.split()
                    if fields1[3].split("-")[0] == f_name:
                        pass
                    else:    
                        f_chr = fields1[0]
                        f_strt = eval(fields1[1])
                        f_end = eval(fields1[2])
                        f_name = fields1[3].split("-")[0]
                        g.seek(0)
                        for line2 in g:
                            if line2.startswith("chr"):
                                fields2 = line2.split()
                                g_chr = fields2[0]
                                g_strt = eval(fields2[1])
                                g_end = eval(fields2[2])
                                g_val = fields2[4]
                                if g_chr == f_chr:
                                    if (f_strt-1000) > g_end or g_strt > (f_end+1000):
                                        pass
                                    elif (((f_strt-1000) >= g_strt) and ((f_strt-1000) < g_end)) or \
                                        (((f_end+1000) >= g_strt) and ((f_end+1000) < g_end)) or \
                                        ((g_strt >= (f_strt-1000)) and (g_strt < (f_end+1000))):
                                        peaks[geneName].append((g_strt, float(g_val)))
                                        peaks_n[geneName].append(g_name)
                                        print "\t".join([f_chr, str(f_strt), str(f_end), geneName, g_val, str(g_strt), fields1[5]])
                                        #break

                        data.append([geneName, f_chr, f_strt, f_end, f_chr, fields1[5] ])
    f.close()
    g.close()
    goi.close()
    
    ###################################################################################################

f = open("~/beaf32-hippo-chip/hippo_gene_exons.tsv")

fig = plt.figure(figsize=[20,10])
fig.subplots_adjust(
    left = 0.1,
    bottom = 0.01,
    right = 0.99,
    top = 0.93,
    wspace = 0.1,
    hspace = 0.2,
)

n = len(orderedNames)
if n%2 != 0:
    n = n+1

maxRange = 0
for gene in data:
    curRange = gene[3]-gene[2] + 2000
    if maxRange < curRange:
        maxRange = curRange

print

for i, name in enumerate(orderedNames):
    for gene in data:
        if gene[0] == name:
            plt.subplot(n//2,2,i+1)
            rstart = gene[2]-1000
            rend = gene[3]+1000
            chrom = gene[-2]
            rmid = rstart + (rend-rstart)/2
            strand = gene[-1]
        
            if strand == '+':
                pend = gene[2]
                pstart = pend - 500
                gstart = gene[2]
                gend = gene[3]
                plt.axis([rstart,rstart+maxRange,-0.5,1.8])
                l1 = "%s: %d" % (chrom, gstart)
                l2 = "%d" % rend
                labels = [l1,l2]
                plt.xticks([gstart], [l1])
                plt.arrow(gstart,1.5,500,0,color='k',head_width=0.4,head_length=200,width=0.1)
            elif strand == '-':
                pstart = gene[3]
                pend = pstart + 500
                gstart = gene[3]
                gend = gene[2]
                plt.axis([rend,rend-maxRange,-0.5,1.8])
                l1 = "%d" % rstart
                l2 = "%s: %d" % (chrom, gstart)
                labels = [l1,l2]
                plt.xticks([gstart], [l2])
                plt.arrow(gstart,1.5,-500,0,color='k',head_width=0.4,head_length=200,width=0.1)
        
            plt.gca()
            plt.gca().add_patch(plt.Rectangle((rstart,-0.5),(rend-rstart),1.0,facecolor="white",edgecolor='k',zorder=1, lw=2))
            f.seek(0)
            for line in f:
                fields=line.split()
                if line.startswith("#"):
                    pass
                elif name2sym[name] == fields[0].split("-")[0]:
                    bbb = zip(eval(fields[3]), eval(fields[4]))
                    for j, region in enumerate(bbb):
                        plt.gca().add_patch(plt.Rectangle((region[0],-0.5),(region[1]-region[0]),1.0,facecolor='dodgerblue',edgecolor='k',zorder=3, lw=2))
                        if j < (len(bbb)-1):
                            mid = (bbb[j+1][0] - region[1])//2 + region[1]
                            Xs = [region[1], mid, bbb[j+1][0]]
                            Ys = [0.5,1.5,0.5]
                            plt.plot(Xs,Ys,'k-', lw=2)
            plt.plot([gstart,gstart],[0.5,1.5],'k-',linewidth=2)
            reg    = []
            pro    = []
            exon   = []
            intron = []

            #**********************************************#
            
            r = rstart
            ranges = [r]
            while r < rend:
                r += 500
                ranges.append(r)
            
            rcounts = {}
            for r in ranges[:-1]:
                rcounts[r] = 0
            
            for kk, (peak, score) in enumerate(peaks[name]):
                if score < 2.5:
                    plt.plot(peak,0,'d',ms=15,mfc='none', alpha=0.8, zorder=11)
                else:
                    plt.plot(peak,0,'d',ms=15,mfc='red', alpha=0.8, zorder=11)
                if rstart <= peak <= rend:
                    if peaks_n[name][kk] not in reg:
                        reg.append(peaks_n[name][kk])
                if pstart <= peak <= pend:
                    if peaks_n[name][kk] not in pro:
                        pro.append(peaks_n[name][kk])
                for j, region in enumerate(bbb):
                    if region[0] < peak < region[1]:
                        if peaks_n[name][kk] not in exon:
                            exon.append(peaks_n[name][kk])
                    if j < (len(bbb)-1):
                        if region[1] <= peak <= bbb[j+1]:
                            if peaks_n[name][kk] not in intron:
                                intron.append(peaks_n[name][kk])
                for j, r in enumerate(ranges):
                    if peak > r:
                        pass
                    elif j == 0:
                        pass
                    else:
                        rcounts[ranges[j-1]] += 1
                        break
            reg = str(len(reg))
            pro = str(len(pro))
            exon = str(len(exon))
            intron = str(len(intron))

            #**********************************************#

            plt.gca().axes.get_yaxis().set_ticks([])
            plt.gca().patch.set_alpha(0.0)
            plt.gca().spines['top'].set_visible(False)
            plt.gca().spines['bottom'].set_visible(False)
            plt.gca().spines['left'].set_visible(False)
            plt.gca().spines['right'].set_visible(False)
            plt.gca().axes.get_xaxis().set_visible(False)
            plt.gca().axes.get_yaxis().set_visible(True)


plt.subplot(n//2,2,n-2)
plt.gca()
plt.gca().axis('off')
plt.plot([0,5000],[1,1],'k-', lw=3)
plt.axis([maxRange,0,0.5,1.5])
plt.text(3000,1.1,'5 kb', fontsize = 18, weight='bold')
plt.gca().patch.set_alpha(0.0)

fig.savefig("~/beaf32-hippo-chip/beafFigure.png", dpi=600)
quit()
plt.show()
