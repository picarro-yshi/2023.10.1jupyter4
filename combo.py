## Create Combo Log plot
import numpy as np
import os
import h5py     #need by Windows
import tables   #need by windows
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Arial"

from spectral_logger1 import SpectralLogReader as slog
import GoldenCalibrationUtilities as GCU
COMBOTAG = "../2023.10.2ui4/par1/combo_stop.txt"

def remove_single(x1, x2, y2):  # x1<x2, use x1 as the standard
    n1 = len(x1)
    n2 = len(x2)
    idx = []
    for i in range(n1):
        a = str(x1[i])[:7]
        b = str(x2[i])[:7]
        if a != b:
            idx.append(i)

        if len(idx) == n2 - n1:
            break

    x2 = np.delete(x2, idx)
    y2 = np.delete(y2, idx)

    # for i in range(n1):
    #     print(x1[i], x2[i])
    return x2, y2

def combo_study(fnr, gas, cid, day, r1, r2, r3, r4, peak_percent=10, savefig=False):
    note = ''
    read = slog(os.path.join(fnr, 'ComboResults'), verbose=True)
    _, _, max_row = read.get_spectra_row('broadband', 100, pull_results=True)
    print('max_row: ', max_row)

    if r1 > max_row or r2 > max_row or r3 > max_row or r4 > max_row:
        note = 'Row, Range numbers must be\nless than max_row: %s' % max_row
    else:
        ### RMSE vs concentration plot
        if r4 == 0:
            r4 = max_row
        x = list(range(r3, r4))
        ycombo = [0] * len(x)
        yconc = [0] * len(x)

        for j, row in enumerate(x):
            with open(COMBOTAG, "r") as f:
                temp = int(f.read())
            if temp:
                break

            a, b, _ = read.get_spectra_row('broadband', row, pull_results=True) #spectrum (arrays), results (data key), max_row
            # print(b.keys())
            # exit()
            residuals = a['residuals']
            partial_fit = a['partial_fit']

            n = len(residuals)
            f = 0
            r = 0
            for i in range(n):
                f += partial_fit[i] ** 2
                r += residuals[i] ** 2
            ycombo[j] = ((f / n) ** 0.5) / ((r / n) ** 0.5)
            yconc[j] = b['gasConcs_'+str(cid)] * 1000000
        print('finish loading data')

        col1 = 'steelblue'
        col2 = 'dimgrey'
        F, [ax1, ax2] = plt.subplots(2, 1, dpi=150)
        ax1.set_facecolor('#Ededee')
        ax2.set_facecolor('#Ededee')
        ax1.grid(c='white')
        ax2.grid(c='white')
        F.tight_layout()

        ax1.set_title('%s Combo Log: %s'%(day, gas))
        ax1.set_xlabel('row#')
        ax1.plot(x, ycombo, color=col1, lw=1)
        ax1.set_ylabel('Phi/RMSE(residue)', color=col1, fontsize=12)
        ax1.spines['right'].set_color(col2)

        ax1a = ax1.twinx()
        ax1a.plot(x, yconc, color=col2, lw=1)
        ax1a.set_ylabel('Concentration, ppm', color=col2, fontsize=12)
        ax1a.spines['left'].set_color(col1)

        ### partial fit of 2 combo log and their difference
        # option1: r1, r2
        if r1:
            idx1 = r1
            idx2 = r2
            peak = yconc[r1]
            peak2 = yconc[r2]
            label1 = str(r1)
            label2 = str(r2)

        # option2:peak and peak * percentage
        # find the peak index:
        else:
            peak = max(yconc)
            idx1_y = yconc.index(peak)
            idx1 = idx1_y + r3
            print(idx1, peak)
            print(len(yconc))

            peak2 = peak * peak_percent / 100
            a = min(enumerate(yconc[idx1_y:]), key=lambda xx: abs(xx[1] - peak2))
            idx2 = a[0] + idx1
            print(idx2, peak2)

            label1 = "peak"
            label2 = 'peak %s%%' % peak_percent

        ax1a.plot(idx1, peak, marker="o", markersize=10, markerfacecolor="red", alpha=0.5, markeredgecolor='none')
        ax1a.plot(idx2, peak2, marker="o", markersize=10, markerfacecolor="red", alpha=0.5, markeredgecolor='none')
        # ax1a.annotate('row %s: %s ppm' % (idx1, round(peak, 2)), (idx1, peak), fontsize=5)
        # ax1a.annotate('row %s: %s ppm' % (idx2, round(peak2, 2)), (idx2, peak2), fontsize=5)

        data = read.get_spectra_row('broadband', idx1, pull_results=True)
        nu1 = data[0]['nu']
        ab1 = data[0]['partial_fit']
        data = read.get_spectra_row('broadband', idx2, pull_results=True)
        nu2 = data[0]['nu']
        ab2 = data[0]['partial_fit']

        n1 = len(nu1)
        n2 = len(nu2)
        print(n1, n2)

        if n1 < n2:
            x1, y1 = nu1, ab1
            x2, y2 = remove_single(nu1, nu2, ab2)
        elif n2 < n1:
            x1, y1 = remove_single(nu2, nu1, ab1)
            x2, y2 = nu2, ab2
        else:
            x1, y1 = nu1, ab1
            x2, y2 = nu2, ab2

        mm = max(y1)
        y10 = y1/mm
        mm = max(y2)
        y20 = y2/mm

        ax2.plot(x1, y10, col1, lw=1, label=label1)
        ax2.legend()
        ax2.set_ylabel('Scaled Partial Fit', color=col1, fontsize=12)
        ax2.set_xlabel('nu (cm-1)')

        ax2.plot(x2, y20, col1, lw=1, linestyle='dashed', label=label2)
        ax2.legend()  # need this line to show legend
        ax2.spines['right'].set_color(col2)

        # spectrum of difference
        ax2a = ax2.twinx()
        ax2a.plot(x1, y10 - y20, col2, lw=1)
        ax2a.set_ylabel('row difference %s vs %s' % (idx1, idx2), color=col2, fontsize=12)
        ax2a.spines['left'].set_color(col1)

        if savefig:
            F.savefig(os.path.join(fnr, day + ' ComboLog.png'), bbox_inches='tight')
        plt.show()

    return note


if __name__ == "__main__":
    # basepath = r'/mnt/r/crd_G9000/AVXxx/3610-NUV1022/R&D/Calibration/'          ##Linux
    basepath = '/Volumes/Data/crd_G9000/AVXxx/3610-NUV1022/R&D/Calibration'       ##Mac
    # basepath = 'R:\crd_G9000\AVXxx\\3610-NUV1022\R&D\Calibration'               ## Windows

    gas = '9233 - Norbornane'
    date = '20230908test'
    peak_percent = 10  # % of peak height

    fnr = os.path.join(basepath, gas, date)
    print(fnr)
    fnrp = os.path.join(fnr, 'par')

    ## check data drive is attached or not
    if not os.path.exists(fnr):
        print('Error, did not find data. Please check if data exist or attach the data/R drive.')
    elif not os.path.exists(fnrp):
        print('Error, did not find experiment parameters.')
    else:
        f = open(os.path.join(fnrp, 't1.txt'), 'r')
        temp = f.read().splitlines()
        day = temp[0][:8]

        f = open(os.path.join(fnrp, 'cid.txt'), 'r')
        temp = f.read().splitlines()
        cid = int(temp[0])

        r1 = 100     # compare 2 rows
        r2 = 200
        r3 = 0      # 0 for minimum row#
        r4 = 1000      # 0 for maximum row#
        # combo_study(fnr, gas, cid, day, r1, r2, r3, r4)
        combo_study(fnr, gas, cid, day, r1, r2, r3, r4, savefig=True)

