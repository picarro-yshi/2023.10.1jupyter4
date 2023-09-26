## calibration factor for gas tank experiment, adapted from jupyter notebook
import os
import numpy as np
import time
import datetime
import h5py     #need by Windows
import tables   #need by windows
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Arial"

from spectral_logger1 import SpectralLogReader as slog
import GoldenCalibrationUtilities as GCU
from loadprivate import loadprivate

def calibration_gastank(fnr, gas, cid, tank_conc, t1, t2, t3, row=500, showgraph=False, savefig=False):
    gas_name = 'broadband_gasConcs_' + str(cid)  ## broadband_gasConcs_176
    cal_name = 'broadband_eCompoundOutputs_'+ str(cid) +'_calibration' ## broadband_eCompoundOutputs_176_calibration

    ## epoch_time = datetime.datetime(2021,11,24,8,0).timestamp()
    ta = datetime.datetime(int(t1[:4]), int(t1[4:6]), int(t1[6:8]), int(t1[8:10]), int(t1[10:])).timestamp()
    tc = datetime.datetime(int(t3[:4]), int(t3[4:6]), int(t3[6:8]), int(t3[8:10]), int(t3[10:])).timestamp()
    date = t1[:8]

    def h(name):
        j = ht.index(name)  # finds the index of the list ht that corresponds to 'name'
        return dat[:, j]
    ht, dat = loadprivate(fnr)

    x1 = h('time')
    y1 = h(gas_name) * 1e6

    idx = (x1 > ta+300) & (x1 < ta + 1560)  # zero1, 20 min
    x0 = x1[idx]
    start0 = x0[0]  # zero start epoch
    x0 = (x0 - x0[0]) / 60  # min, start at 0
    y0 = y1[idx]
    # print(len(x0))

    idx = (x1 > tc-1800) & (x1 < tc-300)       # IPA, time during which values are constant
    x = x1[idx]
    start = x[0]         # IPA start epoch
    x = (x - x[0]) / 60  # min, start at 0
    y = y1[idx]
    # print(len(x))
    # exit()

    gap = 5  # in minutes, look at data in packages of 5 min
    set_conc = []
    data_std = []
    yt = []  # ranged y in min
    xt = []  # ranged x in min

    # zero
    x_pt = int((x0[-1] - x0[0]) / gap)  # x axis points
    for j in range(x_pt):
        indices = (x0 >= j * gap) & (x0 <= (j + 1) * gap)
        if True in indices:
            set_conc.append(0)
            data_std.append(np.std(y0[indices]))
            yt.append(np.mean(y0[indices]))
            xt.append(x0[indices][-1])

    # gas
    x_pt = int((x[-1] - x[0]) / gap)  # x axis point
    for j in range(x_pt):
        indices = (x >= j * gap) & (x <= (j + 1) * gap)
        if True in indices:
            set_conc.append(tank_conc)
            data_std.append(np.std(y[indices]))
            yt.append(np.mean(y[indices]))
            xt.append(x[indices][-1] + int((start - start0) / 60))

    set_conc = np.array(set_conc)
    data_std = np.array(data_std)
    yt = np.array(yt)
    xt = np.array(xt)

    extra_bit = 0.1 / 1000.0
    zero_error = data_std[np.argmin(yt)]
    units = 'ppm'

    xsim = np.arange(0, np.nanmax(set_conc), extra_bit)
    a = np.polyfit(set_conc, yt, 1.0)
    ysim = np.polyval(a, xsim)
    resids = yt - np.polyval(a, set_conc)
    kf = 3

    if np.all(data_std == 0):
        print('std 0')
        MDL = GCU.calcStuffWeighted(set_conc, yt, zero_error, resids, units=units, k=kf)
    else:
        print('std no 0')
        MDL = GCU.calcStuffWeighted_points(set_conc, yt, data_std, zero_error, units=units, k=kf)
        # (X, Y, point by point uncertainty, uncertainty at 0, unit, sigma)

    set_cal = h(cal_name)[100] * 1e6  ## all value same
    fit = set_conc * MDL[0] + MDL[1]
    print(MDL[-2])
    slope = a[0]
    intercept = a[1]
    cal = set_cal / slope    #cal
    print('initial concentration/cal in library: %.4f ppm' % set_cal)
    if intercept >= 0:
        print('calibration: %.4f [actual] + %.4f' % (slope, intercept))
    else:
        print('calibration: %.4f [actual] - %.4f' % (slope, intercept))
    print('new concentration estimate: %.4f' % cal)

    fncal = os.path.join(fnr, 'par', 'calibration_factor.txt')
    with open(fncal, 'w') as f:
        f.write(str(cal))
    # exit()

    ########### plots #########
    ## Raw data
    F, A = plt.subplots(dpi=150)   #figsize=(6.5, 4.0)
    A.set_facecolor('#Ededee')
    A.grid(c='white')
    idx = (x1 > ta) & (x1 < tc)  # zero1, 30 min
    xtruc = list(x1[idx])
    ytruc = list(y1[idx])
    x_t = [xtruc[0]]    # xtruc data that will be marked
    xmak = [time.strftime('%H:%M', time.localtime(xtruc[0]))]

    n = len(xtruc)
    for i in range(1, n):
        clock0 = time.strftime('%M:%S', time.localtime(xtruc[i-1]))
        clock = time.strftime('%M:%S', time.localtime(xtruc[i]))
        if (clock0[1] == '9' and clock[1]=='0'):
            x_t.append(xtruc[i])
            xmak.append(time.strftime('%H:%M', time.localtime(xtruc[i])))

    A.plot(xtruc, ytruc, label='Raw')
    A.set_xticks(x_t, xmak)
    A.legend()
    A.set_xlabel('Clock time', fontsize=12)
    A.set_ylabel('Conc. (ppm)', fontsize=12)
    A.set_title('%s Raw Data\n%s'%(date, gas))
    F.autofmt_xdate()    #x-label will not overlap
    plt.show(block=False)
    F1 = F
    if savefig:
        plt.savefig(os.path.join(fnr, date +' Raw.png'), bbox_inches='tight')
    # exit()

    ## Fitting plot
    F, A = plt.subplots(dpi=150)   #figsize=(6.5, 4.0)
    A.set_facecolor('#Ededee')
    A.grid(c='white')
    A.plot(xt, yt, label='Measurement')
    A.plot(xt, set_conc, '*', label='Tank Conc. from Airgas')
    A.legend(loc=0, fontsize=10)
    A.set_xlabel('Time (min)', fontsize=12)
    A.set_ylabel('Conc. (ppm)', fontsize=12)
    A.set_title('Gas Tank Calibration: cal = %.3f\n%s' %(cal, gas), fontsize=14)
    plt.show(block=False)
    F2 = F
    if savefig:
        plt.savefig(os.path.join(fnr, date +' measured.png'), bbox_inches='tight')
    # exit()

    ### check if compound has the correct spectrum and RMSE
    read = slog(os.path.join(fnr, 'ComboResults'), verbose=True)
    data, _, max_row = read.get_spectra_row('broadband', row, pull_results=True)
    nu = data[0]['nu']
    k = data[0]['absorbance']
    residuals = data[0]['residuals']
    partial_fit = data[0]['partial_fit']
    model = data[0]['model']

    F, [A1, A2] = plt.subplots(2, 1, dpi=150)   #figsize=(10, 8)
    A1.set_facecolor('#Ededee')
    A2.set_facecolor('#Ededee')
    A1.grid(c='white', which='major', lw=1)
    A1.grid(c='white', which='minor', lw=0.5)
    A2.grid(c='white', which='major', lw=1)
    A2.grid(c='white', which='minor', lw=0.5)
    A1.minorticks_on()
    A2.minorticks_on()
    A1.plot(nu, k, label='data', color='#DC61D0', lw=1) #, alpha=0.75)
    A1.plot(nu, model, label='model', color='#00ffff', lw=1) #, alpha=0.4)
    A1.legend(loc=0, fontsize=10)
    A1.set_ylabel('Absorbance', fontsize=12)
    A1.set_title('%s Combo Log\n%s'%(date, gas))

    A2.plot(nu, residuals, label='residuals', color='#7068bb', lw=1)   #, alpha=0.75
    A2.set_ylabel('Residuals', color='#7068bb', fontsize=12)
    A4 = A2.twinx()
    A4.plot(nu, partial_fit, label='partial fit', color='#A6ce6a', lw=1 ) #, alpha=0.75)
    A2.set_xlabel('nu (cm-1)', fontsize=12)
    A4.set_ylabel('Partial fit', color='#A6ce6a', fontsize=12)
    A2.legend(loc=0, fontsize=10)
    F3 = F
    F4 = F
    if savefig:
        plt.savefig(os.path.join(fnr, date + ' RawSpectraFit.png'), bbox_inches='tight')
    plt.show(block=False)

    if showgraph:
        plt.waitforbuttonpress()  # press any key to close all plots
        plt.close()

    return F1, F2, F3, F4, max_row


if __name__ == "__main__":
    # basepath = r'/mnt/r/crd_G9000/AVXxx/3610-NUV1022/R&D/Calibration/'          ##Linux
    basepath = '/Volumes/Data/crd_G9000/AVXxx/3610-NUV1022/R&D/Calibration'       ##Mac
    # basepath = 'R:\crd_G9000\AVXxx\\3610-NUV1022\R&D\Calibration'               ## Windows

    gas = '702 - Ethanol'
    date = '20230209t3test'
    row1 = 1200

    fnr = os.path.join(basepath, gas, date)
    print(fnr)
    fnrp = os.path.join(fnr, 'par')

    if not os.path.exists(fnr):
        print('Error, did not find data. Please check if data exist or attach the data/R drive.')
    elif not os.path.exists(fnrp):
        print('Error, did not find experiment parameters.')
    else:
        f = open(os.path.join(fnrp, 't1.txt'), 'r')
        temp = f.read().splitlines()
        t1 = temp[0] + temp[1] + temp[2]

        f = open(os.path.join(fnrp, 't2.txt'), 'r')
        temp = f.read().splitlines()
        t2 = temp[0] + temp[1] + temp[2]

        f = open(os.path.join(fnrp, 't3.txt'), 'r')
        temp = f.read().splitlines()
        t3 = temp[0] + temp[1] + temp[2]
        print(t1, t2, t3)

        f = open(os.path.join(fnrp, 'cid.txt'), 'r')
        temp = f.read().splitlines()
        cid = int(temp[0])

        f = open(os.path.join(fnrp, 'tankconc.txt'), 'r')
        temp = f.read().splitlines()
        tank_conc = float(temp[0])

        calibration_gastank(fnr, gas, cid, tank_conc, t1, t2, t3, row1, showgraph=True)
        # calibration_gastank(fnr, gas, cid, tank_conc, t1, t2, t3, row1, showgraph=True, savefig=True)


