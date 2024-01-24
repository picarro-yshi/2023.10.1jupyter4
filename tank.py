## calibration factor for gas tank experiment, adapted from jupyter notebook
import os
import time
import platform
import numpy as np
import h5py     #need by Windows
import tables   #need by windows
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Arial"

from spectral_logger1 import SpectralLogReader as slog
import GoldenCalibrationUtilities as GCU
from loadprivate import loadprivate

def calibration_gastank(fnr, gas, cid, tank_conc, t1, t2, t3, row=100, showgraph=False, savefig=False):
    gas_name = 'broadband_gasConcs_' + str(cid)
    cal_name = 'broadband_eCompoundOutputs_'+ str(cid) +'_calibration'
    date = t1[:8]

    # t1 format: 202312110941
    epo1 = int(time.mktime(time.strptime(t1, "%Y%m%d%H%M")))
    epo3 = int(time.mktime(time.strptime(t3, "%Y%m%d%H%M")))

    def h(name):
        j = ht.index(name)  # finds the index of the list ht that corresponds to 'name'
        return dat[:, j]
    ht, dat = loadprivate(fnr)

    x1 = h('time')
    y1 = h(gas_name) * 1e6
    idx = (x1 > epo1 + 300) & (x1 < epo1 + 1560)  # zero1, 20 min
    x0 = x1[idx]
    start0 = x0[0]  # zero start epoch
    x0 = (x0 - x0[0]) / 60  # min, start at 0
    y0 = y1[idx]

    idx = (x1 > epo3 - 1800) & (x1 < epo3 - 300) # time during which values are constant
    x2 = x1[idx]
    start = x2[0]  # start epoch
    x2 = (x2 - x2[0]) / 60  # min, start at 0
    y = y1[idx]

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
    x_pt = int((x2[-1] - x2[0]) / gap)  # x axis point
    for j in range(x_pt):
        indices = (x2 >= j * gap) & (x2 <= (j + 1) * gap)
        if True in indices:
            set_conc.append(tank_conc)
            data_std.append(np.std(y[indices]))
            yt.append(np.mean(y[indices]))
            xt.append(x2[indices][-1] + int((start - start0) / 60))

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
    cal = set_cal / slope

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

    ########### Generate x axis tick marks: every 10 mins #########
    # main axis, bottom
    idx = (x1 > epo1) & (x1 < epo3) # time during which values are constant
    x = x1[idx]
    y = y1[idx]
    x_t = [int(x[0])]  # epoch time
    xmak = [time.strftime('%H:%M', time.localtime(x[0]))]  # mark

    n= len(x)
    for i in range(1, n):
        clock0 = time.strftime('%M:%S', time.localtime(x[i-1]))
        clock =  time.strftime('%M:%S', time.localtime(x[i]))
        if (clock0[1] == '9' and clock[1] == '0'):
            x_t.append(int(x[i]))
            xmak.append(time.strftime('%H:%M', time.localtime(x[i])))

    x_t.append(x[-1])
    xmak.append(time.strftime('%H:%M', time.localtime(x[-1])))

    # secondary axis, top
    x_t2 = [int(x[0])]
    xmak2 = ["0"]
    counter = 1
    for i in range(1, n):
        t = x[i] - x[0]  # s
        if t > counter * 600 > (x[i - 1] - x[0]):  # 10 min
            x_t2.append(x[i])
            xmak2.append(str(int(t/60)))
            counter += 1

    x_t2.append(x[-1])
    xmak2.append(str(int((x[-1] - x[0])/60)))

    ########### plots #########
    #1 Raw data
    F, A = plt.subplots(dpi=150)   #figsize=(6.5, 4.0)
    A.set_facecolor('#Ededee')
    A.grid(c='white')

    A.plot(x, y, label='Raw')
    A.legend()
    A.set_xlabel('Clock time', fontsize=12)
    A.set_ylabel('Conc. (ppm)', fontsize=12)
    A.set_title('%s Raw Data\n%s'%(date, gas))

    A.set_xticks(x_t, xmak, fontsize=10, rotation=40, ha='right')
    A.tick_params(pad=0)

    # top axis for min
    ax2 = A.twiny()
    ax2.set_xticks(x_t2, xmak2, fontsize=8)
    ax2.tick_params(length=2, grid_alpha=0.5, pad=-1)
    ax2.plot(x, y, linewidth= 0)
    ax2.set_xlabel("min:", fontsize=8)
    ax2.xaxis.set_label_coords(-0.01, 1.01)

    # F.autofmt_xdate()    #x-label will not overlap
    plt.show(block=False)
    F1 = F
    if savefig:
        plt.savefig(os.path.join(fnr, date +' Raw.png'), bbox_inches='tight')
    # exit()

    #2 Fitting plot
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

    #3 Combo plot, check if compound has the correct spectrum and RMSE
    read = slog(os.path.join(fnr, 'ComboResults'), verbose=True)
    data, _, max_row = read.get_spectra_row('broadband', row, pull_results=True)
    nu = data['nu']
    k = data['absorbance']
    residuals = data['residuals']
    partial_fit = data['partial_fit']
    model = data['model']

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
        # press q to close figures 1 by 1
        plt.show()

    return F1, F2, F3, F4, max_row


if __name__ == "__main__":
    if platform.system() == "Linux":
        basepath = r'/mnt/r/crd_G9000/AVXxx/3610-NUV1022/R&D/Calibration/'
    elif platform.system() == "Windows":
        basepath = 'R:/crd_G9000/AVXxx/3610-NUV1022/R&D/Calibration'
    else:
        basepath = '/Volumes/Data/crd_G9000/AVXxx/3610-NUV1022/R&D/Calibration'

    # for test
    gas = '702 - Ethanol'
    date = '20230209t3test'
    row1 = 1200

    # gas = '6574 - 1,1,2-Trichloroethane'
    # date = '20231026t2bad'
    # row1 = 30

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


