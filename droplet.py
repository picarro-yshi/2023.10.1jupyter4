## calibration factor for droplet experiment, adapted from jupyter notebook
import os
import numpy as np
import time
import datetime
import h5py  # need by Windows
import tables  # need by windows
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Arial"

from spectral_logger1 import SpectralLogReader as slog
from loadprivate import loadprivate

## all parameters must be valid before sending to below functions; no sanity check
## need .h5 files in PrvateData/broadband and ComboResults folder
def calibration_droplet(fnr, gas, cid, weight, MW, t1, t2, t3, pct = 4, row=500, showgraph=False, savefig=False):
    gas_name = 'broadband_gasConcs_' + str(cid)  ## broadband_gasConcs_176
    cal_name = 'broadband_eCompoundOutputs_'+ str(cid) +'_calibration' ## broadband_eCompoundOutputs_176_calibration

    ## epoch_time = datetime.datetime(2021,11,24,8,0).timestamp()
    ta = datetime.datetime(int(t1[:4]), int(t1[4:6]), int(t1[6:8]), int(t1[8:10]), int(t1[10:])).timestamp()
    tb = datetime.datetime(int(t2[:4]), int(t2[4:6]), int(t2[6:8]), int(t2[8:10]), int(t2[10:])).timestamp()
    tc = datetime.datetime(int(t3[:4]), int(t3[4:6]), int(t3[6:8]), int(t3[8:10]), int(t3[10:])).timestamp()
    date = t1[:8]

    def h(name):
        j = ht.index(name)  # finds the index of the list ht that corresponds to 'name'
        return dat[:, j]
    ht, dat = loadprivate(fnr)

    x1 = h('time')
    y1 = h(gas_name) * 1e6
    MFC11 = h('MFC1_flow')*1000.0
    MFC21 = h('MFC2_flow')

    idx = (x1 > ta) & (x1 < tc)      # write in 2 step cause error
    x = x1[idx]              #time
    y = y1[idx]              #concentration
    MFC1 = MFC11[idx]
    MFC2 = MFC21[idx]

    set_cal = h(cal_name)[0]         ## all value same
    print('cal value in lib')
    print(set_cal*1e6)

    idx2 = (x > ta) & (x < ta+1500)  # 25min baseline
    baseline1 = y[idx2]
    zero1 = np.mean(baseline1)
    std1 = np.std(baseline1)
    print('zero1, std1:', zero1, std1)
    # print(zero1)
    # print(std1)

    idxtb = np.where(abs(x-tb)<10)
    tbb = idxtb[0][-1] - 3*12  ## index 3 min before recorded sample adding time
    i = tbb
    tcc = 0
    i += 500       # 40min after add sample, about 12 points/min
    ## Find end time: zero 2 < zero1 + pct*sigma
    while i < len(x):
        zero2 = np.mean(y[i-250: i])
        if zero2 < zero1 + std1*pct:  ## baseline end here
            tcc = i                   ## experiment end point
            break
        i += 200

    if tcc:
        print('calculation end time')
        print(time.ctime(x[tcc]))
    else:
        print('baseline still above mean+sigma')

    y0 = y[tbb:i]-zero1        ##truncated y, zero-ed
    xtruc = x[tbb:i]           ## truncated x in s
    x0 = (xtruc-xtruc[0])/60   ## in min,start from 0
    x3 = (x-x[0])/60
    # plt.plot(x0, y0)
    # plt.plot(x3, y)    ## raw plot
    # plt.show()
    # exit()
    MFC10=MFC1[tbb:i]
    MFC20=MFC2[tbb:i]

    ## integral under curve
    s = 0
    S = [0]
    for j in range(len(y0)-1):
        flow_tot = MFC10[j] + MFC20[j]  # sccm
        mflow_tot = flow_tot / 24455  # 22400 #mol/min #24466
        dt = x0[j + 1] - x0[j]
        dY = 0.5 * (y0[j + 1] + y0[j])
        moles = dY * mflow_tot * dt
        s += moles    # micromoles
        S.append(s)
    # print(S[-1])

    vol_in = float(weight) / float(MW)
    vol_ratio = vol_in / (S[-1]) * 1E6
    print('ratio')
    print(vol_ratio)

    cal = vol_ratio * set_cal * 1e6
    print('calibration factor')
    print(cal)     #### this is!
    fncal = os.path.join(fnr, 'par', 'calibration_factor.txt')
    with open(fncal, 'w') as f:
        f.write(str(cal))
    # exit()

    ########### plots #########
    ## flow control, use all data
    F, A1 = plt.subplots(dpi=150)   # figsize=(6.5, 4.0)
    A1.plot(x3, MFC1, label='Dilution', color='#03A609')  # dkgreen
    A1.set_xlabel('Time (minutes)', fontsize=12)
    A1.set_ylabel('Dilution, ZA (sccm)', fontsize=14, color='#03A609')
    A1.set_facecolor('#Ededee')
    A1.grid(c='white')
    A3 = A1.twinx()
    A3.plot(x3, MFC2, label='Beaker', color='#000088')   #dkblue
    A3.set_ylabel('Bubbler (sccm)', fontsize=12, color='#000088')
    # A1.set_xlim(T[cal_index[0]], T[cal_index[-1]])
    A1.set_title(date + ' Flow Rate')
    plt.show(block=False)
    F1 = F
    if savefig:
        plt.savefig(os.path.join(fnr, date +' Flowrates.png'), bbox_inches='tight')
    # exit()

    #Integrated droplet
    F, [A1, A2] = plt.subplots(2, 1, dpi=150, figsize=(10, 8))   #figsize=(6.5, 4.0)
    A1.set_facecolor('#Ededee')
    A2.set_facecolor('#Ededee')
    A1.grid(c='white')
    A2.grid(c='white')
    A1.plot(x0, y0, c='#B5292f', lw=1, label='Time series')
    A2.plot(x0, S, c='#4188BA', lw=1, label=('Integrated, total: %.2f' % S[-1]))
    A1.set_title('Droplet Calibration: cal = %.3f\n%s' %(cal, gas), fontsize=16)
    A1.legend(loc='best', fontsize=10)
    A2.legend(loc='best', fontsize=10)
    A1.set_xticklabels([])
    A2.set_xlabel('Time (minutes)', fontsize=12)
    A1.set_ylabel('Sample (ppm)', fontsize=12)
    A2.set_ylabel('Sample (µ moles)', fontsize=12)
    # F.tight_layout()
    plt.show(block=False)
    F2 = F
    if savefig:
        plt.savefig(os.path.join(fnr, date +' Integration.png'), bbox_inches='tight')
    # exit()

    #Raw data
    F, A = plt.subplots(dpi=150)   #figsize=(6.5, 4.0)
    A.set_facecolor('#Ededee')
    A.grid(c='white')
    F.tight_layout()
    x_t = [x[0]]     # x data that will be marked
    xmak = [time.strftime('%H:%M', time.localtime(x[0]))]

    n= len(x)
    for i in range(1, n):
        clock0 = time.strftime('%M:%S', time.localtime(x[i-1]))
        clock =  time.strftime('%M:%S', time.localtime(x[i]))
        # if (clock0[:2] == '29' and clock[:2]=='30') or (clock0[:2] == '59' and clock[:2]=='00'):
        if (clock0[:2] == '59' and clock[:2]=='00'):
            x_t.append(x[i])
            xmak.append(time.strftime('%H:%M', time.localtime(x[i])))

    A.plot(x, y, label=' %s (ppm)' % gas)
    A.set_xticks(x_t, xmak)
    A.set_ylabel('Conc. (ppm)', fontsize=12)
    A.set_xlabel('Clock time', fontsize=12)
    A.set_title('%s Raw Data\n%s'%(date, gas))
    F.autofmt_xdate()    # x-label will not overlap
    plt.show(block=False)
    F3 = F
    if savefig:
        plt.savefig(os.path.join(fnr, date +' Raw.png'), bbox_inches='tight')
    # exit()

    ### check if compound has the correct spectrum and RMSE
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
    A1.grid(c='white', which='major',lw=1)
    A1.grid(c='white', which='minor', lw=0.5)
    A2.grid(c='white', which='major',lw=1)
    A2.grid(c='white', which='minor', lw=0.5)
    A1.minorticks_on()
    A2.minorticks_on()

    A1.plot(nu, k, label='data', color='#DC61D0', lw=1) #, alpha=0.75)
    A1.plot(nu, model, label='model', color='#00ffff', lw=1) #, alpha=0.4)
    A1.legend(loc=0, fontsize=10)
    A1.set_ylabel('Absorbance', fontsize=12)
    A1.set_title('%s Combo Log\n%s'%(date, gas))

    A2.plot(nu, residuals, label='residuals', color='#7068bb', lw=1)  #, alpha=0.75
    A2.set_ylabel('Residuals', color='#7068bb', fontsize=12)
    A4 = A2.twinx()
    A4.plot(nu, partial_fit, label='partial fit', color='#A6ce6a', lw=1)  #, alpha=0.75)
    A2.set_xlabel('nu (cm-1)', fontsize=12)
    A4.set_ylabel('Partial fit', color='#A6ce6a', fontsize=12)
    A2.legend(loc=0, fontsize=10)
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

    pct = 4  # zero 2 < zero1 + pct*sigma, end experiment
    # sample = '7219 - Indene'
    sample = '9233 - Norbornane'
    date = '20230908test'
    row1 = 1200

    sample = '164514 - Methoxyperfluorobutane'
    date = '20231020d3'

    fnr = os.path.join(basepath, sample, date)
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

        f = open(os.path.join(fnrp, 'weight.txt'), 'r')
        temp = f.read().splitlines()
        weight = temp[0]
        f = open(os.path.join(fnrp, 'molecular_weight.txt'), 'r')
        temp = f.read().splitlines()
        MW = temp[0]

        _, _, _, _, max_row = calibration_droplet(fnr, sample, cid, weight, MW, t1, t2, t3, showgraph=True)
        print(max_row)
        # calibration_droplet(fnr, sample, cid, weight, MW, t1, t2, t3, row1, showgraph=True)
        # calibration_droplet(fnr, sample, cid, weight, MW, t1, t2, t3, row1, showgraph=True, savefig=True)

