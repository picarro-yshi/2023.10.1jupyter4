## calibration factor for droplet experiment, adapted from jupyter notebook
import os
import time
import platform
import numpy as np
import h5py  # need by Windows
import tables  # need by windows
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Arial"

from spectral_logger1 import SpectralLogReader as slog
from loadprivate import loadprivate

BASELINE_TIME = 40*60  # s, baseline time


## all parameters must be valid before sending to below functions; no sanity check
## need .h5 files in PrvateData/broadband and ComboResults folder
def calibration_droplet(fnr, gas, cid, weight, MW, t1, t2, t3, row=100, pct=4, showgraph=False, savefig=False):
    gas_name = 'broadband_gasConcs_' + str(cid)
    cal_name = 'broadband_eCompoundOutputs_'+ str(cid) +'_calibration'
    date = t1[:8]

    # t2 format: 202312110941
    epo2 = int(time.mktime(time.strptime(t2, "%Y%m%d%H%M")))
    # get start time of usable data
    epo1 = epo2 - 1800  # 30 min baseline before add sample
    epo3_raw = int(time.mktime(time.strptime(t3, "%Y%m%d%H%M")))

    def h(name):
        j = ht.index(name)  # finds the index of the list ht that corresponds to 'name'
        return dat[:, j]
    ht, dat = loadprivate(fnr)

    # raw data
    x1 = h('time')
    y1 = h(gas_name) * 1e6
    MFC11 = h('MFC1_flow')*1000.0
    MFC21 = h('MFC2_flow')

    # calculate baseline before level
    idx = (x1 > epo1) & (x1 < epo1 + 1500)  # 25min baseline
    baseline = y1[idx]
    zero1 = np.mean(baseline)
    std1 = np.std(baseline)
    print('zero1, std1:', zero1, std1)

    # get end time of usable data
    # start track baseline 40 min after add sample
    t = epo2 + 2400
    epo3 = 0

    # check every 10 min
    while t < epo3_raw:
        idx = (x1 > t) & (x1 < t + BASELINE_TIME)
        baseline = y1[idx]
        zero2 = np.mean(baseline)

        if zero2 < zero1 + std1 * pct:  ## baseline end here
            print('dropped')
            epo3 = t + BASELINE_TIME  ## experiment end point
            break
        t += 600

    if epo3:
        print('calculation end time: ')
        print(time.ctime(epo3))
    else:
        epo3 = epo3_raw
        print('baseline still above mean+sigma')

    # truncate to get usable data
    idx = (x1 > epo1) & (x1 < epo3)
    x = x1[idx]  # epoch time
    y = y1[idx] - zero1  # concentration
    MFC1 = MFC11[idx]
    MFC2 = MFC21[idx]

    # integral under curve
    x0 = (x-x[0])/60  # min, float, time elapsed
    s = 0
    S = [0]
    for j in range(len(y)-1):
        flow_tot = MFC1[j] + MFC2[j]  # sccm
        mflow_tot = flow_tot / 24455  # 22400 #mol/min #24466
        dt = x0[j + 1] - x0[j]
        dY = 0.5 * (y[j + 1] + y[j])
        moles = dY * mflow_tot * dt
        s += moles    # micro moles
        S.append(s)
    # print(S[-1])

    set_cal = h(cal_name)[0]  # all value same
    print('cal value in the library')
    print(set_cal*1e6)

    vol_in = float(weight) / float(MW)
    vol_ratio = vol_in / (S[-1]) * 1E6
    print('ratio')
    print(vol_ratio)

    cal = vol_ratio * set_cal * 1e6
    print('calibration factor')
    print(cal)  # this is!

    fncal = os.path.join(fnr, 'par', 'calibration_factor.txt')
    with open(fncal, 'w') as f:
        f.write(str(cal))
    # exit()

    ########### Generate x axis tick marks #########
    # every 30 mins and no more than 20 marks
    # main axis, bottom
    x_t = []  # epoch time
    xmak = []  # mark

    n= len(x)
    for i in range(1, n):
        clock0 = time.strftime('%M:%S', time.localtime(x[i-1]))
        clock =  time.strftime('%M:%S', time.localtime(x[i]))
        if (clock0[:2] == '29' and clock[:2]=='30') or (clock0[:2] == '59' and clock[:2]=='00'):
            x_t.append(int(x[i]))
            xmak.append(time.strftime('%H:%M', time.localtime(x[i])))

    if len(xmak) > 20:
        # print("kick out 30 min marks")
        if "30" in xmak[0]:
            x_t = x_t[1::2]
            xmak = xmak[1::2]

    # get odd index marks
    while len(xmak) > 20:
        x_t = x_t[::2]
        xmak = xmak[::2]

    # add start time
    x_t.insert(0, int(x[0]))
    xmak.insert(0, time.strftime('%H:%M', time.localtime(x[0])))

    # add end time
    tag = 1  # attach to end
    gap = x_t[-1] - x_t[-2]
    if len(xmak) > 10 and x[-1] - x_t[-1] < gap/2:
        tag = 0  # attach as the end

    if tag:
        x_t.append(x[-1])
        xmak.append(time.strftime('%H:%M', time.localtime(x[-1])))
    else:
        x_t[-1] = x[-1]
        xmak[-1] = time.strftime('%H:%M', time.localtime(x[-1]))

    # secondary axis, top
    x_t2 = []
    xmak2 = []
    counter = 1
    for i in range(1, n):
        t = x[i] - x[0]  # s
        if t > counter * 1800 > (x[i - 1] - x[0]):  # 30 min
            x_t2.append(x[i])
            xmak2.append(str(int(t/60)))
            counter += 1

    if len(xmak2) > 20:
        # print("kick out 30 min marks")
        if "30" in xmak2[0]:
            x_t2 = x_t2[1::2]
            xmak2 = xmak2[1::2]

    # get odd index marks
    while len(xmak2) > 20:
        x_t2 = x_t2[::2]
        xmak2 = xmak2[::2]

    x_t2.insert(0, int(x[0]))
    xmak2.insert(0, "0")

    tag = 1  # attach to end
    gap = x_t2[-1] - x_t2[-2]
    if len(xmak2) > 10 and x[-1] - x_t2[-1] < gap/2:
        print('tag =0 ')
        tag = 0  # attach as the end

    if tag:
        x_t2.append(x[-1])
        xmak2.append(str(int((x[-1] - x[0])/60)))
    else:
        x_t2[-1] = x[-1]
        xmak2[-1] = str(int((x[-1] - x[0])/60))

    ########### plots #########
    #1 flow control, use all data
    F, A1 = plt.subplots(dpi=150)   # figsize=(6.5, 4.0)
    A1.set_facecolor('#Ededee')
    A1.grid(c='white')

    A1.plot(x, MFC1, label='Dilution', color='#03A609')  # dkgreen
    A1.set_xlabel('Clock time', fontsize=12)
    A1.set_ylabel('Dilution, ZA (sccm)', fontsize=14, color='#03A609')
    A1.set_title(date + ' Flow Rate')

    A1.set_xticks(x_t, xmak, fontsize=10, rotation=40, ha='right')
    A1.tick_params(pad=0)

    # top axis for min
    ax2 = A1.twiny()
    ax2.set_xticks(x_t2, xmak2, fontsize=8)
    ax2.tick_params(length=2, grid_alpha=0.5, pad=-1)
    ax2.plot(x, MFC1, linewidth= 0)
    ax2.set_xlabel("min:", fontsize=8)
    ax2.xaxis.set_label_coords(-0.01, 1.01)

    A3 = A1.twinx()
    A3.plot(x, MFC2, label='Bubbler', color='#000088')   #dkblue
    A3.set_ylabel('Bubbler (sccm)', fontsize=12, color='#000088')

    A3.set_xticks(x_t, xmak, fontsize=10, rotation=40, ha='right')
    A3.tick_params(pad=0)

    plt.tight_layout()
    plt.show(block=False)
    F1 = F
    if savefig:
        plt.savefig(os.path.join(fnr, date +' Flowrates.png'), bbox_inches='tight')
    # plt.show()
    # exit()

    #2 Integrated droplet
    F, [A1, A2] = plt.subplots(2, 1, dpi=150, figsize=(10, 8))   #figsize=(6.5, 4.0)
    A1.set_facecolor('#Ededee')
    A2.set_facecolor('#Ededee')
    A1.grid(c='white')
    A2.grid(c='white')

    A1.plot(x, y, c='#B5292f', lw=1, label='Time series')
    A1.set_title('Droplet Calibration: cal = %.3f\n%s' %(cal, gas), fontsize=16)
    A1.set_xticklabels([])
    A1.legend(loc='best', fontsize=10)
    A1.set_ylabel('Sample (ppm)', fontsize=12)

    # top axis for min
    ax2 = A1.twiny()
    ax2.set_xticks(x_t2, xmak2, fontsize=8)
    ax2.tick_params(length=2, grid_alpha=0.5,pad=-1)
    ax2.plot(x, y, linewidth= 0)
    ax2.set_xlabel("min:", fontsize=8)
    ax2.xaxis.set_label_coords(-0.01, 1.01)

    A2.plot(x, S, c='#4188BA', lw=1, label=('Integrated, total: %.2f' % S[-1]))
    A2.legend(loc='best', fontsize=10)
    A2.set_xlabel('Clock time', fontsize=12)
    A2.set_ylabel('Sample (µ moles)', fontsize=12)

    A2.set_xticks(x_t, xmak, fontsize=10, rotation=40, ha='right')
    A2.tick_params(pad=0)

    F.tight_layout()
    plt.show(block=False)
    F2 = F
    if savefig:
        plt.savefig(os.path.join(fnr, date +' Integration.png'), bbox_inches='tight')
    # plt.show()
    # exit()

    #3 Raw data
    F, A = plt.subplots(dpi=150)   #figsize=(6.5, 4.0)
    A.set_facecolor('#Ededee')
    A.grid(c='white')

    A.plot(x, y, label=' %s (ppm)' % gas)
    A.set_ylabel('Conc. (ppm)', fontsize=12)
    A.set_xlabel('Clock time', fontsize=12)
    A.set_title('%s Raw Data\n%s'%(date, gas))
    A.set_xticks(x_t, xmak, fontsize=10, rotation=40, ha='right')
    A.tick_params(pad=0)

    # Shrink current axis's height by 20% on the bottom
    # box = A.get_position()
    # A.set_position(
    #     [box.x0, box.y0 + box.height * 0.05, box.width, box.height * 0.9]
    # )

    ax2 = A.twiny()
    ax2.set_xticks(x_t2, xmak2, fontsize=8)
    # ax2.tick_params(direction='out', length=6, width=2, colors='r',
    #                grid_color='r', grid_alpha=0.5)
    ax2.tick_params(length=2, grid_alpha=0.5, pad=-1)
    ax2.plot(x, y)
    ax2.set_xlabel("min:", fontsize=8)
    ax2.xaxis.set_label_coords(-0.01, 1.01)

    F.tight_layout()
    plt.show(block=False)
    F3 = F
    if savefig:
        plt.savefig(os.path.join(fnr, date +' Raw.png'))
    # plt.show()
    # exit()

    #4 Combo plot: check if compound has the correct spectrum and RMSE
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

    def close_figure(event):
        plt.close('all')

    if showgraph:
        # enable if you want to press one key and close all figures at a time
        # plt.gcf().canvas.mpl_connect('key_press_event', close_figure)

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

    # pct = 4  # zero 2 < zero1 + pct*sigma, end experiment
    # sample = '7219 - Indene'
    sample = '9233 - Norbornane'
    date = '20230908test'
    row1 = 1200

    # sample = '8864 - 1,3-Diethylbenzene'
    # date = '20240105'
    # date = '20231229'

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

        # F1, F2, F3, F4, max_row = calibration_droplet(fnr, sample, cid, weight, MW, t1, t2, t3, showgraph=True)
        # print(max_row)

        calibration_droplet(fnr, sample, cid, weight, MW, t1, t2, t3, row1, showgraph=True)
        # calibration_droplet(fnr, sample, cid, weight, MW, t1, t2, t3, row1, showgraph=True, savefig=True)

