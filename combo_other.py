## plot individual combo log spectrum
import numpy as np
import os
import h5py     #need by Windows
import tables   #need by windows
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Arial"

from spectral_logger1 import SpectralLogReader as slog

#plot the combo data array dimension
def array_dimension(fnr, r1, r2):
    note = ''
    read = slog(os.path.join(fnr, 'ComboResults'), verbose=True)
    _, _, max_row = read.get_spectra_row('broadband', 100, pull_results=True)
    print('max_row: ', max_row)
    if r1 > max_row or r2 > max_row:
        note = 'Row, Range numbers must be\nless than max_row: %s' % max_row
        print(note)
    else:
        x = list(range(r1, r2))
        for j, row in enumerate(x):
            data = read.get_spectra_row('broadband', row, pull_results=True)
            nu = data[0]['nu']
            print(row, len(nu))
    return note


## use animation to automatic scan row numbers
def plot_combo(fnr, gas, combokey = "partial_fit", row=None):
    note = ''
    viewtime = 1000  ## ms, animation interval
    read = slog(os.path.join(fnr, 'ComboResults'), verbose=True)
    _, _, max_row = read.get_spectra_row('broadband', 100, pull_results=True)
    # print('max row: ' + str(max_row))

    if row is None:
        rowrange = np.arange(100, int(max_row/100)*100, 100)

        figure = plt.figure()
        ax = figure.add_subplot(111)

        def gen():
            for i in rowrange:
                yield i

        def animate(i):
            data = read.get_spectra_row('broadband', i, pull_results=True)
            nu = data[0]['nu']
            y = data[0][combokey]
            ax.clear()
            ax.plot(nu, y, linewidth=0.5)
            ax.set_ylabel(combokey)
            ax.set_xlabel('nu, cm-1')
            ax.set_title(gas + ', row: ' + str(i))

        anim = FuncAnimation(figure, animate, frames=gen(), repeat=False, interval=viewtime)
        plt.show()

    else:
        if row > max_row or row < 0:
            note = 'Row number must be between 0 and %s' % max_row
        else:
            
            data = read.get_spectra_row('broadband', row, pull_results=True)
            nu = data[0]['nu']
            y = data[0][combokey]
            # model = data[0]['model']
            # print(model)
            print(len(nu))
    
            fig, ax = plt.subplots()
            ax.plot(nu, y)
            ax.set_xlabel('nu, cm-1')
            ax.set_ylabel(combokey)
            ax.set_title(gas + ', row: '+str(row))
            plt.show()
    return note


if __name__ == "__main__":
    # basepath = r'/mnt/r/crd_G9000/AVXxx/3610-NUV1022/R&D/Calibration/'          ##Linux
    basepath = '/Volumes/Data/crd_G9000/AVXxx/3610-NUV1022/R&D/Calibration'       ##Mac
    # basepath = 'R:\crd_G9000\AVXxx\\3610-NUV1022\R&D\Calibration'               ## Windows
    COMBOKEYS = ['partial_fit', 'absorbance', 'model', 'residuals']
    combokey = COMBOKEYS[0]

    # gas = '9233 - Norbornane'
    # date = '20230908test'
    # row1 = 1200

    gas = '7868 - 2,4,4-trimethyl-1-pentene'
    date = '20230906d1'
    row1 = 1999

    fnr = os.path.join(basepath, gas, date)
    print(fnr)
    fnrp = os.path.join(fnr, 'par')

    if not os.path.exists(fnr):
        print('Error, did not find data. Please check if data exist or attach the data/R drive.')
    elif not os.path.exists(fnrp):
        print('Error, did not find parameters data.')
    else:
        # plot each row of spectrum
        # plot_combo(fnr, gas)            # find row number animation
        # plot_combo(fnr, gas, combokey)            # find row number
        plot_combo(fnr, gas, row=row1)  # row number

        # array list, mode number
        # r1 = 1100
        # r2 = 1150
        # array_dimension(fnr, r1, r2)



