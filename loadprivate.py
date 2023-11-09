# read private log

import numpy as np
import os
import h5py  # need by Windows
import tables  # need by windows


def loadprivate(fnr, verbose=True):  # from Chris's datautility3
    ht = []  # head tags
    dat = None
    fd = os.path.join(fnr, 'PrivateData', 'broadband')
    N = len(os.listdir(fd))
    ct = 1
    for f in os.listdir(fd):
        ext = f[-3:]
        fn = os.path.join(fd, f)
        if ext == 'dat':
            fob = open(fn, 'r')
            header = fob.readline(-1)
            ht = string.split(string.strip(header, ' '))
            mydata = []
            stdlen = -1
            for lines in fob:
                data = string.split(string.strip(lines, ' '))
                row = []
                try:
                    for k in range(len(ht)):
                        try:
                            row.append(float(data[k]))
                        except:
                            row.append(-1.0)
                    if stdlen == -1:
                        stdlen = len(data)
                    if len(data) == stdlen:
                        mydata.append(row)
                except:
                    if verbose: print('!', end=' ')
            fob.close()

        elif ext == '.h5':
            fob5 = tables.open_file(fn)
            D = [g for g in fob5.root]
            if ht == []:
                ht = D[0].colnames

            for tag in ht:
                if tag == ht[0]:
                    datarray = np.array([D[0].col(tag)])
                else:
                    datarray = np.vstack((datarray, D[0].col(tag)))

            mydata = np.transpose(datarray)
            fob5.close()

        try:
            if dat is None:
                dat = mydata
            else:
                dat = np.vstack((dat, mydata))
        except:
            print('dat array dimension changed in file %s. Try cut the first ten minutes and re-run.' % f)

        if verbose:
            print("Loading %d / %d... " % (ct, N))
        ct += 1

    return ht, dat


if __name__ == "__main__":
    rdrive = '/Volumes/Data/crd_G9000/AVXxx/3610-NUV1022/R&D/Calibration'
    compound = '638186 - trans-1,2-Dichloroethylene'
    compound = '6574 - 1,1,2-Trichloroethane'
    date = '20230815'
    date = '20231103'
    fnr = os.path.join(rdrive, compound, date)

    ht, dat = loadprivate(fnr)
    print(ht)  # head tag
    # print(dat)  # head tag

# ## import customized files from other folder
# helperpath = '../code/Rella-Python/'    ## './' in same folder
# sys.path.append(helperpath)
# import DataUtilities3.load_files_rella as LF
