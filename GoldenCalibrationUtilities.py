# -*- coding: utf-8 -*-
"""
Created on Mon Mar 06 11:40:56 2017
​
@author: chris
"""
#  20170629 - added 'units' to calcStuffWeighted
#  20180209 - 

import numpy as np
import scipy.stats

def calcStuffWeighted(X, Y, w0, w1, txt = '', units = 'ppb', k=3): #directly from semi standard
    n = len(X)
    w = 1.0 / (w0**2 + (w1*X)**2)
    W = w / np.sum(w)

    txt += '*** Weighted LS Fit ***' + '\n'
    txt += 'measurements = %d' % n  + '\n'

    X_bar = np.sum(X*W)
    Y_bar = np.sum(Y*W)
    SS_xx = np.sum(W*(X - X_bar)**2)
    SS_xy = np.sum(W*(X - X_bar)*(Y - Y_bar))
    m = SS_xy / SS_xx
    txt += 'slope = %.4f' % m  + '\n'

    b = Y_bar - m*X_bar
    txt += 'intercept = %.4f [%s]' % (b, units)  + '\n'
    Y_fit = m*X + b
    ss = np.sqrt(np.sum((Y - Y_fit)**2)/(n - 2))
    SS_yy = np.sum(W*(Y - Y_bar)**2)
    ss2 = np.sqrt((SS_yy - m*SS_xy)/ (n - 2))
    txt += 'residual std dev = %.4f' % ss + '\n'
    #txt += 'std dev of regression model = %.4f' % ss2 + '\n'

    G = scipy.stats.norm
    ef = 1.0 - G.cdf(k) #one-sided probability for 3 sigma

    T = scipy.stats.t
    tvalue = T.ppf(1.0 - ef, n - 2)
    txt += 'students 3-sigma t-statistic is %.3f' % tvalue + '\n'
    UCL = b + tvalue*np.sqrt(1.0 / max(W) + 1.0/n + X_bar**2 / SS_xx)*ss2
    UCL_x = tvalue*np.sqrt(1.0/W + 1.0/n + (X - X_bar)**2 / SS_xx)*ss2
    myMDL = (UCL - b)/m
    txt += 'MDL = %.3f [%s]' % (myMDL, units) + '\n'
    residuals = Y - Y_fit
    return m, b, UCL_x, residuals, txt, myMDL

def calcStuffWeighted_points(X, Y, E, E0, txt = '', units = 'ppb', k=3): 
    #modified from semi standard
    # E: point-by-point uncertainty
    # E0: uncertainty at zero (required for MDL calc) 
    # -- get it from your uncertainty propagation (I used the average of E for the zero values)
    n = len(X)
    W = 1/E**2
    W0 = 1/E0**2
    W0 = W0/sum(W)
    W = W/sum(W)
    txt += '*** Weighted LS Fit ***' + '\n'
    txt += 'measurements = %d' % n  + '\n'

    X_bar = np.sum(X*W)
    Y_bar = np.sum(Y*W)

    SS_xx = np.sum(W*(X - X_bar)**2)
    SS_xy = np.sum(W*(X - X_bar)*(Y - Y_bar))
    m = SS_xy / SS_xx
    txt += 'slope = %.4f' % m  + '\n'

    b = Y_bar - m*X_bar
    txt += 'intercept = %.4f [%s]' % (b, units)  + '\n'
    Y_fit = m*X + b
    ss = np.sqrt(np.sum((Y - Y_fit)**2)/(n - 2))
    SS_yy = np.sum(W*(Y - Y_bar)**2)
    ss2 = np.sqrt((SS_yy - m*SS_xy)/ (n - 2))
    txt += 'residual std dev = %.4f' % ss + '\n'
    #txt += 'std dev of regression model = %.4f' % ss2 + '\n'

    G = scipy.stats.norm
    ef = 1.0 - G.cdf(k) #one-sided probability for 3 sigma

    T = scipy.stats.t
    tvalue = T.ppf(1.0 - ef, n - 2)
    txt += 'students %d-sigma t-statistic is %.3f' % (k, tvalue) + '\n'
    UCL = b + tvalue*np.sqrt(1.0 / W0 + 1.0/n + X_bar**2 / SS_xx)*ss2
    sigma_data = tvalue*np.sqrt(1.0/W + 1.0/n + (X - X_bar)**2 / SS_xx)*ss2
    myMDL = (UCL - b)/m
    txt += 'MDL = %.3f [%s]' % (myMDL, units) + '\n'
    residuals = Y - Y_fit
    return m, b, sigma_data, residuals, txt, myMDL

def rawDataSaver(cols, colnames, suffix, fn = 'raw_data_outputs_%s.csv'):
    row_list = [','.join(colnames)]
    for i in range(len(cols[0])):
        row = [c[i] for c in cols]
        rowtxt = ','.join(['%.3f' % r for r in row])
        row_list.append(rowtxt)
    alltxt = '\n'.join(row_list)
    fout = open(fn % suffix, 'w')
    fout.write(alltxt)
    fout.close()

def Grouper(ArrayList, width_max, N_max, mode='logmode', param=0):
    newA = []
    for A in ArrayList:
        X = A[:,0]
        if mode == 'fixed':
            width = max(X) - min(X)
            final_width_max = width_max
        elif mode == 'logmode':
            width = max(np.log(X)) - min(np.log(X))
            final_width_max = width_max
        elif mode == 'mixedMode':
            #added by Rella 20180210 to handle negative x values -- based on the non-negative filter
            width = max(X) - min(X)
            final_width_max = 0.5*np.sqrt((2.0*param)**2 + (width_max*X.mean())**2) + 0.5*width_max*X.mean() 
        #print "width = %.4f; final_width_max=%.4f; %.3f, %.3f" % (width, final_width_max, min(X), max(X))
        N = len(X)
        if width <= final_width_max and N <= N_max:
            newA.append(A)
        else:
            maxdiff = -1.0
            icut = -1
            for i in range(N-1):
                if X[i+1] - X[i] > maxdiff:
                    maxdiff = X[i+1] - X[i] 
                    icut = i + 1
            A1 = A[:icut,:]
            A2 = A[icut:,:]
            newA.extend(Grouper([A1, A2], width_max, N_max, mode=mode, param=param))
    return newA

def findIndex(ht, labels):
    if type(labels) == list:
        print ('list found')
        for label in labels:
            try:
                i = ht.index(label)
                label_used = label
                break #break out of loop
            except ValueError:
                print ('%s column not found' % label)
                i = -1
        return i, label_used
    else:
        label_used = labels
        try:
            i = ht.index(labels)
        except ValueError:
            print ('%s column not found' % label)
        return i

def CleanUp(T, C, P = None):
    good = []
    D = []
    for i in range(len(C) - 2):
        D.append(2*C[i+1] - C[i+2] - C[i])
    D = np.array(D)
    medD = np.median(abs(D))
    bad = []
    for i in range(len(D)):
        if abs(D[i]) < 5.0*medD:
            good.append(i + 1)
        else:
            bad.append(i+ 1)

    print ((len(D) - len(good))/float(len(D)), ' fraction removed!', medD, np.mean(D), np.max(abs(D)))
    if P == None:
        return T[good], C[good]
    else:
        return T[good], C[good], P[good]

def exponentialAverage(T, C, tau):
    new = []
    C_ave = 0.0
    for j in range(len(T)):
        if j != len(T) - 1:
            dt = T[j + 1] - T[j]
        else:
            dt = T[j] - T[j-1]
        alph = np.exp(-dt / tau)
        if j == 0:
            C_ave = C[j]
        else:
            C_ave = C_ave*alph + C[j]*(1.0 - alph)
        new.append(C_ave)
    return np.array(new)