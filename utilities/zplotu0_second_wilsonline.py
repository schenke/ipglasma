#!/usr/bin/env python3
import numpy as np
#import collections
#from enum import Enum
import math
import pylab as pl
#from collections import Counter
#import pandas as pd
from numpy import *
from scipy import optimize, special
from scipy import *
from scipy.optimize import curve_fit
#from scipy import asarray as ar,exp 
from sys import argv, exit 
istep = 15
istep = int(istep)
ievent = 15 
NADD = 0
ievent = int(ievent)
aa = np.loadtxt("V_second_1.txt".format(ievent))
length =512+NADD

x = []
y = []
dx = 30/512
T00 = []
temp = aa[:,2]
for ii in range(length):
    mid = temp[int(ii*length):int(ii*length+length)]
    T00.append(mid)
    x.append(-1*length/2.*dx + ii*dx)
    #if ( ii <length):
    #    y.append(-1*length/2.*dx + ii*dx)
T00 = np.array(T00)
x=np.array(x)
#y=np.array(y)
y=x
X, Y = np.meshgrid(x, y)
pl.figure(10020)
levels = np.linspace(-1,1,50)
zxx = [-25,25]
zyy= [0,0]
for kk in range(1):
    print(np.sum(T00))
    cs = pl.contourf(X, Y, T00, levels=levels,colorscale='Electric', fill=False,cmap=pl.cm.jet)#cmap="Blues")
    pl.xlim(-17.5,17.5)
    pl.ylim(-17.5,17.5)
    pl.plot(zxx,zyy,'k--')
    pl.plot(zyy,zxx,'k--')
    pl.xlabel('fm', fontsize=15)# make axis labels
    pl.ylabel('fm', fontsize=15)
    pl.xticks(fontsize=15)
    pl.yticks(fontsize=15)
    pl.colorbar()
    pl.title('Pb-Pb, V1,second shift,Q(0,0)_real, b={}fm'.format(istep), fontsize=13)# give plot a title
    pl.savefig("Png/V1_Q_00_real{}.png".format(istep))
    
T00 = []
temp = aa[:,3]
for ii in range(length):
    mid = temp[int(ii*length):int(ii*length+length)]
    T00.append(mid)
T00 = np.array(T00)
pl.figure(1002013)
#levels = np.linspace(-2.0,2.0,21)
for kk in range(1):
    print(np.sum(T00))
    cs = pl.contourf(X, Y, T00, levels=levels,colorscale='Electric', fill=False,cmap=pl.cm.jet)#cmap="Blues")
    pl.xlim(-17.5,17.5)
    pl.ylim(-17.5,17.5)
    pl.plot(zxx,zyy,'k--')
    pl.plot(zyy,zxx,'k--')
    pl.xlabel('fm', fontsize=15)# make axis labels
    pl.ylabel('fm', fontsize=15)
    pl.xticks(fontsize=15)
    pl.yticks(fontsize=15)
    pl.colorbar(cs)
    pl.title('Pb-Pb, V1, second shift,Q(0,0)_imag, b={}fm'.format(istep), fontsize=13)# give plot a title
    pl.savefig("Png/V1_Q_00_imag{}.png".format(istep))

T00 = []
temp = aa[:,4]
for ii in range(length):
    mid = temp[int(ii*length):int(ii*length+length)]
    T00.append(mid)
T00 = np.array(T00)
pl.figure(1002014)
#levels = np.linspace(-2.0,2.0,21)
for kk in range(1):
    print(np.sum(T00))
    cs = pl.contourf(X, Y, T00, levels=levels,colorscale='Electric', fill=False,cmap=pl.cm.jet)#cmap="Blues")
    pl.xlim(-17.5,17.5)
    pl.ylim(-17.5,17.5)
    pl.plot(zxx,zyy,'k--')
    pl.plot(zyy,zxx,'k--')
    pl.xlabel('fm', fontsize=15)# make axis labels
    pl.ylabel('fm', fontsize=15)
    pl.xticks(fontsize=15)
    pl.yticks(fontsize=15)
    pl.colorbar(cs)
    pl.title('Pb-Pb, V1, second shift,Q(0,0)_imag, b={}fm'.format(istep), fontsize=13)# give plot a title
    pl.savefig("Png/V1_Q_4_imag{}.png".format(istep))

T00 = []
temp = aa[:,5]
for ii in range(length):
    mid = temp[int(ii*length):int(ii*length+length)]
    T00.append(mid)
T00 = np.array(T00)
pl.figure(1002015)
#levels = np.linspace(-2.0,2.0,21)
for kk in range(1):
    print(np.sum(T00))
    cs = pl.contourf(X, Y, T00, levels=levels,colorscale='Electric', fill=False,cmap=pl.cm.jet)#cmap="Blues")
    pl.xlim(-17.5,17.5)
    pl.ylim(-17.5,17.5)
    pl.plot(zxx,zyy,'k--')
    pl.plot(zyy,zxx,'k--')
    pl.xlabel('fm', fontsize=15)# make axis labels
    pl.ylabel('fm', fontsize=15)
    pl.xticks(fontsize=15)
    pl.yticks(fontsize=15)
    pl.colorbar(cs)
    pl.title('Pb-Pb, V1, second shift,Q(0,0)_imag, b={}fm'.format(istep), fontsize=13)# give plot a title
    pl.savefig("Png/V1_Q_5_imag{}.png".format(istep))

T00 = []
temp = aa[:,6]
for ii in range(length):
    mid = temp[int(ii*length):int(ii*length+length)]
    T00.append(mid)
T00 = np.array(T00)
pl.figure(1002016)
#levels = np.linspace(-2.0,2.0,21)
for kk in range(1):
    print(np.sum(T00))
    cs = pl.contourf(X, Y, T00, levels=levels,colorscale='Electric', fill=False,cmap=pl.cm.jet)#cmap="Blues")
    pl.xlim(-17.5,17.5)
    pl.ylim(-17.5,17.5)
    pl.plot(zxx,zyy,'k--')
    pl.plot(zyy,zxx,'k--')
    pl.xlabel('fm', fontsize=15)# make axis labels
    pl.ylabel('fm', fontsize=15)
    pl.xticks(fontsize=15)
    pl.yticks(fontsize=15)
    pl.colorbar(cs)
    pl.title('Pb-Pb, V1, second shift,Q(0,0)_imag, b={}fm'.format(istep), fontsize=13)# give plot a title
    pl.savefig("Png/V1_Q_6_imag{}.png".format(istep))

T00 = []
temp = aa[:,7]
for ii in range(length):
    mid = temp[int(ii*length):int(ii*length+length)]
    T00.append(mid)
T00 = np.array(T00)
pl.figure(1002017)
#levels = np.linspace(-2.0,2.0,21)
for kk in range(1):
    print(np.sum(T00))
    cs = pl.contourf(X, Y, T00, levels=levels,colorscale='Electric', fill=False,cmap=pl.cm.jet)#cmap="Blues")
    pl.xlim(-17.5,17.5)
    pl.ylim(-17.5,17.5)
    pl.plot(zxx,zyy,'k--')
    pl.plot(zyy,zxx,'k--')
    pl.xlabel('fm', fontsize=15)# make axis labels
    pl.ylabel('fm', fontsize=15)
    pl.xticks(fontsize=15)
    pl.yticks(fontsize=15)
    pl.colorbar(cs)
    pl.title('Pb-Pb, V1, second shift,Q(0,0)_imag, b={}fm'.format(istep), fontsize=13)# give plot a title
    pl.savefig("Png/V1_Q_7_imag{}.png".format(istep))

T00 = []
temp = aa[:,8]
for ii in range(length):
    mid = temp[int(ii*length):int(ii*length+length)]
    T00.append(mid)
T00 = np.array(T00)
pl.figure(1002018)
#levels = np.linspace(-2.0,2.0,21)
for kk in range(1):
    print(np.sum(T00))
    cs = pl.contourf(X, Y, T00, levels=levels,colorscale='Electric', fill=False,cmap=pl.cm.jet)#cmap="Blues")
    pl.xlim(-17.5,17.5)
    pl.ylim(-17.5,17.5)
    pl.plot(zxx,zyy,'k--')
    pl.plot(zyy,zxx,'k--')
    pl.xlabel('fm', fontsize=15)# make axis labels
    pl.ylabel('fm', fontsize=15)
    pl.xticks(fontsize=15)
    pl.yticks(fontsize=15)
    pl.colorbar(cs)
    pl.title('Pb-Pb, V1, second shift,Q(0,0)_imag, b={}fm'.format(istep), fontsize=13)# give plot a title
    pl.savefig("Png/V1_Q_8_imag{}.png".format(istep))

T00 = []
temp = aa[:,9]
for ii in range(length):
    mid = temp[int(ii*length):int(ii*length+length)]
    T00.append(mid)
T00 = np.array(T00)
pl.figure(1002019)
#levels = np.linspace(-2.0,2.0,21)
for kk in range(1):
    print(np.sum(T00))
    cs = pl.contourf(X, Y, T00, levels=levels,colorscale='Electric', fill=False,cmap=pl.cm.jet)#cmap="Blues")
    pl.xlim(-17.5,17.5)
    pl.ylim(-17.5,17.5)
    pl.plot(zxx,zyy,'k--')
    pl.plot(zyy,zxx,'k--')
    pl.xlabel('fm', fontsize=15)# make axis labels
    pl.ylabel('fm', fontsize=15)
    pl.xticks(fontsize=15)
    pl.yticks(fontsize=15)
    pl.colorbar(cs)
    pl.title('Pb-Pb, V1, second shift,Q(0,0)_imag, b={}fm'.format(istep), fontsize=13)# give plot a title
    pl.savefig("Png/V1_Q_9_imag{}.png".format(istep))



aa = np.loadtxt("V_second_2.txt".format(ievent))
#length =720

x = []
#dx = 30/length
T00 = []
temp = aa[:,2]
for ii in range(length):
    mid = temp[int(ii*length):int(ii*length+length)]
    T00.append(mid)
    x.append(-1*length/2.*dx + ii*dx)
T00 = np.array(T00)
np.savetxt("T00",T00)
x=np.array(x)
y=x
X, Y = np.meshgrid(x, y)
pl.figure(310020)
#levels = np.linspace(0,200,50)


for kk in range(1):
    print(np.sum(T00))
    cs = pl.contourf(X, Y, T00, levels=levels,colorscale='Electric', fill=False,cmap=pl.cm.jet)#cmap="Blues")
    pl.xlim(-17.5,17.5)
    pl.ylim(-17.5,17.5)
    pl.plot(zxx,zyy,'k--')
    pl.plot(zyy,zxx,'k--')
    pl.xlabel('fm', fontsize=15)# make axis labels
    pl.ylabel('fm', fontsize=15)
    pl.xticks(fontsize=15)
    pl.yticks(fontsize=15)
    pl.colorbar()
    pl.title('Pb-Pb, V2, second shift,Q(0,0)_real, b={}fm'.format(istep), fontsize=13)# give plot a title
    pl.savefig("Png/V2_Q_00_real{}.png".format(istep))

T00 = []
temp = aa[:,3]
for ii in range(length):
    mid = temp[int(ii*length):int(ii*length+length)]
    T00.append(mid)
T00 = np.array(T00)
pl.figure(10022301)
#levels = np.linspace(-2.0,2.0,21)
for kk in range(1):
    print(np.sum(T00))
    cs = pl.contourf(X, Y, T00, levels=levels,colorscale='Electric', fill=False,cmap=pl.cm.jet)#cmap="Blues")
    pl.xlim(-17.5,17.5)
    pl.ylim(-17.5,17.5)
    pl.plot(zxx,zyy,'k--')
    pl.plot(zyy,zxx,'k--')
    pl.xlabel('fm', fontsize=15)# make axis labels
    pl.ylabel('fm', fontsize=15)
    pl.xticks(fontsize=15)
    pl.yticks(fontsize=15)
    pl.colorbar(cs)
    pl.title('Pb-Pb, V2, second shift,Q(0,0)_imag, b={}fm'.format(istep), fontsize=13)# give plot a title
    pl.savefig("Png/V2_Q_00_imag{}.png".format(istep))


T00 = []
temp = aa[:,4]
for ii in range(length):
    mid = temp[int(ii*length):int(ii*length+length)]
    T00.append(mid)
T00 = np.array(T00)
pl.figure(10022304)
#levels = np.linspace(-2.0,2.0,21)
for kk in range(1):
    print(np.sum(T00))
    cs = pl.contourf(X, Y, T00, levels=levels,colorscale='Electric', fill=False,cmap=pl.cm.jet)#cmap="Blues")
    pl.xlim(-17.5,17.5)
    pl.ylim(-17.5,17.5)
    pl.plot(zxx,zyy,'k--')
    pl.plot(zyy,zxx,'k--')
    pl.xlabel('fm', fontsize=15)# make axis labels
    pl.ylabel('fm', fontsize=15)
    pl.xticks(fontsize=15)
    pl.yticks(fontsize=15)
    pl.colorbar(cs)
    pl.title('Pb-Pb, V2, second shift,Q(0,0)_imag, b={}fm'.format(istep), fontsize=13)# give plot a title
    pl.savefig("Png/V2_Q_04_imag{}.png".format(istep))

T00 = []
temp = aa[:,5]
for ii in range(length):
    mid = temp[int(ii*length):int(ii*length+length)]
    T00.append(mid)
T00 = np.array(T00)
pl.figure(10022305)
#levels = np.linspace(-2.0,2.0,21)
for kk in range(1):
    print(np.sum(T00))
    cs = pl.contourf(X, Y, T00, levels=levels,colorscale='Electric', fill=False,cmap=pl.cm.jet)#cmap="Blues")
    pl.xlim(-17.5,17.5)
    pl.ylim(-17.5,17.5)
    pl.plot(zxx,zyy,'k--')
    pl.plot(zyy,zxx,'k--')
    pl.xlabel('fm', fontsize=15)# make axis labels
    pl.ylabel('fm', fontsize=15)
    pl.xticks(fontsize=15)
    pl.yticks(fontsize=15)
    pl.colorbar(cs)
    pl.title('Pb-Pb, V2, second shift,Q(0,0)_imag, b={}fm'.format(istep), fontsize=13)# give plot a title
    pl.savefig("Png/V2_Q_05_imag{}.png".format(istep))

T00 = []
temp = aa[:,6]
for ii in range(length):
    mid = temp[int(ii*length):int(ii*length+length)]
    T00.append(mid)
T00 = np.array(T00)
pl.figure(10022306)
#levels = np.linspace(-2.0,2.0,21)
for kk in range(1):
    print(np.sum(T00))
    cs = pl.contourf(X, Y, T00, levels=levels,colorscale='Electric', fill=False,cmap=pl.cm.jet)#cmap="Blues")
    pl.xlim(-17.5,17.5)
    pl.ylim(-17.5,17.5)
    pl.plot(zxx,zyy,'k--')
    pl.plot(zyy,zxx,'k--')
    pl.xlabel('fm', fontsize=15)# make axis labels
    pl.ylabel('fm', fontsize=15)
    pl.xticks(fontsize=15)
    pl.yticks(fontsize=15)
    pl.colorbar(cs)
    pl.title('Pb-Pb, V2, second shift,Q(0,0)_imag, b={}fm'.format(istep), fontsize=13)# give plot a title
    pl.savefig("Png/V2_Q_06_imag{}.png".format(istep))

T00 = []
temp = aa[:,7]
for ii in range(length):
    mid = temp[int(ii*length):int(ii*length+length)]
    T00.append(mid)
T00 = np.array(T00)
pl.figure(10022307)
#levels = np.linspace(-2.0,2.0,21)
for kk in range(1):
    print(np.sum(T00))
    cs = pl.contourf(X, Y, T00, levels=levels,colorscale='Electric', fill=False,cmap=pl.cm.jet)#cmap="Blues")
    pl.xlim(-17.5,17.5)
    pl.ylim(-17.5,17.5)
    pl.plot(zxx,zyy,'k--')
    pl.plot(zyy,zxx,'k--')
    pl.xlabel('fm', fontsize=15)# make axis labels
    pl.ylabel('fm', fontsize=15)
    pl.xticks(fontsize=15)
    pl.yticks(fontsize=15)
    pl.colorbar(cs)
    pl.title('Pb-Pb, V2, second shift,Q(0,0)_imag, b={}fm'.format(istep), fontsize=13)# give plot a title
    pl.savefig("Png/V2_Q_07_imag{}.png".format(istep))

T00 = []
temp = aa[:,8]
for ii in range(length):
    mid = temp[int(ii*length):int(ii*length+length)]
    T00.append(mid)
T00 = np.array(T00)
pl.figure(10022308)
#levels = np.linspace(-2.0,2.0,21)
for kk in range(1):
    print(np.sum(T00))
    cs = pl.contourf(X, Y, T00, levels=levels,colorscale='Electric', fill=False,cmap=pl.cm.jet)#cmap="Blues")
    pl.xlim(-17.5,17.5)
    pl.ylim(-17.5,17.5)
    pl.plot(zxx,zyy,'k--')
    pl.plot(zyy,zxx,'k--')
    pl.xlabel('fm', fontsize=15)# make axis labels
    pl.ylabel('fm', fontsize=15)
    pl.xticks(fontsize=15)
    pl.yticks(fontsize=15)
    pl.colorbar(cs)
    pl.title('Pb-Pb, V2, second shift,Q(0,0)_imag, b={}fm'.format(istep), fontsize=13)# give plot a title
    pl.savefig("Png/V2_Q_08_imag{}.png".format(istep))

T00 = []
temp = aa[:,9]
for ii in range(length):
    mid = temp[int(ii*length):int(ii*length+length)]
    T00.append(mid)
T00 = np.array(T00)
pl.figure(10022309)
#levels = np.linspace(-2.0,2.0,21)
for kk in range(1):
    print(np.sum(T00))
    cs = pl.contourf(X, Y, T00, levels=levels,colorscale='Electric', fill=False,cmap=pl.cm.jet)#cmap="Blues")
    pl.xlim(-17.5,17.5)
    pl.ylim(-17.5,17.5)
    pl.plot(zxx,zyy,'k--')
    pl.plot(zyy,zxx,'k--')
    pl.xlabel('fm', fontsize=15)# make axis labels
    pl.ylabel('fm', fontsize=15)
    pl.xticks(fontsize=15)
    pl.yticks(fontsize=15)
    pl.colorbar(cs)
    pl.title('Pb-Pb, V2, second shift,Q(0,0)_imag, b={}fm'.format(istep), fontsize=13)# give plot a title
    pl.savefig("Png/V2_Q_09_imag{}.png".format(istep))

