from ibp_boxes import *
bs=[b for b in boxes()]
mb=[b for b in bs if b['tag']=='MB']
# consecutive boxes sorted by X descending? print order
xs=[float(b['X'].mid()) for b in mb]; print(xs[:3], xs[-3:])
import numpy as np
for idx in (1200, 1800, 2100, 2350):
    b0,b1=mb[idx],mb[idx+40]
    dx=float(b1['X'].mid()-b0['X'].mid())
    for n in range(4):
        fd=(float(b1['F'][n].mid())-float(b0['F'][n].mid()))/dx
        mid=(float(b1['F'][n+1].mid())+float(b0['F'][n+1].mid()))/2
        print(idx,n,'FD',fd,'F^(n+1) formula/hull',mid, 'rad', float(b0['F'][n+1].rad()))
