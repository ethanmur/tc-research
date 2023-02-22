#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 11:25:00 2021
https://skywatch.colorado.edu

@author: zhwa2432
"""

import urllib
import matplotlib.pyplot as plt
import numpy as np
import glob 
''' 
url='https://sundowner.colorado.edu/weather/atoc8/wxobs20160807.txt'
lines  = urllib.request.urlopen(url).readlines()
for line in lines[3:]:
    entries = line.decode("utf-8").split()
    print(entries[0:3],entries[1][0:2])
    #for a in entries[0:3]:
    #    print(a)


url='https://seb.noaa.gov/pub/flight/ASPEN_Data/20210816H1/D20210816_085906QC.frd'
lines  = urllib.request.urlopen(url).readlines()
for line in lines[0:10]:
    entries = line.decode("utf-8").split()
    print(entries)
'''    
def get_frd(folder):
    outdir='/Users/zhwa2432/data/CRL/frd/'
    outdir='/Users/zhwa2432/data/CRL/frd/'+folder
    #url='https://seb.noaa.gov/pub/flight/ASPEN_Data/'+folder
    #url='https://seb.noaa.gov/pub/flight/ASPEN_Data/20210816H1/'
    url='https://seb.omao.noaa.gov/pub/acdata/2021/AVAPS/'+folder+'ASPEN_DATA/'
    print(url)
    lines  = urllib.request.urlopen(url).readlines()
    for line in lines[1:]:
        entries = line.decode("utf-8").split()
        if len(entries) > 1:
            a=entries[1]
            b1=a.find('"') #; print(b1)
            if b1 < 0: continue
            #b1=a.index('"')
            b2=a[b1+1:-1].index('"')
            filename=a[b1+1:b2+b1+1]
            #print(b1,b2, filename[-3:])
            if filename[-3:]=='frd':
                print(filename)
                text_files = glob.glob(outdir + filename, recursive = True)
                if len(text_files) < 1 :
                    #url='https://seb.noaa.gov/pub/flight/ASPEN_Data/20210816H1/'+filename
                    #url1='https://seb.noaa.gov/pub/flight/ASPEN_Data/'+folder +filename
                    url1='https://seb.omao.noaa.gov/pub/acdata/2021/AVAPS/'+folder+'ASPEN_DATA/' +filename
                    print(url1)
                    #lines  = urllib.request.urlopen(url).readlines()
                    urllib.request.urlretrieve(url1, outdir + filename)
               
folder='20221008H1/'
#get_frd(folder)
outdir='/Users/zhwa2432/data/CRL/'+folder

import ssl
ssl._create_default_https_context = ssl._create_unverified_context

url0='https://seb.omao.noaa.gov/pub/acdata/2022/CRL/'+folder
lines0  = urllib.request.urlopen(url0).readlines()
print(lines0)
num=0
dir_all=[]
for line0 in lines0[0:-1]:
    entries = line0.decode("utf-8").split()
    print(entries,len(entries))
    if len(entries) > 1:
            a=entries[1]
            b1=a.find('"') #; print(b1)
            if b1 < 0: continue
            #b1=a.index('"')
            b2=a[b1+1:-1].index('"')
            if b2 < 22: continue
            filename=a[b1+1:b2+b1+1]
            dir_all.append(filename)
            print(' here:  ',num,'  ',filename)
            num=num+1
for i in range(0,num):#(298,304):
    filename=dir_all[i]
    print(i, '  ',filename )
    #get_frd(folder)
    url1=url0+filename
    print(url1)
    ##lines  = urllib.request.urlopen(url).readlines()
    urllib.request.urlretrieve(url1, outdir + filename)


'''
## find a file   
outdir='/Users/zhwa2432/data/CRL/frd/'
fname='D20210925_221941_PQC.frd'
text_files = glob.glob(outdir + fname, recursive = True)
print(len(text_files))
'''

'''
num=0
url0='https://seb.omao.noaa.gov/pub/flight/ASPEN_Data/'
url0='https://seb.omao.noaa.gov/pub/acdata/2021/AVAPS/20210816H1/ASPEN_DATA/'
lines0  = urllib.request.urlopen(url0).readlines()
for line0 in lines0[20:]:
    entries = line0.decode("utf-8").split()
    #print(entries,len(entries))
    if len(entries)<= 1 : continue  
    a=entries[1]
    b1=a.find('"')
    if b1 < 0: continue
    b2=a[b1+1:-1].index('"')
    folder=a[b1+1:b2+b1+1]
    print(b1,b2, folder, len(folder))
    if len(folder) != 11 : continue
    print(b1,b2, folder, len(folder))
    num=num+1
    #get_frd(folder)
print(num)
'''