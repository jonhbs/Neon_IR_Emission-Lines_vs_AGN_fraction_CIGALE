#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 11:56:50 2023

@author: Jonhatan Bernal Salinas
jobernals@unal.edu.co
"""

import numpy as np
import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import norm


def bootstrapCorr(n,fract,data,var1,var2,pdfhisto=False,kdehisto=False,conf=95):
    '''
    This functions do the bootstraping for the correlation coefficients
    Inputs:
    n: Number of samples
    fract: Fraction of data from the intial 
           dataframe to make the samples
    data: Data
    var1, var2: Variables to be correlated
    pdfhisto: If you want to print the histogram with a pdf
    snshisto_ If you want to print the histogram with the kde
    conf: Confidence Interval
    '''

    df = data.to_pandas()

    corr_bootstrap = []
    for i in range(n):
        sample = df.sample(frac=fract,replace=True, random_state=i)
        corr = sample[var1].corr(sample[var2])
        corr_bootstrap.append(corr)
     
    
    med = np.median(corr_bootstrap) #median
    mean = np.mean(corr_bootstrap) #mean
    stdv = np.std(corr_bootstrap) #standard deviation
    
    #Percentiles:
    per1=np.percentile(corr_bootstrap,(100-conf)/2)
    per2=np.percentile(corr_bootstrap,100-((100-conf)/2))
    
    print('Bootstraps results for the correlation coefficient between',var1,
          'y',var2,':')
    print('Samples lenght:',len(sample),'of',len(df))
    print('Number of samples:',n)
    print('Median:',np.median(corr_bootstrap))
    print('Mean:',mean)
    print('Variance:',np.var(corr_bootstrap))
    print('Standard Desviation:',stdv)
    print(conf,'% confidence interval:',np.percentile(corr_bootstrap,5),'-',
          np.percentile(corr_bootstrap,95))
    
    bins = np.linspace(min(corr_bootstrap), max(corr_bootstrap),25)
    
    if pdfhisto == True:
        plt.figure(figsize=(6,3.75))
        plt.hist(corr_bootstrap, bins=bins, density=True, histtype='step',label='Corr Bootstrap')
        plt.axvline(x=mean,ls='--',color='gray',label='PCC Mean: %.2f' %mean)
        plt.axvline(x=med,ls='--',color='blue',label='PCC Med: %.2f' %med)
        plt.axvline(x=per1,ls='--',color='yellow',label='CI')
        plt.axvline(x=per2,ls='--',color='yellow')
        plt.xlabel(r'Correlation Coefficient',fontsize=14)
        plt.ylabel('Density',fontsize=14)
        plt.grid(color='k', linestyle='--', linewidth=0.1)
        # Plot the PDF.
        xmin, xmax = plt.xlim()
        x = np.linspace(xmin, xmax, 100)
        mu, std = norm.fit(corr_bootstrap)# mean and standard deviation
        p = norm.pdf(x, mu, std)
        plt.plot(x, p, 'k', linewidth=0.5)
        plt.legend()
        plt.show()
    
    if kdehisto == True:
        plt.figure(figsize=(6,3.75))
        sns.histplot(corr_bootstrap,stat='density',kde=True,bins=bins,
                     color='k',label='Corr Bootstrap')
        plt.axvline(x=mean,ls='--',color='gray',label='PCC Mean: %.2f' %mean)
        plt.axvline(x=med,ls='--',color='blue',label='PCC Med: %.2f' %med)
        plt.axvline(x=per1,ls='--',color='yellow',label='CI')
        plt.axvline(x=per2,ls='--',color='yellow')
        plt.xlabel(r'Correlation Coefficient',fontsize=14)
        plt.ylabel('Density',fontsize=14)
        plt.grid(color='k', linestyle='--', linewidth=0.1)
        plt.legend()
        plt.show()
    
    return corr_bootstrap,mean,stdv,med,per1,per2