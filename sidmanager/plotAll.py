"""
05-05-2013
@author: Nicolas Moreno Ch.
e-mail: nicolas.morenochaparro@kaust.edu.sa

It centralize different kind of plots, to be instansiated from any script, in order to simplify the plotting functions. In addition bring some homogeneity in all the plot I generate.

"""

import sys, os
import numpy as np
import pylab as plt
#import sketchPlot as sp
import itertools
import matplotlib
from matplotlib.transforms import offset_copy
matplotlib.rcParams['legend.fancybox'] = True


"""
This module plots and save all the type of postprocessing plots I made. The input are the initial and the final index for the the simulations

"""

def plotHistogram(array, bins, xlabel, ylabel, name, normed, legend=None):
    plt.figure() 
    plot = plt.hist(array, bins, normed=normed, alpha=0.8)
    #plt.hist(rog0, 10, normed=False, color='r', alpha=0.3)
    #plt.xlim(11,14)
    if legend != None:
        plt.figtext(0.815, 0.013, legend, color='black', weight='roman',
           size='x-small')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(name+'.eps')

def plotTrajectories(array, xdata, xlabel, ylabel, name, maxX, maxY, grid, marker, linestyle,legend,stds='a',exten='.eps'):
    p = plt.figure()
    #   plt.plot(xdata, array, c='b', marker=marker, ls='-')
    #plt.plot(xdata, array, c='b')
    plt.plot(xdata, array, marker=marker, linestyle=linestyle)
    if isinstance(stds, str)==False:
        plt.fill_between(xdata, array-stds, array+stds,
                         alpha=0.05, facecolor='#E6F2FC', linewidth=1, linestyle=':')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if legend != None:
        plt.legend(legend)
    plt.grid(grid)
    if maxX[1]>0:
        plt.xlim(maxX)
    if maxY[1]>0:
        plt.ylim(maxY)
    plt.savefig("%s%s"%(name, exten))

def plotComparison(ydata, xdata, xlabel, ylabel, name, maxX, maxY, grid, marker, linestyle,legend, sd = 0, position=[0,0],exten='.eps',
                   title ='', scale = 'normal', deb = False, annotate = 1, inset=False, insetxlim = [0, 50], insetylim = [0.4, 0.7],
                   mksize=6, offset=[0.3, 0.3, 10, 10], freqData = 1, legTitle ="", xshift = False, filledM = True, fittedC = 0, linewidth=1.5):

    markers = itertools.cycle(['o','D','v', '^', '*', 's', 'x', '>', '.', '<', 'D'])
    #marker = ['o','+','v', '^', '*', 's', 'x', '>', '.', '<']
    
    fig = plt.figure(figsize=(8., 7.2), dpi=100)
    
    if annotate!=1:
        weig, hei =0.65, 0.7
        weig2, hei2 = 0.25, 0.25
        
        ax = fig.add_axes([0.15, 0.15, weig, hei])
        ax.set_ylabel(ylabel, {'fontsize'   : 30 })
        ax.set_xlabel(xlabel, {'fontsize'   : 30 })
        if inset:
            inset_ax = fig.add_axes([0.15+weig-weig2, 0.15+hei-hei2, weig2, hei2]) # X, Y, width, height
    else:
        ax = plt.subplot(1,1,1)
        transOffset = offset_copy(ax.transData, fig=fig,
                            x = offset[0], y=offset[1], units='inches')
        plt.ylabel(ylabel, {'fontsize'   : 30 })
        plt.xlabel(xlabel, {'fontsize'   : 30 })

    if title!='':
        plt.title(title)
    ncurves  = len(ydata)
    loc = 0
    q = []
    frac = 0.1
   
    for p in range(ncurves):
        if filledM == False:
            fm = 'none'
        else:
            fm = color=plt.cm.jet(frac)
        #A modified fractionation of the curves ids so I can change the contrast when only two curves are plotted.
        if ncurves == 2:
            frac = frac+p*0.8
        else:
            frac = (1.*p)/ncurves  
            
        loc+=(len(xdata[p])/ncurves) #I modified this addition in order to deal with plots with many points so the labels are not placed too close
        
        #loc +=1
        pts = np.size(xdata[p])
        freq = int(pts/(pts*freqData))
        if scale =='log':
            plt.loglog(xdata[p], ydata[p], marker=markers.next(), markersize=mksize, linestyle=linestyle,  color=plt.cm.jet(frac), markevery=freq, mfc=fm, mec=plt.cm.jet(frac), mew=linewidth)
        else:
           mnext =markers.next()
           if isinstance(sd, int):
               q = (ax.plot(xdata[p], ydata[p], marker=mnext, markersize=mksize,linestyle=linestyle, color=plt.cm.jet(frac), markevery=freq, mfc=fm, mec=plt.cm.jet(frac), mew=linewidth))
           else: 
               ##xshift variable can take values of 0 or 1 if the data want to be shifted to make easy their visualizaiton.
               xd = [x+sh for (x,sh) in zip(xdata[p], np.linspace(0,1,np.size(xdata[p]))*xshift)]
           #print xd
               ax.errorbar(xd, ydata[p], sd[p],  marker=mnext, markersize=mksize,linestyle=linestyle, color=plt.cm.jet(frac),  markevery=freq, mfc=fm, mec=plt.cm.jet(frac), mew=linewidth)

        if isinstance(fittedC, int)==False:
            ax.plot(fittedC[0], fittedC[1], linewidth=linewidth,linestyle='-', color=plt.cm.jet(frac))
            
        if loc+2 >= len(xdata[p]): #This works as a cycle to change the posistion of the label dynacmically
            loc = 1
        midx, midy = len(xdata[p])-loc, len(ydata[p])-loc

        
        
        if position != [0, 0]:
            midx, midy = position[0], position[1]

        if annotate == 1:
            #plt.plot([xdata[p][midx],xdata[p][midx]+0.1], [ydata[p][midy], ydata[p][midy]+0.1], '-k', lw = 0.5,  color=plt.cm.jet(1.*p/ncurves))
            #plt.text(xdata[p][midx], ydata[p][midy], '%s' % (legend[p]), transform=transOffset, color=plt.cm.jet(1.*p/ncurves))            
            # plt.annotate(legend[p], (xdata[p][midx], ydata[p][midy]),
            #        backgroundcolor='w',
            #        color='b',
            #        va='center',
            #        ha='center',
            #   bbox=dict(boxstyle="round", alpha=0.5))
            ix = (xdata[p][np.isfinite(xdata[p])]).mean()/offset[2]  ##Included this offset value as a third parameter when the lables are not properly placed in the plot
            iy = (ydata[p][np.isfinite(ydata[p])]).mean()/offset[3]
            #print ix, iy, 'ixiy'
            ax.annotate('%s' % (legend[p]),
            xy=(xdata[p][midx], ydata[p][midy]), xycoords='data',
            xytext=(xdata[p][midx]+ix, ydata[p][midy]+iy), textcoords='data', fontsize = 20,
            arrowprops=dict(arrowstyle="->", linewidth=3,#linestyle="dashed",
                            color=plt.cm.jet(frac),
                            patchB=None,
                            shrinkB=0,
                            connectionstyle="arc3",
                            ),
            )

    if annotate == 0:
            pepe = 1
            ax.legend(legend,bbox_to_anchor=(1.05, 1), loc = 2, borderaxespad=0., shadow = True, title=legTitle, numpoints=1, fontsize = 15)
        
    
    if maxX[1]>0:
        ax.set_xlim(maxX)
    if maxY[1]>0:
        ax.set_ylim(maxY)
    ax.grid(grid)
    if deb: plt.show()


    if inset:
        inset_ax.set_ylabel(ylabel, {'fontsize'   : 15 })
        inset_ax.set_xlabel(xlabel, {'fontsize'   : 15 })
        markers = itertools.cycle(['o','D','v', '^', '*', 's', 'x', '>', '.', '<', 'D'])
        for p in range(ncurves):
            inset = True
            if inset:
                #axIn = fig.add_axes([0.2, 0.5, 0.4, 0.3]) # inset axes
                #axIn = plt.subplot()
                #axIn.tight_layout()
                # inset


                inset_ax.plot(xdata[p], ydata[p], marker=markers.next(), markersize=mksize , linestyle=linestyle, color=plt.cm.jet(1.*p/ncurves), markevery=freq, mfc=fm)
                #inset_ax.set_title('zoom near origin')
                # set axis range
                inset_ax.set_xlim(insetxlim)
                inset_ax.set_ylim(insetylim)
                # set axis tick locations
                #inset_ax.set_yticks([0, 0.005, 0.01])
            #inset_ax.set_xticks([-0.1,0,.1]);

        
    plt.savefig(name+exten)

def plotSketch(array, xdata, xlabel, ylabel, name, maxX, dataLabel, title):
    ax = plt.axes()
    ax.plot(xdata, array, '-k', lw = 0.5)
    if title!=0:
        ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    i = 10   ##Point index of the data to use as starting point of the label line
    multx =10
    multy =10
    inf  = array[i]
    sup  = np.abs(array[i+1]-array[i])*mult+array[i]
    supX = np.abs(xdata[i+1]-xdata[i])*mult+xdata[i]

    ax.text(supX, sup, dataLabel)
    ax.plot([xdata[10],supX],[inf, sup], '-k', lw = 0.5 )


    if maxX>0:
        ax.set_xlim(0, maxLim)
    #ax.set_ylim(0,1)
    #XKCDify the axes -- this operates in-place
    ax = XKCDify(ax, xaxis_loc=0.0, yaxis_loc=0.0,
            xaxis_arrow='+-', yaxis_arrow='+-',
                        expand_axes=True)

    plt.savefig(name+'.eps')

    
#plotAll(init,fin,index)
