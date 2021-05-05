#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Updated: 4/1020
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, Delaunay
import os, sys, random, os.path
import utilities as ut
import shutil
from InputClass import InputDatabase
class Centers:
    def __init__(self, origin, nbot, radii, list_length):
        self.nbot = nbot
        self.list_length = list_length
        self.origin = origin
        self.radii = radii
        self.two_centers = []


    # create a list of non-overlapping center of mass in a circle
    # list_length is center to center distance bw small and large sphere
    def genRandomCircle(self, circle_radius):
        self.circle_radius = circle_radius
        self.type = 'CIRCLE'
        origin = self.origin
        radii = self.radii # sphere radii
        nbot = self.nbot

        list_centers = []
        bot_radius = []

        for botId in range(nbot):
            # consider each bot as a big disk with this radius
            bot_radius.append(0.5 * (self.list_length[botId] + radii[0] + radii[1]))

        cycle = 0.0
        botId = 0  # current number of bot
        while botId < nbot:
            cycle += 1
            if cycle > 1e6:
                print("only found %d (requested %d) within %0.0f cycle" %
                      (botId, nbot, cycle))
                break

            ro = random.uniform(0.0, circle_radius - bot_radius[botId])
            phi = random.uniform(0.0, 2 * np.pi)
            x = ro * np.cos(phi) + origin[0]
            y = ro * np.sin(phi) + origin[1]

            candicate = [x, y]
            addto = True

            if botId > 0:
                for i in range(len(list_centers)):
                    x = list_centers[i]
                    if ut.distance2d(candicate, x) < (bot_radius[botId] + bot_radius[i]):
                        addto = False

            if addto:
                print("found bot ", botId + 1)
                botId += 1
                list_centers.append(candicate)

        two_centers = []

        for botId in range(len(list_centers)):
            c = list_centers[botId]
            # print(" Confined nBot ", bot)
            # print("c = ", c)
            # generate a random direction
            phi = random.uniform(0.0, 2 * np.pi)
            c1x = c[0] + (bot_radius[botId] - radii[0]) * np.cos(phi)
            c1y = c[1] + (bot_radius[botId] - radii[0]) * np.sin(phi)
            c2x = c[0] - (bot_radius[botId] - radii[1]) * np.cos(phi)
            c2y = c[1] - (bot_radius[botId] - radii[1]) * np.sin(phi)
            two_centers.append([[c1x, c1y], [c2x, c2y]])

        self.two_centers = two_centers
        return
    #i = 0
    #for b in two_centers:
    #    print(" bot %d  %0.6f  %0.6f  %0.6f  %0.6f" %
    #          (i, b[0][0], b[0][1], b[1][0], b[1][1]))
    #    i += 1
    #
    def genUniformRectangle(self, box, nrow, ncol, randupdown):
        if max(self.list_length) != min(self.list_length):
            print('Centers::genUniformRectangle::Error::This only works with identical bots')
            print('please make a change to bot lengths')
            quit()

        if nrow*ncol != self.nbot:
            print('Centers::genUniformRectangle::Error::Invalid nrow, ncol')
            quit()

        self.type = 'RECTANGLE'
        self.box = box
        origin = self.origin
        xlo = -0.5 * box[0] + origin[0]
        xhi = 0.5 * box[0] + origin[0]
        ylo = -0.5 * box[1] + origin[1]
        yhi = 0.5 * box[1] + origin[1]
        cc_dist = self.list_length[0]

        ax = box[0] / ncol
        ay = box[1] / nrow
        #print('ax = %g, ay = %g, nx = %g, ny = %g' % (ax, ay, nx, ny))
        two_centers = []
        for i in range(ncol):
            for j in range(nrow):
                xo = xlo + (i + 0.5) * ax
                yo = ylo + (j + 0.5) * ay
                x1 = x2 = xo
                y1 = yo + 0.5 * cc_dist
                y2 = yo - 0.5 * cc_dist
                c = [[x1, y1], [x2, y2]]
                # if we want a random orientation up and down
                if randupdown:
                    random.shuffle(c)
                two_centers.append(c)
        self.two_centers = two_centers

    # box [Lx, Ly] center around the origin
    def genRandomRectangle(self, box):
        self.box = box
        self.type = 'RECTANGLE'
        radii = self.radii
        nbot = self.nbot
        origin = self.origin

        print('genRandomRectangle: box = [%g, %g]' % (box[0], box[1]))

        xlo = -0.5*box[0] + origin[0]
        xhi =  0.5*box[0] + origin[0]
        ylo = -0.5*box[1] + origin[1]
        yhi =  0.5*box[1] + origin[1]

        # twocenters = []
        # if both nx and ny are non-zero
        # we create aa grid nx x ny swimmers
        # each bot is enclosed in a disk whose centers are specified in coms
        # and whose radius is bot_radius
        coms = []
        bot_radius = []
        for botId in range(self.nbot):
            bot_radius.append(0.5 * (self.list_length[botId] + radii[0] + radii[1]))

        cycle = 0
        botId = 0  # current number of bot
        while botId < nbot:
            cycle += 1
            if cycle > 1e6:
                print("only found %d (requested %d) within %0.0f cycle" %
                      (botId, nbot, cycle))
                break

            x = random.uniform(xlo + bot_radius[botId], xhi - bot_radius[botId])
            y = random.uniform(ylo + bot_radius[botId], yhi - bot_radius[botId])
            candicate = [x, y]
            addto = True

            if botId > 0:  # check to see if this bot overlaps with previous bots
                for i in range(len(coms)):
                    x = coms[i]
                    if ut.distance2d(candicate, x) < (bot_radius[botId] + bot_radius[i]):
                        addto = False

            if addto:
                print(" found bot ", botId)
                botId += 1
                coms.append(candicate)

        # after all bot coms were found
        two_centers = []

        for botId in range(len(coms)):
            c = coms[botId]
            # generate a random direction
            phi = random.uniform(0.0, 2 * np.pi)
            c1x = c[0] + (bot_radius[botId] - radii[0]) * np.cos(phi)
            c1y = c[1] + (bot_radius[botId] - radii[0]) * np.sin(phi)
            c2x = c[0] - (bot_radius[botId] - radii[1]) * np.cos(phi)
            c2y = c[1] - (bot_radius[botId] - radii[1]) * np.sin(phi)

            two_centers.append([[c1x, c1y], [c2x, c2y]])

        self.two_centers = two_centers

    def plot(self):
        type = self.type
        R = self.radii[0]
        r = self.radii[1]
        o = self.origin
        #nbot = len(self.two_centers)

        fig, ax = plt.subplots()

        if type == 'CIRCLE':
            rc = self.circle_radius
            ax.set(xlim=(-rc + o[0], rc + o[0]), ylim=(-rc + o[1], rc + o[1]))
            # plot the circular boundary
            ax.add_artist(plt.Circle((o[0], o[1]), rc, color='r', fill=False))
        elif type == 'RECTANGLE':
            box = self.box
            Lx = box[0]
            Ly = box[1]
            ax.set(xlim=(-0.5*Lx + o[0], 0.5*Lx + o[0]), ylim=(-0.5*Ly + o[1], 0.5*Ly + o[1]))
            ax.add_artist(plt.Rectangle((-0.5*Lx + o[0], -0.5*Ly + o[1]), Lx, Ly, color='black', fill=False))
        else:
            print('Centers::Plot::Error::Unknown plot type')
            quit()

        i = 0
        for x in self.two_centers:
            c = 'C%d' % (i%10)
            cl = [x[0][0], x[0][1]]
            cs = [x[1][0], x[1][1]]
            ax.add_artist(plt.Circle(cl, R, color=c, fill=False))
            ax.add_artist(plt.Circle(cs, r, color=c, fill=False))
            i += 1
        ax.set_aspect('equal')
        plt.show()
        return

    def getTwoCenters(self):
        return self.two_centers

    def getNbot(self):
        return self.nbot

"""
rcircle = 10
list_length = [4, 4, 4, 4, 2]
a = Centers([5, 5], 5, [2, 1], list_length)
a.genRandomCircle(rcircle)
#a.genRandomRectangle([20,20])
a.plot()
"""
"""
list_length = [4, 4, 4, 4]
b = Centers([5,4],4,[2,1],list_length)
b.genUniformRectangle([30,30],1,4,False)
b.plot()
"""


