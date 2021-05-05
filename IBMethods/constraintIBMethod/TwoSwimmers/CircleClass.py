import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, Delaunay
import os, sys, random, os.path
import utilities as ut
import shutil
from InputClass import InputDatabase


class Circle:
    def __init__(self, mesh_method, origin, radius, ds, num_circe):
        self.origin = origin
        self.radius = radius
        self.ds = ds
        self.ncircle = num_circe
        self.vertices = []
        self.mesh_method = mesh_method

        if mesh_method == 'POLAR':
            self.genMeshPolar()
        elif mesh_method == 'SQUARE':
            self.genMeshSquare()
        else:
            print('Circle::Init::Error::mesh method not found')
            quit()

    def genMeshPolar(self):
        circle = []
        origin = self.origin
        radius = self.radius
        ds = self.ds
        PI = np.pi

        for c in range(self.ncircle):
            X = []
            Y = []
            r = radius + c * ds
            nPoints = int(2 * PI * r / ds)
            dTheta = 2*PI/nPoints

            for i in range(nPoints):
                theta = i * dTheta
                x = r * np.cos(theta) + origin[0]
                y = r * np.sin(theta) + origin[1]
                circle.append([x, y])
        self.vertices = circle
        return

    def genMeshSquare(self):
        #print('FUNCTION: gen_square_circle')
        quarters = [[], [], [], []]
        radius = self.radius
        ds = self.ds
        origin = self.origin

        nx = int(np.ceil(radius/ds)) + self.ncircle
        rmax = nx * ds

        for i in range(nx, -1, -1):  # Nx,... 0 inclusive
            for j in range(nx + 1):
                x = i * ds
                y = j * ds
                r = np.sqrt(x * x + y * y)
                if r >= radius and r < rmax:
                    quarters[0].append([x + origin[0], y + origin[1]])  # 1
                    quarters[2].append([-x + origin[0], -y + origin[1]])  # 3
                    quarters[1].append([-x + origin[0], y + origin[1]])  # 2
                    quarters[3].append([x + origin[0], -y + origin[1]])  # 4

        # end for
        # sort the vertices
        m = len(quarters[0])
        circle = []
        for i in range(4):
            for j in range(m):
                index = j
                if i % 2 == 1:
                    index = m - j - 1
                circle.append(quarters[i][index])
        self.vertices = circle
        return

    def print(self):
        print('Circle::Print')
        print('number of vertices %d' % (len(self.vertices)))
        print('circle radius %g' %(self.radius))
        print('number of circles %d' %(self.ncircle))
        #print('='*60)


    def saveVertices(self, dirname, filename):
        ut.make_directory(dirname)
        f = open(dirname + '/' + filename + '.vertex', 'w')
        f.write('%d  # number of circles %d, mesh %s\n'
                % (len(self.vertices), self.ncircle, self.mesh_method))
        for v in self.vertices:
            f.write('%g %g\n' % (v[0], v[1]))

        #print('total points on circle : ', m)
        print("OutFile: ", filename + '.vertex')
        print('=' * 60)

    def getVertices(self):
        return self.vertices

    def getRadius(self):
        return self.radius

    def getNumVertices(self):
        return len(self.vertices)

    def getDx(self):
        return self.ds

    def plot(self):
        fig, ax = plt.subplots()
        X = []
        Y = []
        for pt in self.vertices:
            X.append(pt[0])
            Y.append(pt[1])
            # end for
            # colors = np.random.rand(len(circle))
            # area = np.pi*(0.5*d_dx)**2

        ax.scatter(X, Y, c='b', s=1, alpha=0.5)
        ax.axis('equal')
        ax.grid(True)
        plt.show()
        return

"""
a =Circle('POLAR', [1, 1], 5, 0.2, 2)
a.print()
a.plot()
a.saveVertices('2d','circe.vertex')
"""