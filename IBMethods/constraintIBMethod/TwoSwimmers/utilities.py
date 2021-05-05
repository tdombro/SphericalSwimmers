# remove comments
import numpy as np
import os



# rotate a vector by an angle theta around unit vector u
# u = []
def rotate3D(theta, unitVec, vector):
    ux = unitVec[0]
    uy = unitVec[1]
    uz = unitVec[2]
    c = np.cos(theta)
    s = np.sin(theta)
    # make rotation matrix
    # see https://en.wikipedia.org/wiki/Rotation_matrix
    R11 = c + ux*ux*(1 - c)
    R12 = ux*uy*(1 - c) - uz*s
    R13 = ux*uz*(1 - c) + uy*s

    R21 = uy*ux*(1 - c) + uz*s
    R22 = c + uy*uy*(1 - c)
    R23 = uy*uz*(1 - c) - ux*s

    R31 = uz*ux*(1 - c) - uy*s
    R32 = uz*uy*(1 - c) + ux*s
    R33 = c + uz*uz*(1 - c)

    x = vector[0]
    y = vector[1]
    z = vector[2]

    xPrime = x*R11 + y*R12 + z*R13
    yPrime = x*R21 + y*R22 + z*R23
    zPrime = x*R31 + y*R32 + z*R33
    vector[0] = xPrime
    vector[1] = yPrime
    vector[2] = zPrime
    return

# return u cross v
def cross3D(u, v):
    wx = u[1]*v[2] - u[2]*v[1]
    wy = -(u[0]*v[2] - v[0]*u[2])
    wz = u[0]*v[1] - v[0]*u[1]
    return [wx, wy, wz]

# angle between 2 vector u and v
def angleRad(u, v):
    ndim = len(u)
    uv = 0
    # dot
    for d in range(ndim):
        uv += u[d]*v[d]
    # magnitude u and v
    magU = 0
    magV = 0
    for d in range(ndim):
        magU += u[d]*u[d]
        magV += v[d]*v[d]
    if magU*magV != 0:
        costheta = uv/(np.sqrt(magU)*np.sqrt(magV))
    return np.arccos(costheta)
# normalize a vector
# work for 2D and 3D
def norm(u):
    ndim = len(u)
    mag = 0.0
    for d in range(ndim):
        mag += u[d]*u[d]
    if mag != 0:
        mag = np.sqrt(mag)
        for d in range(ndim):
            u[d] /= mag
    return

# work for 2D and 3D
def distance(u, v):
    distSq = 0
    ndim = len(u)
    for d in range(ndim):
        distSq += (u[d] - v[d])*(u[d] - v[d])
    return np.sqrt(distSq)

# box = [[xlo, xhi],[ylo, yhi]]
def checkOutside(point, box):
    dx = dy = 0
    x = point[0]
    y = point[1]
    xlo = box[0][0]
    xhi = box[0][1]
    ylo = box[1][0]
    yhi = box[1][1]
    if x < xlo:
        dx = 1
    if x > xhi:
        dx = -1
    if y < ylo:
        dy = 1
    if y > yhi:
        dy = -1
    out = (dx != 0 or dy != 0)
    return out, [dx, dy]

#get the first number in the fist line in file
def get_fist_num_infile(s_infile):
    h = open(s_infile, 'r')
    fist_num = int(h.read().split()[0])
    h.close()
    return fist_num

def remove_comment(txt):
  cc = ['#', '!', '//']
  for c in cc:
      n = txt.find(c)
      if n > -1: return txt[0:n]
  return txt


# af = array of float
# f = float
# b = bool
# i = int
# s = string
# aw = array of words
# ARG = filename, key, type OR
#       filename, session, key, type
def read_from_input(*argc):
    narg = len(argc)
    file_name = argc[0]
    f = open(file_name, 'r')

    if narg == 3:
        key  = argc[1]
        type = argc[2]
        lines = f.readlines()
    elif narg == 4:
        session = argc[1]
        key     = argc[2]
        type    = argc[3]
        allLines = f.readlines()
        lines = []
        found = False
        for aline in allLines:
            if found:
                lines.append(aline)
            words = remove_comment(aline).split()
            if (len(words) > 0) and words[0] == session:
                found = True
            for w in words:
                if w == '}': found = False

    else:
        print('read_from_input: INVALID ARGUMENTS')
        quit(1)

    #print(lines)

    for aline in lines:
        words = remove_comment(aline).split()

        if (len(words) > 2) and words[0] == key:
            if type == 'af' or type == 'ArrayFloat':
                af = []
                for c in words[2:]:
                    cc = c.replace(",", "")  # remove trailing "," if any
                    af.append(float(cc))
                return af
            elif type == 'f' or type == 'Float':
                return float(words[2])
            elif type == 'b' or type == 'Bool':
                s = words[2].lower()
                if s[0] == '"' or s[0] == "'":
                    w = s[1:-1]
                else:
                    w = s
                if w == 'y' or w == 'yes' or 'true':
                    return True
                else:
                    return False
            elif type == 'i' or type == 'Int':
                return int(words[2])
            elif type == 's' or type == 'String':
                # remove leading and tailing " if any
                s = ' '.join(words[2:])
                if s[0] == '"' or s[0] == "'":
                    return s[1:-1]
                else:
                    return s
            elif type == 'aw' or 'ArrayWord':
                # print(words[2:])
                # print(words[2:-1])
                arr = []
                for w in words[2:]:
                    wcut = ''
                    for c in w:
                        if c == '"' or c == "'" or c == ",":
                            continue
                        wcut += c
                    if wcut != '':
                        arr.append(wcut)
                return arr
            else:
                return words[2:]

    print('read_from_input: KEY %s NOT FOUND in %s' % (key, file_name))
    quit(1)
    return

def read_from_input2(s_file_name, s_key, s_subkey, s_type):
    f = open(s_file_name, 'r')
    lines = f.readlines()
    outer_dict = {}
    print('FUNCTION: read_from_input')
    i = 0
    while i < len(lines):
        words = remove_comment(lines[i]).split()
        l = len(words)
        i += 1

        if l == 0: continue
        if l == 1:
            key = words[0]
            # if the next line not starting by '{'
            # now i already points to the next line
            if lines[i][0] != '{':
                print('Error: invalid input file at line %d' % (i + 1))
                quit()
            else:
                i += 1
                inter_dict = {}
                while lines[i][0] != '}':
                    wo = remove_comment(lines[i]).split()
                    inter_dict[wo[0]] = wo[2:]
                    i += 1
                i += 1
            outer_dict[key] = inter_dict
        if l >= 3:
            key = words[0]
            outer_dict[key] = words[2:]



    if s_key == '':
        ret = outer_dict[s_subkey]
    else:
        ret = outer_dict[s_key][s_subkey]

    if s_type == 'f':
        print('=' * 30)
        return float(ret)
    elif s_type == 's':
        print('=' * 30)
        return ret
    elif s_type == 'lf':
        l_ret = []
        for x in ret:
            l_ret.append(float(x))
        return l_ret
    else:
        print('Error: unsupported type: ', s_type)
        print('=' * 30)
        quit()

# ******************************************
# Distance between 2 identified vertices
# vertice should has from [x, y] as a list
# ******************************************
def distance2d(u, v):
    # print(u,v)
    return np.sqrt((u[1] - v[1]) * (u[1] - v[1]) + (u[0] - v[0]) * (u[0] - v[0]))


# domain = [Lx, Ly]
def distance2dPBC(u, v, domain, pbc):
    dx = u[0] - v[0]
    dy = u[1] - v[1]
    Lx = domain[0]
    Ly = domain[1]
    if pbc:
        if dx > Lx / 2: dx -= Lx
        if dx < -Lx / 2: dx += Lx
        if dy > Ly / 2: dy -= Ly
        if dy < -Ly / 2: dy += Ly
    return np.sqrt(dx*dx + dy*dy)

# depth first search
def dfs(graph, start):
    visited, stack = set(), [start]
    while stack:
        vertex = stack.pop()
        if vertex not in visited:
            visited.add(vertex)
            stack.extend(graph[vertex] - visited)
    return visited

# ******************************************
# displace set of points by a vector
# ******************************************
def displace(in_array, vector):
    if len(in_array) == 0:
        return
    m = len(in_array[0])
    for i in range(len(in_array)):
        for j in range(m):
            in_array[i][j] += vector[j]

# ******************************************
# rotate set of Points by an angle theta radian
# ******************************************
def rotate2d(in_array_2d, angle):
    # rotation maxtrix
    # | cos -sin|
    # | sin  cos|
    cos = np.cos(angle)
    sin = np.sin(angle)
    for i in range(len(in_array_2d)):
        x = in_array_2d[i][0]
        y = in_array_2d[i][1]
        in_array_2d[i][0] = x * cos - y * sin
        in_array_2d[i][1] = sin * x + cos * y

def make_directory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
    return

# angle btw 2 vectors
def angle2d(first, middle, end):
    u = [first[0] - middle[0], first[1] - middle[1]]
    v = [end[0] - middle[0], end[1] - middle[1]]
    uv = u[0]*v[0] + u[1]*v[1]
    au = np.sqrt(u[0]*u[0] + u[1]*u[1])
    av = np.sqrt(v[0]*v[0] + v[1]*v[1])
    costheta = uv/(au*av)
    return np.arccos(costheta)*180.0/np.pi

# work for list of [], that is [[],[]]
def shift_array2(in_array, n):
    ret_array = []
    for x in in_array:
        ret_array.append([x[0] + n, x[1] + n])
    return ret_array

def save_to_file(list, filename, dirname):
    if len(list) == 0:
        return

    make_directory(dirname)
    f = open(dirname + '/' + filename, 'w')
    f.write('# %s\n' % (os.path.basename(__file__)))
    f.write('# %d lines\n' % (len(list)))
    for row in list:
        for x in row:
            f.write('%g ' % x)
        f.write('\n')

    f.close()
    print('OutFile: %s' % filename)
    return


def plot_list(data, plottype, scale, title, legend, label, figname, dirname, showplot):
    import matplotlib.pyplot as plt
    # scale = [x_scale, y_scale]
    if len(data) == 0:
        return

    make_directory(dirname)
    fig, ax = plt.subplots()
    x = []
    n = len(data[0])
    # print('n = %d' % n)
    # 1st col = x
    # 2nd col = y
    # nth col = yn
    y = []
    for i in range(n - 1):
        y.append([])

    for v in data:
        x.append(v[0]/scale[0])
        for i in range(n - 1):
            y[i].append(v[i + 1]/scale[i + 1])

    lines = []
    if plottype == 'SCATTER':
        for i in range(n - 1):
            ax.scatter(x, y[i], s=1, alpha=0.5)
            #ax.add_artist(legend[i])
    elif plottype == 'LINE':
        for i in range(n - 1):
            lines += ax.plot(x, y[i], marker='.', linestyle='-')
    else:
        print('plot_list::error::invalid plottype')
        quit()

    ax.legend(lines[:], legend, loc='upper right', frameon=False)
    ylo = []
    yhi = []
    for i in range(n - 1):
        ylo.append(min(y[i]))
        yhi.append(max(y[i]))

    ax.set(ylim=(1.1*min(ylo), 1.1*max(yhi)))
    #ax.set_aspect('equal')
    plt.title(title, color='green')
    plt.grid(True)
    plt.xlabel(label[0])
    plt.ylabel(label[1])
    plt.savefig(dirname + '/' + figname)

    if showplot:
        plt.show()
    plt.close()
    return


def cast_type(datType, dat):
    if datType == 'float':
        return float(dat)
    if datType == 'double':
        return float(dat)
    if datType == 'int':
        return int(dat)
    if datType == 'string':
        return str(dat)
    return dat


def save_points_file(points, s_format, text, dirname, filename):
    if len(points) == 0:
        return
    make_directory(dirname)
    fout = open(dirname + '/' + filename, 'w')
    m = len(points)
    if len(points) == 0: return
    num_column = len(points[0])

    fformat = []
    for i in range(num_column):
        fformat.append('%' + s_format[i] + '  ')

    fout.write('%i %s\n' % (m, text))
    for i in range(m):
        for col in range(num_column):
            fout.write(fformat[col] % (points[i][col]))
        fout.write('\n')
    fout.close()


def plot_vertices_and_bonds(l_vertices, l_bonds, dirname, filename):
    import matplotlib.pyplot as plt
    #print('-'*30)
    #print("FUNCTION: plot_vertices_and_bonds: ")
    fig, ax = plt.subplots()

    listx = []
    listy = []
    for v in l_vertices:
        listx.append(v[0])
        listy.append(v[1])

    for b in l_bonds:
        v1 = b[0]
        v2 = b[1]
        x = [l_vertices[v1][0], l_vertices[v2][0]]
        y = [l_vertices[v1][1], l_vertices[v2][1]]
        ax.plot(x, y, linestyle ="solid", color = 'g', linewidth = 1)

    #scale = 1e-2*len(l_vertices)
    ax.scatter(listx, listy, color='b', alpha=0.9, s=5e-2 * len(l_vertices))
    ax.set_aspect('equal')
    plt.grid(True)
    if filename != '':
        make_directory(dirname)
        plt.savefig(dirname + '/' + filename + '.png')
    else:
        plt.show()
    plt.close()
    #print('=' * 30)
    return  # // plot_sphere()


def make_input2d(botnames, anchornames, outfile):
    f = open(outfile, 'w')
    f.write("\tposn_shift = 0.0, 0.0\n")
    f.write("\tmax_levels = MAX_LEVELS\n\n")
    f.write('\tstructure_names = ')
    f.write( '"' + botnames[0] + '"')
    
    for name in (botnames[1:] + anchornames):
        f.write(', "' + name + '"')
    
    f.write('\n')
    
    for name in (botnames + anchornames):
        f.write("\n\t\t" + name + " { " + "\n")
        f.write("\t\t level_number = MAX_LEVELS - 1 " + "\n" +
                "\t\t uniform_spring_stiffness = 0.0\n")
        f.write("\t\t}\n")
    f.close()
def read_xy_file(filename):
    lines = open(filename, 'r').readlines()
    vertices =[]
    for line in lines:
        w = remove_comment(line).split()
        if len(w) <= 1:
            continue
        vertices.append([float(w[0]), float(w[1])])

    return vertices

def plot_vertices_list(list_vertices, point_sizes):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()

    for i in range(len(list_vertices)):
        vert = list_vertices[i]
        x = []
        y = []
        color = 'C%d' % ((i + 1) % 10)
        #color = 'C%d' % (i + 1)
        for v in vert:
            x.append(v[0])
            y.append(v[1])

        #scale = 1e-2*len(l_vertices)
        ax.scatter(x, y, color=color, alpha=0.9, s=point_sizes)
        ax.set_aspect('equal')
        plt.grid(True)

    plt.show()
    plt.close()
    #print('=' * 30)
    return  # // plot_sphere()

def triangle(in_vert):
    from scipy.spatial import Delaunay
    bonds = []

    tri = Delaunay(in_vert)
    indptr = tri.vertex_neighbor_vertices[0]
    vertices = tri.vertex_neighbor_vertices[1]
    simplices = tri.simplices
    for k in range(len(in_vert)):
        ver = vertices[indptr[k]:indptr[k + 1]]
        for kn in ver:
            if kn > k:
                bonds.append([k, kn])
    return simplices, bonds

def polygon_area(points):
    from scipy.spatial import ConvexHull
    hull = ConvexHull(points)
    return hull.volume


def constraintib_database(*argc):
    outfile = argc[0]
    structure_names = argc[1]
    trans_mom = argc[2]
    rot_mom = argc[3]
    lag_pos_update = argc[4]
    centers = argc[5]

    all_names = [name for name in structure_names]

    if len(argc) > 6:
        bdry_name = argc[6]
        all_names.append(argc[6])

    # IBStandardInitializer
    f = open(outfile, 'w')
    f.write('IBStandardInitializer {\n')
    f.write("\t max_levels      = MAX_LEVELS\n")
    f.write("\t structure_names = ")
    for name in all_names[0:-1]:
        f.write('"' + name + '"' + ', ')
    f.write('"' + all_names[-1] + '"' + '\n')

    for name in all_names:
        f.write("\t " + name + " {\n")
        f.write("\t\tlevel_number = MAX_LEVELS - 1 \n" +
                "\t }\n")
    f.write("}\n\n")

    # ConstraintIBKinematics
    f.write('ConstraintIBKinematics {\n')
    id = 0
    for name in structure_names:
        f.write("\t " + name + " {\n")
        f.write('\t\tstructure_names                  = ' + '"' + name + '"\n')
        f.write("\t\tstructure_levels                 = MAX_LEVELS - 1\n")
        f.write("\t\tcalculate_translational_momentum = %d,%d,%d\n" % (trans_mom[0], trans_mom[1], trans_mom[2]))
        f.write("\t\tcalculate_rotational_momentum    = %d,%d,%d\n" % (rot_mom[0], rot_mom[1], rot_mom[2]))
        f.write('\t\tlag_position_update_method       = ' + '"' + lag_pos_update + '"\n')
        f.write('\t\ttagged_pt_identifier             = MAX_LEVELS - 1, 0\n')
        # f.write('\t\tbody_mesh                        = BODY_MESH\n')
        f.write('\t\tradius_0                         = RADIUS_LARGE\n')
        f.write('\t\tradius_1                         = RADIUS_SMALL\n')
        f.write('\t\tamplitude                        = AMPLITUDE\n')
        f.write('\t\tfrequency                        = FREQUENCY\n')
        cent = centers[id]
        for j in range(len(cent)):
            f.write('\t\tcenter_%d                         = %g, %g\n' % (j, cent[j][0], cent[j][1]))
        f.write("\t }\n")
        id += 1

    if len(argc) > 6:
        f.write("\t " + bdry_name + " {\n")
        f.write('\t\tstructure_names                  = ' + '"' + bdry_name + '"\n')
        f.write("\t\tstructure_levels                 = MAX_LEVELS - 1\n")
        f.write("\t\tcalculate_translational_momentum = 0, 0, 0\n")
        f.write("\t\tcalculate_rotational_momentum    = 0, 0, 0\n")
        f.write('\t\tlag_position_update_method       = \"CONSTRAINT_VELOCITY\"\n')
        f.write('\t\ttagged_pt_identifier             = MAX_LEVELS - 1, 0\n')
        f.write("\t }\n")

    f.write("}\n")
    f.close()
    print('copy %s into input2d' % outfile)
    return