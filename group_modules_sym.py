""" Group theory modules for testing with space groups """

import sys
import re
import math
import requests
import operator as op
import numpy as np
from numpy import linalg as LA
from fractions import Fraction
from itertools import zip_longest
from os import listdir


class Group:
    def __init__(self, num, matrix, index=None):
        # num: ITA number
        # matrix: transformation matrix
        self.num = num
        self.matrix = matrix
        self.index = index

class Subgroup:
    def __init__(self, num, normal=True, bieberbach=True):
        # num: ITA number
        # normal: boolean for normality
        self.num = num
        self.normal = normal
        self.bieberbach = bieberbach
    def __eq__(self, other):
        return isinstance(other, Subgroup) \
               and self.num == other.num \
               and self.normal == other.normal
    def __hash__(self):
        return hash((self.num, self.normal))

class Decomp:
    def __init__(self, gnum, snum, matrix, bgroups):
        self.gnum = gnum
        self.snum = snum
        self.matrix = matrix
        self.bgroups = bgroups


class TableEntry:
    def __init__(self, decomps):
        self.decomps = decomps
        self.names = self.getNames()
        self.printEntry()
    def getNames(self):
        with open('tex/sgroup_labels.tex') as file:
            return [line.strip() for line in file]
    def printEntry(self):
        with open('tex/table_entries/non_semidirect/entry_{}'.format(sys.argv[1]), 'w') as file:
            for decomp in self.decomps:
                self.printDecomp(decomp, file)
    def printDecomp(self, d, file):
        gname = self.names[int(d.gnum)-1]
        sname = self.names[int(d.snum)-1]
        bnames = self.getBNames(d.bgroups)
        normalB = self.getNormalB(d.bgroups)
        bieberbach = self.getBieberbach(d.bgroups)
        file.write('${}$ & ${}$ & ${}$ & ${}$ & {} & {} & {} \\\\\n'.format(d.gnum, gname, sname, bnames, bieberbach, 'Y', normalB))
    def getBNames(self, bgroups):
        bnames = self.names[int(bgroups[0].num)-1]
        if not bgroups[0].bieberbach: bnames = '{{\color{{red}}{}}}'.format(bnames)
        for i in range(1,len(bgroups)):
            name = self.names[int(bgroups[i].num)-1]
            if not bgroups[i].bieberbach: name = '{{\color{{red}}{}}}'.format(name)
            bnames += ',' + name
        return bnames
    def getNormalB(self, bgroups):
        if bgroups[0].normal: normalB = 'Y'
        else: normalB = 'N'
        for i in range(1,len(bgroups)):
            if bgroups[i].normal: normalB += ',' + 'Y'
            else: normalB += ',' + 'N'
        return normalB
    def getBieberbach(self, bgroups):
        if bgroups[0].bieberbach: bieberbach = 'Y'
        else: bieberbach = '{{\color{{red}}N}}'
        for i in range(1, len(bgroups)):
            if bgroups[i].bieberbach: bieberbach += ',' + 'Y'
            else: bieberbach += ',' + '{{\color{{red}}N}}'
        return bieberbach



        

# helper function for loadGroup
# sets the chosen column of the row depending on
# the string containing the coefficient of x,y,z
def setRowVal(row, ind, a): 
    if a == "":
        row[ind] = 1
    elif a == '-':
        row[ind] = -1
    else:
        row[ind] = float(Fraction(a).limit_denominator())


# helper function for loadGroup
# returns a row of a matrix from the string that
# describes the transformed coordinate x,y or z
def getRow(str):
    a = ""
    row = np.zeros(4)
    for i in str:
        if i == 'x':
            setRowVal(row, 0, a)
            a = ""
            continue
        if i == 'y':
            setRowVal(row, 1, a)
            a = ""
            continue
        if i == 'z':
            setRowVal(row, 2, a)
            a = ""
            continue
        if i == '+':
            continue
        a += i
    if a != "":
        row[3] = float(Fraction(a).limit_denominator())
    return row


# loads a group as a 3d array of 4x4 matrices
def loadGroup(filename):
    with open(filename) as file:
        A = [line.split() for line in file]
    G = np.zeros((4, 4, len(A)))
    for i in range(len(A)):
        G[0,:,i] = getRow(A[i][0])
        G[1,:,i] = getRow(A[i][1])
        G[2,:,i] = getRow(A[i][2])
        G[:,3,i] = G[:,3,i] % 1
        G[3,:,i] = [0, 0, 0, 1]
    return G


def loadCosetRep(rep):
    g = np.zeros((4,4))
    g[0,:] = getRow(rep[0])
    g[1,:] = getRow(rep[1])
    g[2,:] = getRow(rep[2])
    g[:,3] = g[:,3] % 1
    g[3,:] = [0, 0, 0, 1]
    return g


def getGenerator(rep):
    return rep[0] + "," + rep[1] + "," + rep[2] + "\n"


def getBiebElems(A):
    biebElems = [[],[]]
    for reps in A:
        matrices = []
        gens = []
        for rep in reps:
            g = loadCosetRep(rep)
            if not isSgroupElem(g):
                matrices.append([g])
                gens.append(getGenerator(rep))
        biebElems[0].append(matrices)
        biebElems[1].append(gens)
    return biebElems


# symGroups = []
biebGroups = []
matCombs = []
repCombs = []
Bgroups = []
# Sgroups = []


def combineElems(matTerms, repTerms, matAccum, repAccum):
    last = (len(matTerms) == 1)
    n = len(matTerms[0])
    for i in range(n):
        matItem = matAccum + matTerms[0][i]
        repItem = repAccum + repTerms[0][i]
        if last:
            matCombs.append(matItem)
            repCombs.append(repItem)
        else:
            combineElems(matTerms[1:], repTerms[1:], matItem, repItem)       


# loads a group as general positions (actions on R^3)
def loadGenPos(filename):
    with open(filename) as file:
        A = [line.replace(" ",",") for line in file]
    return A


# loads a matrix
def load_matrix(M):
    X = np.zeros((4,4))
    X[0] = [float(Fraction(M[0]).limit_denominator()), float(Fraction(M[1]).limit_denominator()), float(Fraction(M[2]).limit_denominator()), float(Fraction(M[3]).limit_denominator())]
    X[1] = [float(Fraction(M[4]).limit_denominator()), float(Fraction(M[5]).limit_denominator()), float(Fraction(M[6]).limit_denominator()), float(Fraction(M[7]).limit_denominator())] 
    X[2] = [float(Fraction(M[8]).limit_denominator()), float(Fraction(M[9]).limit_denominator()), float(Fraction(M[10]).limit_denominator()), float(Fraction(M[11]).limit_denominator())] 
    X[3] = [0, 0, 0, 1]
    return X


def id_matrix():
    return ['1','0','0','0','0','1','0','0','0','0','1','0']


# load file with all {Gamma_B < (Gamma)^Z} for a given Gamma
# returns a list of Group objects
# def get_subgroups(filename):
#     global symGroups
#     symGroups = loadSymGroups()
#     subgroups = []
#     with open(filename, 'r') as file:
#         while True:
#             try:
#                 num, index = file.readline().split()
#             except ValueError:
#                 break
#             row1 = file.readline()
#             row2 = file.readline()
#             row3 = file.readline()
#             mat = row1 + row2 + row3
#             subgroups.append(Group(num, mat.split(), index))
#             file.readline()
#     return subgroups


def get_fnames(gnum):
    fnames = listdir('s_matrices/normal_split/s_{}_matrices/'.format(gnum))
    fnames = [fname[:-4] for fname in fnames]
    fnames.sort(key=lambda x: int(x))
    return [fname + '.dat' for fname in fnames]
    


def get_subgroups(filename):
    global biebGroups
    biebGroups = loadBiebGroups()
    subgroups = []
    with open(filename, 'r') as file:
        while True:
            try:
                num, index = file.readline().split()
            except ValueError:
                break
            row1 = file.readline()
            row2 = file.readline()
            row3 = file.readline()
            mat = row1 + row2 + row3
            subgroups.append(Group(num, mat.split(), index))
            file.readline()
    return subgroups


def get_subgroup_text(filename):
    strings = []
    text = ''
    with open(filename, 'r') as file:
        while True:
            try:
                num, index = file.readline().split()
            except ValueError:
                break
            for i in range(4): text += file.readline()
            strings.append('{} {}\n{}'.format(num, index, text))
            text = ''
    return strings


def grouper(n, iterable, fillvalue=None):
    args = [iter(iterable)] * n
    return zip_longest(fillvalue=fillvalue, *args)


def split_matrices(gnum):
    with open('s_matrices/normal/s_{}_normal.dat'.format(gnum)) as f:
        for i, g in enumerate(grouper(5, f, fillvalue=''), 1):
            with open('s_matrices/normal_split/s_{}_matrices/{}.dat'.format(gnum, i), 'w') as fout:
                fout.writelines(g)




# unique matrices
def unique_matrices(filename):
    strings = get_subgroup_text(filename)
    uniq_strings = list(set(strings))
    with open(filename + '_uniq', 'w') as file:
        [file.write(uniq_string) for uniq_string in uniq_strings]


def classify_b_subgroups(gnum):
    filename = 'b_matrices/b_matrices_{}.dat'.format(gnum)
    strings = get_subgroup_text(filename)
    subgroups = get_subgroups(filename)
    bgroups = get_b_groups(subgroups)
    normal_file = open('b_matrices/normal/b_{}_normal.dat'.format(gnum), 'w')
    not_normal_file = open('b_matrices/not_normal/b_{}_not_normal.dat'.format(gnum), 'w')
    for subgroup, B, text in zip(subgroups, bgroups, strings):
        G = getFundDom(gnum, subgroup.matrix)
        if isNormal(B, G):
            normal_file.write(text)
        else:
            not_normal_file.write(text)
    normal_file.close()
    not_normal_file.close()

def classify_s_subgroups(gnum):
    filename = 's_matrices/s_matrices_{}.dat'.format(gnum)
    strings = get_subgroup_text(filename)
    subgroups = get_subgroups(filename)
    sgroups = get_s_groups(subgroups)
    normal_file = open('s_matrices/normal/s_{}_normal.dat'.format(gnum), 'w')
    not_normal_file = open('s_matrices/not_normal/s_{}_not_normal.dat'.format(gnum), 'w')
    for subgroup, S, text in zip(subgroups, sgroups, strings):
        G = getFundDom(gnum, subgroup.matrix)
        if isNormal(S, G):
            normal_file.write(text)
        else:
            not_normal_file.write(text)
    normal_file.close()
    not_normal_file.close()

# mod group element by P1
def modP1(X):
    X[:,3] = [X[0,3]%1, X[1,3]%1, X[2,3]%1, 1]
    return X

  
# determine if two group elements are equal (mod P1)
def gEq(X1, X2):
    threshold = 0.00000008
    if LA.norm(X1 - X2) < threshold:
    # if LA.norm(modP1(X1) - modP1(X2)) < threshold:
    # if LA.norm(modP1(modP1(X1) - modP1(X2))) < threshold:
        return True
    return False


# returns true if matrix is contained in 
# the group, otherwise returns false
def gContains(group, X):
    threshold = 0.00000008
    for i in range(group.shape[2]):
        # equal if norm of difference is very small
        if LA.norm(X - group[:,:,i]) < threshold:
        # if LA.norm(modP1(X) - modP1(group[:,:,i])) < threshold:
            return True
    return False


# returns true if B is a subset of the group
# used for finding B groups within groups
def gInside(group, B):
    for i in range(B.shape[2]):
        if not gContains(group, B[:,:,i]):
            return False
    return True


def loadSymGroups():
    with open('dat_files/sgroups.dat') as file:
        S = [line.strip() for line in file]
    return S


def loadBiebGroups():
    with open('dat_files/bgroups.dat') as file:
        B = [line.strip() for line in file]
    return B

# group operation: matrix mult with t mod Z3 
def mult(X1, X2):
    p = np.dot(X1, X2)
    p[:,3] = [p[0,3]%1, p[1,3]%1, p[2,3]%1, 1]
    return p


# returns the product BS of two groups consisting
# of all pairs b o s for b in B and s in S
def groupMult(B, S):
    Bsize = B.shape[2]
    Ssize = S.shape[2]
    BS = np.zeros((4, 4, Bsize*Ssize))
    count = 0
    for b in range(Bsize):
        for s in range(Ssize):
            BS[:,:,count] = mult(B[:,:,b], S[:,:,s])
            count += 1
    return BS


# returns the group inverse of the given matrix
def ginv(X):
    I = LA.inv(X)
    # t mod Z3
    I[:,3] = [I[0,3]%1, I[1,3]%1, I[2,3]%1, 1]
    return I


def affInv(X):
    return LA.inv(X)


# returns true if group is a subgroup, otherwise false
def subgroup(S):
    # check that identity is contained
    if not gContains(S, np.eye(4)):
        #print('no identity')
        return False
    # check closure under group operation
    for i in range(1, S.shape[2]):
        for j in range(1, S.shape[2]):
            if not gContains(S, mult(S[:,:,i], S[:,:,j])):
                return False
    # check closure under inverses
    for i in range(1, S.shape[2]):
        if not gContains(S, ginv(S[:,:,i])):
            return False
    return True


# determines if two groups are the same
def gEquals(G1, G2):
    # check if G1 is a subset of G2
    for i in range(G1.shape[2]):
        if not gContains(G2, G1[:,:,i]):
            return False
    # check if G2 is a subset of G1
    for i in range(G2.shape[2]):
        if not gContains(G1, G2[:,:,i]):
            return False
    return True


def cardinality(G):
    count = 0
    threshold = 0.0000008
    for i in range(G.shape[2]):
        for j in range(G.shape[2]):
            if i != j and LA.norm(G[:,:,i]-G[:,:,j]) > threshold:
                count += 1

# Check that Lagrange's Thm holds
def lagrangeCheck(G, B, S):
    return G.shape[2] == (B.shape[2] * S.shape[2])


# Determines if a group has no duplicate elements
def gUnique(G):
    threshold = 0.00000008
    for i in range(G.shape[2]):
        for j in range(G.shape[2]):
            if i != j and LA.norm(G[:,:,i]-G[:,:,j]) < threshold:
                return False
    return True


# determines if B is normal in G
def isNormal(B, G):
    for b in range(B.shape[2]):
        for g in range(G.shape[2]):
            # check if g o b o g^-1 is in B for all g in G, b in B
            if not gContains(B, mult(mult(G[:,:,g],B[:,:,b]), ginv(G[:,:,g]))):
                return False
    return True


# conjugates a group by a transformation
# X * g * X^-1  for all g in G
def conjugate1(G, X):
    Gconj = np.zeros((4, 4, G.shape[2]))
    for i in range(G.shape[2]):
        Gconj[:,:,i] = modP1(np.dot(np.dot(X, G[:,:,i]), LA.inv(X)))
    return Gconj


# conjugates a group by a transformation
# X * g * X^-1  for all g in G
def conjugate2(G, X):
    Gconj = np.zeros((4, 4, G.shape[2]))
    for i in range(G.shape[2]):
        Gconj[:,:,i] = modP1(np.dot(np.dot(X, G[:,:,i]), modP1(LA.inv(X))))
    return Gconj


# conjugates a group by a transformation
# X * g * X^-1  for all g in G
def conjugate3(G, X):
    Gconj = np.zeros((4, 4, G.shape[2]))
    for i in range(G.shape[2]):
        Gconj[:,:,i] = mult(mult(X, G[:,:,i]), ginv(X))
    return Gconj


def conjugate4(G, X):
    Gconj = np.zeros((4, 4, G.shape[2]))
    for i in range(G.shape[2]):
        Gconj[:,:,i] = modP1(mult(mult(X, G[:,:,i]), ginv(X)))
    return Gconj


# conjugates a group element by a transformation
def conj(X, Y):
    return modP1(np.dot(np.dot(Y, X), LA.inv(Y)))


# n choose k
def nCk(n, k):
    k = min(k,n-k)
    if k == 0: return 1
    numer = reduce(op.mul, xrange(n, n-k, -1))
    denom = reduce(op.mul, xrange(1, k+1))
    return numer//denom



def isSgroupElem(X):
    Z = equivTranslates(X)
    for x in Z:
        if sgroupElemTest(x):
            return True
    return False


# determines if a group element is symmorphic (ie. an element
# of a symmorphic space group) 
def sgroupElemTest(X):
    X_ = np.copy(X)
    for i in range(6):
        X_ = np.dot(X_,X)
        if gEq(X,X_):
            return True
    return False


def equivTranslates(X):
    Z = [np.copy(X) for i in range(8)]
    a = np.copy(X[0,3])
    b = np.copy(X[1,3])
    c = np.copy(X[2,3])
    # Z[0][:,3] = [a, b, c, 1]
    Z[1][:,3] = [a-1, b, c, 1]
    Z[2][:,3] = [a, b-1, c, 1]
    Z[3][:,3] = [a, b, c-1, 1]
    Z[4][:,3] = [a-1, b-1, c, 1]
    Z[5][:,3] = [a-1, b, c-1, 1]
    Z[6][:,3] = [a, b-1, c-1, 1]
    Z[7][:,3] = [a-1, b-1, c-1, 1]
    return Z


def getSgroupElems(G):
    S = []
    for i in range(1,G.shape[2]):
        if isSgroupElem(G[:,:,i]):
            S.append(i)
    return S



def isTrans(X):
    is_trans = False
    t1,t2,t3 = X[0,3], X[1,3], X[2,3]
   










        
def classifyGroup(G):
    for i in range(G.shape[2]):
        print("Coset " + str(i+1) + ":")
        print(isSgroupElem(G[:,:,i]))


##############################################################################
##############################################################################

def printGenerators(gens):
    gens = '(' + gens.replace(',', ', ')
    gens = gens.replace('\n', ')\n(')
    print(gens[:-1])


# print the matrix entries of a group
def printGroup2(G):
    for i in range(G.shape[2]):
        print(i+1)
        print(G[:,:,i])


# prints the Cosets of a group
def printGroup(G, f):
    for i in range(G.shape[2]):
        print(printMatrix(G, i), file=f)


# helper function for printGroup
def printMatrix(G, matNum):
    x = printRow(G, 0, matNum)
    y = printRow(G, 1, matNum)
    z = printRow(G, 2, matNum)
    return x + ", " + y + ", " + z


# helper function for printGroup
def printRow(G, rowNum, matNum):
    xstr = printStr(G, 0, rowNum, matNum, "x")
    ystr = printStr(G, 1, rowNum, matNum, "y")
    zstr = printStr(G, 2, rowNum, matNum, "z")
    tstr = printStr(G, 3, rowNum, matNum, "")
    if ystr != "" and ystr[0] != '-':
        ystr = '+' + ystr;
    if zstr != "" and zstr[0] != '-':
        zstr = '+' + zstr;
    if tstr != "" and tstr[0] != '-':
        tstr = '+' + tstr;
    string = xstr + ystr + zstr + tstr
    if string[0] == '+':
        string = string[1:]
    return string


# helper function for printGroup
def printStr(G, colNum, rowNum, matNum, coord):
    if G[rowNum, colNum, matNum] == 0:
        string = ""
    elif G[rowNum, colNum, matNum] == 1:
        string = coord
    elif G[rowNum, colNum, matNum] == -1 and coord != "":
        string = "-" + coord
    else:
        string = str(Fraction(G[rowNum,colNum,matNum]).limit_denominator()) + coord
    return string


def printMat(Z, f):
    print('transformation matrix:', file=f)
    print(Z[0] + "    " + Z[1] + "    " + Z[2] + "    " + Z[3], file=f)
    print(Z[4] + "    " + Z[5] + "    " + Z[6] + "    " + Z[7], file=f)
    print(Z[8] + "    " + Z[9] + "    " + Z[10] + "    " + Z[11], file=f)


###########################################################################
################## Interface with Bilbao server ###########################
###########################################################################


def getCosets(sup, sub, filename, M):
    data = {'super':sup, 'sub':sub, 'x1':M[0], 'x2':M[1], 'x3':M[2], 'x4':M[3], 'y1':M[4], 'y2':M[5], 'y3':M[6], 'y4':M[7], 'z1':M[8], 'z2':M[9], 'z3':M[10], 'z4':M[11], 'what':'left'}
    r = requests.post('http://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-cosets', data)
    start = r.text.index('Coset 1') 
    end = r.text.index('<', start)
    res = r.text[start:end]
    res = re.sub(r"Coset \d+:\n\n",'',res)
    res = re.sub(r"\(|\)",'',re.sub(r"\(\n",'',res))
    res = res.replace(',',' ')
    res = res.replace('\n\n','\n')
    with open(filename,'w') as fid:
        fid.write(res)

def identifyGroup(generators, fname, f_err):
    data = {'tipog':'gesp', 'generators':generators}
    t = requests.post('http://www.cryst.ehu.es/cgi-bin/cryst/programs/checkgr.pl', data).text
    start = t.index('grupo') + 14
    end = t.index('\n', start) - 2
    M = t[start:end].split(',')
    num = M[1]
    start = t.index('pre') + 4
    end = t.index('<', start)
    try:
        bmatrix = LA.inv(load_matrix(t[start:end].split()))
    except ValueError:
        print(t[start:end].split(), f_err)
        print(generators, f_err)
        print('group:' + sys.argv[1], f_err)
        print('file: ' + fname, f_err)
        return 0, None, None
    bmatrix = [str(Fraction(b).limit_denominator()) for b in bmatrix[0:3,:].flatten()]
    B = loadBiebGroups()
    if num in B: 
        return num, bmatrix, True
    return num, bmatrix, False


###########################################################################
###########################################################################


def removeNonGroups():
    ind2remove = []
    for i in range(len(Bgroups)):
        if not subgroup(Bgroups[i]): ind2remove.append(i)
    for i in sorted(ind2remove, reverse=True):
        del Bgroups[i]
        del repCombs[i]


def bgroup_test(gnum, snum, mat, k):
    getCosets(gnum, snum, 'dat_files/g_s_cosets', mat)
    with open('dat_files/g_s_cosets') as file:
        A = [line.split() for line in file]
    A = [A[i:i+k] for i in range(0,len(A),k)]
    A.pop(0)
    return A


def getBgroups(gnum, snum, mat, k):
    getCosets(gnum, snum, 'dat_files/g_s_cosets', mat)
    with open('dat_files/g_s_cosets') as file:
        A = [line.split() for line in file]
    A = [A[i:i+k] for i in range(0,len(A),k)]
    A.pop(0)
    # symElems = getSymElems(A)
    biebElems = getBiebElems(A)
    del matCombs[:]
    del repCombs[:]
    global Bgroups
    del Bgroups[:]
    combineElems(biebElems[0], biebElems[1], [], '')
    Bgroups = [np.insert(np.dstack(grp),0,np.eye(4,4),axis=2) for grp in matCombs]
    # removeNonGroups()


def get_b_groups(subgroups):
    bnum = '0'
    bgroups = []
    for gamma_B in subgroups:
        if bnum != gamma_B.num:
            bnum = gamma_B.num
            B = loadGroup('dat_files/bieberbach_cosets/b{}.dat'.format(bnum))
            bgroups.append(B)
        else:
            bgroups.append(B)
    return bgroups


def get_s_groups(subgroups):
    snum = '0'
    sgroups = []
    for gamma_S in subgroups:
        if snum != gamma_S.num:
            snum = gamma_S.num
            S = loadGroup('dat_files/symmorphic_cosets/s{}.dat'.format(snum))
            sgroups.append(S)
        else:
            sgroups.append(S)
    return sgroups
 

def getFundDom(gnum, mat):
    getCosets(gnum, '1', 'dat_files/coset_file', mat)
    return loadGroup('dat_files/coset_file')


def decomp(gnum, subgroups, sgroups, fname, f_err):
    f = open(fname, 'w')
    decomps = []
    for subgroup, S in zip(subgroups, sgroups):
        snum = subgroup.num
        Z = subgroup.matrix
        G = getFundDom(gnum, Z)
        # getSgroups(gnum, bnum, Z, B.shape[2])
        getBgroups(gnum, snum, Z, S.shape[2])
        # S = semidirectTest(G, B, bnum, Z, subgroup.index)
        B = semidirectTest(G, S, snum, Z, subgroup.index, f, fname, f_err)
        # if S != []: decomps.append(Decomp(gnum, bnum, Z, S))
        if B != []: decomps.append(Decomp(gnum, snum, Z, B))
    f.close()
    # TableEntry(decomps)
    return decomps

def semidirectTest(G, S, snum, Z, index, f, fname, f_err):
    bgrps = []
    for i in range(len(Bgroups)):
        B = Bgroups[i]
        if gInside(G, S) and gInside(G, B) and gEquals(G, groupMult(B, S)):
            generators = 'x,y,z\n' + repCombs[i]
            bnum, bmatrix, is_bieberbach = identifyGroup(generators, fname, f_err)
            if bnum == 0: 
                continue
            printMat(Z, f)
            print('Gamma: {}'.format(sys.argv[1]), file=f)
            printGroup(G, f)
            print('Gamma_S: {} normal, index: {}'.format(snum, index), file=f)
            printGroup(S, f)
            # S is normal
            if is_bieberbach:
                if isNormal(B, G):
                    print('Gamma_B: {} normal'.format(bnum), file=f)
                    bgrps.append(Subgroup(bnum, True, True))
                else:
                    print('Gamma_B: {} not normal'.format(bnum), file=f)
                    bgrps.append(Subgroup(bnum, False, True))
                printMat(bmatrix, f)
                print(generators.replace(',', ', '), file=f)
            else:
                if isNormal(B, G):
                    print('Gamma/Gamma_S: {} normal'.format(bnum), file=f)
                    bgrps.append(Subgroup(bnum, True, False))
                else:
                    print('Gamma/Gamma_S: {} not normal'.format(bnum), file=f)
                    bgrps.append(Subgroup(bnum, False, False))
                printMat(bmatrix, f)
                print(generators.replace(',', ', '), file=f)
    return list(set(bgrps))



# def semidirectTest(G, S, snum, Z, index):
#     sgrps = []
#     for i in range(len(Sgroups)):
#         S = Sgroups[i]
#         if gInside(G, S) and gInside(G, B) and gEquals(G, groupMult(B, S)):
#             generators = 'x,y,z\n' + repCombs[i]
#             snum, smatrix, is_symmorphic = identifyGroup(generators)
#             printMat(Z)
#             print('Gamma: {}'.format(sys.argv[1]))
#             printGroup(G)
#             # B is normal
#             print('Gamma_B: {} normal, index: {}'.format(bnum, index))
#             printGroup(B)
#             if is_symmorphic:
#                 if isNormal(S, G):
#                     print('Gamma_S: {} normal'.format(snum))
#                     sgrps.append(Subgroup(snum, True, True))
#                 else:
#                     print('Gamma_S: {} not normal'.format(snum))
#                     sgrps.append(Subgroup(snum, False, True))
#                 printMat(smatrix)
#                 print(generators.replace(',', ', '))
#             else:
#                 if isNormal(S, G):
#                     print('Gamma/Gamma_B: {} normal'.format(snum))
#                     sgrps.append(Subgroup(snum, True, False))
#                 else:
#                     print('Gamma/Gamma_B: {} not normal'.format(snum))
#                     sgrps.append(Subgroup(snum, False, False))
#                 printMat(smatrix)
#                 print(generators.replace(',', ', '))
#     return list(set(sgrps))





