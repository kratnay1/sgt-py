""" Group theory modules for testing with space groups """

import sys
import re
import math
import os
import itertools
import requests
import asyncio
import async_timeout
import aiohttp
import aiofiles
from bs4 import BeautifulSoup
import operator as op
import numpy as np
from functools import reduce
from numpy import linalg as LA
from fractions import Fraction


class Group:
    def __init__(self, num, matrix, index=None, supergroup=None, lin_rep=None):
        # num: ITA number
        # matrix: transformation matrix
        self.num = num
        self.matrix = matrix
        self.index = index
        self.supergroup = supergroup
        self.lin_rep = lin_rep

class Subgroup:
    def __init__(self, num, normal=True, symmorphic=True):
        # num: ITA number
        # normal: boolean for normality
        self.num = num
        self.normal = normal
        self.symmorphic = symmorphic
    def __eq__(self, other):
        return isinstance(other, Subgroup) \
               and self.num == other.num \
               and self.normal == other.normal
    def __hash__(self):
        return hash((self.num, self.normal))

class Decomp:
    def __init__(self, gnum, bnum, matrix, sgroups):
        self.gnum = gnum
        self.bnum = bnum
        self.matrix = matrix
        self.sgroups = sgroups


class TableEntry:
    def __init__(self, decomps):
        self.decomps = decomps
        self.names = self.getNames()
        self.printEntry()
    def getNames(self):
        with open('tex/sgroup_labels.tex') as file:
            return [line.strip() for line in file]
    def printEntry(self):
        with open('tex/table_entries/entry_{}'.format(sys.argv[1]), 'w') as file:
            for decomp in self.decomps:
                self.printDecomp(decomp, file)
    def printDecomp(self, d, file):
        gname = self.names[int(d.gnum)-1]
        bname = self.names[int(d.bnum)-1]
        snames = self.getSNames(d.sgroups)
        normalS = self.getNormalS(d.sgroups)
        symmorphic = self.getSymmorphic(d.sgroups)
        file.write('${}$ & ${}$ & ${}$ & ${}$ & {} & {} & {} \\\\\n'.format(d.gnum, gname, bname, snames, symmorphic, 'Y', normalS))
    def getSNames(self, sgroups):
        snames = self.names[int(sgroups[0].num)-1]
        if not sgroups[0].symmorphic: snames = '{{\color{{red}}{}}}'.format(snames)
        for i in range(1,len(sgroups)):
            name = self.names[int(sgroups[i].num)-1]
            if not sgroups[i].symmorphic: name = '{{\color{{red}}{}}}'.format(name)
            snames += ',' + name
        return snames
    def getNormalS(self, sgroups):
        if sgroups[0].normal: normalS = 'Y'
        else: normalS = 'N'
        for i in range(1,len(sgroups)):
            if sgroups[i].normal: normalS += ',' + 'Y'
            else: normalS += ',' + 'N'
        return normalS
    def getSymmorphic(self, sgroups):
        if sgroups[0].symmorphic: symmorphic = 'Y'
        else: symmorphic = '{{\color{{red}}N}}'
        for i in range(1, len(sgroups)):
            if sgroups[i].symmorphic: symmorphic += ',' + 'Y'
            else: symmorphic += ',' + '{{\color{{red}}N}}'
        return symmorphic


class LinRep:
    """A wrapper class for a linear representation that allows you to use the multiplication operator (*) to multiply two groups."""

    def __init__(self, lin_rep, cosets):
        self.lin_rep = lin_rep
        self.cosets = cosets

    def __str__(self):
        return self.cosets

    def __mul__(self, other):
        new_lin_rep = groupMult(self.lin_rep, other.lin_rep)
        new_cosets = cosets_from_lin_rep(new_lin_rep)
        return LinRep(new_lin_rep, new_cosets)


class SpaceGroup:
    """A class to represent a space group by both a linear representation of :math:`\\frac{\\Gamma}{P1}` and the general positions (i.e. the actions on :math:`\\mathbb{R}^3`).

        :param num: The ITA number of the space group.
        :param lin_rep: A numpy.ndarray containing the homogenous transformation matrices of the group.
        :param cosets: A string containing the general positions.
        :param matrix: An affine transformation matrix used to conjugate the space group.
        """

    def __init__(self, num, lin_rep, cosets, matrix=None):
        """Constructor"""

        self.num = num
        self.lin_rep = LinRep(lin_rep, cosets[:-1])
        if not matrix:
            self.matrix = id_matrix()
        else:
            self.matrix = matrix

    def write_to_file(self, filename=None):
        """Writes the general positions to a specified file.

        :param filename: A string containing the file name.
        """
        # TODO
        return None

class SpaceGroupPair:
    """A class to store a supgergroup-subgroup pair of space groups.

        :param subgroup:
        :param supergroup:
        :param matrix:
        :param index:
        """
    def __init__(self, subgroup, supergroup, matrix, index):
        """Constructor"""

        self.subgroup = subgroup
        self.supergroup = supergroup
        self.matrix = matrix
        self.index = index

    def is_normal(self):
        """is_normal docstring"""
        return None
    

def get_space_group(gnum, matrix=None):
    """Returns a SpaceGroup object given the ITA number and transformation matrix that relates the group and subgroup (:math:`P1`) basis.

    :param gnum: The group ITA number.
    :param matrix: The transformation matrix.
    """
    return SpaceGroup(gnum, *cosetDecomp(gnum, '1', matrix), matrix)


def get_space_subgroups(supergroup, subgroup, index=None):
    """Returns a list of SpaceGroupPair objects given a space group ITA number and a space subgroup ITA number for a given index or for all possible indices if no index is provided (stopping right before the index of the lattice translation group in the supergroup).

    :param supergroup: The supergroup ITA number.
    :param subgroup: The subgroup ITA number.
    :param index: The index of the subgroup in the supergroup.
    """
    urls, indices = get_chain_urls(supergroup, subgroup, index)
    matrix_urls = get_matrix_urls(supergroup, urls)
    if matrix_urls == []:
        if index:
            raise ValueError('Bilbao Error: The group {} can not be related to the subgroup {} with index {}.'.format(supergroup, subgroup, index))
        raise ValueError('Bilbao Error: This is not possible.')
    get_matrices(supergroup, matrix_urls)
    return get_subgroups('matrices/matrices_{}.dat'.format(supergroup), supergroup)


def get_chain_urls(supergroup, subgroup, index=None):
    if index:
        return ['http://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-subggraph?super={}&sub={}&index={}'.format(supergroup, subgroup, index)], [index]
    urls, indices = [], []
    with open('dat_files/ind_T_gamma.dat', 'r') as file:
        ind_T = int([line.strip() for i,line in enumerate(file) if i==int(supergroup)-1][0])
    index = get_start_index(subgroup, ind_T)
    while index < ind_T:
        if ind_T % index:
            index += 1
            continue
        urls.append('http://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-subggraph?super={}&sub={}&index={}'.format(supergroup, subgroup, index))
        indices.append(index)
        index += 1
    return urls, indices

def get_start_index(subgroup, ind_T):
    S = standardGenPos(subgroup)
    start_index = ind_T / S.shape[2]
    if start_index == 1: return 2
    if int(start_index) == 0: return 2
    return int(start_index)

def ensure_all_futures(supergroup, urls, loop):
    with aiohttp.ClientSession(loop=loop) as session:
        return [asyncio.ensure_future(download_coroutine(supergroup, session, loop, url)) for url in urls]

async def download_coroutine(supergroup, session, loop, url):
    """download_coroutine test

    :param supergroup:
    :param session:
    :param loop:
    :param url:
    """
    with aiohttp.ClientSession(loop=loop) as session:
        with async_timeout.timeout(None):
            async with session.get(url) as response:
                text = await response.text()
                soup = BeautifulSoup(text, 'html5lib')
                matrix_links = soup.find_all('a')[4:-3]
                if matrix_links:
                    async with aiofiles.open('subgroup_chains/matrix_urls_{}.dat'.format(supergroup), 'a') as fd:
                        # [await fd.write('http://www.cryst.ehu.es/cgi-bin/cryst/programs/' + link['href'] +'\n') for link in matrix_links]
                        for link in matrix_links:
                            await fd.write('http://www.cryst.ehu.es/cgi-bin/cryst/programs/' + link['href'] +'\n')

def get_matrix_urls(supergroup, urls):
    dirname = 'subgroup_chains'
    filename = 'subgroup_chains/matrix_urls_{}.dat'.format(supergroup)
    if os.path.isdir(dirname):
        for fname in os.listdir(dirname):
            file_path = os.path.join(dirname, fname)
            if os.path.isfile(file_path):
                os.unlink(file_path)
        os.removedirs(dirname)
    os.makedirs(dirname)
    open(filename, 'w').close()
    loop = asyncio.get_event_loop()
    loop.run_until_complete(asyncio.wait(ensure_all_futures(supergroup, urls, loop)))
    with open(filename, 'r') as file:
        matrix_urls = [line.strip() for line in file]
    os.remove(filename)
    os.removedirs(dirname)
    return matrix_urls

def ensure_matrix_futures(supergroup, loop, matrix_urls):
    with aiohttp.ClientSession(loop=loop) as session:
        return [asyncio.ensure_future(download_matrices(supergroup, session, loop, url)) for url in matrix_urls]

async def download_matrices(supergroup, session, loop, url):
    with async_timeout.timeout(None):
        async with aiohttp.ClientSession(loop=loop) as session:
            async with session.get(url) as response:
                text = await response.text()
                matrices = format_matrix_text(text)
                async with aiofiles.open('matrices/matrices_{}.dat'.format(supergroup), 'a') as fd:
                    await fd.write(matrices)

def format_matrix_text(text):
    soup = BeautifulSoup(text, 'html5lib')
    text = soup.find('center').text.replace('\n\n', '\n')
    text = re.sub('\n', lambda m,c=itertools.count():m.group() if next(c)%3 else '\n\n', text)
    chain = text[:text.find('\n\n')]
    indices = chain[chain.find('[')+1: -1]
    index = str(int(reduce((lambda x,y: int(x)*int(y)), indices.split())))
    snum = str(int(chain[:chain.find('[')].split()[-1]))
    text = text.replace('\n\n', '\n\n{} {}\n'.format(snum, index))
    return text[text.find('\n')+2: text.rfind('\n\n')+2]

def get_matrices(supergroup, matrix_urls):
    filename = 'matrices/matrices_{}.dat'.format(supergroup)
    open(filename, 'w').close()
    loop = asyncio.get_event_loop()
    loop.run_until_complete(asyncio.wait(ensure_matrix_futures(supergroup, loop, matrix_urls)))


# helper function for loadGroup
# sets the chosen column of the row depending on
# the string containing the coefficient of x,y,z
def setRowVal(row, ind, a): 
    if a == "":
        row[ind] = 1
    elif a == '-':
        row[ind] = -1
    else:
        row[ind] = float(Fraction(a))


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
        row[3] = float(Fraction(a))
    return row



def cosetDecomp(supergroup, subgroup, mat=None, filename=None):
    cosets = getCosets(supergroup, subgroup, mat, filename)
    return loadGroup(filename, cosets), cosets

def standardGenPos(group, mat=None, filename=None):
    cosets = getCosets(group, 1, mat, filename)
    return loadGroup(filename, cosets)


# loads a group as a 3d array of 4x4 matrices

def loadGroup(filename=None, coset_text=None):
    """Loads an n-element finite group as an numpy.ndarray with shape (4,4,n).  Each group element is represented as a 4x4 homogenous matrix of the form 
        :math:`{\\cal H}(A, {\\bf a}) = \\left(\\begin{array}{ccc}
        A && {\\bf a} \\\ \\\
                {\\bf 0}^t && 1 \\end{array}\\right)` where :math:`A \\in GL(3, \\mathbb{R})` and :math:`{\\bf a} \\in \\mathbb{R}^3.`

    :param filename: test
    """
    if filename:
        with open(filename) as file:
            A = [line.split() for line in file]
    if coset_text:
        A = [coset.split() for coset in coset_text.split('\n')[:-1]]
    G = np.zeros((4, 4, len(A)))
    for i in range(len(A)):
        G[0,:,i] = getRow(A[i][0])
        G[1,:,i] = getRow(A[i][1])
        G[2,:,i] = getRow(A[i][2])
        G[:,3,i] = G[:,3,i] % 1
        G[3,:,i] = [0, 0, 0, 1]
    return G


def loadCosetRep(rep):
    """Loads a coset representative as a 4x4 homogenous matrix stored as a numpy.ndarray.

    :param rep: A string containing the action on.
    """
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


def getSymElems(A):
    symElems = [[],[]]
    for reps in A:
        matrices = []
        gens = []
        for rep in reps:
            g = loadCosetRep(rep)
            if isSgroupElem(g):
                matrices.append([g])
                gens.append(getGenerator(rep))
        symElems[0].append(matrices)
        symElems[1].append(gens)
    return symElems


symGroups = []
biebGroups = []
matCombs = []
repCombs = []
Bgroups = []
Sgroups = []


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
    X[0] = [float(Fraction(M[0])), float(Fraction(M[1])), float(Fraction(M[2])), float(Fraction(M[3]))]
    X[1] = [float(Fraction(M[4])), float(Fraction(M[5])), float(Fraction(M[6])), float(Fraction(M[7]))] 
    X[2] = [float(Fraction(M[8])), float(Fraction(M[9])), float(Fraction(M[10])), float(Fraction(M[11]))] 
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


def get_subgroups(filename, supergroup_num):
    # global biebGroups
    # biebGroups = loadBiebGroups()
    space_group_pairs = []
    with open(filename, 'r') as file:
        while True:
            try:
                subgroup_num, index = file.readline().split()
            except ValueError:
                break
            row1 = file.readline()
            row2 = file.readline()
            row3 = file.readline()
            matrix = row1 + row2 + row3
            matrix = matrix.split()
            subgroup = SpaceGroup(subgroup_num, *cosetDecomp(subgroup_num, '1'))
            supergroup = SpaceGroup(supergroup_num, *cosetDecomp(supergroup_num, '1', matrix), matrix)
            space_group_pairs.append(SpaceGroupPair(subgroup, supergroup, matrix, index))
            file.readline()
    return space_group_pairs


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
    filename = 's_matrices/s_matrices_{}.dat_uniq'.format(gnum)
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
    # if LA.norm(X1 - X2) < threshold:
    if LA.norm(modP1(X1) - modP1(X2)) < threshold:
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
    BS = np.unique(BS, axis=2)
    for i in range(BS.shape[2]):
        if np.allclose(np.eye(4), BS[:,:,i]):
            res = np.delete(BS, i, axis=2)
    return np.insert(res, 0, np.eye(4), axis=2)



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
    return count

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


def cosets_from_lin_rep(G):
    cosets = ""
    for i in range(G.shape[2]):
        cosets += printMatrix(G, i) + '\n'
    return cosets[:-1]


# prints the Cosets of a group
def printGroup(G):
    for i in range(G.shape[2]):
        print(printMatrix(G, i))


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


def printMat(Z):
    print('transformation matrix:')
    print(Z[0] + "    " + Z[1] + "    " + Z[2] + "    " + Z[3])
    print(Z[4] + "    " + Z[5] + "    " + Z[6] + "    " + Z[7])
    print(Z[8] + "    " + Z[9] + "    " + Z[10] + "    " + Z[11])


###########################################################################
################## Interface with Bilbao server ###########################
###########################################################################


def getCosets(sup, sub, M=None, filename=None):
    if not M:
        M = ['1','0','0','0','0','1','0','0','0','0','1','0']  
    data = {'super':sup, 'sub':sub, 'x1':M[0], 'x2':M[1], 'x3':M[2], 'x4':M[3], 'y1':M[4], 'y2':M[5], 'y3':M[6], 'y4':M[7], 'z1':M[8], 'z2':M[9], 'z3':M[10], 'z4':M[11], 'what':'left'}
    r = requests.post('http://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-cosets', data)
    try:
        start = r.text.index('Coset 1') 
        end = r.text.index('<', start)
        res = r.text[start:end]
        res = re.sub(r"Coset \d+:\n\n",'',res)
        res = re.sub(r"\(|\)",'',re.sub(r"\(\n",'',res))
        res = res.replace(',',' ')
        res = res.replace('\n\n','\n')
        if filename:
            with open(filename,'w') as fid:
                fid.write(res)
        else:
            return res
    except ValueError as e:
        print("Bilbao Error: Sorry, no decomposition found. Please, check your transformation matrix!")

def identifyGroup(generators):
    data = {'tipog':'gesp', 'generators':generators}
    t = requests.post('http://www.cryst.ehu.es/cgi-bin/cryst/programs/checkgr.pl', data).text
    start = t.index('grupo') + 14
    end = t.index('\n', start) - 2
    M = t[start:end].split(',')
    num = M[1]
    start = t.index('pre') + 4
    end = t.index('<', start)
    smatrix = LA.inv(load_matrix(t[start:end].split()))
    smatrix = [str(Fraction(s)) for s in smatrix[0:3,:].flatten()]
    S = loadSymGroups()
    if num in S: 
        return num, smatrix, True
    return num, smatrix, False


###########################################################################
###########################################################################


def removeNonGroups():
    ind2remove = []
    for i in range(len(Sgroups)):
        if not subgroup(Sgroups[i]): ind2remove.append(i)
    for i in sorted(ind2remove, reverse=True):
        del Sgroups[i]
        del repCombs[i]


def getSgroups(gnum, bnum, mat, k):
    getCosets(gnum, bnum, 'dat_files/g_b_cosets', mat)
    with open('dat_files/g_b_cosets') as file:
        A = [line.split() for line in file]
    A = [A[i:i+k] for i in range(0,len(A),k)]
    A.pop(0)
    symElems = getSymElems(A)
    del matCombs[:]
    del repCombs[:]
    global Sgroups
    del Sgroups[:]
    combineElems(symElems[0], symElems[1], [], '')
    Sgroups = [np.insert(np.dstack(grp),0,np.eye(4,4),axis=2) for grp in matCombs]
    removeNonGroups()


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
    removeNonGroups()


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
    getCosets(gnum, '1', mat, 'dat_files/coset_file')
    return loadGroup('dat_files/coset_file')


def decomp(gnum, subgroups, bgroups):
    decomps = []
    for subgroup, B in zip(subgroups, bgroups):
        bnum = subgroup.num
        Z = subgroup.matrix
        G = getFundDom(gnum, Z)
        getSgroups(gnum, bnum, Z, B.shape[2])
        S = semidirectTest(G, B, bnum, Z, subgroup.index)
        if S != []: decomps.append(Decomp(gnum, bnum, Z, S))
    TableEntry(decomps)
    return decomps


def semidirectTest(G, B, bnum, Z, index):
    sgrps = []
    for i in range(len(Sgroups)):
        S = Sgroups[i]
        if gInside(G, B) and gInside(G, S) and gEquals(G, groupMult(B, S)):
            generators = 'x,y,z\n' + repCombs[i]
            snum, smatrix, is_symmorphic = identifyGroup(generators)
            printMat(Z)
            print('Gamma: {}'.format(sys.argv[1]))
            printGroup(G)
            # B is normal
            print('Gamma_B: {} normal, index: {}'.format(bnum, index))
            printGroup(B)
            if is_symmorphic:
                if isNormal(S, G):
                    print('Gamma_S: {} normal'.format(snum))
                    sgrps.append(Subgroup(snum, True, True))
                else:
                    print('Gamma_S: {} not normal'.format(snum))
                    sgrps.append(Subgroup(snum, False, True))
                printMat(smatrix)
                print(generators.replace(',', ', '))
            else:
                if isNormal(S, G):
                    print('Gamma/Gamma_B: {} normal'.format(snum))
                    sgrps.append(Subgroup(snum, True, False))
                else:
                    print('Gamma/Gamma_B: {} not normal'.format(snum))
                    sgrps.append(Subgroup(snum, False, False))
                printMat(smatrix)
                print(generators.replace(',', ', '))
    return list(set(sgrps))


def SymSemidirectTest(G, S, snum, Z, index):
    sgrps = []
    for i in range(len(Sgroups)):
        S = Sgroups[i]
        if gInside(G, S) and gInside(G, B) and gEquals(G, groupMult(B, S)):
            generators = 'x,y,z\n' + repCombs[i]
            snum, smatrix, is_symmorphic = identifyGroup(generators)
            printMat(Z)
            print('Gamma: {}'.format(sys.argv[1]))
            printGroup(G)
            # B is normal
            print('Gamma_B: {} normal, index: {}'.format(bnum, index))
            printGroup(B)
            if is_symmorphic:
                if isNormal(S, G):
                    print('Gamma_S: {} normal'.format(snum))
                    sgrps.append(Subgroup(snum, True, True))
                else:
                    print('Gamma_S: {} not normal'.format(snum))
                    sgrps.append(Subgroup(snum, False, True))
                printMat(smatrix)
                print(generators.replace(',', ', '))
            else:
                if isNormal(S, G):
                    print('Gamma/Gamma_B: {} normal'.format(snum))
                    sgrps.append(Subgroup(snum, True, False))
                else:
                    print('Gamma/Gamma_B: {} not normal'.format(snum))
                    sgrps.append(Subgroup(snum, False, False))
                printMat(smatrix)
                print(generators.replace(',', ', '))
    return list(set(sgrps))




