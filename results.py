""" NOT FOR USE.

This code was used to process results and generate tables summarizing the decompositions (Chapter 4 of my Master's thesis) as well as a write-up detailing the results (Chapter 5)"""

import group_modules_sym as gm
import itertools
import re
import sys
import os
import asyncio
import async_timeout
import aiohttp
import aiofiles
from bs4 import BeautifulSoup
from functools import reduce
from collections import defaultdict
import numpy as np
from numpy import linalg as LA
from fractions import Fraction


class redirect_output(object):
    """context manager for redirecting stdout/err to files"""
    def __init__(self, stdout='', stderr=''):
        self.stdout = stdout
        self.stderr = stderr

    def __enter__(self):
        self.sys_stdout = sys.stdout
        self.sys_stderr = sys.stderr

        if self.stdout:
            sys.stdout = open(self.stdout, 'w')
        if self.stderr:
            if self.stderr == self.stdout:
                sys.stderr = sys.stdout
            else:
                sys.stderr = open(self.stderr, 'w')

    def __exit__(self, exc_type, exc_value, traceback):
        sys.stdout = self.sys_stdout


class DecompText:
    def __init__(self, matrix, gnum, g_cosets, bnum, b_cosets, snum, symmorphic, normalS, s_matrix, s_cosets):
        self.fd_labels = self.get_fd_labels()
        self.subgroup_labels = [self.get_label(gnum), self.get_label(bnum), self.get_label(snum)]
        self.cosets_list = [self.format_cosets(g_cosets, 0), self.format_cosets(b_cosets, 1), self.format_cosets(s_cosets, 2)]
        self.symmorphic = symmorphic
        self.normalS = normalS
        self.subgroups = [self.format_subgroup(1, matrix), self.format_subgroup(2, s_matrix)]
        self.decomp_text = self.write_decomp()
    def get_label(self, num):
        with open('tex/sgroup_labels.tex', 'r') as file:
            return [line.strip() for i, line in enumerate(file) if i==int(num)-1][0]
    def format_matrix(self, M):
        return r'\begin{pmatrix} ' + M[0] + r' & ' + M[1] + r' & ' + M[2] + r' & ' + M[3] + r' \\ ' + M[4] + r' & ' + M[5] + r' & ' + M[6] + r' & ' + M[7] + r' \\ ' + M[8] + r' & ' + M[9] + r' & ' + M[10] + r' & ' + M[11] + r' \end{pmatrix}'
    def format_subgroup(self, index, matrix):
        if index == 1: 
            sub, mat_label, normal  = r'_B', r'\alpha', r'\ \text{(normal)}'
        if index == 2: 
            mat_label = r'Y'
            if self.symmorphic:
                sub = r'_S'
            else:
                sub = r'/\Gamma_B'
            if self.normalS:
                normal = r'\ \text{(normal)}'
            else:
                normal = ''
        return r'\[ \Gamma' + sub + r' = ' + self.subgroup_labels[index] + normal +  r', \ \ \  ' + mat_label + r' = ' + self.format_matrix(matrix) +  r',\]'
    def get_fd_labels(self):
        return [r'F_{\frac{\Gamma^Z}{P1}}',r'F_{\frac{\Gamma_{B}^{X}}{P1}}',r'F_{\frac{\Gamma^{Y}}{P1}}']
    def format_cosets(self, cosets, fd_index):
        reps = cosets.split('\n')
        reps = list(filter(None, reps))
        reps = ['({})'.format(rep) for rep in reps]
        reps = '; '.join(reps)
        # return r'{\centering $' + self.fd_labels[fd_index] + r' = \{' + reps + r' $\} \par}'
        # return r' ' + self.fd_labels[fd_index] + r' = \{' + reps + r' \}'
        return self.insert_newlines(r' ' + self.fd_labels[fd_index] + r' = \{' + reps + r' \}')
    def insert(self, source_str, insert_str, pos):
        return source_str[:pos]+insert_str+source_str[pos:]
    def insert_newlines(self, text):
        #insert a r'\\' after every 74 chars
        s = text.split('=')
        c = s[-1]
        len_text = len(c) + 7
        if len_text <= 68: return text
        num = int(len_text / 68)
        if num == 1: return s[0] + r' = ' + self.insert(c, r'\\ ', 68)
        # many lines
        pos = 68
        for i in range(num):
            c = self.insert(c, r'\\ ', pos)
            pos += 68
        return s[0] + r' = ' + c
    def write_decomp(self):
        return r'Given the subgroup ' + self.subgroups[0] + r' and setting $Z=\alpha$ and $X=\mathbb{I}$, we find a complementary subgroup ' + self.subgroups[1] + r' resulting in the following decomposition: \begin{gather*}' + self.cosets_list[0] + r'\\ ' + self.cosets_list[1] + r'\\ ' + self.cosets_list[2] + r'. \end{gather*}'


class WriteUp:
    def __init__(self, gnum):
        self.gnum = gnum
        self.g_label = self.get_label()
        self.header = r'\subsection{{$\bf {}$, $\bf {}$}}'.format(gnum, self.g_label)
        # self.subheader = r'\subsubsection{{$\bf {}$, $\bf {}$}}'.format(gnum, self.g_label)
        self.subheader = ''
        self.sym_groups = self.get_sym_groups()
        self.decomps = []
        self.get_decomps()
        self.write_to_file()
    def get_label(self):
        with open('tex/sgroup_labels.tex', 'r') as file:
            return [line.strip() for i, line in enumerate(file) if i==int(self.gnum)-1][0]
    def get_sym_groups(self):
        with open('dat_files/sgroups.dat') as file:
            return [line.strip() for line in file]
    def get_decomps(self):
        filename = 'decomp/one_semi/decomp_{}'.format(self.gnum)
        with open(filename, 'r') as file:
            while True:
                line = file.readline() # transformation matrix: 
                if not line: break
                matrix = ''
                for i in range(3): matrix += file.readline()
                file.readline() # Gamma
                line = ''
                g_cosets = ''
                while True:
                    line = file.readline()
                    if line.startswith('Gamma'):
                        break
                    g_cosets += line
                bnum = line.split()[1]
                b_cosets = ''
                while True:
                    line = file.readline()
                    if line.startswith('Gamma'):
                        break
                    b_cosets += line
                line = line.split()
                snum, normalS = line[1], True
                if snum in self.sym_groups: 
                    symmorphic = True
                else:
                    symmorphic = False
                if line[2] == "not": normalS = False
                file.readline() # transformation matrix:
                s_matrix = ''
                for i in range(3): s_matrix += file.readline()
                s_cosets = ''
                while True:
                    line = file.readline()
                    if not line.strip():
                        break
                    s_cosets += line
                self.decomps.append(DecompText(matrix.split(), self.gnum, g_cosets, bnum, b_cosets, snum, symmorphic, normalS, s_matrix.split(), s_cosets))
    def write_to_file(self):
        with open('tex/decomps/one_semidirect/decomp_{}.tex'.format(self.gnum), 'w') as file:
            file.write(self.header + '\n')
            for decomp in self.decomps:
                file.write(self.subheader + '\n')
                file.write(decomp.decomp_text + '\n')


class BlankTableEntry:
    def __init__(self, gnum, type):
        self.type = type
        self.gnum = gnum
        self.label = self.get_label()
        self.print_entry()
    def get_label(self):
        with open('tex/sgroup_labels.tex', 'r') as file:
            return [line.strip() for i, line in enumerate(file) if i==int(self.gnum)-1][0]
    def print_entry(self):
        if self.type == 0:
            filename = 'tex/table_entries/semidirect/entry_{}'.format(self.gnum)
        else:
            filename = 'tex/table_entries/non_semidirect/entry_{}'.format(self.gnum)
        with open(filename, 'w') as file:
            file.write('${}$ & ${}$ &  &  &  \\\\\n'.format(self.gnum, self.label))


class RevDecompText:
    def __init__(self, matrix, gnum, g_cosets, snum, s_cosets, bnum, bieberbach, normalB, b_matrix, b_cosets):
        self.fd_labels = self.get_fd_labels()
        self.subgroup_labels = [self.get_label(gnum), self.get_label(snum), self.get_label(bnum)]
        self.cosets_list = [self.format_cosets(g_cosets, 0), self.format_cosets(s_cosets, 1), self.format_cosets(b_cosets, 2)]
        self.bieberbach = bieberbach
        self.normalB = normalB
        self.subgroups = [self.format_subgroup(1, matrix), self.format_subgroup(2, b_matrix)]
        self.decomp_text = self.write_decomp()
    def get_label(self, num):
        with open('tex/sgroup_labels.tex', 'r') as file:
            return [line.strip() for i, line in enumerate(file) if i==int(num)-1][0]
    def format_matrix(self, M):
        return r'\begin{pmatrix} ' + M[0] + r' & ' + M[1] + r' & ' + M[2] + r' & ' + M[3] + r' \\ ' + M[4] + r' & ' + M[5] + r' & ' + M[6] + r' & ' + M[7] + r' \\ ' + M[8] + r' & ' + M[9] + r' & ' + M[10] + r' & ' + M[11] + r' \end{pmatrix}'
    def format_subgroup(self, index, matrix):
        if index == 1: 
            sub, mat_label, normal  = r'_S', r'\alpha', r'\ \text{(normal)}'
        if index == 2: 
            mat_label = r'Y'
            if self.bieberbach:
                sub = r'_B'
            else:
                sub = r'/\Gamma_S'
            if self.normalB:
                normal = r'\ \text{(normal)}'
            else:
                normal = ''
        return r'\[ \Gamma' + sub + r' = ' + self.subgroup_labels[index] + normal +  r', \ \ \  ' + mat_label + r' = ' + self.format_matrix(matrix) +  r',\]'
    def get_fd_labels(self):
        return [r'F_{\frac{\Gamma^Z}{P1}}',r'F_{\frac{\Gamma_{S}^{X}}{P1}}',r'F_{\frac{\Gamma^{Y}}{P1}}']
    def format_cosets(self, cosets, fd_index):
        reps = cosets.split('\n')
        reps = list(filter(None, reps))
        reps = ['({})'.format(rep) for rep in reps]
        reps = '; '.join(reps)
        # return r'{\centering $' + self.fd_labels[fd_index] + r' = \{' + reps + r' $\} \par}'
        return self.insert_newlines(r' ' + self.fd_labels[fd_index] + r' = \{' + reps + r' \}')
    def insert(self, source_str, insert_str, pos):
        return source_str[:pos]+insert_str+source_str[pos:]
    def insert_newlines(self, text):
        #insert a r'\\' after every 74 chars
        s = text.split('=')
        c = s[-1]
        len_text = len(c) + 7
        if len_text <= 68: return text
        num = int(len_text / 68)
        if num == 1: return s[0] + r' = ' + self.insert(c, r'\\ ', 68)
        # many lines
        pos = 68
        for i in range(num):
            c = self.insert(c, r'\\ ', pos)
            pos += 68
        return s[0] + r' = ' + c
    def write_decomp(self):
        return r'Given the subgroup ' + self.subgroups[0] + r' and setting $Z=\alpha$ and $X=\mathbb{I}$, we find a complementary subgroup ' + self.subgroups[1] + r' resulting in the following decomposition: \begin{gather*}' + self.cosets_list[0] + r'\\ ' + self.cosets_list[1] + r'\\ ' + self.cosets_list[2] + r'. \end{gather*}'


class RevWriteUp:
    def __init__(self, gnum):
        self.gnum = gnum
        self.g_label = self.get_label()
        self.header = r'\subsection{{$\bf {}$, $\bf {}$}}'.format(gnum, self.g_label)
        # self.subheader = r'\subsubsection{{$\bf {}$, $\bf {}$}}'.format(gnum, self.g_label)
        self.subheader = ''
        self.bieb_groups = self.get_bieb_groups()
        self.decomps = []
        self.get_decomps()
        self.write_to_file()
    def get_label(self):
        with open('tex/sgroup_labels.tex', 'r') as file:
            return [line.strip() for i, line in enumerate(file) if i==int(self.gnum)-1][0]
    def get_bieb_groups(self):
        with open('dat_files/bgroups.dat') as file:
            return [line.strip() for line in file]
    def get_decomps(self):
        filename = 'decomp/one_non_semi_copy/decomp_{}'.format(self.gnum)
        with open(filename, 'r') as file:
            while True:
                line = file.readline() # transformation matrix: 
                if not line: break
                matrix = ''
                for i in range(3): matrix += file.readline()
                file.readline() # Gamma
                line = ''
                g_cosets = ''
                while True:
                    line = file.readline()
                    if line.startswith('Gamma'):
                        break
                    g_cosets += line
                snum = line.split()[1]
                s_cosets = ''
                while True:
                    line = file.readline()
                    if line.startswith('Gamma'):
                        break
                    s_cosets += line
                line = line.split()
                bnum, normalB = line[1], True
                if bnum in self.bieb_groups: 
                    bieberbach = True
                else:
                    bieberbach = False
                if line[2] == "not": normalB = False
                file.readline() # transformation matrix:
                b_matrix = ''
                for i in range(3): b_matrix += file.readline()
                b_cosets = ''
                while True:
                    line = file.readline()
                    if not line.strip():
                        break
                    b_cosets += line
                self.decomps.append(RevDecompText(matrix.split(), self.gnum, g_cosets, snum, s_cosets, bnum, bieberbach, normalB, b_matrix.split(), b_cosets))
    def write_to_file(self):
        with open('tex/decomps/one_non_semidirect_copy/decomp_{}.tex'.format(self.gnum), 'w') as file:
            file.write(self.header + '\n')
            for decomp in self.decomps:
                file.write(self.subheader + '\n')
                file.write(decomp.decomp_text + '\n')


class CompileReport:
    def __init__(self):
        start_text = r'\documentclass[9pt]{article} \usepackage{amsmath, amsthm, amscd, amsfonts,amssymb, graphicx} \usepackage{accents} \usepackage{xcolor} \usepackage{hyperref} \usepackage{booktabs} \usepackage{ltxtable} \usepackage{array} \usepackage[fleqn]{mathtools} \newcolumntype{L}{>{\raggedright \arraybackslash}} \renewcommand{\arraystretch}{1.25} \begin{document}  \section{\bf Semidecomposition Table} \LTXtable{\textwidth}{t1.tex} \section{\bf Reverse semidecomposition Table} \LTXtable{\textwidth}{t2.tex} \section{\bf Examples of decomposition results} '
        end_text = r' \end{document}'
        self.semi_files, self.non_semi_files = self.sort_files()
        self.text = start_text + self.get_text() + end_text
        self.write_text()
    def sort_files(self):
        semi_files = os.listdir('tex/decomps/one_semidirect/')
        semi_files.sort(key=lambda x: int(x.split('_')[-1].split('.')[0]))
        non_semi_files = os.listdir('tex/decomps/one_non_semidirect_copy/')
        non_semi_files.sort(key=lambda x: int(x.split('_')[-1].split('.')[0]))
        return semi_files, non_semi_files
    def get_text(self):
        text = ''
        for semi_file in self.semi_files:
            with open('tex/decomps/one_semidirect/' + semi_file, 'r') as file:
                text += file.read()
        text += r' \section{\bf Examples of reverse decomposition results} '
        for non_semi_file in self.non_semi_files:
            with open('tex/decomps/one_non_semidirect_copy/' + non_semi_file, 'r') as file:
                text += file.read()
        return text
    def write_text(self):
        with open('tex/one_report_final.tex','w') as file:
            file.write(self.text)


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


class GetDecomps:
    def __init__(self, gnum):
        self.gnum = gnum
        self.base_dir = 'decomp/non_semi_sorted/decomp_{}/'.format(gnum)
        self.reduced_dir = self.get_reduced_dir()
        self.bieb_groups = self.get_bieb_groups()
        self.fnames = []
        self.get_fnames()
        self.decomps = []
        self.decomp = []
        self.texts = []
        self.text = []
        self.get_decomps()
    def get_bieb_groups(self):
        if 'non' in self.base_dir:
            filename = 'dat_files/bgroups.dat'
        else:
            filename = 'dat_files/sgroups.dat'
        with open(filename, 'r') as file:
            return [line.strip() for line in file]
    def get_reduced_dir(self):
        if 'non' in self.base_dir:
            return 'decomp/non_semi_reduced/decomp_{}/'.format(self.gnum)
        else:
            return 'decomp/semi_reduced/decomp_{}/'.format(self.gnum)
    def make_base_dir(self):
        if not os.path.exists(self.base_dir):
            os.mkdir(self.base_dir)
    def make_reduced_dir(self):
        if not os.path.exists(self.reduced_dir):
            os.mkdir(self.reduced_dir)
    def get_decomps(self, fname=None):
        if 'non' in self.base_dir:
            filename = 'decomp/non_semidirect/decomp_{}'.format(self.gnum)
        else:
            filename = 'decomp/semidirect/decomp_{}'.format(self.gnum)
        if fname:
            filename = fname
        with open(filename, 'r') as file:
            if fname:
                self.decomp,self.text = [],[]
            else:
                self.decomps,self.texts = [],[]
            while True:
                line = file.readline() # transformation matrix:
                if not line: break
                text = line
                matrix = ''
                for i in range(3): matrix += file.readline()
                text += matrix
                gamma = file.readline() # Gamma
                text += gamma
                line = ''
                g_cosets = ''
                while True:
                    line = file.readline()
                    if line.startswith('Gamma'):
                        break
                    g_cosets += line
                text += g_cosets
                text += line
                line = line.split()
                snum, ind = line[1], line[-1]
                s_cosets = ''
                while True:
                    line = file.readline()
                    if line.startswith('Gamma'):
                        break
                    s_cosets += line
                text += s_cosets
                text += line
                line = line.split()
                bnum, normalB = line[1], True
                if bnum in self.bieb_groups:
                    bieberbach = True
                else:
                    bieberbach = False
                if line[2] == "not": normalB = False
                line = file.readline() # transformation matrix
                text += line
                b_matrix = ''
                for i in range(3): b_matrix += file.readline()
                text += b_matrix
                b_cosets = ''
                while True:
                    line = file.readline()
                    if not line.strip():
                        break
                    b_cosets += line
                text += b_cosets
                text += line
                if fname:
                    self.text.append(text)
                    self.decomp.append(Decomp(self.gnum, snum, matrix, Subgroup(bnum, normalB, bieberbach)))
                else:
                    self.texts.append(text)
                    self.decomps.append(Decomp(self.gnum, snum, matrix, Subgroup(bnum, normalB, bieberbach)))
    def count_indices(self):
        self.get_fnames()
        for fname in self.fnames:
            indices = []
            self.count_uniq_indices(fname, indices)
            uniq_indices = list(set(indices))
            num_indices = len(uniq_indices)
            if num_indices != 1:
                print(fname)
                print(num_indices)
    def count_uniq_indices(self, fname, indices):
        with open(fname, 'r') as file:
            while True:
                line = file.readline() # transformation matrix:
                if not line: break
                matrix = ''
                for i in range(3): matrix += file.readline()
                gamma = file.readline() # Gamma
                line = ''
                g_cosets = ''
                while True:
                    line = file.readline()
                    if line.startswith('Gamma'):
                        break
                    g_cosets += line
                line = line.split()
                snum, ind = line[1], line[-1]
                s_cosets = ''
                while True:
                    line = file.readline()
                    if line.startswith('Gamma'):
                        break
                    s_cosets += line
                line = line.split()
                bnum, normalB = line[1], True
                if bnum in self.bieb_groups:
                    bieberbach = True
                else:
                    bieberbach = False
                if line[2] == "not": normalB = False
                line = file.readline() # transformation matrix
                b_matrix = ''
                for i in range(3): b_matrix += file.readline()
                b_cosets = ''
                while True:
                    line = file.readline()
                    if not line.strip():
                        break
                    b_cosets += line
                indices.append(ind)
    def get_matrices_from_file(self, fname, matrices):
        with open(fname, 'r') as file:
            while True:
                line = file.readline() # transformation matrix:
                if not line: break
                matrix = ''
                for i in range(3): matrix += file.readline()
                gamma = file.readline() # Gamma
                line = ''
                g_cosets = ''
                while True:
                    line = file.readline()
                    if line.startswith('Gamma'):
                        break
                    g_cosets += line
                line = line.split()
                snum, ind = line[1], line[-1]
                s_cosets = ''
                while True:
                    line = file.readline()
                    if line.startswith('Gamma'):
                        break
                    s_cosets += line
                line = line.split()
                bnum, normalB = line[1], True
                if bnum in self.bieb_groups:
                    bieberbach = True
                else:
                    bieberbach = False
                if line[2] == "not": normalB = False
                line = file.readline() # transformation matrix
                b_matrix = ''
                for i in range(3): b_matrix += file.readline()
                b_cosets = ''
                while True:
                    line = file.readline()
                    if not line.strip():
                        break
                    b_cosets += line
                matrices.append([gm.load_matrix(matrix.split()), gm.load_matrix(b_matrix.split())])
    def get_reduced_matrices_from_file(self, fname):
        m = []
        self.get_matrices_from_file(fname, m)
        return [[LA.qr(mat_pair[0])[0], LA.qr(mat_pair[1])[0]] for mat_pair in m]
    def uniq_qr_decomps(self):
        for fname in self.fnames:
            self.uniq_qr_decomps_from_file(fname)
    def uniq_qr_decomps_from_file(self, fname):
        matrices = d.get_reduced_matrices_from_file(fname)
        alphas = [m[0] for m in matrices]
        comps = []
        for i,alpha in enumerate(alphas):
            if i == len(alphas)-1:
                break
            comps.append(alpha.all() == alphas[i+1].all())
        if not all(comps): 
            print(fname)
            print('NON UNIQ QR DECOMP FOUND')
    def sort_files(self):
        self.create_files()
        self.write_files()
    def create_files(self):
        snums = [d.snum for d in self.decomps]
        bnums = [d.bgroups.num for d in self.decomps]
        fnames = list(set((list(zip(snums, bnums)))))
        for fname in fnames:
            f = open(self.base_dir + 'decomp_{}_{}'.format(fname[0], fname[1]), 'w')
            f.close()
    def write_files(self):
        for ind,decomp in enumerate(self.decomps):
            snum = decomp.snum
            bnum = decomp.bgroups.num
            fname = self.base_dir + 'decomp_{}_{}'.format(snum, bnum)
            with open(fname, 'a') as f:
                f.write(self.texts[ind])
    def get_fnames(self):
        fnames = os.listdir(self.base_dir)
        nums = [fname[7:].split('_') for fname in fnames]
        nums.sort()
        buckets = defaultdict(list)
        for pair in nums: 
            buckets[pair[0]].append(pair[1])
        nums = []
        for key in buckets.keys():
            l = buckets[key]
            l.sort(key=lambda x: int(x))
            nums.append([key,l])
        for pair in nums:
            [self.fnames.append(self.base_dir + 'decomp_{}_{}'.format(pair[0], num)) for num in pair[1]]
    def fnames_to_rm(self):
        for fname in self.fnames:
            nums = fname.split('_')[-2:]
            if self.gnum in nums:
                print(fname)
    def combine_files(self):
        if 'non' in self.base_dir:
            filename = 'decomp/non_semi_combined/decomp_{}'.format(self.gnum)
        else:
            filename = 'decomp/semi_combined/decomp_{}'.format(self.gnum)
        with open(filename, 'w') as outfile: 
            for fname in self.fnames:
                with open(self.reduced_dir + fname.split('/')[-1]) as infile:
                    for line in infile:
                        outfile.write(line)
    def write_one_file(self):
        if 'non' in self.base_dir:
            in_dir_name = 'decomp/non_semi_reduced/decomp_{}/'.format(self.gnum)
            out_fname = 'decomp/one_non_semi_copy/decomp_{}'.format(self.gnum)
        else:
            in_dir_name = 'decomp/semi_reduced/decomp_{}/'.format(self.gnum)
            out_fname = 'decomp/one_semi_copy/decomp_{}'.format(self.gnum)
        in_fname = in_dir_name + self.fnames[0].split('/')[-1]
        with open(out_fname, 'w') as outfile:
            with open(in_fname, 'r') as infile:
                for line in infile:
                    outfile.write(line)
    def write_uniq_file(self, fname, ind):
        out_fname = self.reduced_dir + fname.split('/')[-1]
        with open(out_fname, 'w') as outfile:
            outfile.write(self.text[ind])
    def get_ind_to_keep(self):
        index = 0
        for ind,decomp in enumerate(self.decomp):
            if decomp.bgroups.normal:
                index = ind
                break
        return index
    def make_uniq_decomps(self):
        self.make_reduced_dir()
        for fname in self.fnames:
            self.get_decomps(fname)
            ind = self.get_ind_to_keep()
            self.write_uniq_file(fname, ind)


class TableEntry:
    def __init__(self, type):
        self.type = type
        self.gnums = self.get_gnums()
        self.nums = self.get_sym_or_bieb_nums()
        self.dir = self.get_dir()
        self.table_dir = self.get_table_dir()
        self.names = self.get_names()
    def get_gnums(self):
        if self.type == 0:
            fname = 'decomp/stats/semi_decomp_stats.txt'
        else:
            fname = 'decomp/stats/non_semi_decomp_stats.txt'
        with open(fname, 'r') as file:
            return [line.strip() for line in file]
    def get_dir(self):
        if self.type == 0:
            return 'decomp/semi_reduced/'
        else:
            return 'decomp/non_semi_reduced/'
    def get_table_dir(self):
        if self.type == 0:
            return 'tex/table_entries/semidirect/'
        else:
            return 'tex/table_entries/non_semidirect/'
    def get_names(self):
        with open('tex/sgroup_labels.tex') as file:
            return [line.strip() for line in file]
    def get_sym_or_bieb_nums(self):
        if self.type == 0:
            filename = 'dat_files/sgroups.dat'
        else:
            filename = 'dat_files/bgroups.dat'
        with open(filename, 'r') as file:
            return [line.strip() for line in file]
    def print_entries(self):
        for gnum in self.gnums:
            self.print_entry(gnum)
    def get_info(self, dir_name):
        fnames = os.listdir(dir_name)
        nums = [fname[7:].split('_') for fname in fnames]
        nums.sort()
        buckets = defaultdict(list)
        for pair in nums: 
            buckets[pair[0]].append(pair[1])
        nums = []
        for key in buckets.keys():
            l = buckets[key]
            l.sort(key=lambda x: int(x))
            nums.append([key,l])
        return nums
    def get_fnames(self, entry, dir_name):
        return [dir_name + 'decomp_{}_{}'.format(entry[0], bnum) for bnum in entry[1]]
    def get_normal(self, entry, dir_name):
        fnames = self.get_fnames(entry, dir_name)
        normal = []
        for fname in fnames:
            with open(fname, 'r') as file:
                if 'not' in file.read(): 
                    normal.append(False) 
                else:
                    normal.append(True)
        if normal[0]: normal_text = 'Y'
        else: normal_text = 'N'
        for i in range(1,len(normal)):
            if normal[i]: normal_text += ',Y'
            else: normal_text += ',N'
        return normal_text
    def get_bnames(self, bnums):
        bnames = self.names[int(bnums[0])-1]
        if bnums[0] not in self.nums: bnames = '{{\color{{red}}{}}}'.format(bnames)
        for i in range(1,len(bnums)):
            name = self.names[int(bnums[i])-1]
            if bnums[i] not in self.nums: name = '{{\color{{red}}{}}}'.format(name)
            bnames += ',' + name
        return bnames
    def print_entry(self, gnum):
        dir_name = self.dir + 'decomp_{}/'.format(gnum)
        entry_fname = self.table_dir + 'entry_{}'.format(gnum)
        entries = self.get_info(dir_name)
        gname = self.names[int(gnum)-1]
        snames,bnames,normal = [],[],[]
        for entry in entries:
            snames.append(self.names[int(entry[0])-1])
            bnames.append(self.get_bnames(entry[1]))
            normal.append(self.get_normal(entry, dir_name))
        with open(entry_fname, 'w') as outfile:
            for sname,bname,normal_text in zip(snames, bnames, normal):
                outfile.write('${}$ & ${}$ & ${}$ & ${}$ & {} \\\\\n'.format(gnum, gname, sname, bname, normal_text))


def get_gnums(type):
    fname = ''
    if type == 0:
        fname = 'decomp/stats/semi_decomp_stats.txt'
    elif type == 1:
        fname = 'decomp/stats/non_semi_decomp_stats.txt'
    else:
        fname = 'decomp/stats/no_table_entries.txt'
    with open(fname, 'r') as file:
        return [line.strip() for line in file]

