""" NOT FOR USE.

The classes GetBieberbachSubgroups and GetBiebMatrices were used to download Bieberbach space-subgroups of each space group along with the affine matrices used to conjugate the space (super)group. The results are stored in the folder 'b_matrices/' in separate folders for normal and not normal subgroups. Likewise, the classes GetSymSubgroups and GetSymMatrices were used to download symmorphic space-subgroups for each space group and affine matrices used to conjugate the space (super)group.  Results are stored in 's_matrices/' separated by normal and not normal subgroups."""


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


class GetBieberbachSubgroups:
    def __init__(self, gnum):
        self.gnum = gnum
        self.bnums = [4,7,9,19,29,33,76,78,144,145,169,170]
    async def download_coroutine(self, session, loop, bnum):
        with async_timeout.timeout(None):
            async with aiohttp.ClientSession(loop=loop) as session:
                async with session.get('http://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-subggraph?super={}&sub={}&index='.format(self.gnum, bnum)) as response:
                    async with aiofiles.open('b_subgroups/b_subgroups_{}_{}.dat'.format(self.gnum, bnum), 'wb') as fd:
                        while True:
                            chunk = await response.content.read(10240)
                            if not chunk:
                                break
                            await fd.write(chunk)
                    return await response.release()
    def ensure_all_futures(self, loop):
        with aiohttp.ClientSession(loop=loop) as session:
            return [asyncio.ensure_future(self.download_coroutine(session, loop, bnum)) for bnum in self.bnums]
    def get_B_subgroups(self):
        loop = asyncio.get_event_loop()
        loop.run_until_complete(asyncio.wait(self.ensure_all_futures(loop)))
    def get_urls_from_response(self):
        with open('b_subgroups/b_subgroups_{}.dat'.format(self.gnum), 'w') as fd:
            for bnum in self.bnums:
                with open('b_subgroups/b_subgroups_{}_{}.dat'.format(self.gnum, bnum), 'r') as file:
                    text = file.read()
                soup = BeautifulSoup(text, 'html5lib')
                if soup.find_all('a')[3:-2]:
                    fd.write(str(bnum)+'\n')


class GetBiebMatrices:
    def __init__(self, gnum):
        self.gnum = gnum
        self.bnums = self.get_B_subgroups_from_file()
        self.ind_T = self.get_ind_T_Gamma()
        self.urls = self.get_chain_urls()
        self.base_url = 'http://www.cryst.ehu.es/cgi-bin/cryst/programs/'
        # self.matrix_urls = self.get_matrix_urls()
    def get_B_subgroups_from_file(self):
        with open('b_subgroups/b_subgroups_{}.dat'.format(self.gnum), 'r') as file:
            return [line.strip() for line in file]
    def get_ind_T_Gamma(self):
        with open('dat_files/ind_T_gamma.dat', 'r') as file:
            return int([line.strip() for i,line in enumerate(file) if i==int(self.gnum)-1][0])
    def get_chain_urls(self):
        max_ind = self.ind_T / 2
        urls = []
        for bnum in self.bnums:
            ind = self.get_start_index(bnum)
            print(ind)
            if ind == 0:
                print(self.gnum)
            while ind <= max_ind:
                if self.ind_T % ind:
                    ind += 1
                    continue
                urls.append('http://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-subggraph?super={}&sub={}&index={}'.format(self.gnum, bnum, ind))
                ind += 1
        return urls
    def get_start_index(self, bnum):
        B = gm.loadGroup('dat_files/bieberbach_cosets/b{}.dat'.format(bnum))
        ind = self.ind_T / B.shape[2]
        if ind == 1: return 2
        if int(ind) == 0: return 2
        return int(ind)
    def ensure_all_futures(self, loop):
        with aiohttp.ClientSession(loop=loop) as session:
            return [asyncio.ensure_future(self.download_coroutine(session, loop, url)) for url in self.urls]
    async def download_coroutine(self, session, loop, url):
        with async_timeout.timeout(None):
            async with aiohttp.ClientSession(loop=loop) as session:
                async with session.get(url) as response:
                    text = await response.text()
                    soup = BeautifulSoup(text, 'html5lib')
                    matrix_links = soup.find_all('a')[4:-3]
                    if matrix_links:
                        async with aiofiles.open('subgroup_chains/b_matrix_urls_{}.dat'.format(self.gnum), 'a') as fd:
                            [await fd.write(self.base_url + link['href'] +'\n') for link in matrix_links]
    def get_matrix_urls(self):
        filename = 'subgroup_chains/b_matrix_urls_{}.dat'.format(self.gnum)
        open(filename, 'w').close()
        loop = asyncio.get_event_loop()
        loop.run_until_complete(asyncio.wait(self.ensure_all_futures(loop)))
        with open(filename, 'r') as file:
            return [line.strip() for line in file]
    def ensure_matrix_futures(self, loop):
        with aiohttp.ClientSession(loop=loop) as session:
            return [asyncio.ensure_future(self.download_matrices(session, loop, url)) for url in self.matrix_urls]
    async def download_matrices(self, session, loop, url):
        with async_timeout.timeout(None):
            async with aiohttp.ClientSession(loop=loop) as session:
                async with session.get(url) as response:
                    text = await response.text()
                    matrices = self.format_matrix_text(text)
                    async with aiofiles.open('b_matrices/b_matrices_{}.dat'.format(self.gnum), 'a') as fd:
                        await fd.write(matrices)
    def format_matrix_text(self, text):
        soup = BeautifulSoup(text, 'html5lib')
        text = soup.find('center').text.replace('\n\n', '\n')
        text = re.sub('\n', lambda m,c=itertools.count():m.group() if next(c)%3 else '\n\n', text)
        chain = text[:text.find('\n\n')]
        indices = chain[chain.find('[')+1: -1]
        index = str(int(reduce((lambda x,y: int(x)*int(y)), indices.split())))
        bnum = str(int(chain[:chain.find('[')].split()[-1]))
        text = text.replace('\n\n', '\n\n{} {}\n'.format(bnum, index))
        return text[text.find('\n')+2: text.rfind('\n\n')+2]
    def get_matrices(self):
        # filename = 'b_matrices_{}.dat'.format(self.gnum)
        # open(filename, 'w').close()
        loop = asyncio.get_event_loop()
        loop.run_until_complete(asyncio.wait(self.ensure_matrix_futures(loop)))


class GetSymSubgroups:
    def __init__(self, gnum):
        self.gnum = gnum
        self.snums = self.get_sym_groups()
    def get_sym_groups(self):
        with open('dat_files/sgroups.dat', 'r') as file:
            return [line.strip() for line in file]
    async def download_coroutine(self, session, loop, snum):
        with async_timeout.timeout(None):
            async with aiohttp.ClientSession(loop=loop) as session:
                async with session.get('http://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-subggraph?super={}&sub={}&index='.format(self.gnum, snum)) as response:
                    async with aiofiles.open('s_subgroups/s_subgroups_{}_{}.dat'.format(self.gnum, snum), 'wb') as fd:
                        while True:
                            chunk = await response.content.read(10240)
                            if not chunk:
                                break
                            await fd.write(chunk)
                    return await response.release()
    def ensure_all_futures(self, loop):
        with aiohttp.ClientSession(loop=loop) as session:
            return [asyncio.ensure_future(self.download_coroutine(session, loop, snum)) for snum in self.snums]
    def get_S_subgroups(self):
        loop = asyncio.get_event_loop()
        loop.run_until_complete(asyncio.wait(self.ensure_all_futures(loop)))
    def get_urls_from_response(self):
        with open('s_subgroups/s_subgroups_{}.dat'.format(self.gnum), 'w') as fd:
            for snum in self.snums:
                with open('s_subgroups/s_subgroups_{}_{}.dat'.format(self.gnum, snum), 'r') as file:
                    text = file.read()
                soup = BeautifulSoup(text, 'html5lib')
                if soup.find_all('a')[3:-2]:
                    fd.write(str(snum)+'\n')


class GetSymMatrices:
    def __init__(self, gnum):
        self.gnum = gnum
        self.snums = self.get_S_subgroups_from_file()
        self.ind_T = self.get_ind_T_Gamma()
        self.urls = self.get_chain_urls()
        self.base_url = 'http://www.cryst.ehu.es/cgi-bin/cryst/programs/'
        self.matrix_urls = self.get_matrix_urls()
    def get_S_subgroups_from_file(self):
        with open('s_subgroups/s_subgroups_{}.dat'.format(self.gnum), 'r') as file:
            return [line.strip() for line in file]
    def get_ind_T_Gamma(self):
        with open('dat_files/ind_T_gamma.dat', 'r') as file:
            return int([line.strip() for i,line in enumerate(file) if i==int(self.gnum)-1][0])
    def get_chain_urls(self):
        max_ind = self.ind_T / 2
        urls = []
        for snum in self.snums:
            ind = self.get_start_index(snum)
            while ind <= max_ind:
                if self.ind_T % ind:
                    ind += 1
                    continue
                urls.append('http://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-subggraph?super={}&sub={}&index={}'.format(self.gnum, snum, ind))
                ind += 1
        return urls
    def get_start_index(self, snum):
        S = gm.loadGroup('dat_files/symmorphic_cosets/s{}.dat'.format(snum))
        ind = self.ind_T / S.shape[2]
        if ind == 1: return 2
        if int(ind) == 0: return 2
        return int(ind)
    def ensure_all_futures(self, loop):
        with aiohttp.ClientSession(loop=loop) as session:
            return [asyncio.ensure_future(self.download_coroutine(session, loop, url)) for url in self.urls]
    async def download_coroutine(self, session, loop, url):
        with async_timeout.timeout(None):
            async with aiohttp.ClientSession(loop=loop) as session:
                async with session.get(url) as response:
                    text = await response.text()
                    soup = BeautifulSoup(text, 'html5lib')
                    matrix_links = soup.find_all('a')[4:-3]
                    if matrix_links:
                        async with aiofiles.open('subgroup_chains/s_matrix_urls_{}.dat'.format(self.gnum), 'a') as fd:
                            [await fd.write(self.base_url + link['href'] +'\n') for link in matrix_links]
    def get_matrix_urls(self):
        filename = 'subgroup_chains/s_matrix_urls_{}.dat'.format(self.gnum)
        open(filename, 'w').close()
        loop = asyncio.get_event_loop()
        loop.run_until_complete(asyncio.wait(self.ensure_all_futures(loop)))
        with open(filename, 'r') as file:
            return [line.strip() for line in file]
    def ensure_matrix_futures(self, loop):
        with aiohttp.ClientSession(loop=loop) as session:
            return [asyncio.ensure_future(self.download_matrices(session, loop, url)) for url in self.matrix_urls]
    async def download_matrices(self, ssession, loop, url):
        with async_timeout.timeout(None):
            async with aiohttp.ClientSession(loop=loop) as session:
                async with session.get(url) as response:
                    text = await response.text()
                    matrices = self.format_matrix_text(text)
                    async with aiofiles.open('s_matrices/s_matrices_{}.dat'.format(self.gnum), 'a') as fd:
                        await fd.write(matrices)
    def format_matrix_text(self, text):
        soup = BeautifulSoup(text, 'html5lib')
        text = soup.find('center').text.replace('\n\n', '\n')
        text = re.sub('\n', lambda m,c=itertools.count():m.group() if next(c)%3 else '\n\n', text)
        chain = text[:text.find('\n\n')]
        indices = chain[chain.find('[')+1: -1]
        index = str(int(reduce((lambda x,y: int(x)*int(y)), indices.split())))
        snum = str(int(chain[:chain.find('[')].split()[-1]))
        text = text.replace('\n\n', '\n\n{} {}\n'.format(snum, index))
        return text[text.find('\n')+2: text.rfind('\n\n')+2]
    def get_matrices(self):
        filename = 's_matrices/s_matrices_{}.dat'.format(self.gnum)
        open(filename, 'w').close()
        loop = asyncio.get_event_loop()
        loop.run_until_complete(asyncio.wait(self.ensure_matrix_futures(loop)))

