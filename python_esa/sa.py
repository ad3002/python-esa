#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 02.04.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com

import numpy
from PyExp.experiments.abstract_experiment import Timer

class ESA(object):
    ''' Class for Enhance Suffix Array construction in Python.'''


    def __init__(self, corpus, delimiter="$", lcp_file="lcp_file.data"):
        ''' Corpus is a list of text, 
        lcp_file is an output file for computed lcp classes.'''

        self.delimiter = delimiter
        self.corpus = corpus
        self.N = len(corpus)
        self.lcp_file = lcp_file
        self.lcp_result = []

    def compute_esa(self):
        ''' Compute SA, LCP, DOCS, and LCP classes tables. '''
        with Timer("SA construction..."):
            self.compute_suffix_array()
        with Timer("LCP construction..."):
            self.compute_lcp()
        with Timer("DOCS construction..."):
            self.compute_docs()
        with Timer("LCP classes construction..."):
            self.compute_classes()

    def _sa_cmp(self, x, y):
        ''' Comparison function. '''
        while x != self.N - 1 \
                    and y != self.N - 1 \
                    and self.corpus[x] == self.corpus[y] \
                    and self.corpus[x] != self.delimiter \
                    and self.corpus[y] != self.delimiter:
            x += 1
            y += 1
        return cmp(self.corpus[x], self.corpus[y])

    def compute_suffix_array(self):

        with open("sa.temp", "wb") as fh:
            numpy.save(fh, numpy.arange(0, self.N))
        with open("sa.temp", "rb") as fh:
            self.sa = numpy.load(fh)
        with Timer("Sort SA"):
            self.sa = sorted(self.sa, cmp=self._sa_cmp)

    def _get_lcp(self, x, y):

        k = 0
        while x != self.N - 1 \
                    and y != self.N - 1 \
                    and self.corpus[x] == self.corpus[y] \
                    and self.corpus[x] != self.delimiter \
                    and self.corpus[y] != self.delimiter:
            x += 1
            y += 1
            k += 1
        return k

    def compute_lcp(self):
        ''' Compute LCP table.'''
        self.lcp = numpy.zeros(self.N + 1)
        for i in xrange(1, self.N):
            self.lcp[i] = self._get_lcp(self.sa[i - 1], self.sa[i])

    def compute_docs(self):
        ''' Compute DOCS table.'''

        print "Docs init"
        self.docs = numpy.zeros(self.N)
        print "Docs compute"
        docid = 0
        for x in xrange(0, self.N):
            if self.corpus[x] == self.delimiter:
                docid += 1
            self.docs[x] = docid
        self.docN = docid + 1

    def _print_string(self, start, n):
        start = int(start)
        n = int(n)
        result = []
        for x in range(start, start + n):
            word = self.corpus[x]
            if word == self.delimiter:
                break
            result.append(word)
        if not result:
            result = ["---"]
        return " ".join(result)

    def print_sa(self):
        for i in xrange(0, self.N):
            print i, self.sa[i], self.lcp[i], self._print_string(self.sa[i], self.lcp[i]), self.docs[self.sa[i]]

    def _dec_df(self, docid, sp, doclink, stack_i, stack_df):
        beg = 0
        end = sp
        mid = sp / 2

        while beg != mid:
            if doclink[docid] > stack_i[mid]:
                beg = mid
            else:
                end = mid
            mid = (beg + end) / 2
        stack_df[mid] -= 1

    def _output(self, i, j, k, df):

        lbp = max(self.lcp[i], self.lcp[j + 1])

        if i == j:
            return

        text = self._print_string(self.sa[k], self.lcp[k])
        length = self.lcp[k]
        tf = j - i + 1

        data = "%s\t%s\t%s\t%s\n" % (text, tf, df, length)

        if self.lcp_file:
            self.lcp_result.append(data)


    def _doc_id(self, i):
        return int(self.docs[self.sa[i]])

    def compute_classes(self):
        ''' Compute LCP classes. 
            It is a copy of C implementation.
            TODO: make it more pythonic.
        '''

        doclink = [-1] * self.docN
        stack_i = [0] * 100000
        stack_k = [0] * 100000
        stack_df = [0] * 100000
        stack_df[0] = 1
        sp = 1

        for j in xrange(0, self.N):

            print j * 100.0 / self.N, "\r",

            self._output(j, j, 0, 1)

            if doclink[self._doc_id(j)] != -1:
                self._dec_df(self._doc_id(j), sp, doclink, stack_i, stack_df)
            doclink[self._doc_id(j)] = j
            df = 1

            while self.lcp[j + 1] < self.lcp[stack_k[sp - 1]]:
                df = stack_df[sp - 1] + df
                self._output(stack_i[sp - 1], j, stack_k[sp - 1], df)
                sp -= 1

            if self.lcp[stack_k[sp - 1]] == self.lcp[j + 1]:
                stack_k[sp - 1] = j + 1
                stack_df[sp - 1] += df
            else:
                stack_i[sp] = stack_k[sp - 1]
                stack_k[sp] = j + 1
                stack_df[sp] = df
                sp += 1

        if self.lcp_file:
            with open(self.lcp_file, "w") as fh:
                fh.writelines(self.lcp_result)

def sc_get_lcp_for_given_corpus(corpus, lcp_file):
    ''' Shortcut for ESA computation.'''

    esa = ESA(corpus, lcp_file=lcp_file)
    esa.compute_esa()
    data = []
    with open(esa.lcp_file) as fh:
        for line in fh:
            items = line.split("\t")
            data.append((items[2], items[0]))
    data.sort(reverse=True)
    return [x[1] for x in data]





