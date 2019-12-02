import numpy as np
from Qs import Alignment

class BandedDP(Alignment):
    def __init__(self, *args):
        super(BandedDP, self).__init__(*args)
    
    def init_domain(self, s1, s2, d, k):
        self.s1, self.s2, self.n, self.m = s1, s2, len(s1), len(s2)
        
        # Bandwidth
        self.k_h = k//2
        
        # pushing d away from the corners
        if d < self.k_h - self.m and d > self.n - self.k_h:
            # do nothing, high bandwidth, maybe run quad mem local alignment
            d = d
        elif d < self.k_h - self.m: 
            d = self.k_h - self.m
        elif d > self.n - self.k_h:
            d = self.n - self.k_h
        else:
            d = d
        # Origin I've defined, top right corner of band in matrix coords
        (self.io, self.jo) = (d - self.k_h, 0) if d - self.k_h >= 0 else (0, self.k_h - d)
        # The i, j values we iterate over
        self.J, self.I = self.find_domain()
        
        self.dl, self.du = d+self.k_h, d-self.k_h
        
        self.local_score = float("-inf")
        
    def find_domain(self):
        if self.jo == 0:# Case 1
#            print('CASE1')
            i = lambda j: range(1 + self.io + j, min(self.io + j + 2*self.k_h, self.n) + 1 )
            return range(1, min(self.n - self.io, self.m) + 1 ), i
        # i == 0 now
        if self.jo <= 2*self.k_h:
#            print('CASE2')
            i = lambda j: range(1 + max(0, j - self.jo), min(j - self.jo + 2*self.k_h, self.n) + 1 )
            return range(1, self.m + 1), i
        if self.jo <= self.m:
#            print('CASE3')
            i = lambda j: range(1 + max(0, j - self.jo), j - self.jo + 2*self.k_h + 1)
            return range(1 + self.jo - 2*self.k_h, self.m + 1), i
#        print('CASE4')
        return range(1 + self.jo - 2*self.k_h, self.m + 1), lambda j: range(1, j - self.jo + 2*self.k_h + 1)
    
#    def score(self, l1, l2):
#        return self.scoring_matrix[self.codes[l1]][self.codes[l2]]

    def d(self, i, j):
        ii = i-j+self.jo-self.io
        return ii, j  -self.J[0]+1# correction
    
    def di_to_i(self, di, j):
        return di+j-self.jo+self.io

    def align(self, s1, s2, d, k):
        self.init_domain(s1, s2, d, k)
        
        diag_data = np.zeros((1+2*self.k_h+1, 1+len(list(self.J))), dtype=int)
        diag_ptrs = np.zeros((1+2*self.k_h+1, 1+len(list(self.J))), dtype=int)

        loc_mx_crds = None
        for j in self.J:# go through each column
            for i in self.I(j):# go across band
                up   = diag_data[self.d(i-1, j)] + self.delete(self.s1[i-1]) if i-j >= self.du else float("-inf")
                left = diag_data[self.d(i, j-1)] + self.insert(self.s2[j-1]) if i-j <= self.dl else float("-inf")
                
                l = [0, diag_data[self.d(i-1, j-1)] + self.match(self.s1[i-1], self.s2[j-1]), up, left]
                # 0:, 1:\, 2:|, 3:_
                diag_data[self.d(i, j)] = max(l)
                diag_ptrs[self.d(i, j)] = np.argmax(l)

            _, j_ = self.d(i, j)
            mx = np.max(diag_data[:,j_])
            if mx > self.local_score:
                self.local_score = mx
                loc_mx_crds = [self.di_to_i(np.argmax(diag_data[:,j_]), j), j]

        # # find alignment:
        [i, j] = loc_mx_crds
        alignment = []
        dir_, prev_dir = -1, -1
        while dir_!=0:
            dir_ = diag_ptrs[self.d(i, j)]
            if prev_dir == 1:
                alignment.append( (i, j) )
            if dir_ == 1:
                i -= 1
                j -= 1
            if dir_ == 2: i -= 1
            if dir_ == 3: j -= 1
            prev_dir = dir_
        [p1, p2] = list(map(list, zip(*alignment)))
        p1.reverse()
        p2.reverse()
        return int(self.local_score), p1, p2
   