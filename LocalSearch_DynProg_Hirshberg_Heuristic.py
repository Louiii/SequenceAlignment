import numpy as np

'''
Part 1
Q1
'''
class Alignment(object):
    def __init__(self, alphabet, scoring_matrix):
        self.codes = { l:i for i, l in enumerate(alphabet)}
        self.codes.update({'_':len(alphabet)})
        self.scoring_matrix = scoring_matrix

        self.sc    = lambda a, b: self.scoring_matrix[self.codes[a]][self.codes[b]]
        self.match = lambda a, b: self.sc(a, b)
        self.delete = lambda a : self.sc(a, '_')
        self.insert = lambda a : self.sc('_', a)

    def total_score(self, aligned_seq_a, aligned_seq_b):
        return sum([self.sc(a, b) for a, b in zip(aligned_seq_a, aligned_seq_b)])

class QuadSpLocal(Alignment):
    def __init__(self, *args):
        super(QuadSpLocal, self).__init__(*args)
        self.local_score = float("-inf")
        self.ptrs = None
        
    def backtrack(self, loc_mx_crds):
        [i, j] = loc_mx_crds
        alignment = []
        dir_, prev_dir = -1, -1
        while dir_!=0:
            dir_ = self.ptrs[i, j]
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
        return [int(self.local_score), p1, p2]

    def align(self, s1, s2):
        n, m = len(s1), len(s2)
        matrix = np.zeros( (n+1, m+1) )
        self.ptrs   = np.zeros( (n+1, m+1), dtype=int)
        col = np.zeros(n+1)
        loc_mx_crds = None
        
        for j in range(1, m+1):# go through each column
            top =  0
            new_col = np.zeros(n+1)
            new_col[0] = top
            pointers = np.zeros(n+1, dtype=int)
            for i in range(1, n+1):
                l = [0, col[i-1] + self.match(s1[i-1], s2[j-1]), 
                     top + self.delete(s1[i-1]),  col[i] + self.insert(s2[j-1])]
                # 0:, 1:\, 2:|, 3:_
                pointers[i] = np.argmax(l)
                new_col[i] = l[pointers[i]]
                # print('i, j =',str((i,j)),', ', str(l), ', pointer: ', str(pointers[i]))
                top = new_col[i]
            col = new_col
            matrix[:, j] = col
            self.ptrs[:, j] = pointers
            mx = np.max(col)
            if mx > self.local_score:
                self.local_score = mx
                loc_mx_crds = [np.argmax(col), j]
        return self.backtrack(loc_mx_crds)
        
'''
Part 1
Q2
'''
class QuadTime(Alignment):
    def __init__(self, *args):
        super(QuadTime, self).__init__(*args)

    def score(self, l1, l2):
        return self.scoring_matrix[self.codes[l1]][self.codes[l2]]

    def iterate_col(self, j):
        top =  0
        new_col = np.zeros(self.n+1, dtype=int)
        for i in range(1,self.n+1):
            l = [0,
                 self.col[i-1] + self.match(self.s1[i-1], self.s2[j-1]), 
                 top + self.delete(self.s1[i-1]), 
                 self.col[i] + self.insert(self.s2[j-1])]
            # 0:.,1:\, 2:|, 3:_
            new_col[i] = max(l)
            top = new_col[i]
        self.col = new_col

    def local_max(self, s1, s2):
        self.s1, self.s2 = s1, s2
        self.n,  self.m = len(self.s1), len(self.s2)
        self.col = np.zeros(self.n+1, dtype=int)
        loc_mx_crds = None
        local_mx = 0
        for j in range(1, self.m+1):# go through each column
            self.iterate_col(j)
            mx = np.max(self.col)
            if mx > local_mx:
                local_mx = mx
                loc_mx_crds = [np.argmax(self.col), j]
        return local_mx, loc_mx_crds

class Needleman(Alignment):
    def __init__(self, *args):
        super(Needleman, self).__init__(*args)

    def compute_matrix(self, a, b, an, bn):
        self.matrix[:, 0] = np.cumsum([0]+[self.delete(l) for l in a])
        self.matrix[0, :] = np.cumsum([0]+[self.insert(l) for l in b])
        for i in range(1, an + 1):# rows
            for j in range(1, bn + 1):# cols
                self.matrix[i, j] = max(self.matrix[i-1, j-1] + self.match(a[i-1], b[j-1]), 
                                        self.matrix[i-1, j]   + self.delete(a[i-1]), 
                                        self.matrix[i, j-1]   + self.insert(b[j-1]))

    def backtrack(self, a, b, an, bn):
        aligned_a, aligned_b = [], []
        i, j = an, bn
        while i > 0 or j > 0:

            if j > 0 and self.matrix[i, j] == self.matrix[i, j - 1] + self.insert(b[j - 1]):
                aligned_a.insert(0, '_' * len(b[j - 1]))
                aligned_b.insert(0, b[j - 1])
                j -= 1

            elif i > 0 and self.matrix[i, j] == self.matrix[i - 1, j] + self.delete(a[i - 1]):
                aligned_a.insert(0, a[i - 1])
                aligned_b.insert(0, '_' * len(a[i - 1]))
                i -= 1

            elif i > 0 and j > 0 and self.matrix[i, j] == self.matrix[i - 1, j - 1] + self.match(a[i - 1], b[j - 1]):
                aligned_a.insert(0, a[i - 1])
                aligned_b.insert(0, b[j - 1])
                i -= 1
                j -= 1

            else:
                print()
        return aligned_a, aligned_b

    def align(self, a, b, semi_global=True):
        an, bn = len(a), len(b)
        self.matrix = np.zeros((1+an, 1+bn))
        self.compute_matrix(a, b, an, bn)
        return self.backtrack(a, b, an, bn)

class Hirschberg(Alignment):
    def __init__(self, *args):
        super(Hirschberg, self).__init__(*args)
        self.needleman = Needleman(*args)
        self.path_a, self.path_b = [], []

    def lin_space(self, a, b, reverse=False):
        # compute matrix but only store one row
        # we only need the final row as it contains the global scores to get to that point
        an, bn = len(a), len(b)
        cur_row = np.zeros(bn + 1, dtype=int)
        pre_row = np.cumsum([0]+[self.insert(b[j - 1]) for j in range(1, bn+1)])

        for i in range(1, an + 1):
            cur_row[0] = self.delete(a[i - 1]) + pre_row[0]
            for j in range(1, bn + 1):
                cur_row[j] = max(pre_row[j-1] + self.match(a[i-1], b[j-1]), 
                                 pre_row[j]   + self.delete(a[i-1]), 
                                 cur_row[j-1] + self.insert(b[j-1]))
            pre_row = cur_row
            cur_row = np.zeros(bn + 1, dtype=int)
        
        if reverse: return pre_row[::-1]
        return pre_row

    def recurse(self, a, b, i1, j1, i2, j2):
        aligned_a, aligned_b = [], []
        an, bn = len(a), len(b)

        if an == 0:# a has nothing in, but is being aligned with b
            for i in range(bn):
                aligned_a.append('_' * len(b[i]))
                aligned_b.append(b[i])
        elif bn == 0:# similarly to above
            for i in range(an):
                aligned_a.append(a[i])
                aligned_b.append('_' * len(a[i]))
        elif an == 1:
            # global alignments for small problem
            aligned_a, aligned_b = self.needleman.align(a, b)
            # append to the paths for the alignment
            self.path_a.append(i1+len(aligned_a)-aligned_a.count('_')-1)
            self.path_b.append(j1+len(aligned_b)-aligned_b.count('_')-1)
        else:
            mid_a = an // 2

            # global align the all possible sub problems (for all j in b) 
            # (the array returned is a row with the global scores for each point)
            left_part  = self.lin_space(a[:mid_a], b)
            right_part = self.lin_space(a[mid_a:][::-1], b[::-1], reverse=True)
            
            # maximise over the two sub problems
            mid_b = np.argmax(left_part + right_part)

            aligned_a_left,  aligned_b_left  = self.recurse(a[:mid_a], b[:mid_b], i1, j1, i2+mid_a, j2+mid_b)
            aligned_a_right, aligned_b_right = self.recurse(a[mid_a:], b[mid_b:], i1+mid_a, j1+mid_b, i2, j2)
            
            # reconstruct alignment
            aligned_a = aligned_a_left + aligned_a_right
            aligned_b = aligned_b_left + aligned_b_right

        return aligned_a, aligned_b

    def align(self, a, b):
        self.seq_a, self.seq_b = a, b
        self.len_a, self.len_b = len(a), len(b)
        return self.recurse(self.seq_a, self.seq_b, 0, 0, self.len_a-1, self.len_b-1)

'''
Part 1
Q3
'''

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

class BLAST:
    def __init__(self, alphabet, scoring_matrix, seq1, seq2, word_length, bandwidth):
        self.codes = { l:i for i, l in enumerate(alphabet)}
        self.codes.update({'_':len(alphabet)})
        self.alp, self.sc_mx = alphabet, scoring_matrix
        self.seq1, self.seq2 = seq1, seq2
        self.n, self.m = len(seq1), len(seq2)
        self.k, self.bandwidth = word_length, bandwidth
        
    def score(self, s1, s2):
        ''' Assumes len(s1)==len(s2) . Need to change this! '''
        return sum([self.sc_mx[self.codes[a]][self.codes[b]] for a, b in zip(s1, s2)])
    
    def FASTA(self, diags):
        ''' 
        order diagonals sensibly and do BandedDP search and return highest scoring alignment
        '''
        b = BandedDP(self.alp, self.sc_mx)
        mxsc, mp1, mp2 = -1, [], []
        for d in diags:
            [sc, p1, p2] = b.align(self.seq1, self.seq2, d, self.bandwidth)
            if sc > mxsc:
                mxsc, mp1, mp2 = sc, p1, p2
        return [mxsc, mp1, mp2]
    
    def makeSeq2Lookup(self, k):
        self.seq1Lookup = {}
        for i in range(len(self.seq1)+1-k):
            word = self.seq1[i:i+k]
            if word in self.seq1Lookup:
                self.seq1Lookup[word].append(i)
            else:
                self.seq1Lookup[word] = []
    
    def makeSeq1Lookup(self, k):
        self.seq2Lookup = {}
        for i in range(len(self.seq2)+1-k):
            word = self.seq2[i:i+k]
            if word in self.seq2Lookup:
                self.seq2Lookup[word].append(i)
            else:
                self.seq2Lookup[word] = []
                
    def make_neighbours_database(self):
        self.neighbours = {}
        for i in range(len(self.seq1)+1-self.k):
            word = self.seq1[i:i+self.k]
            scores = [(self.score(word, w), w) for w in self.dictionary]
            scores = sorted(scores, key=lambda x: x[0])[-30:]# goes from low to high
            self.neighbours[(i, word)] = list(zip(*scores[::-1]))[1]
            
    def find_seeds(self):
        self.seeds = {diag:[] for diag in range(-self.n, self.m)}
        for (i, word), neighbours in self.neighbours.items():# linear in seq1
            for neighbour in neighbours:# constant time, length of neighbours is 50
                if neighbour in self.seq2Lookup:
                    for j in self.seq2Lookup[neighbour]:# sub linear
                        self.seeds[i - j].append( (i, j) )
                        
    def prune_seeds(self):
        lengths = [(diag, len(seeds)) for diag, seeds in self.seeds.items()]
        lengths = sorted(lengths, key=lambda x: x[1])[-20:]# top 20 diagonals
        lengths = [d for (d, _) in lengths]
        selected_seeds = {}
        for diag in lengths:
            selected_seeds[diag] = self.seeds[diag]
        self.seeds = selected_seeds
        
    def highScPs(self):
        # extend each neighbour
        self.high_scoring_pairs = []
        diags = {diag:0 for diag in self.seeds.keys()}
        for diag, sds in self.seeds.items():
            if len(sds) > 0:
                final_i = sds[0]
                for i in sds:#, j in sds:
                    if i > final_i:# dont run over seeds that have already been expanded
                        j = i - diag
                        # extend seed
                        seed_start_i = i
                        seed_end_i   = i + self.k
                        
                        query_start_i = j
                        query_end_i   = j + self.k
                        
                        # seed is in seq1
                        (i1, j1) = (seed_start_i, query_start_i)
                        sc = self.score(self.seq1[seed_start_i:seed_end_i], self.seq2[query_start_i:query_end_i])
                        pr_sc = sc
                        while i1 > 0 and j1 > 0:
                            # go up
                            i1 -= 1
                            j1 -= 1
                            sc = self.score(self.seq1[i1:seed_end_i], self.seq2[j1:query_end_i])
                            delta_score = sc - pr_sc
                            if delta_score <= 0: 
                                sc = pr_sc
                                (i1, j1) = (i1+1, j1+1)
                                break
                            pr_sc = sc
                        
                        (i2, j2) = (seed_end_i, query_end_i)
                        sc = self.score(self.seq1[i1:seed_end_i], self.seq2[j1:query_end_i])
                        pr_sc = sc
                        while i2 < len(self.seq1)-1 and j2 < len(self.seq2)-1:
                            # go down
                            i2 += 1
                            j2 += 1
                            sc = self.score(self.seq1[i1:i2], self.seq2[j1:j2])
                            delta_score = sc - pr_sc
                            if delta_score <= 0: 
                                sc = pr_sc
                                (i2, j2) = i2-1, j2-1
                                break
                            pr_sc = sc
                            
                            
                        diags[diag] += sc
                        final_i = i2
        return [(diag, score) for diag, score in diags.items()]
            
    def find_seeds_heuristic(self, threshold=5000):
        '''
        This fn computes the diagonals with the most seeds on.
        To keep the whole program linear time and space it samples above
        a certain threshold, which can be thought of as a resolution for
        large sequences (~> 20000).
        '''
        diagonals = {diag:[] for diag in range(-self.n, self.m)}
        count = {diag:0 for diag in range(-self.n, self.m)}
        for poss_word in self.dictionary:# constant time
            if poss_word in self.seq1Lookup and poss_word in self.seq2Lookup:
                l1, l2 = len(self.seq1Lookup[poss_word]), len(self.seq2Lookup[poss_word])
                if l1*l2 > threshold:# sample
                    seen = set()
                    for _ in range(threshold):
                        x, y = np.random.randint(l1), np.random.randint(l2)
                        if (x, y) not in seen:
                            seen.add((x, y))
                            diagonals[x-y].append( x )
                            count[x-y] += 1
                else:
                    for i in self.seq1Lookup[poss_word]:
                        for j in self.seq2Lookup[poss_word]:
                            diagonals[i-j].append( i ) 
                            count[i-j] += 1
        counts = [(d, c) for d, c in count.items()]
        counts = sorted(counts, key=lambda x: x[1])[-30:]
        diags = [d for (d, c) in counts]
        self.seeds = {}
        for d in diags:
            xs = diagonals[d]
            xs = sorted(xs)
            self.seeds[d] = xs
            
    def align(self):
        self.makeSeq2Lookup(self.k)
        self.dictionary = product(self.k, self.alp, [''])
        self.makeSeq1Lookup(self.k)
        self.find_seeds_heuristic()
        self.prune_seeds()
        diags = self.highScPs()
        diags = sorted(diags, key=lambda x: x[1])[-10:]
        diags = [diag for (diag, _) in diags]

        return self.FASTA(diags)

    def time_break_down_align(self):
        import time
        now = time.time()
        self.times = [now]
        self.makeSeq2Lookup(self.k)# time index 1
        now = time.time()
        self.times.append(now)
        self.dictionary = product(self.k, self.alp, [''])# time index 2
        now = time.time()
        self.times.append(now)
        self.makeSeq1Lookup(self.k)# time index 1
        now = time.time()
        self.times.append(now)
        self.find_seeds_heuristic()
#        self.find_seeds()# time index 4      <------------
        now = time.time()
        self.times.append(now)
        self.prune_seeds()# time index 5      <------------
        now = time.time()
        self.times.append(now)
        diags = self.highScPs()# time index 6      <------------
        now = time.time()
        self.times.append(now)
        diags = sorted(diags, key=lambda x: x[1])[-10:]
        diags = [diag for (diag, _) in diags]
        
        thing = self.FASTA(diags)# time index 7      <------------
        now = time.time()
        self.times.append(now)
        t0 = self.times[0]
        self.times = [t-t0 for t in self.times]
#        print(self.times)
        return thing
      
def product(k, alphabet, current):
    if k==0:return current
    return product(k-1, alphabet, [c+a for c in current for a in alphabet])

''' 
Wrappers:
'''

def dynprog(alphabet, scoring_matrix, sequence1, sequence2):
    ql = QuadSpLocal(alphabet, scoring_matrix)
    return ql.align(sequence1, sequence2)

def dynproglin(alphabet, scoring_matrix, sequence1, sequence2):
    
    def localMaxCoords(sequence1, sequence2, alphabet, scoring_matrix):
        q = QuadTime(alphabet, scoring_matrix)
        return q.local_max(sequence1, sequence2)
    
    # reduce to global problem:
    mx, [i_end, j_end] = localMaxCoords(sequence1, sequence2, alphabet, scoring_matrix)
    sequence1, sequence2 = sequence1[:i_end][::-1], sequence2[:j_end][::-1]
    mx, [i_start_r, j_start_r] = localMaxCoords(sequence1, sequence2, alphabet, scoring_matrix)
    
    subsequence1, subsequence2 = sequence1[:i_start_r][::-1], sequence2[:j_start_r][::-1]
    
    # global alignment:
    algo = Hirschberg(alphabet, scoring_matrix)
    a, b = algo.align(subsequence1, subsequence2)
    p1, p2 = algo.path_a, algo.path_b
    
    p1 = list(i_end - i_start_r + np.array(p1))
    p2 = list(j_end - j_start_r + np.array(p2))
    return [mx, p1, p2]

def heuralign(alphabet, scoring_matrix, sequence1, sequence2, word_length=3, bandwidth=18):
    if len(sequence1) >= len(sequence2):
        b = BLAST(alphabet, scoring_matrix, sequence1, sequence2, word_length, bandwidth)
        return b.align()
    else:
        b = BLAST(alphabet, scoring_matrix, sequence2, sequence1, word_length, bandwidth)
        [score, path1, path2] = b.align()
        return [score, path2, path1]

