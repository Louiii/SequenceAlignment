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
