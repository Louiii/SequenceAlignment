import numpy as np

"""
Part 2
"""

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

class MyCostFn(Alignment):
    def __init__(self):
        ''' Assume alphabet of ABC. '''
        scoring_matrix = [[2,-2,-2,-3],[-2,2,-2,-3],[-2,-2,2,-3],[-3,-3,-3,0]]
        super(MyCostFn, self).__init__("ABC", scoring_matrix)
        self.local_score = float("-inf")
        self.ptrs = None
    
    def backtrack(self, loc_mx_crds):
        [i, j] = loc_mx_crds
        alignment = []
        dir_, prev_dir = -1, -1
#        print('backtracing...')
        while dir_ != 2:# 2==' '
            dir_ = self.ptrs[i, j]
            if prev_dir == 3:
                alignment.append( (i, j) )
            if dir_ == 3:
                i -= 1
                j -= 1
            if dir_ == 0: i -= 1
            if dir_ == 1: j -= 1
            prev_dir = dir_
        [p1, p2] = list(map(list, zip(*alignment)))
        p1.reverse()
        p2.reverse()
        return [int(self.local_score), p1, p2]
    
    def traceback(self, i, j, up=True):
        ''' return the length traced back in the indel directions '''
        count = 0
        direction = 0 if up else 1# left
        dir_ = -1
        while dir_ != 2:
            dir_ = self.ptrs[i, j]
            if dir_ == direction: 
                if up:
                    i -= 1
                else:
                    j -= 1
                count += 1
            else:
                return count

    def align(self, s1, s2):
        n, m = len(s1), len(s2)
        self.matrix =  np.zeros( (n+1, m+1), dtype=int)
        self.ptrs   = 2*np.ones( (n+1, m+1), dtype=int)
        loc_mx_crds = None
        
        for j in range(1, m+1):# go through each column
            top =  0
            for i in range(1, n+1):
                ''' 
                        for the delete option check how long the 'up' route is
                        for the insert option check how long the 'right' route is
                        
                        use the normal score function with an added option to add
                        the length of the route, append the score function so it
                        returns the appropriate score
                        
                        (from scoring matrix: 
                              _ = -3 )
                            
                            _ _ = -2
                          _ _ _ = -1
                        _ _ _ _ =  0
                '''
                # for case0 we need to check the up path
                count = self.traceback(i-1, j, up=True)
                pen = self.delete(s1[i-1]) if count==0 else -2 if count==1 else -1 if count==2 else 0
                case0 = top + pen
                # for case1 we need to check the left path
                count = self.traceback(i, j-1, up=False)
                pen = self.insert(s2[j-1]) if count==0 else -2 if count==1 else -1 if count==2 else 0
                case1 = self.matrix[i, j-1] + pen
                
                l = [case0, case1, 0, self.matrix[i-1, j-1] + self.match(s1[i-1], s2[j-1])]
                # 0:|, 1:_, 2:, 3:\
                self.ptrs[i, j] = np.argmax(l)
                self.matrix[i, j] = l[self.ptrs[i, j]]
                top = self.matrix[i, j]

            mx = np.max(self.matrix[:, j])
            if mx > self.local_score:
                self.local_score = mx
                loc_mx_crds = [np.argmax(self.matrix[:, j]), j]
                
        return self.backtrack(loc_mx_crds)
    
#a = 'CAABCCCCCCBACC'
#b = 'CAABBACC'
a = 'ABCABCABCCCCCCCBBBBBCCCAAACACABBBBBBBBB'
b = 'ABCABCBBBBBB'

m = MyCostFn()
[s, p1, p2] = m.align(a, b)


from tests import check_indices
print('a = ',str(a),'\nb = ',str(b),'\n')
check_indices("ABC", [[2,-2,-2,-3],[-2,2,-2,-3],[-2,-2,2,-3],[-3,-3,-3,0]], a, b, p1, p2, return_alignment=True)






