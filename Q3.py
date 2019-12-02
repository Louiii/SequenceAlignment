import numpy as np
from BDP import BandedDP

class Word:
    def __init__(self, w, i):
        self.w = w
        self.ind = i
        
class Point:
    def __init__(self, i, j):
        self.i = i
        self.j = j

class Hsp:
    def __init__(self, score, word1, word2, start, end):
        self.score = score
        self.word1 = word1
        self.word2 = word2
        self.start = start
        self.end = end
    
def visu(hs):
    import matplotlib.pyplot as plt
    
    crds = lambda x, y, l: [(x+i, y+i) for i in range(l)]
    all_crds = []
    for h in hs: all_crds += crds(*h)
    xs, ys = zip(*all_crds)
    m = np.zeros((max(xs)+1, max(ys)+1))
    for (x, y) in all_crds: m[x, y] += 1
    
    plt.imshow(m)
    plt.colorbar();


'''
[20 Marks] A Heuristic procedure that runs in sub-quadratic time (similar to FASTA and BLAST)

will be assessed on trade-offs between running time and quality 
of output. A typical input will come with a planted alignment, which consists 
of segments of matches of different lengths separated by some random stuff 
(so that I will know that there is an alignment, which is at least as good as
the planted one). 
'''

class BLAST:
    def __init__(self, alphabet, scoring_matrix, seq1, seq2, word_length, bandwidth):
        self.alp, self.sc_mx = alphabet, scoring_matrix
        self.seq1, self.seq2 = seq1, seq2
        self.n, self.m = len(seq1), len(seq2)
        self.k, self.bandwidth = word_length, bandwidth
        self.database_seq = self.make_table(seq1, word_length)
    
    def score(self, s1, s2):
        ''' Assumes len(s1)==len(s2) . Need to change this! '''
        return sum([3 if s1[i]==s2[i] else -4 for i in range(len(s1))])
    
    def make_table(self, seq, k):
        ''' list<Word> '''
        return [Word(seq[i:i+k], i) for i in range(len(seq)+1-k)]
    
    def _50_high_scoring_k_tuples(self, query_word):
        ''' (ordered) top 50 scoring k tuples in database associated with query word '''
        scores = [(self.score(word.w, query_word), word.ind) for word in self.database_seq]
        scores = sorted(scores, key=lambda x: x[0])[-50:]
        return list(zip(*scores))[1]# indices in seq1 of top scoring k-tuples
    
    def HSP(self):
        # extend each neighbour
        high_scoring_pairs = []
        # high-scoring-pairs (use some datastructure to store HSPs)
        for query, neighbours in self.query_neighbourhood_words.items():
            query_start_i = query.ind# query_database[query]
            query_end_i   = query_start_i + self.k
#            print('Query word: ', query)
            for seed in neighbours:# indices of words in seq1
                # extend seed
                seed_start_i = seed# database_seq[seed]
                seed_end_i   = seed_start_i + self.k
                
                # seed is in seq1
                (i1, j1) = (seed_start_i, query_start_i)
                sc = self.score(self.seq1[seed_start_i:seed_end_i], self.seq2[query_start_i:query_end_i])
                pr_sc = sc
                while i1 > 0 and j1 > 0:
                    # go up
                    i1 -= 1
                    j1 -= 1
                    sc = self.score(self.seq1[i1:seed_end_i], self.seq2[j1:query_end_i])
#                    print('up sc, s1, s2: ',str((sc, self.seq1[i1:seed_end_i], self.seq2[j1:query_end_i])))
                    delta_score = sc - pr_sc
        #            print(delta_score)
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
#                    print('dn sc, s1, s2: ',str((sc, self.seq1[i1:i2], self.seq2[j1:j2])))
                    delta_score = sc - pr_sc
        #            print(delta_score)
                    if delta_score <= 0: 
                        sc = pr_sc
                        (i2, j2) = i2-1, j2-1
                        break
                    pr_sc = sc
                high_scoring_pairs.append( (sc, self.seq1[i1:i2], self.seq2[j1:j2], i1, j1) )
#                print('done seed: ',seed,', ',str(high_scoring_pairs[-1]))
            
        high_scoring_pairs = list(set(high_scoring_pairs))
        high_scoring_pairs = sorted(high_scoring_pairs, key=lambda x: x[0])
        
        return high_scoring_pairs
    
    def FASTA(self, hsps, num_to_search=10):
        ''' 
        order diagonals sensibly and do BandedDP search and return highest scoring alignment
        '''
    #    # sort diagonals (i-j), had high_scoring_pairs
    #    diag = dict()
    #    for h in hsps:
    #        k = h[3]-h[4]# i-j
    #        if k in diag:
    #            diag[k].append( h )
    #        else:
    #            diag[k] = [h]
    #    print('\nDiagonals\n')
    #    print(diag)
        
    #    while True:
    #        curr_best = hsps.pop()
    #        
    #        # find possible connections
    #        check_sw = lambda p1, p2: p1.end.j >= p2.start.j and p1.end.i >= p2.start.i
    #        options = []
    #        for h in hsps:
    #            if check_sw(curr_best, h):
        diagonals = dict()#np.zeros(self.n+self.m, dtype=int)
        for h in hsps:
            ''' there may be lots of hsp on one diagonal because they are all from the same section 
            may want to find a way to avoid this
            '''
            d = h.start.i - h.start.j
            if d in diagonals:
                diagonals[d] += h.score
            else:
                diagonals[d] = h.score
        diagonals = [(d, s) for d, s in diagonals.items()]
        
        # highest to lowest scoring diagonals
        ''' could do based on how different too! '''
        diagonals = sorted(diagonals, key=lambda x: x[1])[::-1][:num_to_search]
        
        diagonals_to_search = tuple(zip(*diagonals))[0]
        print('diagonals_to_search: ',str(diagonals_to_search))
        
    
        b = BandedDP(self.alp, self.sc_mx)
        mxsc, mp1, mp2 = -1, [], []
        for d in diagonals_to_search:
            print('searching diagonal...',str(d))
#            b = BandedDP(self.alp, self.sc_mx, self.seq1, self.seq2, d, self.bandwidth)
#            [sc, p1, p2] = b.align()
            [sc, p1, p2] = b.align(self.seq1, self.seq2, d, self.bandwidth)
            
            print('Got score: ', str(sc))
            if sc > mxsc:
                mxsc, mp1, mp2 = sc, p1, p2
        return [mxsc, mp1, mp2]
        
    def align(self):
        self.query_database = self.make_table(self.seq2, self.k)
        self.query_neighbourhood_words = {query:self._50_high_scoring_k_tuples(query.w) for query in self.query_database}
        
        
        self.high_scoring_pairs = self.HSP()
        
        pts = []
        for h in reversed(self.high_scoring_pairs[-50:]):# change the -50
            pts.append( Hsp(h[0], h[1], h[2], Point(h[3], h[4]), Point(h[3]+len(h[1]), h[4]+len(h[2]))) )
        
#        print(self.high_scoring_pairs)
#        print([(h[3],h[4]) for h in self.high_scoring_pairs])
        visu([(h[3],h[4], len(h[2])) for h in self.high_scoring_pairs])
        
        return self.FASTA(pts)

#seq1 = 'ACBAAAAACCAAAAB'#   'ABCDAAABCD'
#seq2 = 'DAAAAADAAAACCCD'#   'DAAABDD'
#ap = 'ABCD'
#m = [[ 1,-5,-5,-5,-1],[-5, 1,-5,-5,-1],[-5,-5, 5,-5,-4],[-5,-5,-5, 6,-4],[-1,-1,-4,-4,-9]]
#k, bandwidth = 3, 8
#
#b = BLAST(ap, m, seq1, seq2, k, bandwidth)
#[mxsc, mp1, mp2] = b.align()
#
#
#from tests import score_alignment_local
#
#st1, st2, sc = score_alignment_local(ap, m, seq1, seq2, mp1, mp2)
#
#print('score: ',str(sc))
#print('alignment:\n',str(st1),'\n',str(st2))
#
