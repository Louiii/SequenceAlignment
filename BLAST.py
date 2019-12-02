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

def score(s1, s2):
    ''' Assumes len(s1)==len(s2) . Need to change this! '''
    return sum([3 if s1[i]==s2[i] else -4 for i in range(len(s1))])

def make_table(seq, k):
    ''' list<Word> '''
    return [Word(seq[i:i+k], i) for i in range(len(seq)+1-k)]

def _50_high_scoring_k_tuples(query_word):
    ''' (ordered) top 50 scoring k tuples in database associated with query word '''
    scores = [(score(word.w, query_word), word.ind) for word in database_seq]
    scores = sorted(scores, key=lambda x: x[0])[-50:]
    return list(zip(*scores))[1]# indices in seq1 of top scoring k-tuples

def HSP():
    # extend each neighbour
    high_scoring_pairs = []
    # high-scoring-pairs (use some datastructure to store HSPs)
    for query, neighbours in query_neighbourhood_words.items():
        query_start_i = query.ind# query_database[query]
        query_end_i   = query_start_i + k
        print('Query word: ', query)
        for seed in neighbours:# indices of words in seq1
            # extend seed
            seed_start_i = seed# database_seq[seed]
            seed_end_i   = seed_start_i + k
            
            # seed is in seq1
            (i1, j1) = (seed_start_i, query_start_i)
            sc = score(seq1[seed_start_i:seed_end_i], seq2[query_start_i:query_end_i])
            pr_sc = sc
            while i1 > 0 and j1 > 0:
                # go up
                i1 -= 1
                j1 -= 1
                sc = score(seq1[i1:seed_end_i], seq2[j1:query_end_i])
                print('up sc, s1, s2: ',str((sc, seq1[i1:seed_end_i], seq2[j1:query_end_i])))
                delta_score = sc - pr_sc
    #            print(delta_score)
                if delta_score <= 0: 
                    sc = pr_sc
                    (i1, j1) = (i1+1, j1+1)
                    break
                pr_sc = sc
            
            (i2, j2) = (seed_end_i, query_end_i)
            sc = score(seq1[i1:seed_end_i], seq2[j1:query_end_i])
            pr_sc = sc
            while i2 < len(seq1)-1 and j2 < len(seq2)-1:
                # go down
                i2 += 1
                j2 += 1
                sc = score(seq1[i1:i2], seq2[j1:j2])
                print('dn sc, s1, s2: ',str((sc, seq1[i1:i2], seq2[j1:j2])))
                delta_score = sc - pr_sc
    #            print(delta_score)
                if delta_score <= 0: 
                    sc = pr_sc
                    (i2, j2) = i2-1, j2-1
                    break
                pr_sc = sc
            high_scoring_pairs.append( (sc, seq1[i1:i2], seq2[j1:j2], i1, j1) )
            print('done seed: ',seed,', ',str(high_scoring_pairs[-1]))
        
    high_scoring_pairs = list(set(high_scoring_pairs))
    high_scoring_pairs = sorted(high_scoring_pairs, key=lambda x: x[0])
    
    return high_scoring_pairs

def FASTA(hsps):
    ''' 
        score hsps/gapless-alignments using eg. BLOSUM
        keep the best
        connect them using simple affine gap
    '''
    
    '''
        look for diagonals with many matching words on
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
#                
#            
    return
    
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

if __name__=="__main__":
    # should be AAAAACCAAAA
    #           AAAAAD_AAAA <- FASTA test
    seq1 = 'ACBAAAAACCAAAAB'#   'ABCDAAABCD'
    seq2 = 'DAAAAADAAAACCCD'#   'DAAABDD'
    k = 3
    database_seq = make_table(seq1, k)
    
    
    query_database = make_table(seq2, k)
    query_neighbourhood_words = {query:_50_high_scoring_k_tuples(query.w) for query in query_database}
    
    
    high_scoring_pairs = HSP()
    
    pts = []
    for h in reversed(high_scoring_pairs[-50:]):# change the -50
        pts.append( Hsp(h[0], h[1], h[2], Point(h[3], h[4]), Point(h[3]+len(h[1]), h[4]+len(h[2]))) )
    
    print(high_scoring_pairs)
    print([(h[3],h[4]) for h in high_scoring_pairs])
    visu([(h[3],h[4], len(h[2])) for h in high_scoring_pairs])
    
    x = FASTA(pts)