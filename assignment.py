from Qs import *
from Q3 import BLAST

from tests import *

def dynprog(alphabet, scoring_matrix, sequence1, sequence2):
    ''' 
    Input:
        alphabet:       set of letters in our alphabet.
        scoring_matrix: scores for each pair of letters in the alphabet.
        sequence1:      first sequence to align.
        sequence2:      second sequence to align.

    Output:
        list[
            int:        the score of the best local alignment found by the algorithm 
            list[int]:  one for each input sequences, that realise the matches/mismatches in the alignment.
            list[int]: 
        ]
    '''
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



print("Quadratic Memory\n")
for args, out in tests:
    [score, matches1, matches2] = dynprog(*args)
    check_indices(*(args+[matches1, matches2, score]))
#    print((score, matches1, matches2))
    print("Success score!") if out[0]==score else print('Failed score!')
    print("Success p1!") if out[1]==matches1 else print('Failed p1!')
    print("Success p2!\n\n") if out[2]==matches2 else print('Failed p2!\n')

print("\n\nLinear Memory\n")
for args, out in tests:
    [score, matches1, matches2] = dynproglin(*args)
    check_indices(*(args+[matches1, matches2, score]))
    print("Success score!") if out[0]==score else print('Failed score!\nCorrect: ',str(out[0]),'\nMine:   ',str(score))
    print("Success p1!") if out[1]==matches1 else print('Failed p1!\nCorrect:',str(out[1]),'\nMine:   ',str(matches1))
    print("Success p2!\n\n") if out[2]==matches2 else print('Failed p2!\nCorrect:',str(out[2]),'\nMine:   ',str(matches2)+'\n\n')
#
#print("\n\nBLAST\n")
#for args, out in tests:
#    [score, matches1, matches2] = heuralign(*args)
#    check_indices(*(args+[matches1, matches2, score]))
#    print("Success score!") if out[0]==score else print('Failed score!\nCorrect: ',str(out[0]),'\nMine:   ',str(score))
#    print("Success p1!") if out[1]==matches1 else print('Failed p1!\nCorrect:',str(out[1]),'\nMine:   ',str(matches1))
#    print("Success p2!\n\n") if out[2]==matches2 else print('Failed p2!\nCorrect:',str(out[2]),'\nMine:   ',str(matches2)+'\n\n')
