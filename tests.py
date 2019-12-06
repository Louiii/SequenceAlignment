import numpy as np
from LocalSearch_DynProg_Hirshberg_Heuristic import *

a1, m1 = "ABC", [[1,-1,-2,-1],[-1,2,-4,-1],[-2,-4,3,-2],[-1,-1,-2,0]]
a2, m2 = "ABCD",[[ 1,-5,-5,-5,-1],[-5, 1,-5,-5,-1],[-5,-5, 5,-5,-4],[-5,-5,-5, 6,-4],[-1,-1,-4,-4,-9]]
# print('m; \n'+str(np.array(m1)))
tests = [([a1, m1, "AABBAACA", "CBACCCBA"], (5, [3,5,6], [1,2,3]) ),
         # ([a1, m1, "ABCACA", "BAACB"], (None,  [1,3,4,5], [0,1,3,4]) ),
         ([a2, m2, "AAAAACCDDCCDDAAAAACC", "CCAAADDAAAACCAAADDCCAAAA"], (39, [5, 6, 7, 8, 9, 10, 11, 12, 18, 19], [0, 1, 5, 6, 11, 12, 16, 17, 18, 19])),
         ([a2, m2, "AACAAADAAAACAADAADAAA", "CDCDDD"], (17, [2, 6, 11, 14, 17], [0, 1, 2, 3, 4])),
         ([a2, m2, "DDCDDCCCDCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCDDDCDADCDCDCDCD", "DDCDDCCCDCBCCCCDDDCDBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBDCDCDCDCD"], (81, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 40, 41, 42, 43, 44, 45, 46, 47, 48, 50, 51, 52, 53, 54, 55, 56, 57, 58], [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 61, 62, 63, 64, 65, 66, 67, 68, 69])),
         ([a2, m2, "DDCDDCCCDCBCCCCDDDCDBBBBBBBBBBBBBBBBBBBBAAAAAAADDDDDDCCCACADCAAAAAAAADCADDDCCCCCCCCCAAAACDDDAAAADDCADCADBBBBBBBBBBBBBBBBBBBBBDCDCDCDCD", "DDCDDCCCDCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCDDDCDADCDCDCDCDAAAAACDCCCCCCDDDDDDDAAAAAAACCCDDDDADCADCCCDAAAADCADCCCCCCDAAAADCDCDADDDDCCCADCADCACCCCC"], (150, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 57, 59, 65, 66, 67, 68, 69, 70, 71, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 87, 88, 89, 90, 91, 96, 97, 98, 99, 100, 101, 102, 103], [48, 50, 51, 52, 54, 55, 57, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 88, 89, 90, 91, 92, 94, 95, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 121, 122, 123, 124, 125, 127, 128, 129, 130, 133, 134, 135, 136, 137, 138])),
         ]

def generateDNA(alphabet, length):
    # generate two random sequences from the alphabet
    a, b = np.random.choice(range(len(alphabet)), length), np.random.choice(range(len(alphabet)), length)
    a, b = [alphabet[int(i)] for i in a], [alphabet[int(i)] for i in b]
    
    # find places to copy one into the other
    subsets = np.sort(list(set(list(np.random.choice(list(range(length)), length//8)))))
    starts, ends = [s for i, s in enumerate(subsets) if i%2==0], [s for i, s in enumerate(subsets) if i%2==1]
    starts, ends = starts[:min(len(starts), len(ends))], ends[:min(len(starts), len(ends))]
    
    # copy across the subsections
    for i1, i2 in zip(starts, ends):
        b[i1:i2] = a[i1:i2]
    return ''.join(a), ''.join(b)

def check_indices(alphabet, scoring_matrix, sequence1, sequence2, path1, path2):
    scoring_matrix = np.array(scoring_matrix)
    index = {i: alphabet.index(i) for i in alphabet}
    index['_'] = len(alphabet)

    score = 0
    x, y = path1[0], path2[0]

    for c in range(len(path1)):
        i = path1[c]
        j = path2[c]
        while x < i:
            score += scoring_matrix[index[sequence1[x]], index['_']]
            x += 1
        while y < j:
            score += scoring_matrix[index[sequence2[y]], index['_']]
            y += 1
        score += scoring_matrix[index[sequence1[i]], index[sequence2[j]]]
        
        x += 1
        y += 1
        
    return score



import time, sys

data = dict()# dp, ln, blast, checkdp, checkln, lndpcheck, checkblast, %blast
ap, sm = "ABCD", [[ 1,-5,-5,-5,-1],[-5, 1,-5,-5,-1],[-5,-5, 5,-5,-4],[-5,-5,-5, 6,-4],[-1,-1,-4,-4,-9]]
for length in [10]+[500*i for i in range(1, 20)]:#[20, 50, 100, 150, 200, 250, 400, 500, 1000, 2000, 3000, 10000][:-4]:
    a, b = generateDNA(ap, length)
#    print((a, b))
    print("SEQUENCE LENGTHS: ",str(length))
    
    sys.stdout.write('Running quad space...')
    start = time.time()
    scdp, p1, p2 = dynprog(ap, sm, a, b)
    end = time.time()
    dptime = end - start
    sys.stdout.write(str('Done in '+str(dptime)+' seconds\nRunning lin space...'))

    start = time.time()
    scln, pl1, pl2 = dynproglin(ap, sm, a, b)
    end = time.time()
    lntime = end - start
    sys.stdout.write(str('Done in '+str(lntime)+' seconds\nChecking indices...'))
    cdp = check_indices(ap, sm, a, b, p1, p2)
    cln = check_indices(ap, sm, a, b, pl1, pl2)
    print((scdp, scln, cdp, cln))
    toPrint1 = 'is correct.' if cdp==scdp else 'Failed!'
    toPrint2 = 'is correct.' if cdp==scln else 'Failed!'
    print('QuadTime: '+toPrint1+'\nLinSpace: '+toPrint2)
    
    sys.stdout.write('Running BLAST...')
    start = time.time()
    schr, ph1, ph2 = heuralign(ap, sm, a, b)
    end = time.time()
    bstime = end - start
    sys.stdout.write(str('Done in '+str(bstime)+' seconds\n'))
    ch = check_indices(ap, sm, a, b, ph1, ph2)
    print('BLAST paths are valid.') if ch==schr else print('BLAST paths are invalid!')
    print('BLAST score: ',str(schr), ', % of optimal: ', str(round(100*schr/scdp)),'%\n\n')

    dpcheck = True if cdp==scdp else False
    lncheck = True if cln==scln else False
    lndpcheck = True if cdp==scln else False
    blastcheck = True if ch==schr else False
    
    data[length] = (dptime, lntime, bstime, dpcheck, lncheck, lndpcheck, blastcheck, 100*schr/scdp)


import matplotlib.pyplot as plt
x = np.sort(np.array(list(data.keys())))
Y = np.array([[data[key][0], data[key][1], data[key][2]] for key in x])
plt.plot(x, Y[:, 0], 'rx--', label='Quad space')
plt.plot(x, Y[:, 1], 'gx--', label='Linear space')
plt.plot(x, Y[:, 2], 'bx--', label='BLAST')
plt.xlabel('Length of sequences')
plt.ylabel('Time of algorithms')
plt.title('Time profiling')
plt.legend()
plt.savefig('TimeComplexity', dpi=400)
plt.show()


print('checks')
print([[v[3], v[4], v[5], v[6]] for k, v in data.items()])

print('perc optimality')
print([v[7] for k, v in data.items()])




