import numpy as np

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

def check_indices(alphabet, scoring_matrix, sequence1, sequence2, path1, path2, score_to_check=None):
    scoring_matrix = np.array(scoring_matrix)
    index = {i: alphabet.index(i) for i in alphabet}
    index['_'] = len(alphabet)

    score = 0
    x, y = path1[0], path2[0]
    
#    aligned1, aligned2 = sequence1[x], sequence2[y]
    for c in range(len(path1)):
        i = path1[c]
        j = path2[c]
        
#        aligned1 += sequence1[x]
#        aligned2 += sequence2[y]

        while x < i:
            score += scoring_matrix[index[sequence1[x]], index['_']]
            x += 1
            
#            if x == i: break
#            aligned1 += sequence1[x]
#            aligned2 += '_'

        while y < j:
            score += scoring_matrix[index[sequence2[y]], index['_']]
            y += 1
            
#            if y == j: break
#            aligned2 += sequence2[y]
#            aligned1 += '_'

        score += scoring_matrix[index[sequence1[i]], index[sequence2[j]]]
        
        x += 1
        y += 1
        
    
#    print('Aligned sequences:\n',str(aligned1),'\n', str(aligned2))
    if score_to_check is not None:
        print(str(score), ' is the correct score!' if score == score_to_check else ' is NOT the correct score!')