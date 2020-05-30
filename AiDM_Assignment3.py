import numpy as np
import scipy.sparse as sps

data = np.load('/Users/admin/Dropbox/A3_data/user_movie.npy')

num_users = max(data[:, 0])
num_movies = max(data[:, 1])
ones = np.repeat(1, len(data))
q = sps.csc_matrix((ones, (data[:, 1] - 1, data[:, 0] - 1)), shape=(num_movies, num_users))
np.random.seed(1)
def intersect(a, b):
    """ return the intersection of two lists """
    return list(set(a) & set(b))

def union(a, b):
    """ return the union of two lists """
    return list(set(a) | set(b))
# Calculates Jaccard Similarity. x and y are sets. example: set([1,2,3])
def jsim(x, y):
    return len(intersect(x,y))/1.0*len(union(x,y))#sum(x.intersection(y)) / sum(x.union(y))




def minhash(x, p):

    # finnd minhash for data in row according to permuation p.
    # permute rows of matrix
    y = x[p, :]
    y.sort_indices()
    #print ((y.indices[y.indptr[:-1]]))
    # find position of first i in each col
    return y.indices[y.indptr[:-1]]

#
def minhash_signature(x, k):
    # signatures  matrix of length K in matrix x
    [n_row, n_col] = x.shape
    # allocate memory
    y = np.zeros((k, n_col))
    for i in range(k):
        p = np.random.permutation(n_row)
        y[i, :] = minhash(x, p)
    return y

#
# V is the signature matrix
v = minhash_signature(q, 100)
print jsim(v[:,0],v[:1])
# print v.shape
# for j in np.arange(len(v[0,:])):
#    print jsim(v[:,j],v[:,j+1])
#
# LSH implementation
sim_list = [[] for x in range(num_users)]
band = 5
row = 5
buckets = 5000
sim_count = 0

for i in range(band):
    table = [[] for x in range(buckets)]
    # Sum the signatures of each user in the band and use that as the bucket number
    for j in range(5):
        sigs = 0
        for r in range(row):
            sigs += int(v[j][r + i * row])

            bucket = sigs
            table[bucket].append(v[j][r + i * row])
    # For each bucket we check the similarities of all the pairs within that bucket
    for b in range(buckets):
        lb = len(table[b])
        if lb > 1:
            for c, d in itertools.combinations(table[b], 2):
                print jsim(v[:, c], v[:, d])

#real_count = 0


# i

# # Function which tests the similarity of the pair of signatures and if it's >0.5 we store it in list sim_list
# def similarity(sig1, sig2, u1, u2):
#     com = np.sum(sig1 == sig2)
#     sim = float(com) / (2 * number_hash - com)
#     if sim > 0.5:
#         sim_list[u1].append(u2)


dd
real_count = 0

# Checking the real similarity and outputing the results
f1 = open('results.txt', 'w')
for i in range(len(sim_list)):
    if len(sim_list[i]) > 0:
        sim_list[i] = list(set(sim_list[i]))
        for j in sim_list[i]:
            col1 = q[:, i].nonzero()
            col2 = q[:, j].nonzero()
            intersect = len(np.intersect1d(col1, col2))
            total = intersect + len(np.setxor1d(col1, col2))

            real_sim = float(intersect) / total
            if real_sim > 0.5:
                real_count += 1
                f1.write(str(i + 1) + ", " + str(j + 1) + "\n")
f1.close()
