# -*- coding: utf-8 -*-
import numpy as np
from time import time
import residue_module as rem

# Generates all unique permutations of the rows
rows = rem.generate_residue_rows()
rperms = []
for r in rows:
    rp = rem.gen_perms(r)
    rperms.append(rp)

rperms_alt = []
for r in rows:
    for rp in rem.gen_perms(r):
        rperms_alt.append(rp)

def generate_residues():
    results = []
    # The last four entries in "rows" are excluded from being the top row, since they are all divisible by sqrt(2).
    # This prevents any matrices from being created that are completely divisble by sqrt(2).
    for i1 in range(len(rows)-4):
        print(f'Progress: {i1} / {(len(rows)-4)}')
        r1 = rows[i1]
        # row 2 can be any entry in "rows" that don't come before the previous row.
        # This prevents repeats and provides a major speed-up.
        for i2 in range(i1,len(rows)):
            # row 2 can be any (unique) permutation of rows[i2], as drawn from "rperms".
            # Note that "rperms" only needs to be generated once, which saves significant time.
            for r2 in rperms[i2]:
                # Checks if row 2 is orthogonal to row 1.
                if not rem.orthogonal(r1,r2):
                    continue
                # The same mechanism for generating row 2 is repeated for row 3, and so on.
                for i3 in range(i2,len(rows)):
                    for r3 in rperms[i3]:
                        if (not rem.orthogonal(r1,r3)) or (not rem.orthogonal(r2,r3)):
                            continue
                        for i4 in range(i3,len(rows)):
                            for r4 in rperms[i4]:
                                if (not rem.orthogonal(r1,r4)) or (not rem.orthogonal(r2,r4)) or (not rem.orthogonal(r3,r4)):
                                    continue
                                for i5 in range(i4,len(rows)):
                                    for r5 in rperms[i5]:
                                        if (not rem.orthogonal(r1,r5)) or (not rem.orthogonal(r2,r5)) or (not rem.orthogonal(r3,r5)) or (not rem.orthogonal(r4,r5)):
                                            continue
                                        res = np.array([r1,r2,r3,r4,r5])
                                        # The 6th row is fixed by the previous rows, so no permutations are required.
                                        # If there is no valid 6th row, then it is skipped.
                                        r6 = np.zeros(6)
                                        skip = False
                                        for i in range(6):
                                            l = rem.normalize(res[:,i])
                                            if l == -1:
                                                skip = True
                                            r6[i] = l
                                        res = np.append(res,[r6],axis=0)
                                        # Checks that the last row is orthonormal, and that all of the columns are orthonormal.
                                        if (not skip) & rem.normal(r6) & rem.new_row_orthogonal(res,r6) & rem.cols_orthogonal(res):
                                            results.append(res)
                                            # Checks that the newly entered residues is unique. If not, it is removed from "results".
                                            for prev in results[:-1]:
                                                if rem.equal(prev,res):
                                                    results.pop()
    return np.array(results)

start_time = time.time()
residues = generate_residues()
# rem.save_residues(residues,'residues.txt')
print('Residues generated: %d' % len(residues))
print('Time taken: %d s' % (time.time()-start_time))

# Checks whether the residues have paired rows.
# Sorts them into three lists: no pairs, one pair, multiple pairs.
# They are then stored as txt files.
no_pair_residues = []
one_pair_residues = []
multi_pair_residues = []
for res in residues:
    num_pairs = 0
    for i in range(5):
        for j in range(i+1,6):
            # Paired zero rows are excluded, as these are not useful for row reduction.
            if (np.allclose(res[i],res[j])) & (not np.allclose(res[i],np.zeros(6))):
                num_pairs += 1
    if num_pairs == 0:
        no_pair_residues.append(res)
    elif num_pairs == 1:
        one_pair_residues.append(res)
    else:
        multi_pair_residues.append(res)

no_pair_residues = np.array(no_pair_residues)
print(f'{len(no_pair_residues)} residues with no row pairings.')
# rem.save_residues(no_pair_residues,'no_pair_residues.txt')
one_pair_residues = np.array(one_pair_residues)
print(f'{len(one_pair_residues)} residues with one row pairing.')
# rem.save_residues(one_pair_residues,'one_pair_residues.txt')
multi_pair_residues = np.array(multi_pair_residues)
print(f'{len(multi_pair_residues)} residues with multiple row pairings.')
# rem.save_residues(multi_pair_residues,'multi_pair_residues.txt')
