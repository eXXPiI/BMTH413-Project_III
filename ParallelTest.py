
"""
Preamble
# Program: Name
# Author: Jonathan Myers
# Date: Mon Jan 27 03:06:29 2020
# Purpose: None.
# Arguments: None.
# Loads: None.
# Calls: None.
# Returns: None.
"""

import numpy as np
import multiprocessing as mp

def howmany_within_range(row, minimum, maximum):
    """Returns how many numbers lie within `maximum` and `minimum` in a given `row`"""
    count = 0
    for n in row:
        if minimum <= n <= maximum:
            count = count + 1
    return count

np.random.RandomState(100)
arr = np.random.randint(0, 10, size=[100, 5])
data = arr.tolist()

pool = mp.Pool(2)

# Step 2: `pool.apply` the `howmany_within_range()`
results = [pool.apply(howmany_within_range, args=(row, 4, 8)) for row in data]

# Step 3: Don't forget to close
pool.close()

print(results[:100])

