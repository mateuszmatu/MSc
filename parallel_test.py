from joblib import Parallel, delayed
import math
import time


def printing(i):
    time.sleep(1)
    print(i)

s = time.perf_counter()
Parallel(n_jobs=8)(delayed(printing)(i) for i in range(100))
e = time.perf_counter()

print(e-s)

s = time.perf_counter()
for i in range(10):
    time.sleep(1)
    print(i)
e = time.perf_counter()
print(e-s)

