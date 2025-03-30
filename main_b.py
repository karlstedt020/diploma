import mne
import ctypes
from datetime import datetime, timedelta
import numpy as np
from sklearn.metrics import mean_absolute_percentage_error as mape, mean_squared_error as mse
from matplotlib import pyplot as plt
from argparse import ArgumentParser
import wfdb

from generate_cheb import generate


def read_edf(name):
    record = mne.io.read_raw_edf(NAME, preload=True)
    info = record.info
    channels = record.ch_names
    record_1, times=record.get_data(return_times=True, picks='sig')
    return record_1

def read_mit(name):
    record = wfdb.rdrecord(name)
    signals = record.p_signal
    fields = record.__dict__
    return np.transpose(signals)

NAME = 'Rh10010.edf'
NAME2 = 'rec_1'
N = 32
M = 3

parser = ArgumentParser()
parser.add_argument("-N", default=32, type=int)
parser.add_argument("-M", default=32, type=int)
args = parser.parse_args()
N = args.N
M = min(N, args.M)
chebyshev_matrix = generate(N)

ls = []
for arr in chebyshev_matrix:
    ls.extend(arr)
chebyshev_matrix = ls
print("Chebyshev matrix ready")

#record_1 = read_mit(NAME2)
record_1 = read_edf(NAME)
f = ctypes.CDLL("./libtest.so").solve_c
f.restype = ctypes.POINTER(ctypes.c_double)
f.argtypes = (ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.c_int,ctypes.c_int)
numbers_array = ctypes.c_double * N
chebyshev_array = ctypes.c_double * (N * N)
print("Total records:", len(record_1[0, :]))

total_time = 0
measures = len(record_1[0, :])
total_err = 0
total_avg_err = 0
deviants = list()
total_mse = 0
total_mape = 0
total_M = 0
total_N = 0
reconstructed = []
for batch_i in range(0, measures, N):
    if batch_i+N > 5000:
        break
    if (batch_i % 10000) // 1000 != ((batch_i - N) % 10000) // 1000:
        print(batch_i)
    mean = record_1[0, batch_i:batch_i+N].mean()
    slice_ = record_1[0, batch_i:batch_i+N].copy()
    min_1 = slice_.min()
    max_1 = slice_.max()
    if abs(min_1 - max_1) > 1.:
        print(max_1, min_1)
        currM = N // 2
        print("UES")
    else:
        currM = M
    total_M += currM
    total_N += N
    if not np.allclose([max_1 - min_1], [0]):
        slice_ = (slice_ - min_1) / (max_1 - min_1)
    numbers = numbers_array(*list(slice_))
    cheb = chebyshev_array(*chebyshev_matrix)
    start = datetime.now()
    res = f(numbers, cheb, ctypes.c_int(N), ctypes.c_int(currM))
    end = datetime.now()
    indeces = []
    for i in range(N):
        indeces.append(res[i])
    indeces = np.array(indeces)
    #indeces = indeces * (max_ - min_) + min_
    indeces = indeces * (max_1 - min_1) + min_1
    reconstructed.extend(indeces)
    print("Batch", batch_i // 32, "MAPE:", mape_res := mape(y_true=record_1[0, batch_i:batch_i+N], y_pred=indeces), "MSE:", mse_res := mse(y_true=record_1[0, batch_i:batch_i+N], y_pred=indeces))
    #if mape_res >= 20:
        #deviants.append({"original": record_1[0, batch_i:batch_i+N], "restored": indeces})
        #fig, ax = plt.subplots()
        #ax.plot(np.arange(N), record_1[0, batch_i:batch_i+N], linestyle='-', label="Original signal")
        #ax.plot(np.arange(N), indeces, linestyle='-', label="Reconstructed signal")
        #ax.legend()
        #plt.tight_layout()
        #plt.show()
    total_mape += mape(y_true=record_1[0, batch_i:batch_i+N], y_pred=indeces)
    total_mse += mse(y_true=record_1[0, batch_i:batch_i+N], y_pred=indeces)
    total_time += (end-start).total_seconds()


avg_time = total_time / (measures // N)
print("Average batch time:", avg_time)
print("Average MAPE:", total_mape / (measures // N))
print("Average MSE:", total_mse / (measures // N))
print("Total N:", total_N, "")
print("Total M:", total_M, "")
fig, ax = plt.subplots()
ln = len(reconstructed)
ax.plot(np.arange(ln), record_1[0, :ln], linestyle='-', label="Original signal")
ax.plot(np.arange(ln), reconstructed, linestyle='-', label="Reconstructed signal")
ax.legend()
plt.tight_layout()
plt.show()