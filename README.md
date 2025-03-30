Shell command for building C++ library: 
```shell
g++ -fPIC -shared -O2 -o libtest.so saved_abc.cpp abc_extern.cpp
```

Shell command for launching script: 
```shell
python3 main_b.py -N 50 -M 25
```

## Integration plan

Application, after receiving full ECG data and putting it in a temporary storage, calls for compressing application (maybe several instances of it). Compressed data is put in another clickhouse database, then standard compression methods are applyed to reduce the size even further.
