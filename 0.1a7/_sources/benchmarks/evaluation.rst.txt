.. _evaluation:

================
Evaluation Study
================

For each run, we compute running times as follows: for each fold of a 10-fold cross-validation experiment, the running time of the kernel corresponds to the running time for the computation of the kernel matrix that performed best on the validation experiment.

We next present our experimental results.
As mentioned above, we evaluate the graph kernels by performing graph classification on unlabeled, node-labeled and node-attributed benchmark datasets.

Node-Labeled Graphs
^^^^^^^^^^^^^^^^^^^





Unabeled Graphs
^^^^^^^^^^^^^^^
The following Table illustrates the average CPU running times of the compared kernels on the 6 datasets that contain unlabeled graphs.


+---------+---------------+---------------+---------------+----------------+---------------+----------------+
|         | IMDB-B        | IMDB-M        | REDDIT-B      | REDDIT-M-5K    | REDDIT-M-12K  | COLLAB         |
+=========+===============+===============+===============+================+===============+================+
| VH      | 0.07s         | 0.15s         | 0.67s         | 2.20s          | 6.37s         | 1.12s          |
+---------+---------------+---------------+---------------+----------------+---------------+----------------+
| RW      | 7m 20.94s     | 13m 40.75s    | TIMEOUT       | TIMEOUT        | TIMEOUT       | 13h 38m 11.49s |
+---------+---------------+---------------+---------------+----------------+---------------+----------------+
| SP      | 11.51s        | 7.92s         | 4h 48m 11.19s | 12h 40m 19.50s | TIMEOUT       | 1h 9m 5.50s    |
+---------+---------------+---------------+---------------+----------------+---------------+----------------+
| GR      | 22m 45.89s    | 21m 44.30s    | 44m 45.42s    | 44m 6.52s      | 53m 14.22s    | 2h 58m 1.14s   |
+---------+---------------+---------------+---------------+----------------+---------------+----------------+
| WL-VH   | 4.49s         | 6.16s         | 16m 2.65s     | OUT-OF-MEM     | OUT-OF-MEM    | 38m 42.24s     |
+---------+---------------+---------------+---------------+----------------+---------------+----------------+
| WL-SP   | 1m 32.66s     | 1m 40.46s     | TIMEOUT       | TIMEOUT        | TIMEOUT       | 10h 27m 41.97s |
+---------+---------------+---------------+---------------+----------------+---------------+----------------+
| NH      | 21.83s        | 26.07s        | 23m 3.42s     | 2h 44m 44.66s  | 9h 11m 23.67s | 35m 49.96s     |
+---------+---------------+---------------+---------------+----------------+---------------+----------------+
| NSPDK   | 4m 18.12s     | 2m 49.45s     | TIMEOUT       | TIMEOUT        | TIMEOUT       | TIMEOUT        |
+---------+---------------+---------------+---------------+----------------+---------------+----------------+
| Lo-θ    | 5h 19m 27.17s | 6h 33m 6.55s  | TIMEOUT       | TIMEOUT        | TIMEOUT       | TIMEOUT        |
+---------+---------------+---------------+---------------+----------------+---------------+----------------+
| SVM-θ   | 39.40s        | 1m 0.57s      | 19m 24.73s    | 23m 14.31s     | 52m 10.36s    | 5m 57.31s      |
+---------+---------------+---------------+---------------+----------------+---------------+----------------+
| ODD-STh | 4.47s         | 4.85s         | 1m 53.50s     | 4m 48.92s      | 8m 20.66s     | 2h 1m 9.55s    |
+---------+---------------+---------------+---------------+----------------+---------------+----------------+
| PM      | 1m 28.02s     | 2m 13.01s     | 10m 9.24s     | 51m 45.10s     | 3h 50m 38.60s | 36m 26.14s     |
+---------+---------------+---------------+---------------+----------------+---------------+----------------+
| GH      | 2m 11.15s     | 2m 3.71s      | TIMEOUT       | TIMEOUT        | TIMEOUT       | 5h 51m 32.27s  |
+---------+---------------+---------------+---------------+----------------+---------------+----------------+
| SM      | TIMEOUT       | TIMEOUT       | OUT-OF-MEM    | OUT-OF-MEM     | OUT-OF-MEM    | TIMEOUT        |
+---------+---------------+---------------+---------------+----------------+---------------+----------------+
| PK      | 7.41s         | 14.26s        | 1m 23.42s     | 5m 49.01s      | 20m 41.73s    | 4m 34.26s      |
+---------+---------------+---------------+---------------+----------------+---------------+----------------+
| ML      | 1h 22m 6.04s  | 1h 41m 13.74s | 8h 21m 18.76s | 47m 51.91s     | OUT-OF-MEM    | 9h 24m 15.22s  |
+---------+---------------+---------------+---------------+----------------+---------------+----------------+
| CORE-WL | 36.74s        | 1m 1.82s      | 45m 1.09s     | OUT-OF-MEM     | OUT-OF-MEM    | OUT-OF-MEM     |
+---------+---------------+---------------+---------------+----------------+---------------+----------------+
| CORE-SP | 3m 58.29s     | 4m 29.55s     | 10h 37m 3.94s | TIMEOUT        | OUT-OF-MEM    | TIMEOUT        |
+---------+---------------+---------------+---------------+----------------+---------------+----------------+

Similar to the labeled case, VH is again the fastest kernel on all datasets. The running time of PK, ODD-STh, and WL-VH is also low compared to the other kernels on most datasets. The SVM-θ, NH, PM, SP and CORE-WL kernels were also competitive in terms of running time. Besides achieving low accuracy levels, the Lo-$\vartheta$ kernel is also very computationally expensive. The RW, NSPDK, CORE-SP, WL-SP and GH are also very expensive in terms of running time. The above $6$ kernels did not manage to calculate any kernel matrix on REDDIT-MULTI-5K and REDDIT-MULTI-12K within one day. The SM kernel failed to compute the kernel matrix on IMDB-BINARY, IMDB-MULTI and COLLAB within one day, while it exceeded the maximum available memory on the remaining three datasets. The WL-VH and CORE-WL kernels also exceeded the maximum available memory on REDDIT-MULTI-5K and REDDIT-MULTI-12K. It should be mentioned that these two datasets contain several thousands of graphs, while the size of the graphs is also large (\ie several hundreds of vertices on average).

Node-Attributed Graphs
^^^^^^^^^^^^^^^^^^^^^^
The Table below illustrates the average CPU running time for kernel matrix computation on the 5 classification datasets containing node-attributed graphs.

+----+------------+---------------+---------------+------------+-------------+
|    | ENZYMES    | PROTEINS_full | SYNTHETICnew  | Synthie    | BZR         | 
+====+============+===============+===============+============+=============+
| SP | TIMEOUT    | TIMEOUT       | TIMEOUT       | TIMEOUT    | TIMEOUT     | 
+----+------------+---------------+---------------+------------+-------------+
| SM | TIMEOUT    | OUT-OF-MEM    | TIMEOUT       | TIMEOUT    | 8h 2m 3.79s |
+----+------------+---------------+---------------+------------+-------------+
| GH | 16m 36.12s | 5h 16m 46.48s | 13m 54.36s    | 24m 20.00s | 4m 24.79s   |
+----+------------+---------------+---------------+------------+-------------+
| PK | 15.85s     | 1m 43.58s     | 13.44s        | 34.68s     | 10.40s      | 
+----+------------+---------------+---------------+------------+-------------+
| ML | 26.05s     | 4h 29m 35.69s | 2h 54m 31.22s | 15m 11.29s | 49m 33.60s  |
+----+------------+---------------+---------------+------------+-------------+

In terms of running time, PK is the most efficient kernel since it handled all datasets in less than two minutes. GH and ML are much slower than PK on all datasets. For instance, the average computation time of ML and GH was greater than 4 hours and 5 hours on PROTEINS_full, respectively. The SP and SM kernels, as already discussed, are very expensive in terms of running time, and hence, their usefulness in real-world problems is limited.