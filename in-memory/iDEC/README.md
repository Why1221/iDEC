# iDEC Codes

## Some Important Issues

+ KD-tree: distance type should be `float`, otherwise, the query speedup would be much slower
+ build index: no need to build again, it will be built during the construction

## Potential

+ Duplication after projection: deduplication might help (apply to in-memory only)
+ Steal version (apply to both in-memory and external-memory)