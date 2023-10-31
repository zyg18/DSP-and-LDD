# Public library:
ortools,eigen-3.3.8
# An example for exhibiting how to run it:
chmod +x demo_pruning.sh  

./demo_pruning.sh

# Other algorithms can be run by the same way as demo_pruning.sh, like:
g++  pruning+frank-wolfe.cpp -o pruning+fw.out

./pruning_fw.out comm-EmailEnron_SNAP.edgelist 0 0
