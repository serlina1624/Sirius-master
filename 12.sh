#!/bin/bash
# 12. RAxML
raxmlHPC -T 10 -q clean_partitions.txt -s clean_UCE.phy -m GTRCAT -p 1234 -n trees -# 500
raxmlHPC -T 10 -q nonclean_partitions.txt -s nonclean_UCE.phy -m GTRCAT -p 1234 -n trees -# 500
raxmlHPC -T 10 -q clean_partitions.txt -s clean_UCE.phy -m GTRCAT -p 1234 -b 1234 -n bootstrap -# 500
raxmlHPC -T 10 -q nonclean_partitions.txt -s nonclean_UCE.phy -m GTRCAT -p 1234 -b 1234 -n bootstrap -# 500
raxmlHPC -T 10 -q clean_partitions.txt -s clean_UCE.phy -m GTRCAT -f b -t Clean_bestTree.trees -z Clean_bootstrap.bootstrap -n Clean_tree
raxmlHPC -T 10 -q nonclean_partitions.txt -s nonclean_UCE.phy -m GTRCAT -f b -t Nonclean_bestTree.trees -z Nonclean_bootstrap.bootstrap -n Nonclean_tree
