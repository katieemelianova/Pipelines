# must be run with Python 3 or higher

import subprocess
path = input('give me a fasta file you want to make a bootstrapped tree from')

mySeq = (path)

subprocess.call(["~/Desktop/Software/raxmlHPC-AVX-v8/raxml -T 2 -m GTRGAMMA -p 12345 -# 20 -s " + mySeq + " -n T13 "], shell=True)

subprocess.call(["~/Desktop/Software/raxmlHPC-AVX-v8/raxml -T 2 -m GTRGAMMA -p 12345 -b 12345 -# 1000 -s " + mySeq + " -n T14 "], shell=True)

subprocess.call(["~/Desktop/Software/raxmlHPC-AVX-v8/raxml -T 2 -m GTRCAT -p 12345 -f b -t RAxML_bestTree.T13 -z RAxML_bootstrap.T14 -n T15"], shell=True)



