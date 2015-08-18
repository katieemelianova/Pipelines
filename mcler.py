# script to call all mcl commands
# make sure inout file ends with .abc

import subprocess


path = input('give me the formatted file path')
file = (path)

print(file)

# neg-log10 log10 transforms the blast scores
# seq.mci file wil contain a list of all edge weights, whch can be viewed on a histogram. If too skewed to very high values (all too close), can try to transform them

subprocess.call(["mcxload -abc " + file + " --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o " + file.rstrip('.abc') + ".mci -write-tab " + file.rstrip('.abc') + ".tab"], shell=True)


# -I sets the inflation parameter. 3 was used after testing using anthocyanin biosynthetic pathway genes to backcheck results. -I dictates the 'granularity' of clusters.

subprocess.call(["mcl " + file.rstrip('.abc') + ".mci -I 0.5"], shell=True)

subprocess.call(["mcxdump -icl out." + file.rstrip('.abc') + ".mci.I05 -tabr "+ file.rstrip('.abc') + ".tab -o dump." + file.rstrip('.abc') + ".mci.I05"], shell=True)

