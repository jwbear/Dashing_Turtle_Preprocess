import sys, re
import pandas as pd

# total arguments
n = len(sys.argv)
print("Total arguments passed:", n)
if n < 2:
    print("Requires full path to signal text file.")
    sys.exit(0)
else:
    signal_file = sys.argv[1]

# Arguments passed
print("\nCounting sequences in signal file:", signal_file)
df = pd.read_csv(signal_file, sep ="\t")
#uncompressed read format:
grp_cols = ['contig', 'read_index']
df = df.groupby(grp_cols, group_keys=True)['contig'].count()
df = df.groupby('contig').count()
print(df)
sys.exit(0)


