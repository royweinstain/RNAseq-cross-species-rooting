#!/usr/bin/env python3
import csv, sys
f = sys.argv[1]
arab = argan = poplar = both_ap = all3 = 0
with open(f) as fh:
    reader = csv.DictReader(fh, delimiter='\t')
    for row in reader:
        a = int(row['arabidopsis']) > 0
        r = int(row['argan']) > 0
        p = int(row['poplar']) > 0
        arab += a; argan += r; poplar += p
        both_ap += (r and p)
        all3 += (a and r and p)
print(f"OGs with arabidopsis: {arab}")
print(f"OGs with argan: {argan}")
print(f"OGs with poplar: {poplar}")
print(f"OGs with argan+poplar: {both_ap}")
print(f"OGs with all 3: {all3}")
