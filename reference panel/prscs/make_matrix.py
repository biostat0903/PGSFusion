import sys
import os
from multiprocessing import Pool

OUT = "/gpfs/chencao/zhenghuili/data/geno/pgsfusion/ldblk_1kg/ldblk_1kg_eur"
chr = sys.argv[1]
regions = int(sys.argv[2])

def task(n):
    print(n)
    id = f'{chr}_{n}'
    if os.path.exists(f'{OUT}/{id}.matrix'):
        return
    if not os.path.exists(f'{OUT}/{id}.ld'):
        return
    with open(f'{OUT}/{id}.snplist', 'r') as f:
        snp_set = set([x.strip() for x in f.readlines()])

    pre = None

    with open(f'{OUT}/{id}.ld', 'r') as f, open(f'{OUT}/{id}.matrix', 'w') as f2:
        f.readline()
        for line in f:
            arr = line.strip().split()
            if arr[2] == pre:
                if arr[5] in snp_set:
                    f2.write(' ')
                    f2.write(arr[6])
            else:
                if arr[5] in snp_set:
                    if pre != None:
                        f2.write('\n')
                    f2.write(arr[6])
                    pre = arr[2]
        f2.write('\n')

if __name__ == '__main__':
    Pool(4).map(task, range(1, regions + 1))
    # for i in range(1, regions + 1):
    #     task(i)