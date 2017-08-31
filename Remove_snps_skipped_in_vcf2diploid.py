from sys import argv

skipped_snp_fp = argv[1]
snp_file_input = argv[2]
snp_file_toKeep = argv[3]
snp_file_skipped = argv[4]


skipped_set=set()
with open(skipped_snp_fp) as f:
    for l in f:
        if l.startswith('Variant overlap'):
            chrom, snppos = l.strip().split(',')[0].split(' ')[-1].split(':')
            skipped_set.add((chrom, snppos))

with open(snp_file_input) as f:
    with open(snp_file_toKeep, 'w') as keep:
        with open(snp_file_skipped, 'w') as skip:
            for l in f:
                chrom, snppos = l.strip().split('\t')[0:2]
                if (chrom, snppos) in skipped_set:
                    skip.write(l)
                else:
                    keep.write(l)



