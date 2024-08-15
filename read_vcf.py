import gzip

fpath = "/home/jose/analyses/g2psol/variants/as.haplotypes.vcf.gz.all"
fhand = gzip.open(fpath, 'rt')

def parse_gt(gt):
    return list(map(int, gt.split("/")))

snps_read = 0
for line in fhand:
    if line.startswith("#"):
        continue
    gts = list(map(parse_gt, line.strip().split("\t")[9:]))
    snps_read += 1
print(f"SNPs read {snps_read}")