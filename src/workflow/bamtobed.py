import argparse
import numpy
import pysam


def m6a_bam_to_bed(fn: str, bed: str):
    with pysam.AlignmentFile(fn, 'rb') as bam:
        refs = list(bam.header.references)
        refs.sort()
        feat_chrom = {r: [] for r in refs}

        for read in bam:
            chrom = read.reference_name
            cig = read.cigartuples
            name = read.query_name
            start = read.reference_start
            tmp_feat = []
            p = start
            for tp, num in cig:
                if tp != 1:
                    p += num  # 1-based right end
                else:
                    tmp_feat.append(p)  # 1-based A location

            feat_chrom[chrom].append((name, start + 1, p, tmp_feat, 0))  # 1-based close left, right
    content = ['track name=m6A description="m6A sites in Reads" useScore=0']
    for chrom, feats in feat_chrom.items():
        for name, start, end, tmp_feat, sc in feats:
            bstart = []
            bsize = [1] * len(tmp_feat)
            if len(tmp_feat) == 0:
                continue
            for p in tmp_feat:
                bstart.append(p - start)
            content.append(
                '{chrom}\t{start}\t{end}\t{name}\t{score}\t.\t{start}\t{end}\t0\t{cblock}\t{bsizes}\t{bstarts}'.format(
                    chrom=chrom, name=name, start=start - 1, end=end, score=0, cblock=len(bsize), # 0-bases close left open right
                    bsizes=','.join(map(str, bsize)), bstarts=','.join(map(str, bstart))
                ))
    with open(bed, 'w') as filep:
        filep.write('\n'.join(content))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam', help='BAM File to be converted',required=True, type=str)
    parser.add_argument('--bed', help='Output BED file', required=True, type=str)
    m6a_bam_to_bed(
        bam=bam,
        bed=bed
    )


if __name__ == '__main__':
    main()