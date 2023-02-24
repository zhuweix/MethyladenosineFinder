#!/usr/bin/env python
import pysam
import argparse
import os


def filter_split_zmw(bamfile: str, coveragecutoff: int, outdir: str):
    """
    Split SMRT reads into individual bamfiles that one zmw per file.
    :param bamfile: Name of Bam File
    :param coveragecutoff: Min Coverage for subreads per ZMW
    :param outdir: Output Dir
    :return:
    """

    print('Analyzing ZMW from bam files')
    zmw_count = {}
    with pysam.AlignmentFile(bamfile, "rb", check_sq=False) as bam:
        for read in bam:
            # Only use primary alignment
            readflag = read.flag
            if readflag != 0 or readflag != 16:
                continue
            # Get ZMW name: Name of ZMW read: {movieName}/{ZMWNumber}/{qStart}_{qEnd}
            holeNumber = read.query_name.split('/')[1]
            if holeNumber:
                zmw_count.setdefault(holeNumber, 0)
                zmw_count[holeNumber] += 1

    print('Counting ZMW holes with >={} coverage (subreads)'.format(coveragecutoff))

    filter_ZMWholes = []

    # Save zmws
    print('Saving ZMW reads')
    read_depth = {}
    zmw_read = {}
    count = 0
    with pysam.AlignmentFile(bamfile, "rb", check_sq=False) as bam:
        for read in bam:
            # Get ZMW name: Name of ZMW read: {movieName}/{ZMWNumber}/{qStart}_{qEnd}
            holeNumber = read.query_name.split('/')[1]
            if holeNumber:
                if zmw_count[holeNumber] < coveragecutoff:
                    continue
                zmw_read.setdefault(holeNumber, [0, []])
                zmw_read[holeNumber][0] +=1
                zmw_read[holeNumber][1].append(read)
                if zmw_read[holeNumber][0] == zmw_count[holeNumber]:
                    tmp_dict = {}
                    for read2 in zmw_read[holeNumber][1]:
                        chrom = read.reference_name
                        start = read.reference_start
                        end = read.reference_end
                        tmp_dict.setdefault(chrom, [])
                        tmp_dict[chrom].append((start, end, read2))
                    if len(tmp_dict) > 1:
                        ct = -1
                        chrom = ''
                        for c, reads in tmp_dict.items():
                            if len(reads) > ct:
                                chrom = c
                                ct = len(reads)
                    tmp_reads = tmp_dict[chrom]
                    if len(tmp_reads) < coveragecutoff:
                        continue
                    filter_ZMWholes.append(holeNumber)
                    with pysam.AlignmentFile('{}/tmp.{}.bam'.format(outdir, holeNumber), 'wb', header=bam.header) as outbam:
                        for (_, _, read2) in tmp_reads:
                            outbam.write(read2)
                            del read2
                        count += 1
                    del zmw_read[holeNumber]
    print('{} ZMWs counted'.format(len(filter_ZMWholes)))
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    # Save read list
    with open('{}/zmw.cov.{}.list.txt'.format(outdir, coveragecutoff), 'w') as filep:
        filep.write('\n'.join(filter_ZMWholes))
    print('Finished. {} BAM files are generated.'.format(count))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bamfile')
    parser.add_argument('-c', '--coveragecutoff')
    parser.add_argument('-o', '--outdir')
    args = parser.parse_args()
    filter_split_zmw(bamfile=args.bamfile, outdir=args.outdir, coveragecutoff=int(args.coveragecutoff))