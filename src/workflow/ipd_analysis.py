import pickle
import subprocess
import argparse
import os
import numpy as np
import pysam


def load_score(score_fn: str):
    scorecutoffdict = {}
    with open(score_fn) as filep:
        next(filep)
        for line in filep:
            ent = line.split(',')
            scorecutoffdict[int(ent[0])] = float(ent[1])
    return scorecutoffdict

def ipd_motif_finder(gff: str, motifpositionlist: np.ndarray, avecov: int, mod: str,
                     scorecutoffdict: dict,
                     start: int, end: int, is_strict=False):
    """
    Filter the modified Motif Location from PacBio IPDSummary
    Args:
        gff (str): name of the gff file from PacBio IPDSummary
        motifpositionlist (np.ndarray): numpy array for A locations 1=A/T
        avecov (int): average coverage for Subreads per ZMW
        mod (str): Type of the modification predicted by PacBio IPDSummary
        scorecutoffdict (dict): Dictionary For {Coverage: p-value threshold}
        start (int): Start location
        end (int): End location
        is_strict (bool): Whether only include m6A in ipdSummary

    Returns:
        list: sorted (1-based) positions of the modified nucleotides
    """


    # fill score dict
    if avecov > max(scorecutoffdict.keys()):
        pvalue = max(scorecutoffdict.values())
    elif avecov < min(scorecutoffdict.keys()):
        pvalue = min(scorecutoffdict.values())
    else:
        pvalue = scorecutoffdict[avecov]

    modpositions = {}
    with open(gff) as filep:
        for line in filep:
            if '##' not in line:
                line = line.strip().split('\t')
                if is_strict:
                    if line[2] != mod:  # Include modification as well
                        continue
                pos = int(line[3])
                if start <= pos <= end and motifpositionlist[pos -1] and int(line[5]) >= pvalue:
                    modpositions[int(line[3])] = True
    modpositions = sorted(list(modpositions.keys()))
    return modpositions


def cigar_writer(modpositions: list, start: int, stop: int):
    """
    Generate the Cigar string for the pseudo read of each ZMW.

    Args:
        modpositions (list): sorted list of identified motif locations (1-based)
        start (int): start position for read alignment
        stop (int): end position for read alignment

    Returns:
        str: Cigar string I for modification M for match (no modification)
    """

    cigartuples = []
    lastposition = start - 1
    for position in modpositions:
        cigartuples.append((0, position - lastposition))
        lastposition = position
        cigartuples.append((1, 1))
    cigartuples.append((0, stop - lastposition))
    return cigartuples


def analyze_ipd_zmw(bamfile: str, output: str, motifpositionfile: str, scorecutoffdict: dict,
                    reference: str, coveragecutoff=2, is_strict=True, timeout=600):
    """
    Generate pseudo reads for individual ZMW with modification predicted by IPDSummary.

    Args:
        bamfile (str): input aligned Bam file
        output (str): output dir
        motifpositionfile (str): Stored nucleotide locations for modification analysis.
        reference (str): name of the reference file (should be indexed for IPDSummary)
        coveragecutoff (int): minimal coverage of ZMW
        scorecutoffdict (dict): Dictionary For {Coverage: p-value threshold}
        is_strict (bool): Whether only include m6A in ipdSummary
        timeout (int): Maximal Process time for one ZMW

    Returns:
        None. Stored in {output}/ipd.tmp.{zmw_name}.gff and {output}/ipd.tmp.{zmw_name}.bam

    """
    zmw_dict = {}
    final_bam = []

    with open(motifpositionfile, 'rb') as filep:
        motifpositiondict = pickle.load(filep)

    if not os.path.isdir(output):
        os.mkdir(output)

    subprocess.run(['samtools', 'index', bamfile],
                   check=True)  # Check samtools output, wait for index to fininsh before proceed

    with pysam.AlignmentFile(bamfile, 'rb', check_sq=False) as inbam:
        header = dict(inbam.header)
        # Edit header for new Bam file
        # label unsorted
        header['HD']['SO'] = 'unsorted'
        count_depth = {}
        for read in inbam:
            chrom = read.reference_name
            start = read.reference_start
            end = read.reference_end
            count_depth.setdefault(chrom, [{}, 0])
            count_depth[chrom][1] += 1
            movie, zmw = read.query_name.split('/')[:2]
            for i in range(start, end + 1):
                count_depth[chrom][0].setdefault(i, 0)
                count_depth[chrom][0][i] += 1
        else:
            read_header = read.header

        if not count_depth:
            return

        if len(count_depth) > 1:
            tmp = [(chrom, e[1]) for chrom, e in count_depth.items()]
            chrom = max(tmp, key=lambda x: x[1])[0]

        count_depth = sorted(count_depth[chrom][0].items())
        fcount_depth = [d for d in count_depth if d[1] >= coveragecutoff]
        if not fcount_depth:
            return
        start = fcount_depth[0][0]
        stop = fcount_depth[-1][0]

        cov_depth = np.zeros(stop - start + 1)
        for d in fcount_depth:
            cov_depth[int(d[0]) - start] = int(d[1])
        avecov = int(np.sum(cov_depth) / len(cov_depth))

        num_read = 1
        # split zero-depth
        # not equal length: pick maximum
        # equal length: randomly select one
        pos_locs = []
        if not np.all(cov_depth):
            tmpmulti = np.where(np.diff(cov_depth > 0) != 0)[0]
            tmpmulti = np.append(tmpmulti, len(cov_depth))
            left = 0
            idx = 0
            for idx, right in enumerate(tmpmulti):
                readstart = start + left + 1
                readstop = start + right + 1
                readcov = int(np.mean(cov_depth[left: right + 1]))
                if readcov < coveragecutoff:
                    continue
                tmp_gff = os.path.join(output, 'ipd.tmp.{}_{}.gff'.format(zmw, idx))
                tmp_bam = os.path.join(output, 'ipd.tmp.{}_{}.bam'.format(zmw, idx))

                # calulcate p-value for motif position using ipdSummary from PacBio
                cmd = ['ipdSummary', bamfile, '--reference', reference,
                       '--gff', tmp_gff, '--identify', 'm6A', '-w', '{}:{}-{}'.format(chrom, readstart, readstop),
                       '--quiet', '-j', '1', '--identifyMinCov', '3', '--pvalue', '0.9']  # --identify 6mA

                proc = subprocess.run(cmd, check=True, timeout=timeout)
                # filter for high-quality motif position for each zmw

                motifpositionlist = motifpositiondict[chrom]

                modpositions = ipd_motif_finder(gff=tmp_gff,
                                                motifpositionlist=motifpositionlist,
                                                avecov=readcov,
                                                mod='m6A',
                                                start=readstart,
                                                end=readstop,
                                                scorecutoffdict=scorecutoffdict,
                                                is_strict=is_strict)
                # Remove tmp file
                cmd = ['rm', tmp_gff]
                proc = subprocess.run(cmd, check=True)

                # generate cigar
                cigar = cigar_writer(modpositions, readstart, readstop)
                # generate synthesized read for mod analysis
                new_read = pysam.AlignedSegment(header=read_header)
                new_read.query_name = '{}/{}/{}'.format(movie, zmw, idx)
                new_read.flag = 0
                new_read.cigartuples = cigar
                new_read.mapping_quality = 255
                new_read.reference_name = chrom
                new_read.reference_start = readstart - 1
                new_read.set_tags([('cv', readcov, 'i'),
                                   ('ln', right - left + 1, 'i')])
                new_read.set_tag('dp', cov_depth[left: right + 1].tolist())
                final_bam.append(new_read)
                with pysam.AlignmentFile(tmp_bam, 'wb', header=header) as outbam:
                    outbam.write(new_read)
                left = right + 1
        else:
            if avecov < coveragecutoff:
                return
            start += 1  # 0 based to 1 based
            stop += 1

            tmp_gff = os.path.join(output, 'ipd.tmp.{}.gff'.format(zmw))
            tmp_bam = os.path.join(output, 'ipd.tmp.{}.bam'.format(zmw))

            # calulcate p-value for motif position using ipdSummary from PacBio
            cmd = ['ipdSummary', bamfile, '--reference', reference,
                   '--gff', tmp_gff, '--identify', 'm6A', '-w', '{}:{}-{}'.format(chrom, start, stop),
                   '--quiet', '-j', '1', '--identifyMinCov', '3', '--pvalue', '0.9']  # --identify 6mA

            proc = subprocess.run(cmd, check=True, timeout=timeout)
            # filter for high-quality motif position for each zmw

            motifpositionlist = motifpositiondict[chrom]

            modpositions = ipd_motif_finder(gff=tmp_gff,
                                            motifpositionlist=motifpositionlist,
                                            avecov=avecov,
                                            mod='m6A',
                                            start=start,
                                            end=stop,
                                            scorecutoffdict=scorecutoffdict,
                                            is_strict=is_strict)
            # Remove tmp file
            # cmd = ['rm', tmp_gff]
            # proc = subprocess.run(cmd, check=True)

            # generate cigar
            cigar = cigar_writer(modpositions, start, stop)
            # generate synthesized read for mod analysis
            new_read = pysam.AlignedSegment(header=read_header)
            new_read.query_name = '{}/{}'.format(movie, zmw)
            new_read.flag = 0
            new_read.cigartuples = cigar
            new_read.mapping_quality = 255
            new_read.reference_name = chrom
            new_read.reference_start = start - 1
            new_read.set_tags([('cv', avecov, 'i'),
                               ('ln', stop - start + 1, 'i')])
            new_read.set_tag('dp', cov_depth.tolist())
            final_bam.append(new_read)
            with pysam.AlignmentFile(tmp_bam, 'wb', header=header) as outbam:
                outbam.write(new_read)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bamfile', required=True)
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('-r', '--reference', required=True)
    parser.add_argument('-m', '--motifpositionfile')
    parser.add_argument('-f', '--is_m6a', default=1)
    parser.add_argument('-c', '--coveragecutoff', default=3)
    parser.add_argument('-s', '--scorefn', required=True)
    parser.add_argument('-t', '--timeout', default=600)

    args = parser.parse_args()
    score_dict = load_score(score_fn=args.scorefn)
    analyze_ipd_zmw(
        bamfile=args.bamfile,
        output=args.output,
        motifpositionfile=args.motifpositionfile,
        reference=args.reference,
        coveragecutoff=int(args.coveragecutoff),
        scorecutoffdict=score_dict,
        is_strict=int(args.is_m6a) > 0,
        timeout=int(args.timeout))