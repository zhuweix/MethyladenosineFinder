#!/usr/bin/env python
import argparse
import os

def generate_sbatch_ipd_local(bamdir: str, swarmfile: str, zmwfile: str, outdir: str, score_fn: str,
                   motifmodfile: str, reference: str, coveragecutoff: int, is_strict: int, timeout: int, is_clean: bool):
    """Generate Swarm file for IPDSummary analysis"""
    # load zmw list
    zmw_list = []
    with open(zmwfile) as filep:
        for line in filep:
            zmw_list.append(line.split()[0])
    script_dir = os.path.dirname(os.path.realpath(__file__))

    # generate command
    content = [r'''#!/bin/bash''']
    script_fn = os.path.abspath(os.path.join(script_dir, 'ipd_analysis.py'))
    for zmw in zmw_list:
        content.append('python '
                       '{} '
                       '-b {}/tmp.{}.bam -o {} -m {} -r {} -c {} -f {} -t {} -s {} --is_clean {}'.format(
                        script_fn, bamdir, zmw, outdir, motifmodfile,
                        reference, coveragecutoff, is_strict, timeout, score_fn, is_clean))

    with open(swarmfile, 'w') as filep:
        filep.write('\n'.join(content))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bamdir', required=True)
    parser.add_argument('-c', '--coveragecutoff', default=3, type=int)
    parser.add_argument('-o', '--outdir', required=True)
    parser.add_argument('-s', '--swarmfile', required=True)
    parser.add_argument('-z', '--zmwfile', required=True)
    parser.add_argument('-m', '--motifmodfile', required=True)
    parser.add_argument('-r', '--reference', required=True)
    parser.add_argument('--scorefn', required=True)
    parser.add_argument('-t', '--timeout', default=600, type=int)
    parser.add_argument('-f', '--is_strict_flag', default=1)
    parser.add_argument('--is_clean', default=True)

    args = parser.parse_args()

    generate_sbatch_ipd_local(
        bamdir=args.bamdir,
        swarmfile=args.swarmfile,
        zmwfile=args.zmwfile,
        outdir=args.outdir,
        motifmodfile=args.motifmodfile,
        reference=args.reference,
        is_strict=args.is_strict_flag,
        is_clean=args.is_clean,
        coveragecutoff=args.coveragecutoff,
        timeout=args.timeout,
        score_fn=args.scorefn)