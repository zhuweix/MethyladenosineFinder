#!/usr/bin/env python
import argparse
import os

def generate_sbatch_ipd(bamdir: str, swarmfile: str, zmwfile: str, outdir: str, batch: int, score_fn: str, log: str, job: str,
                   motifmodfile: str, reference: str, coveragecutoff: int, is_strict: int, timeout: int, is_clean: bool):
    """Generate Swarm file for IPDSummary analysis"""
    # load zmw list
    zmw_list = []
    with open(zmwfile) as filep:
        for line in filep:
            zmw_list.append(line.split()[0])
    script_dir = os.path.dirname(os.path.realpath(__file__))
    num_jobs = len(zmw_list)
    num_subjobs = num_jobs // batch
    if num_jobs % batch != 0:
        num_subjobs += 1
    # generate command
    content = [r'''#!/bin/bash
#SBATCH -o {}
#SBATCH -e {}
#SBATCH --cpus-per-task=1
#SBATCH --array=1-{}%2
module load samtools
module load smrtanalysis
'''.format(
        os.path.join(log, '{}.ipd.out'.format(job)),
        os.path.join(log, '{}.ipd.err'.format(job)),
        num_subjobs)]
    for zmw in zmw_list:
        content.append('python '
                       '{}/ipd_analysis.py '
                       '-b {}/tmp.{}.bam -o {} -m {} -r {} -c {} -f {} -t {} -s {} --is_clean {}'.format(
                        script_dir, bamdir, zmw, outdir, motifmodfile,
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
    parser.add_argument('--job', default='maw')
    parser.add_argument('--log', required=True)
    parser.add_argument('-t', '--timeout', default=600, type=int)
    parser.add_argument('-f', '--is_strict_flag', default=1)
    parser.add_argument('--batch', type=int, default=400)
    parser.add_argument('--is_clean', default=True)

    args = parser.parse_args()

    generate_sbatch_ipd(
        bamdir=args.bamdir,
        swarmfile=args.swarmfile,
        zmwfile=args.zmwfile,
        outdir=args.outdir,
        log=args.log,
        job=args.job,
        motifmodfile=args.motifmodfile,
        reference=args.reference,
        is_strict=args.is_strict_flag,
        is_clean=args.is_clean,
        coveragecutoff=args.coveragecutoff,
        timeout=args.timeout,
        score_fn=args.scorefn,
        batch=args.batch)