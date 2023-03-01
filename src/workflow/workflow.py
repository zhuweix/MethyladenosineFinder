import argparse
import os
from .generate_sbatch import generate_sbatch


def main():
    script_dir = os.path.dirname(os.path.realpath(__file__))
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bamfile', required=True,
                        help='Aligned BAM file, pbmm2/pbalign is recommended.')
    parser.add_argument('-r', '--ref', help='Reference File in FASTA format', required=True)
    parser.add_argument('--index', help='Indexed Reference for Adenine locations', required=True)
    parser.add_argument('-o', '--outputdir', default='./output/',
                        help='Output Folder. Default: ./output/',)
    parser.add_argument('-t', '--tmpdir', default='./tmp/',
                        help='Folder for Tmp files. Default: ./tmp '
                             '(Caution: Large Number of Files will be generated.)')
    parser.add_argument('-p', '--prefix', default='output', help='Prefix of Output. Default: output')
    parser.add_argument('-s', '--sbatchdir', default='./sh/', help='Folder for sbatch scripts. Default: ./sh/')
    parser.add_argument('-l', '--logdir', default='./log/', help='Log Folder. Default: ./log/')
    parser.add_argument('--jobname', default='maw', help='Job prefix. Default: maw')
    parser.add_argument('-c', '--coverage', default=6, help='Minimal Subread Coverage. Default: 6')
    parser.add_argument('-@', '--threads', default=12, help='Threads for samtools. Default: 12')
    parser.add_argument('--ipdbatch', default=200, help='Number of subjobs per batch for ipdsummary'
                                                         'Default: 200')
    parser.add_argument('--mem', default='20g', help='Size of memory per cpu. Default 20g')
    parser.add_argument('--gres', default='lscratch', help='local disk for SLURM gres option. Default: lscratch')
    parser.add_argument('--scorefn',default=os.path.abspath(os.path.join(script_dir, '../asset/default_cov_score.csv')))
    parser.add_argument('-f', '--m6Aonly', default=1, help='1=Only Include m6A sites, 0=All modified As. Default=1')
    parser.add_argument('--timeout', default=600, type=int,
                        help='Maximal Time (s) for single ipdSummary job. Default: 600')
    parser.add_argument('--splittime', default=600, type=int, help='Time limit to Split Reads (min). Default=600')
    parser.add_argument('--ipdtime', default=600, type=int, help='Time limit to predict m6A sites (min). Default=600')
    parser.add_argument('--mergetime', default=100, type=int, help='Time limit to Merge Reads (min). Default=100')
    parser.add_argument('--isclean', default=True, help='Whether to remove tmp files. Default=True',
                        action='store_true')
    args = parser.parse_args()
    isclean = args.isclean
    if isinstance(isclean, str):
        isclean = isclean.lower()
        isclean = isclean.capitalize()
        if isclean[0] in ['T', 'F']:
            isclean = eval(isclean)
        else:
            isclean = int(isclean) > 0
    generate_sbatch(
        bam=os.path.abspath(args.bamfile),
        out=os.path.abspath(args.tmpdir),
        prefix=args.prefix,
        sbatchdir=os.path.abspath(args.sbatchdir),
        logdir=os.path.abspath(args.logdir),
        ref=os.path.abspath(args.ref),
        mem=args.mem,
        jobname=args.jobname,
        mod=os.path.abspath(args.index),
        cov=args.coverage,
        gres=args.gres,
        savedir=os.path.abspath(args.outputdir),
        threads=args.threads,
        split_time=args.splittime,
        ipd_time=args.ipdtime,
        merge_time=args.mergetime,
        is_m6A=args.m6Aonly,
        timeout=args.timeout,
        score_fn=os.path.abspath(args.scorefn),
        batchsize=args.ipdbatch,
        is_clean=isclean
        )


if __name__ == "__main__":
    main()