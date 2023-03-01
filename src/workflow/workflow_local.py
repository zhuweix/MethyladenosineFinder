import argparse
import os
from .generate_sbatch_local import generate_sbatch_local


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
    parser.add_argument('-s', '--shdir', default='./sh/', help='Folder for sh scripts. Default: ./sh/')
    parser.add_argument('-c', '--coverage', default=6, help='Minimal Subread Coverage. Default: 6')
    parser.add_argument('-@', '--threads', default=12, help='Threads for samtools. Default: 12')
    parser.add_argument('--scorefn',default=os.path.abspath(os.path.join(script_dir, '../asset/default_cov_score.csv')))
    parser.add_argument('-f', '--m6Aonly', default=1, help='1=Only Include m6A sites, 0=All modified As. Default=1')
    parser.add_argument('--timeout', default=600, type=int,
                        help='Maximal Time (s) for single ipdSummary job. Default: 600')
    parser.add_argument('--isclean', default=True, help='Whether to remove tmp files. Default=True',)
    args = parser.parse_args()
    if isinstance(args.isclean, str):
        isclean = args.isclean.lower()
        isclean = isclean.capitalize()
        if isclean[0] in ['T', 'F']:
            isclean = eval(isclean)
        else:
            isclean = int(isclean) > 0
    generate_sbatch_local(
        bam=os.path.abspath(args.bamfile),
        out=os.path.abspath(args.tmpdir),
        prefix=args.prefix,
        sbatchdir=os.path.abspath(args.shdir),
        ref=os.path.abspath(args.ref),
        mod=os.path.abspath(args.index),
        cov=args.coverage,
        savedir=os.path.abspath(args.outputdir),
        threads=args.threads,
        is_m6A=args.m6Aonly,
        timeout=args.timeout,
        score_fn=os.path.abspath(args.scorefn),
        is_clean=isclean
        )


if __name__ == "__main__":
    main()
