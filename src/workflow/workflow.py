import argparse
from .generate_sbatch import generate_sbatch


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bamfile', required=True,
                        help='Aligned BAM file, pbmm2/pbalign is recommended.')
    parser.add_argument('-o', '--outputdir', default='./output/',
                        help='Output Folder. Default: ./output/',)
    parser.add_argument('-t', '--tmpoutdir', default='./tmp/',
                        help='Folder for Tmp files. Default: ./tmp '
                             '(Caution: Large Number of Files will be generated.)')
    parser.add_argument('-p', '--prefix', default='output', help='Prefix of Output. Default: output')
    parser.add_argument('-s', '--sbatchdir', default='./sh/', help='Folder for sbatch scripts. Default: ./')
    parser.add_argument('-l', '--logdir', default='./log/', help='Log Folder. Default: ./')
    parser.add_argument('--jobname', default='maw', help='Job prefix. Default: maw')
    parser.add_argument('-r', '--reference', help='Reference File in FASTA format', required=True)
    parser.add_argument('-m', '--Indexed Reference for Adenine locations', required=True)
    parser.add_argument('-c', '--coverage', default=6, help='Minimal Subread Coverage. Default: 6')
    parser.add_argument('-@', '--threads', default=1)
    parser.add_argument('--ipd_batch', default=400, help='Number of subjobs per batch for ipdsummary'
                                                         'Default: 400')
    parser.add_argument('--mem', default='20g', help='Size of memory per cpu. Default 20g')
    parser.add_argument('--gres', default='lscratch', help='Partition for local disk. Default: lscratch')
    parser.add_argument('--score_fn',default='../asset/default_cov_score.csv')
    parser.add_argument('-f', '--m6Aonly', default=0)
    parser.add_argument('--timeout', default=600, type=int,
                        help='Maximal Time (s) for single ipdSummary job. Default: 600')
    args = parser.parse_args()

    generate_sbatch(
        bam=args.bamfile,
        out=args.tmpoutdir,
        prefix=args.prefix,
        sbatchdir=args.sbatchdir,
        logdir=args.logdir,
        ref=args.reference,
        mem=args.mem,
        jobname=args.jobname,
        mod=args.modposfile,
        cov=args.coverage,
        gres=args.gres,
        savedir=args.outputdir,
        threads=args.threads,
        is_m6A=args.m6Aonly,
        timeout=int(args.timeout))


def workflow_csv(parameter_fn: str):
    params = {}
    with open(parameter_fn) as filep:
        for line in filep:
            ent = line.split(',')
            if ent[0].lower() == 'bam':
                params['bam'] = ent[1].strip()


if __name__ == "__main__":
    main()