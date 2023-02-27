import argparse
from .indexreference import index_reference
from .filter_ccs import filter_ccs
from .filter_subread import filter_subread

def main():
    parser = argparse.ArgumentParser(
        prog='Preprocess',
        description="Preprocess Data for m6A identification")
    subparser = parser.add_subparsers(dest='command')
    indexref = subparser.add_parser('index', help='Prepare Reference index for A locations')
    zmwfromccs = subparser.add_parser('zmwfromccs', help='Get ZMWs with AveBaseQuality Filter from CCS Reads')
    filtersubread = subparser.add_parser('filtersubread', help='Filter SubReads by zmw list')

    indexref.add_argument('-r', '--ref', type=str, help='Reference File in FASTA Format', required=True)
    indexref.add_argument('-o', '--output', type=str, help='Name of the index file', required=True)

    zmwfromccs.add_argument('-b', '--bam', type=str, help='CCS Read BAM file', required=True)
    zmwfromccs.add_argument('-q', '--qual', default=90, type=int,
                            help='Minimal Average Base Quality to Filter ZMWs. Default=90')
    zmwfromccs.add_argument('-c', '--cov', default=3, type=int,
                            help='Minimal Subread passes in CCS reads. Default=3')
    zmwfromccs.add_argument('-o', '--output', type=str, default='filter_zmw.txt',
                            help='Output ZMW list. Default=filter_zmw.txt',)

    filtersubread.add_argument('-b', '--bam', type=str, required=True,
                               help='Name of the PacBio Subread')
    filtersubread.add_argument('-z', '--zmw', type=str, default='filter_zmw.txt',
                               help='Name of the ZMW list. Default=filter_zmw.txt')
    filtersubread.add_argument('-o', '--output', type=str, required=True,
                               help='Name of the output Subread file')

    args = parser.parse_args()
    if args.command == 'index':
        ref = args.r
        output = args.o
        print('Making Reference Dictionary for {}'.format(ref))
        index_reference(ref=ref, output=output)
        print('Finished')
    elif args.command == 'zmwfromccs':
        bam = args.bam
        qual = args.qual
        cov = args.cov
        output = args.output
        print('Start making ZMW list with Qual >= {}, Subread Pass >= {}'.
              format(qual, cov))
        filter_ccs(
            bam=bam,
            min_qual=qual,
            min_cov=cov,
            output=output
        )
    elif args.command == 'filtersubread':
        bam = args.bam
        zmw = args.zmw
        output = args.output
        filter_subread(
            bam=bam,
            zmw=zmw,
            output=output
        )
