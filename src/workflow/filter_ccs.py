
import pysam


def filter_ccs(bam: str, min_qual: int, min_cov: int, output: str):
    read_score = {}
    with pysam.AlignmentFile(bam, 'rb') as inbam:
        for read in inbam:
            zmw = read.query_name[:-4]  # remove /ccs
            quals = read.get_forward_qualities()
            if len(quals) > 0:
                qual = sum(quals) / len(quals)
            else:
                continue
            np = read.get_tag('np')
            if zmw not in read_score:
                read_score[zmw] = (np, qual)
            else:
                read_score[zmw] = max((np, qual), read_score[zmw])
    zmw_list = []
    for zmw, (np, qual) in read_score.items():
        zmw = zmw.split('/')[1]
        if np < min_cov:
            continue
        if qual < min_qual:
            continue
        zmw_list.append(zmw)
    with open(output, 'w') as filep:
        filep.write('\n'.join(output))
