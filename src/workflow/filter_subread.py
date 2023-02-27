import pysam


def filter_subread(bam: str, zmw: str, output: str):
    zmw_list = []
    with open(zmw) as filep:
        for line in zmw:
            zmw_list.append(line.strip())
    zmw_list = set(zmw_list)

    with pysam.AlignmentFile(bam, 'rb', check_sq=False) as inbam:
        header = inbam.header
        with pysam.AlignmentFile(output, 'wb', header=header) as outbam:
            for read in inbam:
                zmw = read.query_name.split('/')[1]
                if zmw in zmw_list:
                    outbam.write(read)
