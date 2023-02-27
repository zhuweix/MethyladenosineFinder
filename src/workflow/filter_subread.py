import pysam


def filter_subread(bam: str, zmw: str, output: str):
    zmw_list = []
    with open(zmw) as filep:
        for line in filep:
            zmw_list.append(int(line.strip()))
    print('{} ZMWs in the ZMW list file'.format(len(zmw_list)))
    zmw_list = set(zmw_list)
    read_count = 0
    with pysam.AlignmentFile(bam, 'rb', check_sq=False) as inbam:
        header = inbam.header
        with pysam.AlignmentFile(output, 'wb', header=header) as outbam:
            for read in inbam:
                zmw = int(read.query_name.split('/')[1])
                if zmw in zmw_list:
                    outbam.write(read)
                    read_count += 1
    print('{} Subreads written in {}'.format(read_count, output))
