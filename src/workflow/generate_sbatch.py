import os
import argparse


def generate_sbatch(bam: str, out: str, prefix: str, jobname: str, gres: str, score_fn : str, is_clean: bool,
                    sbatchdir: str, ref: str, mod: str, cov: str, savedir: str, split_time: int, ipd_time: int,
                    merge_time: int, mem: str, logdir: str, threads: int, batchsize: int, is_m6A=0, timeout=600):
    if not os.path.isdir(sbatchdir):
        os.makedirs(sbatchdir)
    if not os.path.isdir(logdir):
        os.makedirs(logdir)
    if not os.path.isdir(savedir):
        os.makedirs(savedir)
    sfn = '{}.split.sh'.format(prefix)
    sfn = os.path.join(sbatchdir, sfn)
    fullprefix = os.path.join(out, prefix)
    pfn = '{}.pipeline'.format(prefix)
    pfn = os.path.join('./', pfn)
    script_dir = os.path.dirname(os.path.realpath(__file__))
    filter_split_py = os.path.join(script_dir, 'filter_split_zmw.py')
    log_out = os.path.join(logdir, '{}.split.out'.format(jobname))
    log_err = os.path.join(logdir, '{}.split.err'.format(jobname))
    with open(sfn, 'w') as filep:
        splitbam = [('''#!/bin/bash
#SBATCH -o {}
#SBATCH -e {}''').format(log_out, log_err),
                    ('python {split_py} \\\n'
                    '\t-b {sort} \\\n\t-c {cov} \\\n\t-o {zmw} ; \\\n').format(
            split_py=filter_split_py, sort=bam, zmw=fullprefix + '_zmw', cov=cov)]
        filep.write('\n'.join(splitbam))
        if not os.path.isdir(fullprefix + '_zmw'):
            os.makedirs(fullprefix + '_zmw')
        if not os.path.isdir(fullprefix + '_ipd'):
            os.makedirs(fullprefix + '_ipd')
        if not os.path.isdir(savedir):
            os.makedirs(savedir)
        log_out = os.path.join(logdir, '{}.gen_ipd.out'.format(jobname))
        log_err = os.path.join(logdir, '{}.gen_ipd.err'.format(jobname))
        gen_ipd_py = os.path.join(script_dir, 'generate_sbatch_ipd.py')
        ipdscript = [('python {gen_py} \\\n'
                      '\t-b {zmw} \\\n\t-c {cov} \\\n'
                      '\t-o {ipd} \\\n\t-s {sh} \\\n'
                      '\t-m {mod} \\\n\t-z {zlist} \\\n'
                      '\t-t {timeout} \\\n\t--scorefn {scorefn} \\\n'
                      '\t--log {log} \\\n\t--is_clean {isclean} \\\n'
                      '\t--job {job} \\\n\t--batch {batch}\\\n'
                      '\t-r {ref} \\\n\t-f {strict_flag}\n\n').format(
            gen_py=gen_ipd_py, log=logdir, job=jobname, scorefn=score_fn,
            zmw=fullprefix + '_zmw', ipd=fullprefix + '_ipd', batch=batchsize,
            mod=mod, ref=ref, sh=os.path.join(sbatchdir, prefix + '.ipd_analysis.sh'),
            zlist=os.path.join(fullprefix + '_zmw', 'zmw.cov.{}.list.txt'.format(cov)),
            cov=cov, strict_flag=is_m6A, timeout=timeout, isclean=is_clean)]
        filep.write('\n'.join(ipdscript))

    p1 = ('# Split Subreads by ZMW: Large Number of Files will be generated!\n'
          'sbatch \\\n'
          '\t--mem {mem} \\\n\t--gres={scratch}:40 \\\n'
          '\t--job-name {job}split \\\n\t--time {stime}  \\\n'
          '\t{sh1} \n\n').format(sh1=sfn, scratch=gres, mem=mem, stime=split_time, job=jobname)
    p2 = ('# Predict m6A sites\n'
          '# {sh2} is generated after {sh1}\n'
          'sbatch \\\n'
          '\t--gres={scratch}:10 \\\n\t--time {itime} \\\n'
          '\t--job-name {job}ipd \\\n'
          '\t{sh2} \n\n').format(
        scratch=gres, sh1=sfn, itime=ipd_time, job=jobname,
        sh2=os.path.join(sbatchdir, prefix + '.ipd_analysis.sh'))
    log_out = os.path.join(logdir, '{}.merge.out'.format(jobname))
    log_err = os.path.join(logdir, '{}.merge.err'.format(jobname))
    mfn = '{}.merge.sh'.format(prefix)
    mfn = os.path.join(sbatchdir, mfn)
    with open(mfn, 'w') as filep:
        mergebam = [(
            r'''#!/bin/bash
#SBATCH -o {}
#SBATCH -e {}
module load samtools''').format(log_out, log_err),
            ("find  {ipd} -name '*.bam' > {bamlist} ; \\\n"
             '\tsamtools cat --threads {thread} --no-PG -o {merge} \\\n'
             '\t-b {bamlist} ; \\\n'
             'samtools sort -o {sort} \\\n\t {merge} ; \\\n'
             'samtools index {sort}\n\n').format(
            ipd=fullprefix + '_ipd',
            bamlist=os.path.join(fullprefix + '_ipd', 'bamfile.list'),
            merge=os.path.join(savedir, '{}.merge.bam'.format(prefix)),
            sort=os.path.join(savedir, '{}.sort.bam'.format(prefix)),
            thread=threads)]
        if is_clean:
            mergebam.append(r"xargs -a {} -d'\n' rm -f".format(os.path.join(fullprefix + '_ipd', 'bamfile.list')))
            mergebam.append(r"rm -f {}".format(os.path.join(savedir, '{}.merge.bam'.format(prefix))))
        filep.write('\n'.join(mergebam))
        filep.write('\n')
    p3 = ('# Merge Reads\n'
          'sbatch \\\n'
          '\t--mem {mem} \\\n\t--gres={scratch}:10 --time {mtime} \\\n'
          '\t--cpus-per-task={thread} \\\n\t--time 400 \\\n\t--job-name {job}mer \\\n \t{sh}\n\n').format(
        sh=mfn, scratch=gres, mem=mem, thread=threads, mtime=merge_time, job=jobname)

    with open(pfn, 'w') as filep:
        filep.write(''.join([p1, p2, p3]))



