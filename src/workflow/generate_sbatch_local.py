import os
import argparse


def generate_sbatch_local(bam: str, out: str, prefix: str, score_fn : str, is_clean: bool,
                    sbatchdir: str, ref: str, mod: str, cov: str, savedir: str,
                    threads: int, is_m6A=1, timeout=600):
    if not os.path.isdir(sbatchdir):
        os.makedirs(sbatchdir)
    if not os.path.isdir(savedir):
        os.makedirs(savedir)
    sfn = '{}.split.sh'.format(prefix)
    sfn = os.path.join(sbatchdir, sfn)
    fullprefix = os.path.join(out, prefix)
    pfn = '{}.pipeline'.format(prefix)
    pfn = os.path.join('./', pfn)
    script_dir = os.path.dirname(os.path.realpath(__file__))
    filter_split_py = os.path.join(script_dir, 'filter_split_zmw.py')
    with open(sfn, 'w') as filep:
        splitbam = ['''#!/bin/bash''',
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

        gen_ipd_py = os.path.join(script_dir, 'generate_sbatch_ipd_local.py')
        ipdscript = [('python {gen_py} \\\n'
                      '\t-b {zmw} \\\n\t-c {cov} \\\n'
                      '\t-o {ipd} \\\n\t-s {sh} \\\n'
                      '\t-m {mod} \\\n\t-z {zlist} \\\n'
                      '\t-t {timeout} \\\n\t--scorefn {scorefn} \\\n'
                      '\t--is_clean {isclean} \\\n'
                      '\t-r {ref} \\\n\t-f {strict_flag}\n\n').format(
            gen_py=gen_ipd_py,  scorefn=score_fn,
            zmw=fullprefix + '_zmw', ipd=fullprefix + '_ipd',
            mod=mod, ref=ref, sh=os.path.join(sbatchdir, prefix + '.ipd_analysis.sh'),
            zlist=os.path.join(fullprefix + '_zmw', 'zmw.cov.{}.list.txt'.format(cov)),
            cov=cov, strict_flag=is_m6A, timeout=timeout, isclean=is_clean)]
        filep.write('\n'.join(ipdscript))

    p1 = ('# Split Subreads by ZMW: Large Number of Files will be generated!\n'
          'bash\t{sh1} \n\n').format(sh1=sfn)
    p2 = ('# Predict m6A sites\n'
          '# {sh2} is generated after {sh1}\n'
          'bash\t{sh2} \n\n').format(sh1=sfn, sh2=os.path.join(sbatchdir, prefix + '.ipd_analysis.sh'))
    mfn = '{}.merge.sh'.format(prefix)
    mfn = os.path.join(sbatchdir, mfn)
    with open(mfn, 'w') as filep:
        mergebam = [r'''#!/bin/bash''',
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
          'bash\t{sh}\n\n').format(
        sh=mfn)

    with open(pfn, 'w') as filep:
        filep.write(''.join([p1, p2, p3]))



