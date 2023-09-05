#!/usr/bin/env python
import argparse
import os

def generate_sbatch_ipd(bamdir: str, swarmfile: str, zmwfile: str, outdir: str, batch: int, score_fn: str, log: str, job: str,
                   motifmodfile: str, reference: str, coveragecutoff: int, is_strict: int, timeout: int, is_clean: bool, ipd_time: str):
    """Generate Swarm file for IPDSummary analysis"""
    # load zmw list
    zmw_list = []
    with open(zmwfile) as filep:
        for line in filep:
            zmw_list.append(line.split()[0])
    script_dir = os.path.dirname(os.path.realpath(__file__))
    num_jobs = len(zmw_list)
    batch = batch * 2
    num_subjobs = num_jobs // batch
    if num_jobs % batch != 0:
        num_subjobs += 1
    # generate command
    env = '''
ipdanalysis="python {}/ipd_analysis.py"
bamdir="{}"
outdir="{}"
motifmodfile="{}"
scorefn="{}"
ref="{}"    
'''.format(script_dir, bamdir, outdir, motifmodfile, score_fn, reference)
    prefix = swarmfile[:-3] + '.tmp'
    top_script = [f'''#!/bin/bash
#SBATCH --job-name={job} 
#SBATCH -o {log}/{job}.top.out
#SBATCH -e {log}/{job}.top.err 
prefix={prefix}
''']
    top_script.append(env)
    top_script.append(f'''
task_ids=$(seq 0 {num_subjobs - 1})
previous_task_id=""
for task_id in $task_ids; do
    job_id=$(sbatch ${{prefix}}.${{task_id}}.sh)
done    
    ''')
    cmds = []
    for zmw in zmw_list:
        if not os.path.isfile('{}/tmp.{}.bam'.format(bamdir, zmw)):
            continue
        cmds.append('pbindex $bamdir/tmp.{}.bam'.format(zmw))
        cmds.append('$ipdanalysis '
                       '-b $bamdir/tmp.{}.bam -o $outdir -m $motifmodfile -r $ref -c {} -f {} -t {} -s $scorefn --is_clean {}'.format(
                        zmw, coveragecutoff, is_strict, timeout, is_clean))

    with open(swarmfile, 'w') as filep:
        filep.write('\n'.join(top_script))
    for i in range(num_subjobs):
        with open(f'{prefix}.{i}.sh', 'w') as filep:
            header = (
f'''#!/bin/bash
#SBATCH --job-name={job}.{i}
#SBATCH -o {log}/tmp.{job}.{i}.out
#SBATCH -e {log}/tmp.{job}.{i}.err 
#SBATCH --cpus-per-task=1
#SBATCH --mem=12g
#SBATCH --time={ipd_time}

module load smrtanalysis
module load samtools
''')            
            filep.write(header)
            filep.write(env)
            filep.write('\n'.join(cmds[i*batch:(i+1)*batch]))
            filep.write('\n')

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
    parser.add_argument('--ipd_time', default="8:00:00", type=str)

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
        batch=args.batch,
        ipd_time=args.ipd_time)