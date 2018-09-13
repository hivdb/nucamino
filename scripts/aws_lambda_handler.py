import os
import json

from subprocess import Popen, PIPE, TimeoutExpired

LAMBDA_TASK_ROOT = os.environ['LAMBDA_TASK_ROOT']
NUCAMINO_BIN = os.path.join(LAMBDA_TASK_ROOT, 'nucamino')


def align(event, context):
    profile = event['profile']
    genes = event['genes']
    fasta = event['fasta']
    timeout = event.get('timeout', 60)
    if isinstance(genes, (tuple, set, list)):
        genes = ','.join(genes)
    proc = Popen([NUCAMINO_BIN, 'align', profile, genes, '-f', 'json'],
                 stdin=PIPE, stdout=PIPE, stderr=PIPE, encoding='UTF-8')
    try:
        out, _ = proc.communicate(input=fasta, timeout=timeout)
    except TimeoutExpired:
        proc.kill()
    return json.loads(out)
