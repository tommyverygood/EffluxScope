#!/usr/bin/env python3
"""Recompute Figure 6 intervals from frozen scaffold sufficient statistics.

This evaluates the resampling stage, not external-model prediction generation.
Requires Python >=3.10 and NumPy. No network, model bundle or private input.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import platform
import numpy as np

HERE=Path(__file__).resolve().parent

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def compute(root=HERE):
    meta=json.loads((root/'provenance.json').read_text())
    data_path=root/'figure6_scaffold_sufficient_statistics.csv'
    expected_path=root/'expected_figure6_intervals.json'
    if sha(data_path)!=meta['data_sha256'] or sha(expected_path)!=meta['expected_sha256']:
        raise ValueError('Frozen input checksum mismatch')
    with data_path.open(newline='') as fh:
        rows=list(csv.DictReader(fh))
    n=len(rows)
    assert n==meta['cluster_count']
    assert [int(r['cluster_index']) for r in rows]==list(range(n))
    names=[k for k in rows[0] if k!='cluster_index']
    matrix=np.array([[int(r[k]) for k in names] for r in rows],dtype=np.int64)
    assert matrix.min()>=0
    totals=matrix.sum(axis=0)
    C={k:matrix[:,i] for i,k in enumerate(names)}
    T={k:totals[i] for i,k in enumerate(names)}
    assert T['in_n']+T['out_n']==meta['source_rows']==2624
    assert T['in_n']==meta['inside_count']==1349
    assert T['in_accepted']+T['out_accepted']==meta['accepted_count']==453
    for side in ['in','out']:
        assert np.array_equal(C[f'{side}_pos']+C[f'{side}_neg'],C[f'{side}_n'])
        assert np.all(C[f'{side}_accepted']<=C[f'{side}_n'])
        for model in meta['models']:
            assert np.all(C[f'{side}_{model}_tp']<=C[f'{side}_pos'])
            assert np.all(C[f'{side}_{model}_tn']<=C[f'{side}_neg'])
    rng=np.random.default_rng(meta['seed'])
    draws=np.empty((meta['replicates'],len(names)),dtype=np.float64)
    for i in range(meta['replicates']):
        draws[i]=matrix[rng.integers(0,n,size=n)].sum(axis=0)
    B={k:draws[:,i] for i,k in enumerate(names)}
    result={'b':{},'band':{}}
    for model in meta['models']:
        for side in ['in','out']:
            valid=(B[f'{side}_pos']>0)&(B[f'{side}_neg']>0)
            vals=.5*(B[f'{side}_{model}_tp'][valid]/B[f'{side}_pos'][valid]+B[f'{side}_{model}_tn'][valid]/B[f'{side}_neg'][valid])
            lo,hi=np.quantile(vals,[.025,.975])
            result['b'][f'{model}|{side}']={
                'n':int(T[f'{side}_n']),'n_pos':int(T[f'{side}_pos']),'n_neg':int(T[f'{side}_neg']),
                'tp':int(T[f'{side}_{model}_tp']),'tn':int(T[f'{side}_{model}_tn']),
                'bacc':float(.5*(T[f'{side}_{model}_tp']/T[f'{side}_pos']+T[f'{side}_{model}_tn']/T[f'{side}_neg'])),
                'lo':float(lo),'hi':float(hi),'valid_draws':int(valid.sum())}
    delta=B['in_accepted']/B['in_n']-B['out_accepted']/B['out_n']
    lo,hi=np.quantile(delta,[.025,.975])
    result['band']={'n_in':int(T['in_n']),'n_out':int(T['out_n']),
        'pos_in':float(T['in_pos']/T['in_n']),'pos_out':float(T['out_pos']/T['out_n']),
        'acc_in':float(T['in_accepted']/T['in_n']),'acc_out':float(T['out_accepted']/T['out_n']),
        'acc_diff':float(T['in_accepted']/T['in_n']-T['out_accepted']/T['out_n']),
        'acc_lo':float(lo),'acc_hi':float(hi),'n_clusters':n}
    expected=json.loads(expected_path.read_text());max_error=0.;checked=0
    for model,d in expected['b'].items():
        for k,value in d.items():
            err=abs(result['b'][model][k]-value);max_error=max(max_error,err);checked+=1
            if err>1e-14:raise AssertionError((model,k,result['b'][model][k],value))
    for k,value in expected['band'].items():
        err=abs(result['band'][k]-value);max_error=max(max_error,err);checked+=1
        if err>1e-14:raise AssertionError(('band',k,result['band'][k],value))
    return {'status':'PASS','scope':meta['scope'],'seed':meta['seed'],'replicates':meta['replicates'],
        'clusters':n,'numeric_fields_checked':checked,'max_absolute_error':max_error,
        'python':platform.python_version(),'numpy':np.__version__,
        'data_sha256':sha(data_path),'script_sha256':sha(Path(__file__)),
        'result':result}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path('figure6_interval_verification.json'))
    args=parser.parse_args();report=compute();args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['status','clusters','replicates','numeric_fields_checked','max_absolute_error']},indent=2))

if __name__=='__main__':main()
