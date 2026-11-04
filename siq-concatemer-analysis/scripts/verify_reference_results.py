#!/usr/bin/env python3
"""Compare regenerated outputs against compact historical reference results.

Only the standard library is required. The expected/ fixtures were extracted
from the previous analysis deliverables, not from this package's new outputs.
This command is for the supplied study data, not for a replacement dataset.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TOLERANCE = 1e-10


def read_csv(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def indexed(rows, keys):
    out = {}
    for row in rows:
        key = tuple(str(row[k]) for k in keys)
        if key in out:
            raise AssertionError(f'Duplicate key {keys}: {key}')
        out[key] = row
    return out


def value(x):
    if x is None or x == '':
        return None
    if not isinstance(x, str):
        return x
    if x in ('True', 'true'):
        return True
    if x in ('False', 'false'):
        return False
    if x.startswith('[') or x.startswith('{'):
        try:
            return json.loads(x)
        except json.JSONDecodeError:
            # A half-open interval such as [-7,-1) is a label, not JSON.
            pass
    if re.fullmatch(r'-?[0-9]+', x):
        return int(x)
    try:
        return float(x)
    except ValueError:
        return x


class Verification:
    def __init__(self):
        self.comparisons = Counter()

    def equal(self, expected, actual, label, section, exact=False):
        expected, actual = value(expected), value(actual)
        if isinstance(expected, (list, tuple)):
            assert isinstance(actual, (list, tuple)) and len(expected) == len(actual), (label, expected, actual)
            for i, (e, a) in enumerate(zip(expected, actual)):
                self.equal(e, a, f'{label}[{i}]', section, exact=exact)
        elif isinstance(expected, bool) or expected is None or isinstance(expected, str):
            assert expected == actual, (label, expected, actual)
            self.comparisons[section] += 1
        elif isinstance(expected, (float, int)):
            assert isinstance(actual, (float, int)) and not isinstance(actual, bool), (label, expected, actual)
            if exact:
                assert expected == actual, (label, expected, actual)
            else:
                assert math.isclose(expected, actual, rel_tol=TOLERANCE, abs_tol=TOLERANCE), (label, expected, actual)
            self.comparisons[section] += 1
        else:
            assert expected == actual, (label, expected, actual)
            self.comparisons[section] += 1

    def fields(self, expected, actual, fields, label, section):
        for field in fields:
            exact = ('count' in field or 'denominator' in field or field.endswith('_reads')
                     or field in {'rank', 'n', 'length', 'mh_length', 'n_top_patterns'})
            self.equal(expected[field], actual[field], f'{label}/{field}', section, exact=exact)


def verify(results, fixtures):
    v = Verification()
    datasets = [d for d in ('modify', 'inhibitor') if (results/d).is_dir()]
    if not datasets:
        raise ValueError('No modify/ or inhibitor/ output directory found')
    summaries = {d: indexed(read_csv(results/d/'summary/sample_metrics.csv'), ['Alias']) for d in datasets}
    for dataset in datasets:
        old = read_csv(fixtures/f'historical_mh_summary_{dataset}.csv')
        current = summaries[dataset]
        assert len(old) == len(current), (dataset, 'sample count')
        distribution = indexed(read_csv(results/dataset/'summary/mh_distribution_by_sample.csv'), ['Alias', 'junction_MH_bp'])
        statuses = indexed(read_csv(results/dataset/'summary/mh_status_by_sample.csv'), ['Alias', 'analysis_status'])
        for row in old:
            alias = row['Alias']
            now = current[alias,]
            v.fields(row, now, ['total_reported_reads','pure_deletion_reads','evaluable_deletion_reads',
                               'unevaluable_deletion_reads','evaluable_fraction_of_deletions'], alias, 'strict_mh_sample_summary')
            v.equal(row['MHge2_reads'], now['mh_ge2_reads'], alias+'/MH>=2 reads', 'strict_mh_sample_summary', exact=True)
            v.equal(row['MHge2_fraction_of_evaluable'], now['mh_ge2_fraction_of_evaluable'], alias+'/MH>=2 fraction', 'strict_mh_sample_summary')
            v.equal(row['mean_junction_MH_bp'], now['mean_mh_bp'], alias+'/mean MH', 'strict_mh_sample_summary')
            mh = {int(h): int(r['read_count']) for (a,h),r in distribution.items() if a==alias}
            bins = {'MH0_reads':mh.get(0,0), 'MH1_reads':mh.get(1,0),
                    'MH2to3_reads':mh.get(2,0)+mh.get(3,0), 'MHge4_reads':sum(n for h,n in mh.items() if h>=4)}
            v.fields(row,bins,list(bins),alias,'strict_mh_sample_summary')
            for key in row:
                if key.startswith('excluded_') and key.endswith('_reads'):
                    status = key[len('excluded_'):-len('_reads')]
                    n = statuses.get((alias,status),{}).get('read_count',0)
                    v.equal(row[key],n,alias+'/'+key,'strict_mh_exclusions',exact=True)

    for row in read_csv(fixtures/'historical_sample_metrics.csv'):
        dataset = row['dataset']
        if dataset not in datasets:
            continue
        now = summaries[dataset][row['Alias'],]
        for old,new,scale in [('mean_del_bp','mean_deletion_bp_strict_mh',1),('mean_mh_bp','mean_mh_bp',1),
                              ('deletion_pct_all','deletion_fraction_of_reported',100),
                              ('mh1_pct','mh_ge1_fraction_of_evaluable',100),('mh2_pct','mh_ge2_fraction_of_evaluable',100),
                              ('mh3_pct','mh_ge3_fraction_of_evaluable',100)]:
            v.equal(float(row[old])/scale,now[new],row['Alias']+'/'+new,'historical_sample_metrics')

    groups = {d:indexed(read_csv(results/d/'summary/group_metrics.csv'),['condition','time']) for d in datasets}
    for row in read_csv(fixtures/'historical_group_metrics.csv'):
        dataset = row['dataset']
        if dataset not in datasets:
            continue
        now = groups[dataset][row['condition'],row['time']]
        for old,new,scale in [('mean_del_bp','mean_deletion_bp_strict_mh',1),('mean_mh_bp','mean_mh_bp',1),
                              ('deletion_pct_all','deletion_fraction_of_reported',100),
                              ('mh1_pct','mh_ge1_fraction_of_evaluable',100),('mh2_pct','mh_ge2_fraction_of_evaluable',100),
                              ('mh3_pct','mh_ge3_fraction_of_evaluable',100)]:
            for os,ns in [('_mean','_mean'),('_std','_sd')]:
                v.equal(float(row[old+os])/scale,now[new+ns],str((dataset,row['condition'],row['time'],new+ns)),
                        'historical_group_metrics')

    exact = {d:indexed(read_csv(results/d/'summary/exact_length_mh_by_sample.csv'),['Alias','deletion_length_bp']) for d in datasets}
    for row in read_csv(fixtures/'historical_exact_length_mh_samples.csv'):
        dataset = row['dataset']
        if dataset not in datasets:
            continue
        now = exact[dataset][row['Alias'],row['delSize']]
        label = row['Alias']+'/'+row['delSize']
        v.equal(row['countEvents'],now['evaluable_reads_at_length'],label+'/reads','exact_length_mh',exact=True)
        for h in (1,2,3):
            v.equal(row[f'mh{h}_reads'],now[f'mh_ge{h}_reads'],label+f'/MH>={h} reads','exact_length_mh',exact=True)
            v.equal(float(row[f'mh{h}_pct'])/100,now[f'mh_ge{h}_fraction'],label+f'/MH>={h} fraction','exact_length_mh')

    all_top = {d:indexed(read_csv(results/d/'patterns/all_top5.csv'),['condition','time','rank']) for d in datasets}
    for dataset in datasets:
        expected_rows = [r for r in read_csv(fixtures/'historical_all_top5.csv') if r['dataset']==dataset]
        assert len(expected_rows)==len(all_top[dataset]), (dataset,'all Top5 row count')
        for row in expected_rows:
            label = (dataset,row['condition'],str(row['time_h'])+'h',row['rank'])
            now = all_top[dataset][label[1:]]
            v.equal(row['legacy_id'],now['legacy_id'],str(label)+'/legacy ID','all_length_top5')
            v.equal(row['length'],now['pattern_length'],str(label)+'/length','all_length_top5',exact=True)
            for old,new in [('mean_fraction_all','mean_fraction_all_events'),('sd_fraction_all','sd_fraction_all_events'),
                            ('mean_fraction_deletion','mean_fraction'),('sd_fraction_deletion','sd_fraction')]:
                v.equal(row[old],now[new],str(label)+'/'+new,'all_length_top5')
            counts = value(now['replicate_counts'])
            v.equal([int(row[f'count_{r}']) for r in (1,2,3)],counts,str(label)+'/counts','all_length_top5',exact=True)
            v.equal(row['count_sum'],sum(counts),str(label)+'/count_sum','all_length_top5',exact=True)
            v.equal(row['n_detected'],sum(c>0 for c in counts),str(label)+'/n_detected','all_length_top5',exact=True)
            v.equal([float(row[f'del_frac_{r}']) for r in (1,2,3)],value(now['replicate_fractions']),
                    str(label)+'/deletion fractions','all_length_top5')
            all_denoms=value(now['all_event_denominators'])
            v.equal([float(row[f'all_frac_{r}']) for r in (1,2,3)], [n/d for n,d in zip(counts,all_denoms)],
                    str(label)+'/all-event fractions','all_length_top5')

    if 'modify' in datasets:
        path = results/'modify/patterns'
        table_specs = [
            ('historical_length_group_patterns.csv','length_group_patterns.csv',['length','condition','time','pattern_id'],True),
            ('historical_length_sample_totals.csv','length_sample_totals.csv',['length','alias'],True),
            ('historical_length_patterns.csv','pattern_definitions.csv',['pattern_id'],False),
            ('historical_length_coverage.csv','length_coverage.csv',['length','condition','time'],True),
        ]
        for old_name,new_name,keys,exact_rows in table_specs:
            rows = read_csv(fixtures/old_name)
            new = indexed(read_csv(path/new_name),keys)
            if exact_rows:
                assert len(rows)==len(new), (new_name,'row count',len(rows),len(new))
            for row in rows:
                key=tuple(row[k] for k in keys)
                now=new[key]
                fields=set(row)&set(now)
                # Old labels are aliases, not calculations; preserve numeric
                # results and IDs while allowing descriptive label changes.
                fields -= {'condition_label'}
                v.fields(row,now,fields,str(key),'length_top5_'+new_name.removesuffix('.csv'))
        oldtop={(r['length'],r['condition'],r['time'],r['pattern_id']) for r in
                read_csv(fixtures/'historical_length_group_patterns.csv') if value(r['top5'])}
        newtop={(r['length'],r['condition'],r['time'],r['pattern_id']) for r in read_csv(path/'length_top5.csv')}
        assert oldtop==newtop, 'Exact-length Top5 membership differs'
        v.comparisons['length_top5_membership'] = len(oldtop)
        new=indexed(read_csv(results/'modify/summary/length_standardized.csv'),
                    ['time','comparison','min_reads_per_length','mh_threshold'])
        rows=read_csv(fixtures/'historical_standardized.csv')
        for row in rows:
            key=tuple(row[k] for k in ['time','comparison','min_reads_per_length','mh_threshold'])
            now=new[key]
            v.fields(row,now,['n_common_lengths','common_lengths'],str(key),'length_standardization')
            mapping=[('cy5_adjusted_pct','focus_adjusted_fraction_mean'),('cy5_sd','focus_adjusted_fraction_sd'),
                     ('other_adjusted_pct','comparison_adjusted_fraction_mean'),('other_sd','comparison_adjusted_fraction_sd'),
                     ('difference_pp','difference_fraction'),('cy5_coverage_pct','focus_coverage_fraction_mean'),
                     ('other_coverage_pct','comparison_coverage_fraction_mean')]
            for old,newfield in mapping:
                v.equal(float(row[old])/100,now[newfield],str(key)+'/'+newfield,'length_standardization')
    return {'status':'passed','datasets':datasets,'numeric_tolerance':TOLERANCE,
            'integer_counts':'exact','reference':'Historical analysis files, independently checked against source SIQ workbooks',
            'checks_by_section':dict(v.comparisons),'total_checked_values':sum(v.comparisons.values())}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results',type=Path,required=True)
    p.add_argument('--expected',type=Path,default=ROOT/'expected')
    p.add_argument('--report',type=Path,help='Optional JSON output path')
    args=p.parse_args()
    try:
        report=verify(args.results,args.expected)
    except (AssertionError,ValueError,KeyError,OSError) as exc:
        p.exit(1,'Verification FAILED: '+str(exc)+'\n')
    text=json.dumps(report,indent=2)+'\n'
    if args.report:
        args.report.parent.mkdir(parents=True,exist_ok=True)
        args.report.write_text(text,encoding='utf-8')
    print(text,end='')


if __name__=='__main__':
    main()
