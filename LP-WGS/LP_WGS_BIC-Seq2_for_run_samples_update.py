import sys
import os
import io
import subprocess
import argparse
import logging
import concurrent.futures
import glob
import re
import csv
import time
import datetime
import resource
#import seaborn as sns
import pandas as pd
import numpy as np
import threading
import psutil
import math
#import matplotlib.pyplot as plt

ref_fasta='/data/home/heewon/ref/HG19/hg19.fa'
cyto = pd.read_csv('/data/home/heewon/ref/cytoBand.csv')

#MEMORY Check
def get_memory_usage(message: str = 'debug'):
    usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    usage_mb = usage / 1024 
    mem=psutil.virtual_memory()
    return usage_mb,mem[2]

parser = argparse.ArgumentParser(description='NBS CNV calling In-house Pipeline')
parser.add_argument('sequencing_id', help='sequencing id (example: 220723_M70544_0051_000000000-G9C6V)')
parser.add_argument('sample_id', help='sample id (example: NB2307001)')
parser.add_argument('-t', dest='threads', type=int, default=20, help='threads (default: %(default)s)')
parser.add_argument('-p', dest='fastp', default='quay.io/biocontainers/fastp:0.23.2--h79da9fb_0', help='docker image (default: %(default)s)')
parser.add_argument('-b', dest='bwa', default='quay.io/biocontainers/bwa:0.7.17--hed695b0_7', help='docker image (default: %(default)s)')
parser.add_argument('-s', dest='samtools', default='quay.io/biocontainers/samtools:1.16.1--h6899075_1', help='docker image (default: %(default)s)')
parser.add_argument('-q',dest='seqtk',default='quay.io/biocontainers/seqtk:1.3--h7132678_4', help='docker image (default: %(default)s)')
parser.add_argument('-e',dest='bicseq2',default='mwyczalkowski/bicseq2:latest', help='docker image (default: %(default)s)')
parser.add_argument('-i',dest='picard',default='quay.io/biocontainers/picard:3.0.0--hdfd78af_1', help='docker image (default: %(default)s)')

args = parser.parse_args()

###########dir info##############
raw_dir = f'/data/raw_data/{args.sequencing_id}'

data_dir = f'/data/analysis/project/230412_nbs_heewon/raw_data/{args.sequencing_id}'
os.makedirs(f'{data_dir}', exist_ok=True)

result_dir = f'/data/analysis/project/230412_nbs_heewon/Result/{args.sequencing_id}'
os.makedirs(f'{result_dir}', exist_ok=True)

table_all_dir=f'{result_dir}/all_sample_result_table'
os.makedirs(f'{table_all_dir}', exist_ok=True)
####################################


#############log info###############
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('[%(levelname)5s - %(asctime)s] %(message)s',"%Y-%m-%d %H:%M:%S")

stream_handler = logging.StreamHandler(sys.stdout)
stream_handler.setFormatter(formatter)
stream_handler.setLevel(logging.INFO)
logger.addHandler(stream_handler)

file_handler = logging.FileHandler(f'{result_dir}/CNV_calling_process.log','w')
file_handler.setFormatter(formatter)
file_handler.setLevel(logging.DEBUG)
logger.addHandler(file_handler)
###################################

def init_env_check():
    logging.info(f'{"-"*20}{"NBS CNV calling In-house Pipeline":^50}{"-"*20}')

    for i in [args.fastp, args.bwa, args.samtools, args.seqtk,args.bicseq2,args.picard]:
        a = subprocess.Popen(f'docker images -q {i}', shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()[0].decode('UTF-8')
        if a == '' :
            b = subprocess.Popen(f'docker pull {i}', shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()[0].decode('UTF-8')
            print(b)
            if b.find('Error'):
                logging.info(b)
                logging.info(f'Docker image file Error === ')
                exit()
        logging.info(f'{"docker image":<15}: {i:<60}{"...checked":>13}')

    logging.info(f'{"cmd":<15}: {" ".join(sys.argv)}')
    logging.info(f'{"log file":<15}: {result_dir}/CNV_calling_process.log')

def log_subprocess_output(step, cmd):
    logging.info(f'{"processing":<15}: {step:<30}')
    logging.debug(f'$ {cmd}')
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    while True:
        line = process.stdout.readline()
        if not line :
            break
        logging.debug('> %s', line.decode('utf-8').strip())
    process.wait()

###################
pattern_1 = r'^(.*)_R1_001\.fastq(\.gz)$'
def data_align():
    sample_id_list = []
    for i in os.listdir(f'{raw_dir}'):
        match_check=re.match(pattern_1, i)
        if match_check is not None:
            sample_id_list.append(match_check.group(1))
    sample_id_list_unique = list(set(sample_id_list))
    return sorted(sample_id_list_unique)

def fastp_run(sample_id,sample_per_threads):
    sample_id_origin = re.match(r'^(.*)_S.*', sample_id).group(1)
    fastp_dir = f'{data_dir}/1.FASTP'
    os.makedirs(fastp_dir, exist_ok=True)
    cmd = f'docker run -u {os.getuid()}:{os.getgid()} --rm -v /data/:/data/ {args.fastp} fastp \
        -i {raw_dir}/{sample_id}_R1_001.fastq.gz \
        -I {raw_dir}/{sample_id}_R2_001.fastq.gz \
        -o {fastp_dir}/{sample_id_origin}_1.fastp.fastq \
        -O {fastp_dir}/{sample_id_origin}_2.fastp.fastq \
        -j {fastp_dir}/{sample_id_origin}.fastp.json \
        -h {fastp_dir}/{sample_id_origin}.fastp.html \
        -w {sample_per_threads} \
        -x --detect_adapter_for_pe \
        --trim_poly_g' 
    log_subprocess_output(f'fastp ({sample_id_origin})', cmd)

def run_bwa(sample_id,sample_per_threads): 
    fastp_dir = f'{data_dir}/1.FASTP'
    bwa_dir = f'{data_dir}/2.BWA'
    os.makedirs(bwa_dir, exist_ok=True)
    cmd = f'docker run -u {os.getuid()}:{os.getgid()} --rm -v /data/:/data/ {args.bwa} bwa mem \
        -t {sample_per_threads} \
        -M {ref_fasta} \
        -o {bwa_dir}/{sample_id}.bwa.sam \
        -T 30\
        {fastp_dir}/{sample_id}_1.fastp.fastq \
        {fastp_dir}/{sample_id}_2.fastp.fastq'
    log_subprocess_output('bwa mem '+sample_id, cmd) 

def run_samtools(sample_id):
    bwa_dir = f'{data_dir}/2.BWA'
    samtools_dir = f'{data_dir}/3.SAMTOOLS'
    os.makedirs(samtools_dir, exist_ok=True)
    cmd = f'docker run -u {os.getuid()}:{os.getgid()} --rm -v /data/:/data/ {args.samtools} samtools view \
        -bt {ref_fasta} \
        -o {samtools_dir}/{sample_id}.bwa.bam \
        {bwa_dir}/{sample_id}.bwa.sam'
    view_resource = log_subprocess_output('samtools view '+sample_id, cmd)

    cmd = f'docker run -u {os.getuid()}:{os.getgid()} --rm -v /data/:/data/ {args.samtools} samtools sort \
        -o {samtools_dir}/{sample_id}.bwa.sorted.bam \
        {samtools_dir}/{sample_id}.bwa.bam'
    sort_resource=log_subprocess_output('samtools sort '+sample_id, cmd)

    cmd=f'docker run -u {os.getuid()}:{os.getgid()} --rm -v /data/:/data/ {args.samtools} samtools index \
        {samtools_dir}/{sample_id}.bwa.sorted.bam'
    index_resource = log_subprocess_output('samtool index '+sample_id, cmd)

def mapping_quality(sample_id):
    samtools_dir = f'{data_dir}/3.SAMTOOLS'
    cmd=f'docker run -u {os.getuid()}:{os.getgid()} --rm -v /data/:/data/ {args.samtools} samtools flagstat \
        {samtools_dir}/{sample_id}.bwa.sorted.bam > {samtools_dir}/{sample_id}_mapping_quality.txt'
    flagstat_resource= log_subprocess_output('samtool flagstat '+sample_id, cmd)


def fastq_total_count(sample_id):
    fastp_dir = f'{data_dir}/1.FASTP'
    os.makedirs(fastp_dir, exist_ok=True)
    count1 = int(os.popen(f'cat "{fastp_dir}/{sample_id}_1.fastp.fastq" | wc -l').read().strip()) // 4
    count2 = int(os.popen(f'cat "{fastp_dir}/{sample_id}_2.fastp.fastq" | wc -l').read().strip()) // 4
    # Write total_count to a text file
    total_count = count1 + count2
    with open(f'{fastp_dir}/{sample_id}_total_count.txt', 'w') as file:
        file.write(str(total_count))
    print("Total Read Count:", total_count)


def insert_size_check(sample_id):
    samtools_dir = f'{data_dir}/3.SAMTOOLS'
    cmd=f'docker run -u {os.getuid()}:{os.getgid()} --rm -v /data/:/data/ {args.picard} java -jar /usr/local/share/picard-3.0.0-1/picard.jar CollectInsertSizeMetrics \
        I={samtools_dir}/{sample_id}.bwa.sorted.bam \
        O={samtools_dir}/{sample_id}.insert_size.txt \
        H={samtools_dir}/{sample_id}.insert_size.pdf \
        M=0.5'
    log_subprocess_output('picard insert size '+sample_id, cmd) 


def make_project_config(sample_id):
    sample_bicseq2_dir = f'{result_dir}/{sample_id}/bicseq2'
    os.makedirs(f'{sample_bicseq2_dir}', exist_ok=True)
    down_agilent_list = ['NB230700{}_12M'.format(x) for x in ['07','08','09','10','11','12']]
    if sample_id in down_agilent_list:
        frag_size = 'FRAG_SIZE=700'
    else: 
        frag_size = 'FRAG_SIZE=450'
    sample_config = f'{sample_bicseq2_dir}/{sample_id}_config.sh'
    lines = [
    f'OUTD_BASE="{sample_bicseq2_dir}"',
    'NORMD="$OUTD_BASE/norm"',
    'SEGD="$OUTD_BASE/segmentation"',
    'ANND="$OUTD_BASE/annotation"',
    'CHRLIST="/data/analysis/project/230412_nbs_heewon/bicseq2/chromosomes.dat"',
    'REF_CHR="/data/home/heewon/ref/HG19_CH"',
    'FA_CHR="${REF_CHR}/%s.fa"',
    'MAPD="/data/analysis/project/230412_nbs_heewon/bicseq2"',
    'MER="hg19.CRC.100mer"',
    'MAP_CHR="$MAPD/$MER/$MER.%s.txt"',
    f'SEQD="{sample_bicseq2_dir}/unique_mapping"',
    'SEQ_CHR="$SEQD/%s_%s.seq"',
    'SEQ_OUT="$SEQD/%s.seq"',
    'READ_LENGTH=150',
    f'{frag_size}',
    'BIN_SIZE=15000',
    'NORM_CHR="$NORMD/%s.%s.norm.bin"',
    'NORM_PDF="$NORMD/%s.GC.pdf"',
    'BICSEQ_NORM="/NBICseq-norm_v0.2.4/NBICseq-norm.pl"'
    ]
    with open(sample_config, 'w') as file:
        file.write('\n'.join(lines))
    file.close()

def get_unique_mapped_read(sample_id):
    samtools_dir = f'{data_dir}/3.SAMTOOLS'
    bicseq2_unique_dir = f'{result_dir}/{sample_id}/bicseq2/unique_mapping'
    print(f'{samtools_dir}/{sample_id}.bwa.sorted.bam')
    os.makedirs(f'{bicseq2_unique_dir}', exist_ok=True)
    for i in ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
              'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']:
        cmd = f'samtools view -@ 10 {samtools_dir}/{sample_id}.bwa.sorted.bam \
        {i} | perl /data/analysis/project/230412_nbs_heewon/bicseq2/samtools-0.1.7a_getUnique-0.1.3/misc/samtools.pl unique - | cut -f 4 > {bicseq2_unique_dir}/{sample_id}_{i}.seq'
        log_subprocess_output(f'samtool unique mapping {i} '+sample_id, cmd)


# 2. Use BICseq2-norm to remove the biases in the data.

def run_bicseq2_norm(sample_id):
    print(f'{result_dir}/{sample_id}/bicseq2/{sample_id}_config.sh')
    cmd = f'docker run -u 1004:1004 --rm -v /data/:/data/ {args.bicseq2} bash /BICSEQ2/src/run_norm.sh \
        {sample_id} {result_dir}/{sample_id}/bicseq2/{sample_id}_config.sh'
    info = log_subprocess_output('bicseq2 norm ' +sample_id , cmd)
    return info

# 3. Use BICseq2-seg to detect CNVs based on the normalized data.

def make_seg_config(sample_id):
    norm_dir = f'{result_dir}/{sample_id}/bicseq2/norm'
    chr_bin_list=[]
    for i in ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']:
        chr_bin_list.append([i,f'{norm_dir}/{sample_id}.{i}.norm.bin'])
    seq_config = pd.DataFrame(chr_bin_list,columns= ['ChromName','binFileNorm'])
    seq_config.to_csv(f'{result_dir}/{sample_id}/bicseq2/{sample_id}_seg_config.txt', sep='\t', index=False)
    print(f'{result_dir}/{sample_id}/bicseq2/{sample_id}_seg_config.txt')

def run_bicseq2_seg(sample_id):
    cmd =f'/data/analysis/project/230412_nbs_heewon/bicseq2/NBICseq-seg_v0.7.2/NBICseq-seg.pl \
        --tmp {result_dir}/{sample_id}/bicseq2/tmp\
        --fig {result_dir}/{sample_id}/bicseq2/{sample_id}_lambda_0_1.png\
        --lambda=0.1\
        {result_dir}/{sample_id}/bicseq2/{sample_id}_seg_config.txt \
        {result_dir}/{sample_id}/bicseq2/{sample_id}_lambda_0_1.cnv'
    info = log_subprocess_output('bicseq2 seg ' +sample_id, cmd)
    return info

def unique_rate(sample_id):
    fastp_dir = f'{data_dir}/1.FASTP'
    samtools_dir = f'{data_dir}/3.SAMTOOLS'
    os.makedirs(samtools_dir, exist_ok=True)
    sample_bicseq2_dir = f'{result_dir}/{sample_id}/bicseq2'
    os.system(f"awk '{{sum += $5}} END {{print sum}}' {sample_bicseq2_dir}/{sample_id}_lambda_0_1.cnv > {samtools_dir}/{sample_id}_unique_mapping.txt")
    with open(f'{samtools_dir}/{sample_id}_unique_mapping.txt', 'r') as f1,\
        open(f'{fastp_dir}/{sample_id}_total_count.txt', 'r') as f2,\
        open(f'{samtools_dir}/{sample_id}_unique_mapping_rate.txt', 'w') as output_file:
        value1 = float(f1.readline().strip())
        value2 = float(f2.readline().strip())
        result = (value1 / value2)*100
        print(f'{sample_id} {result}')
        output_file.write(str(result))
        output_file.close()

##########filitering blacklist & backdata and add gene &dgv info 

###기준 통과한 CNV만 남기기
def bicseq_call_cut(cnv_call,dup_cf,del_cf,sex=None):
    cnv_call.rename(columns={'log2.copyRatio':'log2_copyRatio'}, inplace = True)
    auto_chr_cutoff = (~cnv_call['chrom'].isin(['chrX', 'chrY'])) & ((cnv_call['log2_copyRatio'] >= dup_cf) | (cnv_call['log2_copyRatio'] <= del_cf))
    sex_chr_cutoff = (cnv_call['chrom'].isin(['chrX', 'chrY'])) & (
    (sex == 'M') & ((cnv_call['log2_copyRatio'] >= (dup_cf + (-1))) | (cnv_call['log2_copyRatio'] <= (del_cf + (-1)))) |
    (sex != 'M') & ((cnv_call['log2_copyRatio'] >= dup_cf) | (cnv_call['log2_copyRatio'] <= del_cf))
    )
    cnv_cut = cnv_call[auto_chr_cutoff|sex_chr_cutoff].copy()
    if len(cnv_cut)==0:
        return (cnv_cut)
    cnv_cut['pos'] = cnv_cut['start'].astype(str) + '-' + cnv_cut['end'].astype(str)
    cnv_cut['segment_size'] = cnv_cut['end'] - cnv_cut['start']
    cnv_cut['chrom'] = cnv_cut['chrom'].str.replace("chr", "")#.replace({'X': 23, 'Y': 24}).astype(int)
    cnv_cut['cnv_type'] = ['dup' if x >= dup_cf else 'del' for x in cnv_cut['log2_copyRatio']] 
    if sex =='M':
        cnv_cut.loc[(cnv_cut['chrom'].isin(['X', 'Y'])) & (cnv_cut['log2_copyRatio'] >= (dup_cf + (-1))),'cnv_type'] ='dup'
        cnv_cut.loc[(cnv_cut['chrom'].isin(['X', 'Y'])) & (cnv_cut['log2_copyRatio'] <= (del_cf + (-1))),'cnv_type'] ='del'
    cnv_cut_fin = cnv_cut[cnv_cut.columns.difference(['binNum','observed','expeceted'])]
    return(cnv_cut_fin)

###overlap check 
def calculate_overlap_cnv(start1, end1, start2, end2, cnv_type1, cnv_type2):
    if cnv_type1 != cnv_type2:
        overlap = 0
    else:
        overlap = max(0, min(end1, end2) - max(start1, start2))
    return overlap

def calculate_overlap(start1, end1, start2, end2):
    overlap = max(0, min(end1, end2) - max(start1, start2))
    return overlap

def overlap_percent(df1,df2):
    if len(df2)==0:
        return(df2)
    overlaps = []
    cyto_info =[]
    for row1 in df1.itertuples(index=False):
        for row2 in df2.itertuples(index=False):
            if row1.Chr != row2.chrom :
                continue
            overlap = calculate_overlap(row1.Start, row1.End, row2.start, row2.end)
            if overlap != 0:
                overlaps.append(row2)
                cyto_info.append(row1.Cytoband)
    overlaps_df = pd.DataFrame(overlaps)
    overlaps_df['cyto'] = cyto_info
    col_list = overlaps_df.columns.difference(['cyto']).to_list()
    cyto_df = overlaps_df.groupby(col_list)['cyto'].apply(','.join).reset_index()
    return(cyto_df)

##blacklist overlap 추가
def blacklist_overlap_info(blacklist,cnv_list):
    if len(cnv_list)==0:
        return(cnv_list)
    no_overlaps = []
    for bl in blacklist.itertuples(index=False):
        for cnv in cnv_list.itertuples(index=False):
            if (bl.chrom == cnv.chrom):
                overlap = calculate_overlap(bl.start, bl.end, cnv.start, cnv.end)
                if (overlap !=0):
                    no_overlaps.append({'chrom':cnv.chrom, 'start':cnv.start, 'end':cnv.end, 'log2_copyRatio':cnv.log2_copyRatio, 'pos':cnv.pos,
                        'cnv_type':cnv.cnv_type, 'segment_size':cnv.segment_size, 'blacklist_ovelap_p':(overlap/cnv.segment_size)*100,'blacklist_overlap':overlap,'blacklist_info':bl.region_info})
                elif (overlap ==0):
                    no_overlaps.append({'chrom':cnv.chrom, 'start':cnv.start, 'end':cnv.end, 'log2_copyRatio':cnv.log2_copyRatio, 'pos':cnv.pos,
                        'cnv_type':cnv.cnv_type, 'segment_size':cnv.segment_size, 'blacklist_ovelap_p':0,'blacklist_overlap':0,'blacklist_info':'-'})
    no_overlaps_df = pd.DataFrame(no_overlaps)
    no_overlaps_df.drop_duplicates(inplace=True)
    return(no_overlaps_df)

##normal overlap 추가
def normal_overlap_info(normal,cnv_list):
    if len(cnv_list)==0:
        return(cnv_list)
    no_overlaps = []
    for nor in normal.itertuples(index=False):
        for cnv in cnv_list.itertuples(index=False):
            if cnv.chrom not in normal['chrom']:
                no_overlaps.append({'chrom':cnv.chrom, 'start':cnv.start, 'end':cnv.end, 'log2_copyRatio':cnv.log2_copyRatio, 'pos':cnv.pos,
                        'cnv_type':cnv.cnv_type, 'segment_size':cnv.segment_size,'blacklist_info':cnv.blacklist_info,'blacklist_overlap_p':cnv.blacklist_ovelap_p,
                        'blacklist_overlap':cnv.blacklist_overlap,'backdata_info':'-','backdata_percent':'-','backdata_intersect_p':0,
                        'backdata_dgv_info':'-'})
            if (nor.chrom == cnv.chrom):
                overlap = calculate_overlap_cnv(nor.start, nor.end, cnv.start, cnv.end, nor.cnv_type, cnv.cnv_type)
                if (overlap !=0):
                    no_overlaps.append({'chrom':cnv.chrom, 'start':cnv.start, 'end':cnv.end, 'log2_copyRatio':cnv.log2_copyRatio, 'pos':cnv.pos,
                        'cnv_type':cnv.cnv_type, 'segment_size':cnv.segment_size,'blacklist_info':cnv.blacklist_info, 'blacklist_overlap_p':cnv.blacklist_ovelap_p,
                        'blacklist_overlap':cnv.blacklist_overlap,'backdata_info':nor.backdata_info,'backdata_percent':nor.normal_percent,
                        'backdata_intersect_p':(overlap/cnv.segment_size)*100,
                        'backdata_dgv_info':nor.dgv_info_all})
                elif (overlap ==0):
                    no_overlaps.append({'chrom':cnv.chrom, 'start':cnv.start, 'end':cnv.end, 'log2_copyRatio':cnv.log2_copyRatio, 'pos':cnv.pos, 
                        'cnv_type':cnv.cnv_type, 'segment_size':cnv.segment_size,'blacklist_info':cnv.blacklist_info,'blacklist_overlap_p':cnv.blacklist_ovelap_p,
                        'blacklist_overlap':cnv.blacklist_overlap,'backdata_info':'-','backdata_percent':'-','backdata_intersect_p':0,
                        'backdata_dgv_info':'-'})
    no_overlaps_df = pd.DataFrame(no_overlaps)
    no_overlaps_df.drop_duplicates(inplace=True)
    return(no_overlaps_df)



###calling cnv - overlap dgv frequency 추가

def max_frequency_check(df1, df2):
    df1['backdata_index'] = df1.index
    # Merge the two dataframes based on 'chrom' and 'variant_sub_type'
    merged_df = pd.merge(df1, df2, left_on=['chrom', 'cnv_type'], right_on=['chrom', 'variant_sub_type'], how='left')
    merged_df = merged_df[(~merged_df['dgv_start'].isnull().values) | (~merged_df['dgv_end'].isnull().values)]
    # Calculate overlap for each row
    merged_df['intersect'] = merged_df.apply(lambda row: calculate_overlap(row['start'], row['end'], row['dgv_start'], row['dgv_end']), axis=1)
    

    # Calculate intersect percent for each row
    merged_df['dgv_intersect_percent'] = (merged_df['intersect'] / merged_df['segment_size'])

    # Filter rows with intersect_percent >= 80
    filtered_df = merged_df.loc[merged_df['dgv_intersect_percent'] >= 0.8]
    # Group by 'backdata_info' and keep the row with the maximum 'frequency'
    high_intersect_backdata_info = filtered_df.loc[filtered_df.groupby('backdata_index')['frequency'].idxmax()]
    high_intersect_backdata_info['dgv_cutoff_pass'] = 'Y'
    # add backdata row with intersect percent <80
    row_intersect_backdata_info = merged_df[~merged_df['backdata_index'].isin(high_intersect_backdata_info['backdata_index'].tolist())]
    row_intersect_backdata_fin = row_intersect_backdata_info.loc[row_intersect_backdata_info.groupby('backdata_index')['dgv_intersect_percent'].idxmax()]
    row_intersect_backdata_fin.loc[row_intersect_backdata_fin['intersect']==0,'frequency']=0
    row_intersect_backdata_fin['dgv_cutoff_pass'] ='N'

    result_df = pd.concat([high_intersect_backdata_info,row_intersect_backdata_fin], sort=False)
    result_df['dgv_info'] = 'chr' +result_df['chrom']+': '+ result_df['dgv_start'].astype(str)+'-' +result_df['dgv_end'].astype(str)
    result_df.sort_values(by=['backdata_index'],inplace=True)


    return result_df

###gene 정보 추가 
def gene_overlap_info(gene_info, cnv_list):
    if cnv_list.empty:
        return cnv_list

    no_overlaps = []

    for cnv in cnv_list.itertuples(index=False):
        matching_chrom = gene_info[gene_info['chrom'] == cnv.chrom]

        if matching_chrom.empty:
            no_overlaps.append(create_overlap_entry(cnv))
        else:
            for gene in matching_chrom.itertuples(index=False):
                overlap = calculate_overlap(gene.start, gene.end, cnv.start, cnv.end)

                if overlap != 0:
                    no_overlaps.append(create_overlap_entry(cnv, gene, overlap))
                else:
                    no_overlaps.append(create_overlap_entry(cnv))

    no_overlaps_df = pd.DataFrame(no_overlaps)
    no_overlaps_df.drop_duplicates(inplace=True)
    return no_overlaps_df

def create_overlap_entry(cnv, gene=None, overlap=0):
    entry = {
        'chrom': cnv.chrom,
        'start': cnv.start,
        'end': cnv.end,
        'log2_copyRatio': cnv.log2_copyRatio,
        'pos': cnv.pos,
        'cnv_type': cnv.cnv_type,
        'segment_size': cnv.segment_size,
        'cyto': cnv.cyto,
        'blacklist_info': cnv.blacklist_info,
        'blacklist_overlap_p': cnv.blacklist_overlap_p,
        'blacklist_overlap': cnv.blacklist_overlap,
        'backdata_info': cnv.backdata_info,
        'backdata_percent': cnv.backdata_percent,
        'backdata_intersect_p': cnv.backdata_intersect_p,
        'backdata_dgv_info': cnv.backdata_dgv_info,
        'Gene_name':'-',
        'Gene_info':'-'
    }

    if gene:
        entry['Gene_name'] = gene.Gene_name
        entry['Gene_info'] = gene.Gene_info

    return entry

#######연결되는 CNV 합치기 
def can_merge_rows(row1, row2):
    return (
        row2['chrom'] == row1['chrom']
        and row2['start'] == row1['end'] + 1
        and row2['cnv_type'] == row1['cnv_type']
    )

def merge_rows(row1, row2):
    new_row = row1.copy()
    new_row['end'] = row2['end']
    new_row['log2_copyRatio'] = ((row1['end'] - row1['start']) / (row2['end'] - row1['start'])) * row1['log2_copyRatio'] + \
                                ((row2['end'] - row2['start']) / (row2['end'] - row1['start'])) * row2['log2_copyRatio']
    new_row['segment_size'] = row1['segment_size'] +  row2['segment_size'] + 1
    return new_row

def continuous_row_merge(bicseq_result):
    bicseq_fin = bicseq_result.sort_values(['chrom','start'])
    bicseq_fin.reset_index(drop=True,inplace=True)
    for i in range(1,len(bicseq_fin)):
        prev_row = bicseq_fin.loc[i-1]
        curr_row = bicseq_fin.loc[i]
        if can_merge_rows(prev_row, curr_row):
            merged_row = merge_rows(prev_row, curr_row)
            bicseq_fin.drop(index=[i-1,i], inplace=True)
            bicseq_fin.loc[i]= merged_row
        else:
            i += 1
    bicseq_fin2 = bicseq_fin.sort_index()
#    bicseq_cut = bicseq_fin[(bicseq_fin['cnv_type'].isin(['dup','del']))]
    return bicseq_fin2
###############


dgv_hg19 = pd.read_csv('/data/analysis/project/230412_nbs_heewon/ref/DGV_hg19_essential_info.csv')

def blacklist_n_back_filitering(sample_id):
    sample_id_origin = sample_id
    #sample sex info
    sample_info = pd.read_csv(f'/data/analysis/project/230412_nbs_heewon/raw_data/230727_A01980_0012_BH5LNGDRX3/sample_info.csv')
    
    lambda_value = '0_1'
    #back data 
    bd_dgv = pd.read_excel(f'/data/analysis/project/230412_nbs_heewon/Result/230727_A01980_0012_BH5LNGDRX3/backdata_based_0_8_add_dgv.xlsx')
    bd_dgv['chrom'] = bd_dgv['chrom'].astype(str)
    #black list 
    blacklist = pd.read_csv('/data/home/heewon/ref/hg19-blacklist.v2.bed', index_col = 0) 
    blacklist['chrom'] = blacklist['chr'].str.replace("chr", "")
    blacklist['region_info'] = 'chr'+blacklist['chrom']+':'+blacklist['start'].astype(str)+'-'+blacklist['end'].astype(str)+' ('+blacklist['info']+')'
    
    #whole chromosome list
    chrom_list2 = list(map(str,[*range(1,23)]+['X','Y']))
    #cut off list
    black_cutoff= 35
    back_cutoff = 15

    #dgv hg19 info

    
    #cytoband
    cytoband = pd.read_csv("/data/home/heewon/ref/cytoBand.csv")
    cytoband['Chr'] = cytoband['Chr'].astype(str)

    #gene info 
    gene_info = pd.read_csv('/data/analysis/project/230412_nbs_heewon/ref/bio_mart_source_info_exist_gene.csv')
    gene_info['chrom'] = gene_info['chrom'].astype(str)
    sample_sex_info =sample_info[sample_info['sample_id'].str.contains(sample_id)]['sex'].values[0]

    bic_result = pd.read_csv(f'{result_dir}/{sample_id_origin}/bicseq2/{sample_id_origin}_lambda_{lambda_value}.cnv',sep='\t')
    sample_cnv = bicseq_call_cut(bic_result,math.log2(2.6/2),math.log2(1.4/2),sex=sample_sex_info)
    dup_cnv = sample_cnv[sample_cnv['cnv_type']=='dup']
    del_cnv = sample_cnv[sample_cnv['cnv_type']=='del']
    dup_fin = continuous_row_merge(dup_cnv)
    del_fin = continuous_row_merge(del_cnv)
    merged_cnv = pd.concat([dup_fin,del_fin])
    
    sample_cnv_size_cut = merged_cnv[merged_cnv['segment_size']>=200000]
    no_bl_cut = blacklist_overlap_info(blacklist,sample_cnv_size_cut)
    bl_with_non_overlap = no_bl_cut[(no_bl_cut[['chrom', 'start', 'end', 'log2_copyRatio', 'pos', 'segment_size']].duplicated(keep=False)) & (no_bl_cut['blacklist_info']=='-')].index
    no_bl_cut.drop(bl_with_non_overlap,inplace=True)
    no_bl_cut_result = no_bl_cut.groupby(['chrom', 'start', 'end', 'log2_copyRatio', 'pos', 'segment_size', 'cnv_type'])\
        .agg({'blacklist_info': ','.join, 'blacklist_overlap': sum,'blacklist_ovelap_p':sum})\
        .reset_index()
    normal_bl_info_add = normal_overlap_info(bd_dgv,no_bl_cut_result)
    normal_with_non_overlap = normal_bl_info_add[(normal_bl_info_add[['chrom', 'start', 'end', 'log2_copyRatio', 'pos', 'segment_size']].duplicated(keep=False)) & (normal_bl_info_add['backdata_percent']==0)].index
    normal_bl_info_add.drop(normal_with_non_overlap,inplace=True)
    normal_bl_info_add['backdata_percent'] =normal_bl_info_add['backdata_percent'].astype(str)
    normal_bl_info_add['chrom'] = pd.Categorical(normal_bl_info_add['chrom'], 
    categories=chrom_list2, 
    ordered=True)
    n_fin  = overlap_percent(cytoband,normal_bl_info_add)
    n_fin['segment_size_KB'] = n_fin['segment_size'].apply(lambda x: f"{int(x / 1000)}KB" if x % 1000 == 0 else f"{x / 1000:.1f}KB")
    n_fin_add_bd = n_fin.groupby(['chrom', 'start', 'end', 'log2_copyRatio', 'pos', 'segment_size','segment_size_KB', 'cnv_type','cyto','blacklist_info','blacklist_overlap','blacklist_overlap_p'])\
        .agg({'backdata_info':','.join,'backdata_dgv_info':','.join,'backdata_percent': ','.join, 'backdata_intersect_p': sum})\
        .reset_index()
    n_fin_add_bd['backdata_dgv_info'] = n_fin_add_bd['backdata_dgv_info'].str.replace("-,", "")
    n_fin_add_bd['backdata_percent'] = n_fin_add_bd['backdata_percent'].str.replace("0,", "")
    n_fin_add_bd['backdata_info'] = n_fin_add_bd['backdata_info'].str.replace("-,", "")

    n_fin_add_bd['chrom'] = pd.Categorical(n_fin_add_bd['chrom'], 
    categories=chrom_list2, 
    ordered=True)
    n_fin_add_bd.sort_values(by=['chrom','start'],inplace=True)

    n_fin_add_bd_cut = n_fin_add_bd

    if n_fin_add_bd_cut.empty: 
        empty_col_list = n_fin_add_bd_cut.columns.to_list()
        empty_df = pd.DataFrame([pd.Series(['-']*len(empty_col_list),index=empty_col_list)])
        empty_df.to_excel(f'{table_all_dir}/{sample_id_origin}_lambda_0_1_filtering.xlsx',index=None) 
    else:
        g_fin = gene_overlap_info(gene_info, n_fin_add_bd_cut)
        g_fin_col_list = g_fin.columns.difference(['Gene_name','Gene_info']).to_list()
        g_fin_df = g_fin.groupby(g_fin_col_list)\
        .agg({'Gene_name': ','.join, 'Gene_info': ','.join}).reset_index()
        fin_df_col_list = ['chrom','start','end','log2_copyRatio','pos','cnv_type','segment_size','cyto','blacklist_info','blacklist_overlap_p','blacklist_overlap',
           'backdata_info','backdata_percent','backdata_intersect_p','backdata_dgv_info','Gene_name','Gene_info']
        g_fin_df2 = g_fin_df[fin_df_col_list]

        g_fin_df2['Gene_name'] = g_fin_df2['Gene_name'].str.replace("-,", "")
        g_fin_df2['Gene_info'] = g_fin_df2['Gene_info'].str.replace("-,", "")
        dgv_info_add_fin_df = max_frequency_check(g_fin_df2,dgv_hg19)
        dgv_info_add_fin_df.drop_duplicates(inplace=True)
        dgv_info_add_fin_df['frequency'] = dgv_info_add_fin_df['frequency'].astype(str)
        dgv_info_add_fin_df2 = dgv_info_add_fin_df.groupby(fin_df_col_list)\
        .agg({'dgv_info': ','.join, 'frequency': ','.join,'dgv_cutoff_pass':','.join}).reset_index()
        dgv_info_add_fin_df2['chrom'] = pd.Categorical(dgv_info_add_fin_df2['chrom'], 
        categories=chrom_list2, 
        ordered=True)
        dgv_info_add_fin_df2.sort_values(by=['chrom','start'],inplace=True)
        dgv_info_add_fin_df2.rename(columns={'blacklist_overlap_p':'blacklist_intersect_p','blacklist_overlap':'blacklist_intersect_bp','backdata_percent':'backdata_include_sample_p'},inplace=True)
        dgv_info_add_fin_df2.loc[dgv_info_add_fin_df2['backdata_include_sample_p']=='1','backdata_include_sample_p'] = '83/83'
        dgv_info_add_fin_df2.to_excel(f'{table_all_dir}/{sample_id_origin}_lambda_{lambda_value}_filtering.xlsx',index=None) 
        dgv_info_add_fin_df2['sample_id'] = sample_id
        dgv_info_add_fin_df2['sex'] = sample_sex_info
        all_sample_result.append(dgv_info_add_fin_df2)

    print(f'{sample_id_origin} sample filitering & add info finish')

all_sample_result = []
def operation(sample_id, sample_per_threads):
    try:
        sample_id_origin = re.match(r'^(.*)_S.*', sample_id).group(1)
        fastp_run(sample_id, sample_per_threads)
        run_bwa(sample_id_origin, sample_per_threads)
        run_samtools(sample_id_origin)
        mapping_quality(sample_id_origin)
        fastq_total_count(sample_id_origin)
        make_project_config(sample_id_origin)
        get_unique_mapped_read(sample_id_origin)
        run_bicseq2_norm(sample_id_origin)
        make_seg_config(sample_id_origin)
        run_bicseq2_seg(sample_id_origin)
        insert_size_check(sample_id_origin)
        unique_rate(sample_id_origin)
        blacklist_n_back_filitering(sample_id_origin) 
    except: 
        print(f'error {sample_id_origin}')

if __name__ == "__main__":
    operation(args.sample_id, 6)

def main():
    init_env_check()
    sample_id_list = data_align()
    if len(sample_id_list) < args.threads :
        sample_per_threads = int(args.threads / len(sample_id_list))
    else :
        sample_per_threads = 20
    executor = concurrent.futures.ProcessPoolExecutor(args.threads)
    futures = [executor.submit(operation, f'{sample_id}', sample_per_threads = sample_per_threads) for sample_id in sample_id_list]
#230612_A01980_0011_AH5LNJDRX3
main()

###########all_sample_info_xlsx
all_info = pd.concat(all_sample_result)
col1= all_info.columns[-2:].to_list()
col2= all_info.columns[:-2].to_list()
new_col = col1+col2
all_info[new_col].to_excel(f'{table_all_dir}/{args.sequencing_id}_all_sample.xlsx',index=None)