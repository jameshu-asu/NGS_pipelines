import os
import subprocess

from Bio import SeqIO
from itertools import zip_longest

'''by James C. Hu
This scipt contains functions for commonly used CLI tools.
'''


def gunzip(infile: str) -> None:
    '''This function will:
    Unzip the input file.

    infile: filename ending in .gz
    '''
    subprocess.call(['gunzip', infile])
    return None


def fq_2_fa(fq_path:str, fa_path:str) -> None:
    '''This function will generate a fasta file from fastq using biopython.
    fq_path: input path to fastq file.
    fa_path: output path for fasta file.
    '''
    with open(fq_path, 'r') as f1, open(fa_path, 'w') as f2:
        sequences = SeqIO.parse(f1, 'fastq')
        SeqIO.write(sequences, f2, 'fasta')
    return None


def bbduk_adaptTrim(r1: str, r2: str, seq_id: str) -> None:
    '''This function will:
    Perfrom adapter trimming using bbduk

    r1: forward fastq
    r2: reverse fastq
    seq_id: sequence id

    Docs: https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/
    '''
    ref = '/mnt/storage/jupiter-data/databases/contaminants/illumina_dna_prep_abcd.fa'
    out1 = f'./{seq_id}/{seq_id}-adaptTrimQC_R1_2023.fastq'
    out2 = f'./{seq_id}/{seq_id}-adaptTrimQC_R2_2023.fastq'
    log_file = f'./{seq_id}/{seq_id}-adaptTrimQC_2023.log.txt'
    command = f'bbduk.sh in={r1} in2={r2} ref={ref} out={out1} out2={out2} k=19 hdist=1 ktrim=r qtrim=rl mink=11 trimq=30 minlength=75 minavgquality=20 removeifeitherbad=f ottm=t tpe=t overwrite=t 1> {log_file} 2>&1'
    subprocess.call(command, shell=True)
    return None


def bbduk_phix_removal(r1: str, r2: str, seq_id: str) -> None:
    '''This funciton will:
    Perform contaminant removal using bbduk.
        - Removes phix reads

    r1: forward fastq
    r2: reverse fastq
    seq_id: sequence id

    Docs: https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/
    '''
    ref = '/mnt/storage/jupiter-data/databases/contaminants/phix174_ill.ref.fa.gz'
    out1 = f'./{seq_id}/{seq_id}-R1-phixRemoved_2023.fastq'
    out2 = f'./{seq_id}/{seq_id}-R2-phixRemoved_2023.fastq'
    log_file = f'./{seq_id}/{seq_id}phixRemoved_2023.log.txt'
    command = f'bbduk.sh in={r1} in2={r2} ref={ref} out={out1} out2={out2} k=31 hdist=1 overwrite=t 1> {log_file} 2>&1'
    subprocess.call(command, shell=True)
    return None


def bbduk_filter(r1: str, seq_id: str) -> None:
    '''This function will:
    Perform a minimum read length filter on an input fastq files using bbduk.

    r1: fastq file.
    seq_id: sequence id

    Docs: https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/
    '''
    out = f'./{seq_id}/{seq_id}firstDeduplication_filtered_2023.fastq'
    log_file = f'./{seq_id}/{seq_id}firstDeduplication_filtered_2023.log.txt'
    command = f'bbduk.sh in={r1} out={out} minlength=75 overwrite=t 1> {log_file} 2>&1'
    subprocess.call(command, shell=True)
    return None


def bbmap(r1: str, R: str, seq_id: str) -> None:
    '''This function will:
    Perform read mapping to the host genome and partition out the mapped reads.

    r1: fastq
    seq_id: sequence id
    R: Forward or Reverse read.
        - Only options are 'R1' for forward or 'R2' for reverse.
        - INPUTS ARE CASE SENSITIVE.

    Docs: https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/
    '''
    ref = f'/mnt/storage/jupiter-data/databases/Human_GRCh38'
    outu = f'{seq_id}/{seq_id}-dehosted_{R}.fastq'
    outm = f'{seq_id}/{seq_id}-hostmatched_{R}.fastq'
    log_file = f'{seq_id}/{seq_id}dehost_{R}_2023.log.txt'
    command = f'bbmap.sh minid=.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 -Xmx64g path={ref} in={r1} outu={outu} outm={outm} 1> {log_file} 2>&1'
    subprocess.call(command, shell=True)
    return None


def concat_fastq(r1: str, r2: str, seq_id: str, out: str) -> None:
    '''This function will:
    Concatenate fastq files.

    r1: fastq file 1
    r2: fastq file 2
    seq_id: Sequence ID
    out: name of outfile
    '''
    command = f'cat {r1} {r2} > {out}'
    subprocess.call(command, shell=True)
    return None


def bbmerge(r1: str, seq_id: str) -> None:
    '''This funciton will:
    Merge overlapping reads within a fastq file.

    r1: fastq
    seq_id: sequnce id

    Docs: https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmerge-guide/
    '''
    out = f'{seq_id}/{seq_id}dehosted_cat_Merged_2023.fastq'
    outu = f'{seq_id}/{seq_id}dehosted_cat_UnMerged_2023.fastq'
    log_file = f'{seq_id}/{seq_id}dehosted_cat_Merged_2023.log.txt'
    command = f'bbmerge.sh in={r1} out={out} outu={outu} 1> {log_file} 2>&1'
    subprocess.call(command, shell=True)
    return None


def dedupe(r1: str, seq_id: str) -> None:
    '''This funciton will:
    Deduplicate an input fastq files.

    r1: fastq
    seq_id: Sequence ID

    Docs: https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/dedupe-guide/
    '''
    out = f'{seq_id}/{seq_id}firstDeduplication_2023.fastq'
    outd = f'{seq_id}/{seq_id}firstDuplication_2023.fastq'
    log_file = f'{seq_id}/{seq_id}firstDeduplication_2023.log.txt'
    command = f'dedupe.sh ac=f -Xmx180g in={r1} out={out} outd={outd} csf=dedupe.cluster.stats overwrite=t minidentity=100 ac=f 1> {log_file} 2>&1'
    subprocess.call(command, shell=True)
    return None


def stream_editor(r1: str, seq_id: str) -> None:
    '''This funciton will:
    This function converts fastq to fasta using Linux stream editor command.

    r1: fastq
    seq_id: seqquence ID

    Docs: https://www.gnu.org/software/sed/manual/sed.html
    '''
    out = f'{seq_id}/{seq_id}firstDeduplication_filtered_2023.fasta'
    command = f"sed -n '1~4s/^@/>/p;2~4p' {r1} > {out}"
    subprocess.call(command, shell=True)
    return None


def bowtie2(r1: str, seq_id: str) -> None:
    '''This funciton will:
    Generates a SAM file from a fasta file.
        - Could also use fastq file.

    r1: fasta file
    seq_id: sequence id

    Docs: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
    '''
    ref = '/mnt/storage/jupiter-data/databases/RVP_DB/RVP_DB'
    out = f'bowtie2_mapped_2023/{seq_id}_contigsMapped.sam'
    command = f'bowtie2 -x {ref} -f {r1} > {out}'
    subprocess.call(command, shell=True)
    return None


def samtools_view(flag1: str, flag2: str, flag3: str, r1: str, out: str, seq_id: str) -> None:
    '''This funciton will:
    Use samtools to read and extract data from SAM/BAM/CRAM files.

    flag_n: flags used for samtools view
    r1: SAM file
    out: outfile name
    seq: sequence id

    FLAGS:
    0x1         PAIRED          paired-end (or multiple-segment) sequencing technology
    0x2         PROPER_PAIR     each segment properly aligned according to the aligner
    0x4         UNMAP           segment unmapped
    0x8         MUNMAP          next segment in the template unmapped
    0x10        REVERSE         SEQ is reverse complemented
    0x20        MREVERSE        SEQ of the next segment in the template is reverse complemented
    0x40        READ1           the first segment in the template
    0x80        READ2           the last segment in the template
    0x100       SECONDARY       secondary alignment
    0x200       QCFAIL          not passing quality controls
    0x400       DUP             PCR or optical duplicate
    0x800       SUPPLEMENTARY   supplementary alignment

    Docs: http://www.htslib.org/doc/samtools-view.html
    '''
    command = f'samtools view {flag1} {flag2} {flag3} {r1} > {out}'
    subprocess.call(command, shell=True)
    return None


def samtools_sort(r1: str, out: str) -> None:
    '''This funciton will:
    Sort an input SAM/BAM/CRAM file.

    r1: SAM/BAM/CRAM
    Out: output file name.

    Docs: http://www.htslib.org/doc/samtools-sort.html
    '''
    command = f'samtools sort {r1} > {out}'
    subprocess.call(command, shell=True)
    return None


def samtools_index(r1: str) -> None:
    '''This funciton will:
    Index a SAM/BAM/CRAM file

    r1: SAM/BAM/CRAM

    Docs: http://www.htslib.org/doc/samtools-index.html
    '''
    command = f'samtools index {r1}'
    subprocess.call(command, shell=True)
    return None


def samtools_idxstats(r1: str, out: str) -> None:
    '''This funciton will:
    Report alignemtn summary stats:

    r1: SAM/BAM/CRAM
    out: Name of outfile.

    Docs: http://www.htslib.org/doc/samtools-idxstats.html
    '''
    command = f'samtools idxstats {r1} > {out}'
    subprocess.call(command, shell=True)
    return None


def copy_files(file: str, destination: str) -> None:
    '''This funciton will:
    Perform Linux copy funciton.

    file: target file
    destination: file destination
    '''
    command = f'cp {file} {destination}'
    subprocess.call(command, shell=True)
    return None


def move_file(file: str, destination: str) -> None:
    '''This funciton will:
    Perform Linux move command.

    file: target file
    destination: file destination
    '''
    command = f'mv {file} {destination}'
    subprocess.call(command, shell=True)
    return None


seq_id_list = list(set([file[:7] for file in os.listdir('fastq') if file.startswith('I')]))
seq_id_list.sort()

r1_fq_list = [file for file in os.listdir('fastq') if 'R1' in file]
r1_fq_list.sort()
r2_fq_list = [file for file in os.listdir('fastq') if 'R2' in file]
r2_fq_list.sort()


for i, j in zip_longest(r1_fq_list, r2_fq_list):
    if i is None:
        print(f'Missing paired reads for sample {j}')
    if j is None:
        print(f'Missing paired reads for sample {i}')


os.mkdir('logFiles_2023')
os.mkdir('qc_fasta_2023')
os.mkdir('qc_fastq_2023')
os.mkdir('bowtie2_mapped_2023')
os.mkdir('counts_2023')
for seq_id in seq_id_list:
    os.mkdir(seq_id)

fastq_path = './fastq'


for seq_id, r1, r2, in zip(seq_id_list, r1_fq_list, r2_fq_list):
    try:
        gunzip(f'{os.path.join(fastq_path, r1)}')
        gunzip(f'{os.path.join(fastq_path, r2)}')
    except Exception as e:
        print(e)
        pass
    bbduk_adaptTrim(f'{fastq_path}/{r1}', f'{fastq_path}/{r2}', seq_id)
    bbduk_phix_removal(f'{seq_id}/{seq_id}-adaptTrimQC_R1_2023.fastq', f'{seq_id}/{seq_id}-adaptTrimQC_R2_2023.fastq', seq_id)
    bbmap(f'{seq_id}/{seq_id}-R1-phixRemoved_2023.fastq', 'R1', seq_id)
    bbmap(f'{seq_id}/{seq_id}-R2-phixRemoved_2023.fastq', 'R2', seq_id)
    concat_fastq(f'{seq_id}/{seq_id}-dehosted_R1.fastq', f'{seq_id}/{seq_id}-dehosted_R2.fastq', seq_id, f'{seq_id}/{seq_id}-dehosted_cat.fastq')
    bbmerge(f'{seq_id}/{seq_id}-dehosted_cat.fastq', seq_id)
    concat_fastq(f'{seq_id}/{seq_id}dehosted_cat_Merged_2023.fastq', f'{seq_id}/{seq_id}dehosted_cat_UnMerged_2023.fastq', seq_id, f'{seq_id}/{seq_id}dehosted_cat_Merged_UnMerged_2023.fastq')
    dedupe(f'{seq_id}/{seq_id}dehosted_cat_Merged_UnMerged_2023.fastq', seq_id)
    bbduk_filter(f'{seq_id}/{seq_id}firstDeduplication_2023.fastq', seq_id)
    stream_editor(f'{seq_id}/{seq_id}firstDeduplication_filtered_2023.fastq', seq_id)
    bowtie2(f'{seq_id}/{seq_id}firstDeduplication_filtered_2023.fasta', seq_id)
    samtools_view('-h', '-F', '0x900', f'./bowtie2_mapped_2023/{seq_id}_contigsMapped.sam', f'./bowtie2_mapped_2023/{seq_id}_secondaryRemoved.sam', seq_id)
    samtools_view('-h', '-F', '0x4', f'./bowtie2_mapped_2023/{seq_id}_secondaryRemoved.sam', f'./bowtie2_mapped_2023/{seq_id}_secondaryUnMappedRemoved.sam', seq_id)
    samtools_view('-S', '-b', '', f'./bowtie2_mapped_2023/{seq_id}_secondaryUnMappedRemoved.sam', f'./bowtie2_mapped_2023/{seq_id}_secondaryUnMappedRemoved.bam', seq_id)
    samtools_sort(f'bowtie2_mapped_2023/{seq_id}_secondaryUnMappedRemoved.bam', f'bowtie2_mapped_2023/{seq_id}_secondaryUnMappedRemoved_sorted.bam')
    samtools_index(f'bowtie2_mapped_2023/{seq_id}_secondaryUnMappedRemoved_sorted.bam')
    samtools_idxstats(f'bowtie2_mapped_2023/{seq_id}_secondaryUnMappedRemoved_sorted.bam', f'bowtie2_mapped_2023/{seq_id}.txt')


for seq_id in seq_id_list:
    file_copy_list = [file for file in os.listdir(seq_id)]
    for file in file_copy_list:
        if file.endswith('2023.log.txt'):
            copy_files(f'{seq_id}/{file}', './logFiles_2023')
        elif file.endswith('fasta'):
            copy_files(f'{seq_id}/{file}', './qc_fasta_2023')
        elif file.endswith('fastq'):
            copy_files(f'{seq_id}/{file}', './qc_fastq_2023')
        elif file == f'{seq_id}.txt':
            move_file(file, './counts_2023')
