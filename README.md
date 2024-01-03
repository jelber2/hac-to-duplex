# hac-to-duplex

Playing around with [vechat's racon](https://github.com/jelber2/vechat/) for overlap-based error-correction, [Brutal rewrite](https://github.com/natir/br) for kmer-based error-correction, and [peregrine-2021](https://github.com/cschin/peregrine-2021) for overlap-based error-correction on **Human** Ultra-long (N50 > 100Kbp) Nanopore Simplex reads called with dna_r10.4.1_e8.2_400bps_sup_v4.3.0 on dorado v0.5.0.

![plot](https://github.com/jelber2/hac-to-duplex/blob/main/chr20_hac_vs_sup_vs_dup_vs_sup-error-correct.svg)


input is GAF file of HAC reads mapped with minigraph (https://github.com/lh3/minigraph)
to chromosomes of human pangenome of the t2t human assembly as reference and other asssemblies
made with minigraph (https://zenodo.org/records/6983934)

## Initial pre-pipeline steps

### get pangenome
```bash
wget https://zenodo.org/records/6983934/files/chm13-90c.r518.gfa.gz
```

### run minigraph

```bash
~/bin/minigraph/minigraph -t 32 -cx lr ../sandbox3/chm13-90c.r518.gfa.gz /mnt/share/nanopore/Tamara/BIonano-Test/SP62_uHMW_25072023/20230725_1206_P2S-00581-A_PAK86034_9e5850fe/basecalling_dorado/SP62_ULK114.fastq.gz 2> SP62_uHMW_25072023_minigraph.log | pigz -p 12 > SP62_uHMW_25072023.gaf.gz &
```

### Separate the gaf.gz file into separate chromosomes

```bash
mkdir -p SP62_uHMW_25072023_pangenome/gaf
while read i
do
  zcat SP62_uHMW_25072023.gaf.gz |awk -v chromosome="$i" '$6 ~ chromosome":" || $6 == chromosome' OFS='\t' |pigz -p 12 > SP62_uHMW_25072023_pangenome/gaf/${i}.gaf.gz
done < <(cat <(seq 1 22) <(echo -e "M\nX")|perl -pe "s/^/chr/g")
```

### Convert FAST5 files to BLOW5 files then merge them and index them

uses [slow5tools](https://github.com/hasindu2008/slow5tools)

```bash
export HDF5_PLUGIN_PATH=/home/jelber43/.local/hdf5/lib/plugin

~/bin/slow5tools-v1.1.0/slow5tools f2s -d 20230725_1206_P2S-00581-A_PAK86034_9e5850fe -c zstd -s svb-zd /mnt/share/nanopore/Tamara/BIonano-Test/SP62_uHMW_25072023/20230725_1206_P2S-00581-A_PAK86034_9e5850fe/fast5_pass/ > SP62_uHMW_25072023_20230725_1206_P2S-00581-A_PAK86034_9e5850fe.blow5.log 2>&1 &
```

### merge and index blow5 file

```bash
~/bin/slow5tools-v1.1.0/slow5tools merge -t 32 -c zstd -s svb-zd 20230725_1206_P2S-00581-A_PAK86034_9e5850fe/*.blow5 \
--output SP62_uHMW_25072023.blow5 > slow5toolsMerge.log 2>&1 && ~/bin/slow5tools-v1.1.0/slow5tools index SP62_uHMW_25072023.blow5 > slow5toolsIndex.log 2>&1 &
```

### delete the directories containing the blow5 individual files

```bash
rm -r SP62_uHMW_25072023/
```

## Snakemake file

Here it is as a diagram.

![plot](https://github.com/jelber2/hac-to-duplex/blob/main/test.svg)


for the following file, you need to manually set the following values for your system, file names, etc.

`MEMORY_SIZE1="300G"` # maximum memory

`MEMORY_SIZE2="-Xmx300g"` # maximum memory

`BLOW5="../SP62_uHMW_25072023.blow5"` # location of your BLOW5 file

`CWD="/sandbox3/SP62_uHMW_25072023_pangenome"` # current working directory

put your gzipped `gaf` (see above) in `CWD/gaf` in the format `{id}.gaf.gz`, where `{id}` is for example `chr17`

the Snakemake Snakefile:

`/home/jelber43/sandbox3/SP62_uHMW_25072023_pangenome/hac-to-corrected-compared-to-duplex2.smk`

```yaml
#
shell.executable("/bin/bash")

# set input path here
IDS, = glob_wildcards("gaf/{id}.gaf.gz")

# set number of threads here
THREADS=48
THREADS2=range(1,THREADS+1)
THREADS3=expand(["{threads}"], threads=THREADS2)
THREADS4=[str(item).zfill(2) for item in THREADS3]
MEMORY_SIZE1="300G"
MEMORY_SIZE2="-Xmx300g"
BLOW5="../SP62_uHMW_25072023.blow5"
CWD="/sandbox3/SP62_uHMW_25072023_pangenome"

scattergather:
    split=8

# list of rules which are not deployed to slurm
localrules: all, reads1, reads2, reads3

rule all:
        input:
            input1=expand("raft4/finalasm_{id}.bp.hap1.p_ctg.fa", id=IDS),
            input2=expand("raft4/finalasm_{id}.bp.hap2.p_ctg.fa", id=IDS)


# make a list of all reads belonging to a specific chromosome
rule gaf:
        input: "gaf/{id}.gaf.gz"
        output: "fast5/{id}_read_names.txt"
        params:
            threads = THREADS,
            memory = MEMORY_SIZE1
        shell: """zcat {input} |cut -f 1|sort --parallel={params.threads} -S {params.memory} > {output}"""


# use slow5tools get to retrieve reads from a pangenome chromosome
rule slow5tools:
        input:
            reads="fast5/{id}_read_names.txt",
            blow5=BLOW5
        output: "fast5/{id}.blow5"
        params: THREADS
        shell: """
        /home/git/slow5tools-v1.1.0/slow5tools get --skip -t {params} --list {input.reads} --output {output} -c zstd -s svb-zd  {input.blow5}
        """


# use blue-crab to convert blow5 to pod5
rule blue_crab:
        input: "fast5/{id}.blow5"
        output: "fast5/{id}.pod5"
        shell: '''
        eval "$(/home/.local/bin/micromamba shell hook --shell bash)"
        micromamba activate blue-crab
        blue-crab s2p {input} --output {output}
        micromamba deactivate
        '''


# basecall with SUP with dorado 0.5.0 and SUP
rule dorado:
        input: "fast5/{id}.pod5"
        output: "fasta/{id}.bam"
        shell: """/home/git/dorado-0.5.0-linux-x64/bin/dorado basecaller sup {input} > {output}"""


# basecall with duplex  with dorado 0.5.0 and SUP
rule dorado2:
        input: "fast5/{id}.pod5"
        output: "fasta/{id}.duplex.bam"
        params: THREADS
        shell: """/home/git/dorado-0.5.0-linux-x64/bin/dorado duplex sup {input} > {output}"""


# convert BAM to FASTQ
rule bam2fastq:
        input:
            orig = "fasta/{id}.bam",
            duplex = "fasta/{id}.duplex.bam"
        output:
            orig =  "fasta/{id}.fasta.gz",
            duplex = "fasta/{id}.duplex.fasta.gz"
        params: THREADS
        shell: '''
        eval "$(/home/.local/bin/micromamba shell hook --shell bash)"
        micromamba activate samtools
        samtools fasta -@ {params} {input.orig} | pigz -p {params} > {output.orig}
        samtools fasta -T "dx" -@ {params} {input.duplex} | /home/git/seqtk/seqtk seq -l0 -A | paste - - | fgrep "dx:i:1"| tr '\t' '\n' | pigz -p {params} > {output.duplex}
        micromamba deactivate
        '''


# run mm2-fast to overlap reads and filter overlaps with fpa
rule mm2_fast:
        input: "fasta/{id}.fasta.gz"
        output:
            paf = "vechat/{id}.paf.gz",
            mm2_fast = temp(multiext("fasta/{id}.fasta.gz_ava-ont_minimizers_key_value_sorted", "_keys.rmi_PARAMETERS",
                                "_pos_bin", "_val_bin", "_size", "_keys.uint64")),
            mm2_fast2 = temp(directory("{id}_2"))
        params:
            CWD=CWD,
            threads=THREADS
        shell: """
        eval "$(/home/.local/bin/micromamba shell hook --shell bash)"
        micromamba activate git
        mkdir -p {wildcards.id}_2
        cd {wildcards.id}_2
        rm -fr mm2-fast
        git clone --recursive https://github.com/bwa-mem2/mm2-fast.git mm2-fast
        cd mm2-fast
        git checkout 830e8c7
        ./build_rmi.sh {params.CWD}/{input} ava-ont
        make clean && make lhash=1 -j {params.threads}
        cp -r ext/TAL/scripts/ scripts/
        perl -pi -e "s/ext\/build-rmi\/learned-systems-rmi/ext\/TAL\/ext\/build-rmi\/learned-systems-rmi/g" scripts/build-rmi.linear_spline.linear.sh
        perl -pi -e "s/rm /rm -f /g" scripts/build-rmi.linear_spline.linear.sh
        export CXX=/usr/bin/c++
        export RUST_BACKTRACE=full
        ./minimap2 -I 64G -t {params.threads} -x ava-ont {params.CWD}/{input} {params.CWD}/{input} | \
        awk '$11>=500' | tail -n +2 | /home/git/fpa/target/release/fpa drop --same-name --internalmatch - | pigz -p {params.threads} > {params.CWD}/{output.paf}
        micromamba deactivate
        """


# run VECHAT/RACON once
rule vechat:
        input:
            reads="fasta/{id}.fasta.gz",
            paf="vechat/{id}.paf.gz"
        output: "vechat/{id}.vechat.fasta"
        params: THREADS
        shell: """
        /home/git/vechat/build/bin/racon --no-trimming -u -f -p -d 0.2 -s 0.2 -t {params} \
        -b --cudaaligner-batches 1 -c 1 {input.reads} {input.paf} {input.reads} > {output}
        """


# run brutal rewrite
rule brutal_rewrite:
        input: "vechat/{id}.vechat.fasta"
        output: "br/{id}.vechat.br.fasta"
        params: THREADS
        shell: """
        /home/git/br/target/release/br -t {params} -k 19 -i {input} -m graph -o {output}
        """


# run KMER FILTER
rule kmer_filter:
        input: "br/{id}.vechat.br.fasta"
        output: "kmrf/{id}.vechat.br.kmrf.fasta"
        params: THREADS
        shell: """
        /home/git/kmrf/target/release/kmrf -k 17 -i {input} -o {output}
        """


# add parition step
rule partition:
    input: "kmrf/{id}.vechat.br.kmrf.fasta"
    output: scatter.split("partition/{{id}}.vechat.br.kmrf_{scatteritem}.fasta.gz")
    params:
        split = "8",
        memory = MEMORY_SIZE2
    shell: '''
        eval "$(/home/.local/bin/micromamba shell hook --shell bash)"
        micromamba activate bbmap
        partition.sh {params.memory} in={input} out=partition/{wildcards.id}.vechat.br.kmrf_%.fasta.gz ways={params.split}
        for i in `seq 1 8`
        do
          j=$((i-1))
          mv partition/{wildcards.id}.vechat.br.kmrf_${{j}}.fasta.gz partition/{wildcards.id}.vechat.br.kmrf_${{i}}-of-8.fasta.gz
        done
        micromamba deactivate
        '''


# add shredding step
rule shred:
    input: "partition/{id}.vechat.br.kmrf_{scatteritem}.fasta.gz"
    output: "shred/{id}.vechat.br.kmrf.shred_{scatteritem}.fasta.gz"
    params: MEMORY_SIZE2
    shell: '''
        eval "$(/home/.local/bin/micromamba shell hook --shell bash)"
        micromamba activate bbmap
        shred.sh {params} in={input} out={output} median=90000 variance=2500
        micromamba deactivate
        '''

rule cat:
    input: gather.split("shred/{{id}}.vechat.br.kmrf.shred_{scatteritem}.fasta.gz")
    output: "shred/{id}.vechat.br.kmrf.shred.fasta.gz"
    shell: """cat {input} > {output}"""


# now run peregrine-2021

# make a file listing the location of the raw HiFi reads
rule reads1:
        input: "shred/{id}.vechat.br.kmrf.shred.fasta.gz"
        output: "{id}/pg_asm1/reads1"
        shell: """ls `pwd`/{input} > {output}"""

# make a sequence database for the raw HiFi reads
rule pg_build_sdb1:
        input: "{id}/pg_asm1/reads1"
        output:
            database="{id}/pg_asm1/1.seqdb",
            index="{id}/pg_asm1/1.idx"
        params:
            prefix= "{id}/pg_asm1/1"
        shell: """/home/git/peregrine-2021/target/release/pg_build_sdb --input {input} --out_prefix {params.prefix}"""

# make a shimmer index for the raw reads
rule pg_build_idx1:
        input:
            sqdb="{id}/pg_asm1/1.seqdb",
            idx="{id}/pg_asm1/1.idx"
        output: expand(["{id}/pg_asm1/1-0{threads2}-of-0{threads}.dat"], threads2=THREADS4, threads=THREADS, allow_missing=True)
        params:
            k="56",
            r="6",
            w="80",
            prefix= "{id}/pg_asm1/1",
            chunks= THREADS
        shell: """/home/git/peregrine-2021/target/release/pg_build_idx -k {params.k} -r {params.r} -w {params.w} {input.sqdb} {input.idx} {params.prefix} {params.chunks} {params.chunks}"""

# overlap the raw reads
rule pg_ovlp1:
        input:
            sqdb="{id}/pg_asm1/1.seqdb",
            idx="{id}/pg_asm1/1.idx",
            dat=expand(["{id}/pg_asm1/1-0{threads2}-of-0{threads}.dat"], threads2=THREADS4, threads=THREADS, allow_missing=True)
        output: expand(["{id}/pg_asm1/overlap1.{threads2}"], threads2=THREADS4, allow_missing=True)
        params:
            prefix1= "{id}/pg_asm1/1",
            prefix2= "{id}/pg_asm1/overlap1",
            chunks= THREADS
        shell: """/home/git/peregrine-2021/target/release/pg_ovlp {input.sqdb} {input.idx} {params.prefix1} {params.prefix2} {params.chunks} {params.chunks}"""

# error-correct the raw reads
rule pg_ovlp_ec1:
        input:
            sqdb="{id}/pg_asm1/1.seqdb",
            idx="{id}/pg_asm1/1.idx",
            ovlp=expand(["{id}/pg_asm1/overlap1.{threads2}"], threads2=THREADS4, allow_missing=True)
        output: expand(["{id}/pg_asm2/2_{threads2}.fa"], threads2=THREADS4, allow_missing=True)
        params:
            prefix1= "{id}/pg_asm1/overlap1",
            prefix2= "{id}/pg_asm2/2",
            chunks= THREADS
        shell: """/home/git/peregrine-2021/target/release/pg_ovlp_ec {input.sqdb} {input.idx} {params.prefix1} {params.prefix2} {params.chunks} {params.chunks}"""

# compress corrected reads
rule compress_fasta1:
        input: expand(["{id}/pg_asm2/2_{threads2}.fa"], threads2=THREADS4, allow_missing=True)
        output: expand(["{id}/pg_asm2/2_{threads2}.fa.gz"], threads2=THREADS4, allow_missing=True)
        params: THREADS
        shell: """pigz -p {params} {input}"""


# compress overlap1 files
rule compress_overlap1:
        input:
            reads=expand(["{id}/pg_asm2/2_{threads2}.fa"], threads2=THREADS4, allow_missing=True),
            overlap=expand(["{id}/pg_asm1/overlap1.{threads2}"], threads2=THREADS4, allow_missing=True)
        output: expand(["{id}/pg_asm1/overlap1.{threads2}.gz"], threads2=THREADS4, allow_missing=True)
        params: THREADS
        shell: """pigz -p {params} {input.overlap}"""


# make a file listing the locations of the error-corrected reads
rule reads2:
        input:
            reads=expand(["{id}/pg_asm2/2_{threads2}.fa.gz"], threads2=THREADS4, allow_missing=True),
            overlap=expand(["{id}/pg_asm1/overlap1.{threads2}.gz"], threads2=THREADS4, allow_missing=True)
        output: "{id}/pg_asm2/reads2"
        shell: """ls {input.reads} > {output}"""

# make a sequence database for the 1x corrected HiFi reads
rule pg_build_sdb2:
        input: "{id}/pg_asm2/reads2"
        output:
            database="{id}/pg_asm2/2.seqdb",
            index="{id}/pg_asm2/2.idx"
        params:
            prefix= "{id}/pg_asm2/2"
        shell: """/home/git/peregrine-2021/target/release/pg_build_sdb --input {input} --out_prefix {params.prefix}"""

# make a shimmer index for the 1x corrected reads
rule pg_build_idx2:
        input:
            sqdb="{id}/pg_asm2/2.seqdb",
            idx="{id}/pg_asm2/2.idx"
        output:expand(["{id}/pg_asm2/2-0{threads2}-of-0{threads}.dat"], threads2=THREADS4, threads=THREADS, allow_missing=True)
        params:
            k="56",
            r="6",
            w="80",
            prefix= "{id}/pg_asm2/2",
            chunks= THREADS
        shell: """/home/git/peregrine-2021/target/release/pg_build_idx -k {params.k} -r {params.r} -w {params.w} {input.sqdb} {input.idx} {params.prefix} {params.chunks} {params.chunks}"""

# overlap the 1x corrected reads
rule pg_ovlp2:
        input:
            sqdb="{id}/pg_asm2/2.seqdb",
            idx="{id}/pg_asm2/2.idx",
            dat=expand(["{id}/pg_asm2/2-0{threads2}-of-0{threads}.dat"], threads2=THREADS4, threads=THREADS, allow_missing=True)
        output: expand(["{id}/pg_asm2/overlap2.{threads2}"], threads2=THREADS4, allow_missing=True)
        params:
            prefix1= "{id}/pg_asm2/2",
            prefix2= "{id}/pg_asm2/overlap2",
            chunks= THREADS
        shell: """/home/git/peregrine-2021/target/release/pg_ovlp {input.sqdb} {input.idx} {params.prefix1} {params.prefix2} {params.chunks} {params.chunks}"""

# error-correct the 1x corrected reads
rule pg_ovlp_ec2:
        input:
            sqdb="{id}/pg_asm2/2.seqdb",
            idx="{id}/pg_asm2/2.idx",
            ovlp=expand(["{id}/pg_asm2/overlap2.{threads2}"], threads2=THREADS4, allow_missing=True)
        output: expand(["{id}/pg_asm3/3_{threads2}.fa"], threads2=THREADS4, allow_missing=True)
        params:
            prefix1= "{id}/pg_asm2/overlap2",
            prefix2= "{id}/pg_asm3/3",
            chunks= THREADS
        shell: """/home/git/peregrine-2021/target/release/pg_ovlp_ec {input.sqdb} {input.idx} {params.prefix1} {params.prefix2} {params.chunks} {params.chunks}"""

# compress 2x corrected reads
rule compress_fasta2:
        input: expand(["{id}/pg_asm3/3_{threads2}.fa"], threads2=THREADS4, allow_missing=True)
        output: expand(["{id}/pg_asm3/3_{threads2}.fa.gz"], threads2=THREADS4, allow_missing=True)
        params: THREADS
        shell: """pigz -p {params} {input}"""


# compress overlap2 files
rule compress_overlap2:
        input:
            reads=expand(["{id}/pg_asm3/3_{threads2}.fa.gz"], threads2=THREADS4, allow_missing=True),
            overlap=expand(["{id}/pg_asm2/overlap2.{threads2}"], threads2=THREADS4, allow_missing=True)
        output: expand(["{id}/pg_asm2/overlap2.{threads2}.gz"], threads2=THREADS4, allow_missing=True)
        params: THREADS
        shell: """pigz -p {params} {input.overlap}"""


# make a file listing the locations of the 2x corrected reads
rule reads3:
        input:
            reads=expand(["{id}/pg_asm3/3_{threads2}.fa.gz"], threads2=THREADS4, allow_missing=True),
            overlap=expand(["{id}/pg_asm2/overlap2.{threads2}.gz"], threads2=THREADS4, allow_missing=True)
        output: "{id}/pg_asm3/reads3"
        shell: """ls {input.reads} > {output}"""
        
# make a sequence database for the 2x corrected HiFi reads
rule pg_build_sdb3:
        input: "{id}/pg_asm3/reads3"
        output:
            database="{id}/pg_asm3/3.seqdb",
            index="{id}/pg_asm3/3.idx"
        params:
            prefix= "{id}/pg_asm3/3"
        shell: """/home/git/peregrine-2021/target/release/pg_build_sdb --input {input} --out_prefix {params.prefix}"""

# make a shimmer index for the 2x corrected reads
rule pg_build_idx3:
        input:
            sqdb="{id}/pg_asm3/3.seqdb",
            idx="{id}/pg_asm3/3.idx"
        output: expand(["{id}/pg_asm3/3-0{threads2}-of-0{threads}.dat"], threads2=THREADS4, threads=THREADS, allow_missing=True)
        params:
            k="56",
            r="6",
            w="80",
            prefix= "{id}/pg_asm3/3",
            chunks= THREADS
        shell: """/home/git/peregrine-2021/target/release/pg_build_idx -k {params.k} -r {params.r} -w {params.w} {input.sqdb} {input.idx} {params.prefix} {params.chunks} {params.chunks}"""

# overlap the 2x corrected reads
rule pg_ovlp3:
        input:
            sqdb="{id}/pg_asm3/3.seqdb",
            idx="{id}/pg_asm3/3.idx",
            dat=expand(["{id}/pg_asm3/3-0{threads2}-of-0{threads}.dat"], threads2=THREADS4, threads=THREADS, allow_missing=True)
        output: expand(["{id}/pg_asm3/overlap3.{threads2}"], threads2=THREADS4, allow_missing=True)
        params:
            prefix1= "{id}/pg_asm3/3",
            prefix2= "{id}/pg_asm3/overlap3",
            chunks= THREADS
        shell: """/home/git/peregrine-2021/target/release/pg_ovlp {input.sqdb} {input.idx} {params.prefix1} {params.prefix2} {params.chunks} {params.chunks}"""

# error-correct the 2x corrected reads
rule pg_ovlp_ec3:
        input:
            sqdb="{id}/pg_asm3/3.seqdb",
            idx="{id}/pg_asm3/3.idx",
            ovlp=expand(["{id}/pg_asm3/overlap3.{threads2}"], threads2=THREADS4, allow_missing=True)
        output: expand(["{id}/pg_asm4/4_{threads2}.fa"], threads2=THREADS4, allow_missing=True)
        params:
            prefix1= "{id}/pg_asm3/overlap3",
            prefix2= "{id}/pg_asm4/4",
            chunks= THREADS
        shell: """/home/git/peregrine-2021/target/release/pg_ovlp_ec {input.sqdb} {input.idx} {params.prefix1} {params.prefix2} {params.chunks} {params.chunks}"""

# compress 3x corrected reads
rule compress_fasta3:
        input: expand(["{id}/pg_asm4/4_{threads2}.fa"], threads2=THREADS4, allow_missing=True)
        output: expand(["{id}/pg_asm4/4_{threads2}.fa.gz"], threads2=THREADS4, allow_missing=True)
        params: THREADS
        shell: """pigz -p {params} {input}"""

# compress overlap3 files
rule compress_overlap3:
        input:
            reads=expand(["{id}/pg_asm4/4_{threads2}.fa.gz"], threads2=THREADS4, allow_missing=True),
            overlap=expand(["{id}/pg_asm3/overlap3.{threads2}"], threads2=THREADS4, allow_missing=True)
        output: expand(["{id}/pg_asm3/overlap3.{threads2}.gz"], threads2=THREADS4, allow_missing=True)
        params: THREADS
        shell: """pigz -p {params} {input.overlap}"""


# de novo assemble with hifiasm/raft
rule raft1:
        input: expand(["{id}/pg_asm4/4_{threads2}.fa.gz"], threads2=THREADS4, allow_missing=True)
        output:
            multiext("raft1/errorcorrect_{id}", ".ec.bin", ".ec.fa", ".ovlp.source.bin", ".ovlp.reverse.bin", ".bp.r_utg.gfa", ".bp.r_utg.noseq.gfa", 
            ".bp.r_utg.lowQ.bed", ".bp.p_ctg.gfa", ".bp.p_ctg.lowQ.bed", ".bp.p_ctg.noseq.gfa", ".bp.p_utg.gfa", ".bp.p_utg.lowQ.bed", ".bp.p_utg.noseq.gfa",
            ".bp.hap1.p_ctg.gfa", ".bp.hap1.p_ctg.lowQ.bed", ".bp.hap1.p_ctg.noseq.gfa", ".bp.hap2.p_ctg.gfa", ".bp.hap2.p_ctg.lowQ.bed", ".bp.hap2.p_ctg.noseq.gfa",
            ".log")
        params:
            threads= THREADS,
            log = "raft1/errorcorrect_{id}"
        shell: """
        /home/git/hifiasm/hifiasm -o {params.log} -t {params.threads} --write-ec {input} 2> {params.log}.log
        """

# de novo assemble with hifiasm/raft
rule raft2:
        input: "raft1/errorcorrect_{id}.ec.fa"
        output:
            multiext("raft2/getOverlaps_{id}", ".1.ovlp.paf", ".0.ovlp.paf", ".paf")
        params:
            threads = THREADS,
            log = "raft2/getOverlaps_{id}"
        shell: """
        /home/git/hifiasm/hifiasm -o {params.log} -t {params.threads} --dbg-ovec {input}
        cat {params.log}.1.ovlp.paf {params.log}.0.ovlp.paf > {params.log}.paf
        """


# get reference
rule ref:
        input: "raft2/getOverlaps_{id}.paf"
        output: "raft2/{id}.fa.gz"
        shell: """
        wget 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/{wildcards.id}.fa.gz' \
        -O {output}
        """


# get overlap error rates
rule overlap_error_rates:
        input:
            corr = "raft1/errorcorrect_{id}.ec.fa",
            orig = "fasta/{id}.fasta.gz",
            duplex = "fasta/{id}.duplex.fasta.gz",
            ref = "raft2/{id}.fa.gz"
        output:
            corr = "raft2/{id}.corr.paf.gz",
            orig = "raft2/{id}.orig.paf.gz",
            duplex = "raft2/{id}.duplex.paf.gz"
        params: THREADS
        shell: """
        /home/git/mm2-fast/minimap2 --secondary=no -c -I 64G -t {params} -k19 -w13 {input.ref} {input.corr} | cut -f 1-12 | \
        pigz -p {params} > {output.corr}

        /home/git/mm2-fast/minimap2 --secondary=no -c -I 64G -t {params} -x map-ont {input.ref} {input.orig} | cut -f 1-12 | \
        pigz -p {params} > {output.orig}

        /home/git/mm2-fast/minimap2 --secondary=no -c -I 64G -t {params} -k19 -w13 {input.ref} {input.duplex} | cut -f 1-12 | \
        pigz -p {params} > {output.duplex}
        """


# plot overlap error rates
rule overlap_error_rates2:
        input:
            corr = "raft2/{id}.corr.paf.gz",
            orig = "raft2/{id}.orig.paf.gz",
            duplex = "raft2/{id}.duplex.paf.gz"
        output: "raft2/{id}_original_versus_error-corrected_and_duplex_error_rates.svg"
        shell: '''
        eval "$(/home/.local/bin/micromamba shell hook --shell bash)"
        set +u
        micromamba activate r-base

        echo "library(\\"ggplot2\\")" > raft2/{wildcards.id}.R
        echo "test <- read.table(\\"{input.orig}\\", header=F)" >> raft2/{wildcards.id}.R
        echo "test2 <- read.table(\\"{input.corr}\\", header=F)" >> raft2/{wildcards.id}.R
        echo "test3 <- read.table(\\"{input.duplex}\\", header=F)" >> raft2/{wildcards.id}.R
        echo "" >> raft2/{wildcards.id}.R
        echo "test4 <- data.frame(\\"identities\\" = c(test\\$V10/test\\$V11, test2\\$V10/test2\\$V11, test3\\$V10/test3\\$V11)," >> raft2/{wildcards.id}.R
        echo "                    \\"dataset\\" = c(rep(\\"original\\", length(test\\$V10)), rep(\\"corrected\\", length(test2\\$V10)), rep(\\"duplex\\", length(test3\\$V10))))" >> raft2/{wildcards.id}.R
        echo "" >> raft2/{wildcards.id}.R
        echo "test4\\$phred <- -10*log10((1-test4\\$identities))" >> raft2/{wildcards.id}.R
        echo "test4\\$phred[is.infinite(test4\\$phred)] <- 70" >> raft2/{wildcards.id}.R
        echo "" >> raft2/{wildcards.id}.R
        echo "plot4 <- ggplot(test4, aes(x=phred, fill = dataset, colour = dataset)) + geom_density(alpha = 0.1) + theme_classic() +" >> raft2/{wildcards.id}.R
        echo "        theme(axis.title = element_text(size = 17), axis.text = element_text(size =15)," >> raft2/{wildcards.id}.R
        echo "        legend.position=c(0.7, 0.5)) +" >> raft2/{wildcards.id}.R
        echo "    xlab(\\"Number of residue matches/Alignment block length (Read Alignment Identity) Converted to Phred Score, Q70 = no mismatches\\") +" >> raft2/{wildcards.id}.R
        echo "    ylab(\\"Proportion of Reads\\") +" >> raft2/{wildcards.id}.R
        echo "    scale_x_continuous(" >> raft2/{wildcards.id}.R
        echo "        breaks = c(seq(from=0, to=70, by=2))," >> raft2/{wildcards.id}.R
        echo "        limits=c(0,70))" >> raft2/{wildcards.id}.R
        echo "" >> raft2/{wildcards.id}.R
        echo "ggsave(file=\\"{output}\\",plot4,width=16,height=4)" >> raft2/{wildcards.id}.R

        Rscript --vanilla raft2/{wildcards.id}.R

        micromamba deactivate
        '''


# de novo assemble with hifiasm/raft
rule raft3:
        input:
            reads = "raft1/errorcorrect_{id}.ec.fa",
            overlaps = "raft2/getOverlaps_{id}.paf",
            svg = "raft2/{id}_original_versus_error-corrected_and_duplex_error_rates.svg"
        output:
            multiext("raft3/fragmented_{id}", ".reads.fasta", ".long_repeats.txt", ".coverage.txt", ".long_repeats.bed")
        params:
            log = "raft1/errorcorrect_{id}",
            log2 = "raft3/fragmented_{id}"
        shell: """
        COVERAGE=$(grep 'homozygous' {params.log}.log | tail -1 | awk '{{print $6}}')
        /home/git/RAFT/raft -e $COVERAGE -o {params.log2} {input.reads} {input.overlaps}
        """


# if RAFT fails then comment out rule raft3 above and use this raft3 rule instead

# de novo assemble with hifiasm/raft
#rule raft3:
#        input:
#            reads = "raft1/errorcorrect_{id}.ec.fa",
#            overlaps = "raft2/getOverlaps_{id}.paf",
#            plot = "raft2/{id}_original_versus_error-corrected_overlap_error_rates.svg"
#        output:
#            multiext("raft3/fragmented_{id}", ".reads.fasta")
#        params:
#            log = "raft1/errorcorrect_{id}",
#           log2 = "raft3/fragmented_{id}"
#       shell: """
#       cp {input.reads} {params.log2}.reads.fasta
#        """


# de novo assemble with hifiasm/raft
rule raft4:
        input:
            reads = "raft3/fragmented_{id}.reads.fasta",
            ultralong = "fasta/{id}.fasta.gz"
        output:
            multiext("raft4/finalasm_{id}", ".ec.bin", ".ovlp.source.bin", ".ovlp.reverse.bin", ".bp.r_utg.gfa", ".bp.r_utg.noseq.gfa",
            ".bp.r_utg.lowQ.bed", ".bp.p_ctg.gfa", ".bp.p_ctg.lowQ.bed", ".bp.p_ctg.noseq.gfa", ".bp.p_utg.gfa", ".bp.p_utg.lowQ.bed", ".bp.p_utg.noseq.gfa",
            ".bp.hap1.p_ctg.gfa", ".bp.hap1.p_ctg.lowQ.bed", ".bp.hap1.p_ctg.noseq.gfa", ".bp.hap2.p_ctg.gfa", ".bp.hap2.p_ctg.lowQ.bed", ".bp.hap2.p_ctg.noseq.gfa",
            ".uidx.bin", ".ul.ovlp.bin", ".re.uidx.bin", ".re.ul.msk.bin", ".re.uidx.ucr.bin", ".re.ul.ovlp.bin")
        params:
            threads = THREADS,
            log = "raft4/finalasm_{id}"
        shell: """
        /home/git/hifiasm/hifiasm -o {params.log} -t {params.threads} -r2 {input.reads} --ul {input.ultralong}
        """


rule gfatools:
        input:
            hap1 = "raft4/finalasm_{id}.bp.hap1.p_ctg.gfa",
            hap2 = "raft4/finalasm_{id}.bp.hap2.p_ctg.gfa"
        output:
            hap1 = "raft4/finalasm_{id}.bp.hap1.p_ctg.fa",
            hap2 = "raft4/finalasm_{id}.bp.hap2.p_ctg.fa"
        shell: """
        /home/git/gfatools/gfatools gfa2fa {input.hap1} > {output.hap1}
        /home/git/gfatools/gfatools gfa2fa {input.hap2} > {output.hap2}
        """
```

## script to run in the docker container

`/home/jelber43/sandbox3/SP62_uHMW_25072023_pangenome/test.sh`

```bash
#! /bin/bash
cd /sandbox3/SP62_uHMW_25072023_pangenome
# >>> mamba initialize >>>
# !! Contents within this block are managed by 'mamba init' !!
export MAMBA_EXE='/home/.local/bin/micromamba';
export MAMBA_ROOT_PREFIX='/home/micromamba';
__mamba_setup="$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX" 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__mamba_setup"
else
    alias micromamba="$MAMBA_EXE"  # Fallback on help from mamba activate
fi
unset __mamba_setup
# <<< mamba initialize <<<
. "$HOME/.cargo/env"
micromamba activate snakemake
snakemake -j 1 --snakefile hac-to-corrected-compared-to-duplex2.smk --printshellcmds --latency-wait 60 --cores 48 all > SP62_uHMW_25072023_pangenome.log 2>&1
micromamba deactivate
```

## Example output files of interest

One example output file of interest is an SVG file such as

`1) raft2/{id}_original_versus_error-corrected_and_duplex_error_rates.svg`

where `{id}` is the name of the chr, here chr17.

![plot](https://github.com/jelber2/hac-to-duplex/blob/main/chr17_original_versus_error-corrected_and_duplex_error_rates.svg)

`2) raft1/errorcorrect_{id}.ec.fa`

These are the error-corrected reads used by hifiasm for de novo assembly, where `{id}` is the name of the chr, such as "chr17".

`Note`so far, there are no variant calling models trained or with proper parameters for these data, so, so far, I have only used
these reads for de novo assembly with hifiasm or [hifiasm](https://github.com/chhylp123/hifiasm)/[raft](https://github.com/at-cg/RAFT)

I will try to train a [clair3](https://github.com/HKU-BAL/Clair3) model for these data, but so far, I was not getting very good F1 scores with just HG003 and a single chromosome.

`3) raft4/finalasm_{id}.bp.hap1.p_ctg.gfa`
where `{id}` is the name of the chr, such as "chr17".

`4) raft4/finalasm_{id}.bp.hap2.p_ctg.gfa`
where `{id}` is the name of the chr, such as "chr17".

These two files are of interest because they have long phase blocks and are [dual assemblies](https://lh3.github.io/2021/10/10/introducing-dual-assembly).


## commands to run the docker container

```bash
chmod u+x test.sh
nvidia-docker run --gpus all --rm -v /home/jelber43/sandbox3:/sandbox3/ -v /home/jelber43/bin:/git/ -v /var/run/docker.sock:/var/run/docker.sock hac-to-duplex /sandbox3/SP62_uHMW_25072023_pangenome/test.sh &
```


## "dockerfile"

rough steps for how the docker container was made

Ubuntu 20.04 LTS (Focal Fossa) container

```bash
nvidia-docker run --gpus all -it --rm -v /home/jelber43/sandbox3:/sandbox3/ -v /home/jelber43/bin:/git/ -v /var/run/docker.sock:/var/run/docker.sock ubuntu

# Install micromamba
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)

Micromamba binary folder? [~/.local/bin]  /home/.local/bin
Init shell (bash)? [Y/n] y (type)
Configure conda-forge? [Y/n] y (type y)
Prefix location? [~/micromamba] /home/micromamba
source ~/.bashrc

cd /home
micromamba create -n blue-crab -c bioconda -c conda-forge python=3.11.3 -y

micromamba activate blue-crab
python3 -m pip install --upgrade pip

apt-get install zlib1g-dev -y
apt-get install libzstd-dev -y

pip install blue-crab

mkdir -p git
cd /home/git

apt-get install libhdf5-dev zlib1g-dev -y
VERSION=v1.1.0
wget "https://github.com/hasindu2008/slow5tools/releases/download/$VERSION/slow5tools-$VERSION-release.tar.gz" && tar xvf slow5tools-$VERSION-release.tar.gz && cd slow5tools-$VERSION/
./configure
apt install make -y
apt install g++ -y
make zstd=1

curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source ~/.bashrc
# select defaults

cd /home/git

wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.5.0-linux-x64.tar.gz

tar xzf dorado-0.5.0-linux-x64.tar.gz
micromamba deactivate

micromamba create -n git -c bioconda -c conda-forge git -y
micromamba activate git

cd /home/git

git clone https://github.com/cschin/peregrine-2021
cd peregrine-2021
git checkout 6698eb1
cargo build --release

cd /home/git

git clone --recursive https://github.com/bwa-mem2/mm2-fast.git mm2-fast
cd mm2-fast
git checkout 10bde16
make -j 12
apt install pigz

cd /home/git
git clone https://github.com/jelber2/vechat
cd vechat
mkdir -p build
cd build
apt install g++-8 -y
apt install gcc-8 -y
apt install libz-dev -y
apt install cmake -y
apt install libthrust-dev -y


export CC=/usr/bin/gcc-8
export CXX=/usr/bin/g++-8

cmake -DCMAKE_BUILD_TYPE=Release -Dracon_enable_cuda=ON ..

make -j 34

cd /home/git

git clone https://github.com/at-cg/RAFT
cd RAFT
git checkout fb8dcc0
make


cd /home/git

git clone https://github.com/natir/fpa.git
cd fpa
git checkout v0.5.1
cargo build --release

cd /home/git
git clone https://github.com/chhylp123/hifiasm
cd hifiasm
git checkout 1ac574a
make

cd /home/git

git clone https://github.com/natir/kmrf.git
cd kmrf
git checkout 36cad24
cargo build --release

cd /home/git

git clone https://github.com/natir/br.git
cd br
git checkout ad87f92
rm -r rust-toolchain.toml
cargo build --release

micromamba create -n snakemake -c conda-forge -c bioconda snakemake=7.32.4

cd /home/git
git clone https://github.com/lh3/gfatools
cd gfatools/
git checkout 109f927
make

micromamba create -n bbmap -c bioconda bbmap=39.01
micromamba create -n samtools -c bioconda samtools=1.9

apt install libfreetype6 libfreetype6-dev freetype2-demos freetype2-doc

export C_INCLUDE_PATH=/usr/include/freetype2:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=/usr/include/freetype2:$CPLUS_INCLUDE_PATH

micromamba create -n r-base -c bioconda -c conda-forge r-base=4.0.3
micromamba activate r-base
R
install.packages("ggplot2")
1
install.packages("svglite")
apt install rsync

cd /home/git
git clone https://github.com/lh3/seqtk
cd seqtk
git checkout c9458ba
make
```

## Docker container

When it finishes uploading, you can do:

```bash
docker pull jelber2/hac-to-duplex
```
