# If you are starting with pod5 files, you can do this directly

`Get the HAC FASTQ files, then map them to the pangenome with minigraph`

```bash
~/bin/minigraph --version

0.20-r574-dirty

cd ~/sandbox4

cat /mnt/share/nanopore/gridion-data-offload-destination/WGS_HG002_Monarch_08012024/WGS_HG002_Monarch/20240108_0939_P2S-00581-A_PAU53728_d916b678/fastq_*/*.fastq.gz > WGS_HG002_Monarch_08012024.fastq.gz

~/bin/minigraph/minigraph -t 48 -cx lr ../sandbox3/chm13-90c.r518.gfa.gz WGS_HG002_Monarch_08012024.fastq.gz 2> WGS_HG002_Monarch_08012024_minigraph.log | pigz -p 12 > WGS_HG002_Monarch_08012024.gaf.gz &
```



`Separate the gaf.gz file into separate chromosomes`

```bash
mkdir -p WGS_HG002_Monarch_08012024/gaf
while read i
do
  zcat WGS_HG002_Monarch_08012024.gaf.gz |awk -v chromosome="$i" '$6 ~ chromosome":" || $6 == chromosome' OFS='\t' |pigz -p 12 > WGS_HG002_Monarch_08012024_pangenome/gaf/${i}.gaf.gz
done < <(cat <(seq 1 22) <(echo -e "M\nX\nY")|perl -pe "s/^/chr/g"|grep "chr22")
```


`Merge the pod5 files for easier filtering``

Note that `ubuntu.22.04.cuda` is actually `hac-to-duplex2` from 

`docker pull jelber2/hac-to-duplex2``

```bash
nvidia-docker run -it --gpus all --rm -v /home/jelber43/sandbox4:/sandbox4/ -v /mnt/share/nanopore:/nanopore/ -v /var/run/docker.sock:/var/run/docker.sock ubuntu.22.04.cuda

cd /sandbox4

micromamba activate blue-crab

pod5 merge --threads 48 -o 2024-01-09-09-09-08.merged.pod5 --duplicate-ok /nanopore/gridion-data-offload-destination/WGS_HG002_Monarch_08012024/WGS_HG002_Monarch/20240108_0939_P2S-00581-A_PAU53728_d916b678/pod5_*/*.pod5  > pod5merge.log 2>&1 &
```


This snakefile will convert the pod5 to SUP FASTQ and do vechat twice but with the -u option on
the second vechat/racon run

`/home/jelber43/sandbox4/WGS_HG002_Monarch_08012024/hac-to-corrected-with-vechat-twice-with-ultra-long-vechat2-with-u.smk`

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
BLOW5="../2024-01-09-09-09-08.merged.pod5"
CWD="/sandbox4/WGS_HG002_Monarch_08012024_pangenome"

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
rule pod5:
        input:
            reads="fast5/{id}_read_names.txt",
            pod5=BLOW5
        output: "fast5/{id}.pod5"
        params: THREADS
        shell: """
        eval "$(/home/.local/bin/micromamba shell hook --shell bash)"
        micromamba activate blue-crab
        pod5 filter {input.blow5} --output {output} --missing-ok --ids {input.reads}
        """


# basecall with SUP with dorado 0.5.0 and SUP
rule dorado:
        input: "fast5/{id}.pod5"
        output: "fasta/{id}.bam"
        shell: """/home/git/dorado-0.5.0-linux-x64/bin/dorado basecaller sup {input} > {output}"""


# convert BAM to FASTQ
rule bam2fastq:
        input:
            orig = "fasta/{id}.bam"
        output:
            orig =  "fasta/{id}.fasta.gz"
        params: THREADS
        shell: '''
        eval "$(/home/.local/bin/micromamba shell hook --shell bash)"
        micromamba activate samtools
        samtools fasta -@ {params} {input.orig} | pigz -p {params} > {output.orig}
        micromamba deactivate
        '''


# run mm2-fast to overlap reads and filter overlaps with fpa
rule mm2_fast:
        input: "fasta/{id}.fasta.gz"
        output:
            paf = "vechat/{id}.paf.gz",
            mm2_fast = multiext("fasta/{id}.fasta.gz_ava-ont_minimizers_key_value_sorted", "_keys.rmi_PARAMETERS",
                                "_pos_bin", "_val_bin", "_size", "_keys.uint64")
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
        ./minimap2 -I 64G -t {params.threads} --dual=yes -x ava-ont {params.CWD}/{input} {params.CWD}/{input} | \
        awk '$11>=500' | tail -n +2 | /home/git/fpa/target/release/fpa drop --same-name --internalmatch - | pigz -p {params.threads} > {params.CWD}/{output.paf}
        rm -fr {wildcards.id}_2
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
        /home/git/vechat/build/bin/racon -f -p -d 0.2 -s 0.2 -t {params} \
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


# remove duplicates
rule dedupe:
        input: "kmrf/{id}.vechat.br.kmrf.fasta"
        output: "kmrf/{id}.vechat.br.kmrf2.fasta"
        params:
            threads = THREADS,
            memory = MEMORY_SIZE1
        shell: """
        /home/git/seqtk/seqtk seq -Cl0 {input} | paste - - | sort -uk1,1 --parallel={params.threads} -S {params.memory} | tr '\t' '\n' > {output}
        """


# run mm2-fast to overlap reads and filter overlaps with fpa
rule mm2_fast2:
        input: "kmrf/{id}.vechat.br.kmrf2.fasta"
        output:
            paf = "vechat2/{id}.vechat.br.kmrf.paf.gz",
            mm2_fast = multiext("kmrf/{id}.vechat.br.kmrf2.fasta_ava-ont_minimizers_key_value_sorted", "_keys.rmi_PARAMETERS",
                                "_pos_bin", "_val_bin", "_size", "_keys.uint64")
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
        ./minimap2 -I 64G -t {params.threads} --dual=yes -x ava-ont {params.CWD}/{input} {params.CWD}/{input} | \
        awk '$11>=500 && $10/$11>=0.99' | tail -n +2 | /home/git/fpa/target/release/fpa drop --same-name --internalmatch - | pigz -p {params.threads} > {params.CWD}/{output.paf}
        rm -fr {wildcards.id}_2
        micromamba deactivate
        """


# run VECHAT/RACON second time
rule vechat2:
        input:
            reads="kmrf/{id}.vechat.br.kmrf2.fasta",
            paf="vechat2/{id}.vechat.br.kmrf.paf.gz"
        output: "vechat2/{id}.vechat.br.kmrf.vechat.fasta"
        params: THREADS
        shell: """
        /home/git/vechat/build/bin/racon -u -f -p -d 0.2 -s 0.2 -t {params} \
        -b --cudaaligner-batches 1 -c 1 {input.reads} {input.paf} {input.reads} | /home/git/seqtk/seqtk seq -Cl0 > {output}
        """


# add parition step
rule partition:
    input: "vechat2/{id}.vechat.br.kmrf.vechat.fasta"
    output: scatter.split("partition/{{id}}.vechat.br.kmrf.vechat_{scatteritem}.fasta.gz")
    params:
        split = "8",
        memory = MEMORY_SIZE2
    shell: '''
        eval "$(/home/.local/bin/micromamba shell hook --shell bash)"
        micromamba activate bbmap
        partition.sh {params.memory} in={input} out=partition/{wildcards.id}.vechat.br.kmrf.vechat_%.fasta.gz ways={params.split}
        for i in `seq 1 8`
        do
          j=$((i-1))
          mv partition/{wildcards.id}.vechat.br.kmrf.vechat_${{j}}.fasta.gz partition/{wildcards.id}.vechat.br.kmrf.vechat_${{i}}-of-8.fasta.gz
        done
        micromamba deactivate
        '''


# add shredding step
rule shred:
    input: "partition/{id}.vechat.br.kmrf.vechat_{scatteritem}.fasta.gz"
    output: "shred/{id}.vechat.br.kmrf.vechat.shred_{scatteritem}.fasta.gz"
    params: MEMORY_SIZE2
    shell: '''
        eval "$(/home/.local/bin/micromamba shell hook --shell bash)"
        micromamba activate bbmap
        shred.sh {params} in={input} out={output} median=90000 variance=2500
        micromamba deactivate
        '''


rule cat:
    input: gather.split("shred/{{id}}.vechat.br.kmrf.vechat.shred_{scatteritem}.fasta.gz")
    output: "shred/{id}.vechat.br.kmrf.vechat.shred.fasta.gz"
    shell: """cat {input} > {output}"""


# now run peregrine-2021

# make a file listing the location of the raw HiFi reads
rule reads1:
        input: "shred/{id}.vechat.br.kmrf.vechat.shred.fasta.gz"
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
            multiext("raft1/errorcorrect_{id}", ".stitch-fasta.jl", ".raw.fasta", ".stitched.fasta", ".ec.bin", ".ec.fa", ".ovlp.source.bin", ".ovlp.reverse.bin", ".bp.r_utg.gfa", ".bp.r_utg.noseq.gfa", 
            ".bp.r_utg.lowQ.bed", ".bp.p_ctg.gfa", ".bp.p_ctg.lowQ.bed", ".bp.p_ctg.noseq.gfa", ".bp.p_utg.gfa", ".bp.p_utg.lowQ.bed", ".bp.p_utg.noseq.gfa",
            ".bp.hap1.p_ctg.gfa", ".bp.hap1.p_ctg.lowQ.bed", ".bp.hap1.p_ctg.noseq.gfa", ".bp.hap2.p_ctg.gfa", ".bp.hap2.p_ctg.lowQ.bed", ".bp.hap2.p_ctg.noseq.gfa",
            ".log")
        params:
            threads= THREADS,
            log = "raft1/errorcorrect_{id}",
            memory = MEMORY_SIZE1
        shell: """
        zcat {input} | /home/git/seqtk/seqtk seq -Cl0 | paste - - |perl -pe s"/rr_/_/g" | sort --parallel={params.threads} -S {params.memory} -g -k1,1 |tr '\t' '\n' > {params.log}.raw.fasta

        echo "import Pkg; Pkg.add(\\"FASTX\\")" > {params.log}.stitch-fasta.jl
        echo "import Pkg; Pkg.add(\\"BioSequences\\")" >> {params.log}.stitch-fasta.jl
        echo "using BioSequences" >> {params.log}.stitch-fasta.jl
        echo "using FASTX" >> {params.log}.stitch-fasta.jl
        echo "" >> {params.log}.stitch-fasta.jl
        echo "function process_fasta_file(filename::String)" >> {params.log}.stitch-fasta.jl
        echo "    reader = FASTAReader(open(filename))" >> {params.log}.stitch-fasta.jl
        echo "    sequences = Dict{{String, String}}()" >> {params.log}.stitch-fasta.jl
        echo "" >> {params.log}.stitch-fasta.jl
        echo "    for record in reader" >> {params.log}.stitch-fasta.jl
        echo "        # Extract the main identifier before the underscore" >> {params.log}.stitch-fasta.jl
        echo "        main_id = split(identifier(record), '_')[1]" >> {params.log}.stitch-fasta.jl
        echo "        seq = sequence(record)" >> {params.log}.stitch-fasta.jl
        echo "        " >> {params.log}.stitch-fasta.jl
        echo "        # Aggregate sequences by the main identifier" >> {params.log}.stitch-fasta.jl
        echo "        if haskey(sequences, main_id)" >> {params.log}.stitch-fasta.jl
        echo "            sequences[main_id] *= seq" >> {params.log}.stitch-fasta.jl
        echo "        else" >> {params.log}.stitch-fasta.jl
        echo "            sequences[main_id] = seq" >> {params.log}.stitch-fasta.jl
        echo "        end" >> {params.log}.stitch-fasta.jl
        echo "    end" >> {params.log}.stitch-fasta.jl
        echo "" >> {params.log}.stitch-fasta.jl
        echo "    close(reader)" >> {params.log}.stitch-fasta.jl
        echo "" >> {params.log}.stitch-fasta.jl
        echo "    # Print the aggregated sequences" >> {params.log}.stitch-fasta.jl
        echo "    for (id, seq) in sequences" >> {params.log}.stitch-fasta.jl
        echo "        println(\\">\\$id\\")" >> {params.log}.stitch-fasta.jl
        echo "        println(seq)" >> {params.log}.stitch-fasta.jl
        echo "    end" >> {params.log}.stitch-fasta.jl
        echo "end" >> {params.log}.stitch-fasta.jl
        echo "" >> {params.log}.stitch-fasta.jl
        echo "# Check if a filename is passed as a command line argument" >> {params.log}.stitch-fasta.jl
        echo "if length(ARGS) < 1" >> {params.log}.stitch-fasta.jl
        echo "    println(\\"Usage: julia \\", PROGRAM_FILE, \\" <filename.fasta>\\")" >> {params.log}.stitch-fasta.jl
        echo "    return" >> {params.log}.stitch-fasta.jl
        echo "end" >> {params.log}.stitch-fasta.jl
        echo "# Example usage with command line argument" >> {params.log}.stitch-fasta.jl
        echo "process_fasta_file(ARGS[1])" >> {params.log}.stitch-fasta.jl

        /root/.juliaup/bin/julia {params.log}.stitch-fasta.jl {params.log}.raw.fasta > {params.log}.stitched.fasta
        /home/git/hifiasm/hifiasm -o {params.log} -t {params.threads} --write-ec {params.log}.stitched.fasta 2> {params.log}.log
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
            ref = "raft2/{id}.fa.gz"
        output:
            corr = "raft2/{id}.corr.paf.gz",
            orig = "raft2/{id}.orig.paf.gz",
        params: THREADS
        shell: """
        /home/git/mm2-fast/minimap2 --secondary=no -c -I 64G -t {params} -k19 -w13 {input.ref} {input.corr} | cut -f 1-12 | \
        pigz -p {params} > {output.corr}

        /home/git/mm2-fast/minimap2 --secondary=no -c -I 64G -t {params} -x map-ont {input.ref} {input.orig} | cut -f 1-12 | \
        pigz -p {params} > {output.orig}
        """


# plot overlap error rates
rule overlap_error_rates2:
        input:
            corr = "raft2/{id}.corr.paf.gz",
            orig = "raft2/{id}.orig.paf.gz",
        output: "raft2/{id}_original_versus_error-corrected_error_rates.svg"
        shell: '''
        eval "$(/home/.local/bin/micromamba shell hook --shell bash)"
        set +u
        micromamba activate r-base

        echo "library(\\"ggplot2\\")" > raft2/{wildcards.id}.R
        echo "test <- read.table(\\"{input.orig}\\", header=F)" >> raft2/{wildcards.id}.R
        echo "test2 <- read.table(\\"{input.corr}\\", header=F)" >> raft2/{wildcards.id}.R
        echo "" >> raft2/{wildcards.id}.R
        echo "test4 <- data.frame(\\"identities\\" = c(test\\$V10/test\\$V11, test2\\$V10/test2\\$V11)," >> raft2/{wildcards.id}.R
        echo "                    \\"dataset\\" = c(rep(\\"original\\", length(test\\$V10)), rep(\\"corrected\\", length(test2\\$V10))))" >> raft2/{wildcards.id}.R
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
            svg = "raft2/{id}_original_versus_error-corrected_error_rates.svg"
        output:
            multiext("raft3/fragmented_{id}", ".reads.fasta", ".long_repeats.txt", ".coverage.txt", ".long_repeats.bed")
        params:
            log = "raft1/errorcorrect_{id}",
            log2 = "raft3/fragmented_{id}"
        shell: """
        COVERAGE=$(grep 'homozygous' {params.log}.log | tail -1 | awk '{{print $6}}')
        /home/git/RAFT/raft -e $COVERAGE -o {params.log2} {input.reads} {input.overlaps}
        """


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

This script will enable the container and Snakefile to be run non-interactively

`/home/jelber43/sandbox4/WGS_HG002_Monarch_08012024/vechat-twice.sh`

```bash
#! /bin/bash
cd /sandbox4/WGS_HG002_Monarch_08012024
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
snakemake -j 1 --snakefile hac-to-corrected-with-vechat-twice-with-ultra-long-vechat2-with-u.smk --printshellcmds --cores 48 all > WGS_HG002_Monarch_08012024.log 2>&1
micromamba deactivate
```



Run the commands

```bash
nvidia-docker run --gpus all --rm -v /home/jelber43/sandbox4:/sandbox4/ -v /var/run/docker.sock:/var/run/docker.sock ubuntu.22.04.cuda /sandbox4/WGS_HG002_Monarch_08012024/vechat-twice.sh &
```
