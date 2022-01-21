## snakemake --jobs 120 --cluster "sbatch -c {threads} --time 30" --use-envmodules

configfile: "snake_config.yaml"
localrules: all, update_alignments_table, extract_relevant_alignment_info, stage_in, multiqc_pre, sample_names


import os
import re
#import subprocess as sp


# I'll have to type these many times, so simplify them.
# Snakemake recommends not using absolute paths, for portability, but it is causing problems for cluster jobs.
HOME = config["projectHome"]
SCRATCH = config["projectScratch"]

sys.path.append(HOME + '/code')


## Parse sample IDs and names from description file
samples_raw = None
samples_temp = {}
with open(config["sampleinfo"], 'r') as fin:
	next(fin) # skip header line
	# Get sample identifiers as tuples.
	samples_raw = [ (x[config["idcols"]["id"]].rstrip("\n"), \
	                 x[config["idcols"]["name"]].rstrip("\n")) \
	               for x in [line.split("\t") for line in fin] ]
# Convert tuples to dictionary.
def Convert(tup, di):
    di = dict(tup)
    return di
samples = Convert(samples_raw, samples_temp)
#print( samples )

THREADS=len(samples.keys())



# Create file names based on assigned name and original sample ID.
uscore = config["idjoin"]
COMBONAME=[uscore.join([samples[x], x]) for x in samples.keys()]





## One chunk to rule them all
rule all:
	input:
		HOME + "/process/multiqc_pre/multiqc_report.html",
		expand(HOME + "/process/featureCounts/{countsfile}", 
		       countsfile=["exon_se_genecounts.txt", "intron_se_genecounts.txt", "exonsplit_se_genecounts.txt", 
			               "exon_se_counts.txt", "exon_se_counts.txt.jcounts", "intron_se_counts.txt"])
##



# Obtain data from VBCF


rule update_alignments_table:
	output:
		HOME + "/aux/all_alignments_table.txt"
	params:
		url=config["vbcf_alignments_table_url"]
	shell:
		"/usr/bin/wget -c --no-check-certificate --auth-no-challenge --user $(cat $HOME/vbcf_username.txt) --password $(cat $HOME/vbcf_password.txt) -O {output} {params.url}"


rule extract_relevant_alignment_info:
	input:
		rules.update_alignments_table.output
	output:
		HOME + "/aux/sample_vbcf_metadata.txt"
	run:
		r = re.compile('|'.join(samples.keys()))
		with open(str(input), 'r') as fin:
			with open(str(output), 'w') as fout:
				for line in fin:
					l = line.split("\t")
					matched = r.match(l[11])  # 11 sample id, 50 paired, 52 bam path, 54 bam url
					if matched or l[0] == "flowcellId":
						fout.write(line)


rule fetch_alignments:
	input:
		HOME + "/aux/sample_vbcf_metadata.txt"
	output:
		HOME + "/data/{comboname}.aligned.bam"
	params:
		sep=uscore
	shell:
		"/usr/bin/wget -c --no-check-certificate --auth-no-challenge --user $(cat $HOME/vbcf_username.txt) --password $(cat $HOME/vbcf_password.txt) -O {output} $(grep $(perl -e '$a=\"{output}\";$a=~s/.*?{params.sep}(\d+)\.aligned\.bam/$1/;print $a') {input} | cut -f 55)"
##


## Prepare


rule stage_in:
	input:
		HOME + "/data/{comboname}.aligned.bam"
	output:
		SCRATCH + "/data/{comboname}.aligned.bam"
	params:
		project=SCRATCH
	shell:
		"cp {input} {params.project}/data/"


rule fix_chr:
	input:
		SCRATCH + "/data/{comboname}.aligned.bam"
	output:
		SCRATCH + "/data_chrfixed/{comboname}.aligned.bam"
	shell:
		"samtools view -h {input} | sed 's/chrMT/chrM/g' | samtools view -Sb - > {output}"


rule index_bam:
	input:
		SCRATCH + "/data_chrfixed/{comboname}.aligned.bam"
	output:
		SCRATCH + "/data_chrfixed/{comboname}.aligned.bam.bai"
	envmodules:
		"build-env/f2021",
		"samtools/1.12-gcc-10.2.0"
	shell:
		"samtools index {input}"
##

# Already markdup'ed by facility



## Quality control


rule flagstats:
	input:
		SCRATCH + "/data_chrfixed/{comboname}.aligned.bam"
	output:
		SCRATCH + "/qc_pre/{comboname}.flagstats.txt"
	envmodules:
		"build-env/f2021",
		"samtools/1.12-gcc-10.2.0"
	threads: 1 # fastqc is slower than samtools, no point finishing too fast
	shell:
		"samtools flagstat -@ $(({threads} - 1)) {input} > {output}"


rule bamstats:
	input:
		SCRATCH + "/data_chrfixed/{comboname}.aligned.bam"
	output:
		SCRATCH + "/qc_pre/{comboname}.bamstats.txt"
	params:
		ref=config["reference"]
	envmodules:
		"build-env/f2021",
		"samtools/1.12-gcc-10.2.0"
	threads: 1 # fastqc is slower than samtools, no point finishing too fast
	shell:     # The facility's reference has contig names that are not in our reference fasta and it breaks the -r option.
		#"samtools stats -@ $(({threads} - 1)) -p -r {params.ref} {input} > {output}"
		"samtools stats -@ $(({threads} - 1)) -p {input} > {output}"


rule fastqc_mapped:
	input:
		SCRATCH + "/data_chrfixed/{comboname}.aligned.bam"
	output:
		SCRATCH + "/qc_pre/{comboname}.aligned_fastqc.html",
		SCRATCH + "/qc_pre/{comboname}.aligned_fastqc.zip"
	params:
		project=SCRATCH
	threads: 1 # multiple threads not applicable when each input has its own instance
	envmodules:
		"build-env/f2021",
		"fastqc/0.11.9-java-11"
	shell:
		"fastqc -f bam_mapped -t {threads} -o {params.project}/qc_pre {input}"


rule multiqc_pre:
	input:						# don't os.path.dirname here, because it won't check for folder contents
		flagstats=expand(SCRATCH + "/qc_pre/{comboname}.flagstats.txt", comboname=COMBONAME),
		bamstats=expand(SCRATCH + "/qc_pre/{comboname}.bamstats.txt", comboname=COMBONAME),
		fastqc=expand(SCRATCH + "/qc_pre/{comboname}.aligned_fastqc.zip", comboname=COMBONAME)
	output:
		HOME + "/process/multiqc_pre/multiqc_report.html"
	params:
		project=HOME
	envmodules:
		"build-env/f2021",
		"multiqc/1.11-foss-2020b-python-3.8.6"
	shell:
		"multiqc -f -o {params.project}/process/multiqc_pre $(dirname {input.fastqc[0]}) $(dirname {input.flagstats[0]}) $(dirname {input.bamstats[0]})"
##



## Quantify


# --nonOverlap 0 should prevent intron-exon boundary-crossing reads.
# Intron counts are a proxy for pre-mRNA counts. Therefore we count reads spanning exon-intron boundaries as intronic instead of exonic.

# featureCounts and htseq don't allow quantifying paired end data as single end.


rule unpair:
	input:
		bam = SCRATCH + "/data_chrfixed/{comboname}.aligned.bam",
		bai = SCRATCH + "/data_chrfixed/{comboname}.aligned.bam.bai"
	output:
		SCRATCH + "/data_unpaired/{comboname}.aligned.bam"
	params:
		home=HOME,
		scratch=SCRATCH
	shell:
		"{params.home}/code/sequtilities.py T {input.bam} -O {params.scratch}/data_unpaired '' '.bam' --samUnpair"


rule quantifygenes_exons_singleend:
	input:
		expand(SCRATCH + "/data_unpaired/{comboname}.aligned.bam", comboname=COMBONAME)
	output:
		SCRATCH + "/featureCounts/exon_se_genecounts.txt"
	params:
		ann=config["annotation"]["exons"]
	threads: THREADS
	group: "featureCounts"
	envmodules:
		"build-env/f2021",
		"subread/2.0.2-gcc-10.2.0"
	shell:
		"featureCounts -T {threads} -t exon -O --fraction --nonOverlap 0 -a {params.ann} -o {output} {input}"


rule quantifygenes_introns_singleend:
	input:
		expand(SCRATCH + "/data_unpaired/{comboname}.aligned.bam", comboname=COMBONAME)
	output:
		SCRATCH + "/featureCounts/intron_se_genecounts.txt"
	params:
		ann=config["annotation"]["introns"]
	threads: THREADS
	group: "featureCounts"
	envmodules:
		"build-env/f2021",
		"subread/2.0.2-gcc-10.2.0"
	shell:
		"featureCounts -T {threads} -t intron -O --fraction -a {params.ann} -o {output} {input}"


rule quantifygenes_split_singleend:   # Requiring exons to not have overlaps with other features at all (with introns in mind specifically) might also throw out spliced reads spanning different exons.
	input:
		expand(SCRATCH + "/data_unpaired/{comboname}.aligned.bam", comboname=COMBONAME)
	output:
		SCRATCH + "/featureCounts/exonsplit_se_genecounts.txt"
	params:
		ann=config["annotation"]["introns"]
	threads: THREADS
	group: "featureCounts"
	envmodules:
		"build-env/f2021",
		"subread/2.0.2-gcc-10.2.0"
	shell:
		"featureCounts -T {threads} -t exon --splitOnly -O --fraction -a {params.ann} -o {output} {input}"


rule quantifyexons_singleend:   # and exon junctions
	input:
		expand(SCRATCH + "/data_unpaired/{comboname}.aligned.bam", comboname=COMBONAME)
	output:
		exons=SCRATCH + "/featureCounts/exon_se_counts.txt",
		junctions=SCRATCH + "/featureCounts/exon_se_counts.txt.jcounts"
	params:
		ann=config["annotation"]["exons"]
	threads: THREADS
	group: "featureCounts"
	envmodules:
		"build-env/f2021",
		"subread/2.0.2-gcc-10.2.0"
	shell:
		"featureCounts -T {threads} -t exon -f -J -O --fraction --nonOverlap 0 -a {params.ann} -o {output.exons} {input}"


rule quantifyintrons_singleend:
	input:
		expand(SCRATCH + "/data_unpaired/{comboname}.aligned.bam", comboname=COMBONAME)
	output:
		SCRATCH + "/featureCounts/intron_se_counts.txt"
	params:
		ann=config["annotation"]["introns"]
	threads: THREADS
	group: "featureCounts"
	envmodules:
		"build-env/f2021",
		"subread/2.0.2-gcc-10.2.0"
	shell:
		"featureCounts -T {threads} -t intron -f -O --fraction -a {params.ann} -o {output} {input}"


rule sample_names: # featureCounts uses full input paths as sample names, wich sucks for everything downstream.
	input:
		SCRATCH + "/featureCounts/{countsfile}"
		# SCRATCH + "/featureCounts/exon_se_genecounts.txt",
		# SCRATCH + "/featureCounts/intron_se_genecounts.txt",
		# SCRATCH + "/featureCounts/exonsplit_se_genecounts.txt",
		# SCRATCH + "/featureCounts/exon_se_counts.txt",
		# SCRATCH + "/featureCounts/exon_se_counts.txt.jcounts",
		# SCRATCH + "/featureCounts/intron_se_counts.txt"
	output:
		HOME + "/process/featureCounts/{countsfile}"
	params:
		prefix=(SCRATCH + "/data_unpaired/").replace('/','\/')
	shell:  # Skip the top line comment and strip the paths and file extensions from the sample names.
		"""
		head -n2 {input} | tail -n1 | perl -ne 's/.aligned.bam//g;s/{params.prefix}//g;print' > {output} &&
		tail -n +3 {input} >> {output}
		"""
##



## Tidy up


# rule cleanup_home:
# 	input:
# 		HOME + "/process/featureCounts/exon_se_genecounts.txt", # prevent cleanup until project completed
# 		HOME + "/process/featureCounts/intron_se_genecounts.txt",
		# HOME + "/process/featureCounts/exonsplit_se_genecounts.txt",
# 		HOME + "/process/featureCounts/exon_se_counts.txt",
# 		HOME + "/process/featureCounts/exon_se_counts.txt.jcounts",
# 		HOME + "/process/featureCounts/intron_se_counts.txt",
# 		HOME + "/aux/all_alignments_table.txt"
# 	params:
# 		project=HOME
# 	shell:
# 		"rm {params.project}/aux/all_alignments_table.txt"


# rule cleanup_scratch:
# 	input:
# 		HOME + "/process/featureCounts/exon_se_genecounts.txt", # prevent cleanup until project completed
# 		HOME + "/process/featureCounts/intron_se_genecounts.txt",
		# HOME + "/process/featureCounts/exonsplit_se_genecounts.txt",
		# HOME + "/process/featureCounts/exon_se_counts.txt",
		# HOME + "/process/featureCounts/exon_se_counts.txt.jcounts",
		# HOME + "/process/featureCounts/intron_se_counts.txt"
# 	params:
# 		project=SCRATCH
# 	shell:
# 		"rm {params.project}"
##
