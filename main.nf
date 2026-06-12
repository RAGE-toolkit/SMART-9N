// main.nf
nextflow.enable.dsl=2

// -----------------------------------------------------------------------------
// Project-relative includes (robust under EPI2ME instances)
// -----------------------------------------------------------------------------
def MODULES = "${projectDir}/modules"

include { GUPPY_BASECALLER }       from "${MODULES}/guppy_basecaller.nf"
include { GUPPY_BARCODER }         from "${MODULES}/guppy_barcode.nf"
include { GUPPY_PLEX }             from "${MODULES}/guppy_plex.nf"

include { DORADO_BASECALLER }      from "${MODULES}/dorado_basecaller.nf"
include { DORADO_BARCODER }        from "${MODULES}/dorado_barcoder.nf"

include { PLEX_FQ_FILES }          from "${MODULES}/plex_fq_files.nf"
include { PLEX_DIRS }              from "${MODULES}/plex_dirs.nf"
include { NANOSTAT }               from "${MODULES}/nanostat.nf"

include { KRAKEN }                 from "${MODULES}/kraken.nf"
include { MINIMAP2 }               from "${MODULES}/minimap2.nf"
include { DOWNLOAD_REFERENCE }     from "${MODULES}/download_reference.nf"
include { COVERAGE_SUMMARY }       from "${MODULES}/coverage_summary.nf"
include { MEDAKA_SNP_1 }           from "${MODULES}/medaka_snp-1.nf"
include { MEDAKA_CONSENSUS }       from "${MODULES}/medaka_consensus.nf"
include { MASK }                   from "${MODULES}/mask.nf" # use custom mask depth below 20x coverage are discarded with 'N' in the consensus
include { SUMMARY_STATS }          from "${MODULES}/summary_stats.nf" #add N50, coverage, etc. to the final report


// -----------------------------------------------------------------------------
// Param defaults (CLI/GUI can override). No filesystem ops at compile time.
// -----------------------------------------------------------------------------
params.out_dir       = params.output_dir    ?: "${projectDir}/results"
params.fastq_dir     = params.fastq_dir     ?: 'raw_files/fastq'
params.rawfile_dir   = params.rawfile_dir   ?: 'raw_files'          
params.rawfile_type  = params.rawfile_type  ?: 'fastq'              // 'fastq' | 'fast5_pod5'
params.basecaller    = params.basecaller    ?: 'Dorado'             // 'Dorado' | 'Guppy'
params.fq_extension  = params.fq_extension  ?: '.fastq'
params.threads       = (params.threads ?: 5) as int

// Common model/config params should already be set in nextflow.config (as you have)

// -----------------------------------------------------------------------------
// Sample sheet: require GUI/CLI to provide --sample_sheet (CSV with header)
// Expect columns: sampleId, barcode, schema, version
// -----------------------------------------------------------------------------
if( !params.sample_sheet ) {
  exit 1, "No sample sheet provided. Set --sample_sheet via EPI2ME GUI or CLI."
}
def sampleSheet = file(params.sample_sheet)
if( !sampleSheet.exists() ) {
  exit 1, "Sample sheet not found: ${sampleSheet}"
}

// Build CSV channel (light validation)
Channel
  .fromPath(sampleSheet, checkIfExists: true)
  .splitCsv(header: true, sep: ',')
  .map { row ->
      def sid     = (row.sampleId ?: '').toString().trim()
      def barcode = (row.barcode  ?: '').toString().trim()
      def scheme  = (row.schema   ?: '').toString().trim()
      def version = (row.version  ?: '').toString().trim()
      if( !sid || !barcode || !scheme || !version )
          throw new IllegalArgumentException("Sample sheet row missing fields: ${row}")
      tuple(sid, barcode, scheme, version)
  }
  .set { fq_channel }

// -----------------------------------------------------------------------------
// Logs (no mkdirs/prefetching here; processes handle their own outputs)
// -----------------------------------------------------------------------------
log.info ""
log.info "=== RAGE-toolkit/Artic-nf ==="
log.info "Output dir       : ${params.out_dir}"
log.info "Sample sheet     : ${sampleSheet}"
log.info "Run name         : ${params.run_name}"
log.info "Rawfile type     : ${params.rawfile_type}"
log.info "Basecaller       : ${params.basecaller}"
log.info "FASTQ dir        : ${params.fastq_dir}"
log.info "Rawfile dir      : ${params.rawfile_dir}"
log.info "Threads          : ${params.threads}"
log.info "QueueSize        : ${params.queueSize}"
log.info "Mask size        : ${params.mask_depth}"
log.info "Medaka normalise : ${params.medaka_normalise}"	
log.info "Medaka model     : ${params.medaka_model}"
log.info "Sequence length  : ${params.seq_len}" 

if( params.rawfile_type == 'fast5_pod5' ) {
    log.info "---- Basecaller details ----"
    log.info "Basecaller config     : ${params.basecaller_config}"
    log.info "GPU mode              : ${params.gpu_mode}"
    if( params.basecaller_dir ) {
        log.info "Basecaller dir    : ${params.basecaller_dir}"
    }
    if( params.model_dir ) {
        log.info "Model dir         : ${params.model_dir}"
    }
}

log.info "=============================="

def dir_plex_script          = file("${projectDir}/scripts/directory_plex.py", checkIfExists: true)
def align_trim_script	     = Channel.fromPath("${projectDir}/scripts/align_trim.py")
def vcf_merge_script	     = Channel.fromPath("${projectDir}/scripts/vcf_merge.py")
def vcf_filter_script        = Channel.fromPath("${projectDir}/scripts/vcf_filter.py")
def make_depth_mask_script   = Channel.fromPath("${projectDir}/scripts/make_depth_mask.py")
def mask_script              = Channel.fromPath("${projectDir}/scripts/mask.py")
def fasta_header_script      = Channel.fromPath("${projectDir}/scripts/fasta_header.py")
def summary_stats_script     = Channel.fromPath("${projectDir}/scripts/summary_stats.py")
def report_script            = Channel.fromPath("${projectDir}/scripts/report.py")
def concat_script            = Channel.fromPath("${projectDir}/scripts/concat.py")


def medaka_dir 				= Channel.fromPath("${params.out_dir}/medaka")
def summary_stats_dir 		= Channel.fromPath("${params.out_dir}/summary_stats")
def rawfile_dir 			= file("${params.rawfile_dir}")
def basecaller_dir 			= Channel.fromPath("${params.basecaller_dir}")
def basecaller_model 		= Channel.fromPath("${params.model_dir}")

//================output directory========================

def ref_ch = fq_channel.map { sid, item, scheme, version ->
	def ref = file("${params.primer_schema}/${scheme}/${version}/${scheme}.reference.fasta")
	if( !ref.exists() ) throw new IllegalArgumentException("Missing reference for ${scheme}/${version} -> ${ref}")
	tuple(sid, ref)
	}

def bed_ch = fq_channel.map { sid, item, scheme, version ->
	def bed = file("${params.primer_schema}/${scheme}/${version}/${scheme}.scheme.bed")
	if( !bed.exists() ) throw new IllegalArgumentException("Missing bed file for ${scheme}/${version} -> ${bed}")
	tuple(sid, bed)
	}

plex_dirs_channel = fq_channel.map { sample_id, item, scheme, version ->
	tuple(rawfile_dir, dir_plex_script, sample_id, item, scheme, version)
	}


vcf_filter_scr = fq_channel.map { sid, item, scheme, version ->
def script = file("${projectDir}/scripts/vcf_filter.py", checkIfExists: true)
	tuple(sid, script)
	}

make_depth_mask_scr = fq_channel.map { sid, item, scheme, version ->
def script = file("${projectDir}/scripts/make_depth_mask.py", checkIfExists: true)
	tuple(sid, script)
	}

mask_scr = fq_channel.map { sid, item, scheme, version ->
def script = file("${projectDir}/scripts/mask.py", checkIfExists: true)
	tuple(sid, script)
	}

fasta_header_scr = fq_channel.map { sid, item, scheme, version ->
def script = file("${projectDir}/scripts/fasta_header.py", checkIfExists: true)
	tuple(sid, script)
	}

// -----------------------------------------------------------------------------
//  - MINIMAP2 and downstream take the prepared inputs plus fq_channel.
// -----------------------------------------------------------------------------
workflow {

  if( params.rawfile_type == 'fastq' ) {
    
		//running PLEX_DIRS here
		PLEX_DIRS(plex_dirs_channel)

		reads_ch = PLEX_DIRS.out.reads
    .map { file ->
        def base     = file.baseName
        def sampleId = base.replaceFirst(/_barcode.*/, "")
        tuple(sampleId, file)
    }
    .groupTuple()
	
		meta_ch = fq_channel.map { sid, item, scheme, version ->
    tuple(sid, item, scheme, version)
		}

		// MINIMAP channel here
		minimap_channel = reads_ch
    	.join(meta_ch)
    	.join(ref_ch)

		//**********running MINIMAP here
		MINIMAP2(minimap_channel)

		//channel the align_trim script
		align_trim_scr = fq_channel.map { sid, item, scheme, version ->
    	def script = file("${projectDir}/scripts/align_trim.py", checkIfExists: true)
    	tuple(sid, script)
		}

		//ALIGN_TRIM_1 channel
		align_trim_channel = MINIMAP2.out.sorted_bam
				.join(meta_ch)   // → [sid, bam, sid, item, scheme, version]
				.join(bed_ch)    // → [sid, bam, sid, item, scheme, version, sid, bed]
  			.join(align_trim_scr) //[sid, bam, sid, item, scheme, version, sid, bed, align_trim_script]
	}
  else {
    if( params.basecaller == 'Dorado' ) {
      //DORADO_BASECALLER(basecaller_dir, basecaller_model, rawfile_dir)
      //DORADO_BARCODER(fastq_file: DORADO_BASECALLER.out)
      //PLEX_FQ_FILES(DORADO_BARCODER.out, fq_channel)
      //MINIMAP2(input_dir: PLEX_FQ_FILES.out.collect(), fq_channel)
    }
    else {
      // Guppy path
      //GUPPY_BASECALLER(fast5_or_pod5_dir: params.rawfile_dir)
      //GUPPY_BARCODER(fastq_file: GUPPY_BASECALLER.out)
      //PLEX_DIRS(input_dir: GUPPY_BARCODER.out, fq_channel)
      //MINIMAP2(input_dir: PLEX_DIRS.out.collect(), fq_channel)
    }
  }

	//align_trim_channel.view { row ->
  //  "ALIGN_TRIM_INPUT-1 >>> ${row.collect { it instanceof List ? "LIST(${it})" : it }}"
	//}

	//**********running ALIGN_TRIM_1
	ALIGN_TRIM_1(align_trim_channel)

	//ALIGN_TRIM_2 channel
	align_trim_2_channel = MINIMAP2.out.sorted_bam
		.join(meta_ch)   // → [sid, bam, sid, item, scheme, version]
		.join(bed_ch)    // → [sid, bam, sid, item, scheme, version, sid, bed]
		.join(align_trim_scr)

	//**********running ALIGN_TRIM_2
	ALIGN_TRIM_2(align_trim_2_channel)
	
	//MEDAKA-1 channel
	medaka_1_channel = ALIGN_TRIM_1.out.trimmed_bam
		.join(fq_channel)

	//**********running MEDAKA-1
	MEDAKA_1(medaka_1_channel)
	//ALIGN_TRIM_2.out.primertrimmed_bam.view()
	//medaka_1_channel.view()

	//MEDAKA-2 channel
	medaka_2_channel = ALIGN_TRIM_2.out.primertrimmed_bam
		.join(fq_channel)

	//**********running MEDAKA-2
	MEDAKA_2(medaka_2_channel)
	
	//MEDAKA_SNP_1 channel
	medaka_snp_1_channel = MEDAKA_1.out.hdf
		.join(ref_ch)
		.join(fq_channel)

	//**********running MEDAKA-2
	MEDAKA_SNP_1(medaka_snp_1_channel)

	//MEDAKA_SNP_2 channel
	medaka_snp_2_channel = MEDAKA_2.out.hdf
		.join(ref_ch)
		.join(fq_channel)
	
	//**********running MEDAKA_SNP_2
	MEDAKA_SNP_2(medaka_snp_2_channel)

	//channel the vcf_merge script
	vcf_merge_scr = fq_channel.map { sid, item, scheme, version ->
	def script = file("${projectDir}/scripts/vcf_merge.py", checkIfExists: true)
		tuple(sid, script)
	}

	//VCF_MERGE channel
	vcf_merge_channel = MEDAKA_SNP_2.out.vcf
		.join(MEDAKA_SNP_1.out.vcf)
		.join(bed_ch)
		.join(vcf_merge_scr)
		.join(fq_channel)

	//**********running VCF_MERGE
	VCF_MERGE(vcf_merge_channel)

	//LONGSHOT channel
	longshot_channel = VCF_MERGE.out.merged_tbi
		.join(ALIGN_TRIM_2.out.primertrimmed_bam)
		.join(ref_ch)
		.join(fq_channel)

	//**********running LONGSHOT
	LONGSHOT(longshot_channel)	


	//VCF_FILTER channel
	vcf_filter_channel = LONGSHOT.out.vcf
		.join(vcf_filter_scr)
		.join(fq_channel)

	//**********running VCF_FILTER
	VCF_FILTER(vcf_filter_channel)
	
	//MAKE_DEPTH_MASK channel
	make_depth_mask_channel = VCF_FILTER.out.pass_vcf
		.join(ALIGN_TRIM_2.out.primertrimmed_bam)
		.join(ref_ch)
		.join(make_depth_mask_scr)
		.join(fq_channel)

	//**********running VCF_FILTER
	MAKE_DEPTH_MASK(make_depth_mask_channel)

	//MASK channel
	mask_channel = MAKE_DEPTH_MASK.out.coverage_mask
		.join(VCF_FILTER.out.fail_vcf)	
		.join(ref_ch)
		.join(mask_scr)
		.join(fq_channel)

	//**********running VCF_FILTER
	MASK(mask_channel)

	//BCFTOOLS_CONSENSUS channel
	bcftools_consensus_channel = MASK.out.preconsensus
		.join(VCF_FILTER.out.pass_vcf)
		.join(MAKE_DEPTH_MASK.out.coverage_mask)
		.join(fq_channel)

	//**********running BCFTOOLS_CONSENSUS
	BCFTOOLS_CONSENSUS(bcftools_consensus_channel)
	

	//FASTA_HEADER channel
	fasta_header_channel = BCFTOOLS_CONSENSUS.out.consensus_fa
		.join(fasta_header_scr)
		.join(fq_channel)

	//**********running FASTA_HEADER
	FASTA_HEADER(fasta_header_channel)


	//CONCAT_FOR_MUSCLE channel
	contat_for_muscle = FASTA_HEADER.out.fasta
		.join(ref_ch)
		.join(fq_channel)

	//**********running CONCAT_FOR_MUSCLE
	CONCAT_FOR_MUSCLE(contat_for_muscle)


	//MUSCLE channel
	muscle_channel = CONCAT_FOR_MUSCLE.out.muscle_fa
		.join(fq_channel)

	//**********running MUSCLE
	MUSCLE(muscle_channel)

	//CONCAT channel
	concat_channel = MUSCLE.out.muscle_op_fasta.collect()

	CONCAT(concat_channel, concat_script)

	//MAFFT channel
	mafft_channel = CONCAT.out.genome_fa
	
	//**********running MAFFT
	MAFFT(mafft_channel)

	//**********running SUMMARY_STATS
	SUMMARY_STATS(MAFFT.out.mafft_fa, medaka_dir, summary_stats_script)

	//**********running REPORT
	REPORT(SUMMARY_STATS.out.summary, medaka_dir, summary_stats_dir, report_script)
}
