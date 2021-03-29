#WORKFLOW GangSTR : Genotype STRs on all autosomes
workflow GangSTRPipeline {
  String pipeline_version = "1.0"
  meta {
	author: "Bharati Jadhav"
    	email: "bharati.jadhav@mssm.edu"
    	description: "Genotype short tandem repeats using GangSTR"	
  }
  File bamFile
  File baiFile
  File RefFasta
  File ref_fasta_index
  File STR_BED
  String sampleName =basename(bamFile, ".cram")
  #String gotc_docker = "gymreklab/str-toolkit:latest"
  #String gotc_docker2 = "bharatij/bcftools:1.1"
  Int preemptible_tries = 3
  Array[String] chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22']
  scatter (chrom in chroms) {
      	call RunGangSTR { 
        input: 
          RefFasta = RefFasta,
          ref_fasta_index = ref_fasta_index,
          bamFile= bamFile,
          baiFile= baiFile,
          sampleName = sampleName,
          STR_BED = STR_BED,
          chr = chrom,
          preemptible_tries = preemptible_tries
    }
    
    call RunTabix { 
        input: 
          VCF=RunGangSTR.rawVCF,
          preemptible_tries = preemptible_tries
    }
    
    call RunDumpSTR { 
        input: 
          VCF=RunTabix.gzrawFiltVCF,
          preemptible_tries = preemptible_tries
    }
    
    call RunVcfToTab { 
        input: 
            VCF=RunDumpSTR.filtVCF,
            preemptible_tries = preemptible_tries
    }
        
  }
  
  
  output {
    Array[File] samplestats = RunGangSTR.samplestats
    Array[File] gzrawVCF = RunTabix.gzrawVCF
    Array[File] idxrawVCF = RunTabix.idxrawVCF
    Array[File] gzrawFiltVCF = RunTabix.gzrawFiltVCF
    Array[File] idxrawFiltVCF = RunTabix.idxrawFiltVCF
    Array[File] gzVCF = RunVcfToTab.gzVCF
    Array[File] idxVCF = RunVcfToTab.idxVCF
    Array[File] filtbed = RunVcfToTab.filtbed      
  } 
  
}
#Task Definitions
task RunGangSTR {
    File RefFasta
    File ref_fasta_index
    File bamFile
    File baiFile
    String sampleName
    File STR_BED
    String chr
    # Runtime parameters
    Int addtional_disk_size = 2
    String machine_mem_size = 4
    Int preemptible_tries
    Float output_vcf_size = size(bamFile, "GB") * 0.30
    Float ref_size = size(RefFasta, "GB") 
    Float refidx_size = size(ref_fasta_index, "GB") 
    Float str_size = size(STR_BED, "GB")
    Int disk_size = ceil(size(bamFile, "GB") + size(baiFile, "GB") + output_vcf_size + ref_size + refidx_size + str_size) + addtional_disk_size
  
    command {
  		GangSTR \
        	--bam  ${bamFile} \
        	--ref ${RefFasta} \
        	--regions ${STR_BED} \
        	--chrom ${chr} \
        	--out ${sampleName}.${chr}
  }
  
  runtime {
    docker: "gymreklab/str-toolkit:latest"
    memory: machine_mem_size + " GB"
    disks: "local-disk " + disk_size + " SSD"
    cpu: "1"
    zones: "us-central1-c us-central1-b"
    preemptible: preemptible_tries   
    continueOnReturnCode: [0,1]
  }
  
  output {
    File rawVCF = "${sampleName}.${chr}.vcf"
    File samplestats = "${sampleName}.${chr}.samplestats.tab"
  }
}

task RunTabix {
    File VCF
    String outVCF = basename(VCF, ".vcf") + ".rawFilt.vcf"
   
    Int addtional_disk_size = 1
    String machine_mem_size = 2
    
    Int preemptible_tries
    Int disk_size = ceil(size(VCF, "GB") + size(VCF, "GB")) + addtional_disk_size
   
    command {
   		bgzip ${VCF}	
   	   	tabix -p vcf ${VCF}.gz
   		zcat  ${VCF}.gz | grep -v -w "REF=0"  > ${outVCF}
        	bgzip ${outVCF}	
      		tabix -p vcf ${outVCF}.gz
    }  
    
    runtime {
    	docker: "bharatij/bcftools:1.1" 
    	memory: machine_mem_size + " GB"
    	disks: "local-disk " + disk_size + " SSD"
    	cpu: "1"
    	zones: "us-east1-c us-east1-b"
    	preemptible: preemptible_tries
    	continueOnReturnCode: [0]
    }
   
    output {
      	File gzrawVCF = "${VCF}.gz"
      	File idxrawVCF = "${VCF}.gz.tbi"
      	File gzrawFiltVCF = "${outVCF}.gz"
      	File idxrawFiltVCF = "${outVCF}.gz.tbi"
    }
}

task RunDumpSTR {
    File VCF
    Int addtional_disk_size = 2
    String machine_mem_size = 2
   
    Int preemptible_tries
    String outVCF = basename(VCF, ".vcf.gz") + ".DumpSTRfilt"
    Int maxDP = 1000
    Int minDP = 10
    Int disk_size = ceil(size(VCF, "GB") + size(VCF, "GB")) + addtional_disk_size
    command {
    	  	dumpSTR \
    		--vcf ${VCF} \
        	--out ${outVCF} \
        	--filter-span-only \
        	--filter-badCI \
        	--max-call-DP ${maxDP} \
        	--min-call-DP ${minDP} \
        	--drop-filtered
    }  
    runtime {
    	docker: "gymreklab/str-toolkit:latest"
    	memory: machine_mem_size + " GB"
    	disks: "local-disk " + disk_size + " SSD"
    	cpu: "1"
    	zones: "us-central1-c us-central1-b"
    	preemptible: preemptible_tries
    	continueOnReturnCode: [0,1]
    }
   
    output {
    	File filtVCF = "${outVCF}.vcf"   
    }
}

task RunVcfToTab {
    File VCF
    String outBed = basename(VCF, ".vcf") + ".bed"
   
    Int addtional_disk_size = 1
    String machine_mem_size = 2
    
    Int preemptible_tries
    Int disk_size = ceil(size(VCF, "GB") + size(VCF, "GB")) + addtional_disk_size
   
    command {
   		bgzip ${VCF}	
   	   	tabix -p vcf ${VCF}.gz
   	        echo -e "Chrom\tStart\tEnd\tRef\tAlt\tPeriod\tRepeatUnit\tRefCopy\tSampleId\tGT\tQ\tDP\tA1A2_Copy\tRD_EnSpFrrFl" > ${outBed}
   	        bcftools query  -f '%CHROM\t%POS\t%END\t%REF\t%ALT\t%PERIOD\t%RU\t[%INFO/REF\t%SAMPLE\t%GT\t%Q\t%DP\t%REPCN\t%RC]\n' ${VCF}.gz  >> ${outBed}
        
    }  
    
    runtime {
    	docker: "bharatij/bcftools:1.1"
    	memory: machine_mem_size + " GB"
    	disks: "local-disk " + disk_size + " SSD"
    	cpu: "1"
    	zones: "us-east1-c us-east1-b"
    	preemptible: preemptible_tries
    	continueOnReturnCode: [0]
    }
   
    output {
      	File gzVCF = "${VCF}.gz"
      	File idxVCF = "${VCF}.gz.tbi"
      	File filtbed = "${outBed}"
    }
}
