#!/usr/bin/env nextflow

// Input of sample names, conditions, and FAST5s path.
Channel
    .fromPath(params.samples)
    .splitCsv(header: true, sep:'\t')
    .map{ row-> tuple(row.SampleName, row.Condition, file(row.DataPath)) }
    .set{samplesFile_gzipCompression}

// FAST5s gzip compression.
process gzipCompression {
    input:
        tuple val(sample), val(condition), file(f5) from samplesFile_gzipCompression

    output:
        tuple val(sample), val(condition), file(f5) into gzipCompression_singleToMulti

    script:
    if(params.gzipCompression)
    """
        mkdir -p ${params.resultsDir}
        mkdir -p ${params.resultsDir}/${condition}
        mkdir -p ${params.resultsDir}/${condition}/${sample}

        fast5=\$(grep ${sample} ${params.samples} | cut -f3)
        
        compress_fast5 -i \$fast5 -s ${params.resultsDir}/${condition}/${sample}/FAST5/ -c gzip -t ${task.cpus} --recursive
    """
    else
    """
        echo "Skipped."
    """
}

// From single- to multi-reads FAST5s.
process singleToMulti {
    input:
        tuple val(sample), val(condition), file(f5) from gzipCompression_singleToMulti

    output:
        tuple val(sample), val(condition) into singleToMulti_fastqExtraction

    script:
    if(params.singleToMulti)
    """
        for dir in /${params.resultsDir}/${condition}/${sample}/FAST5/*/
        do
            nameTmp=\$(basename \$dir)
            single_to_multi_fast5 -i \$dir -s /${params.resultsDir}/${condition}/${sample}/FAST5/ -f \$nameTmp -t 1
            
            find \$dir -name '*.fast5' -type f | parallel -j ${task.cpus} rm --

            rm -r \$dir
        done
    """
    else
    """
        echo "Skipped"
    """
}

// From FAST5s to FASTQs.
process fastqExtraction {
    input:
        tuple val(sample), val(condition) from singleToMulti_fastqExtraction

    output:
        tuple val(sample), val(condition) into fastqExtraction_fastqFiltering

    script:
    if(params.fastqExtraction)
    """
        mkdir -p ${params.resultsDir}/${condition}/${sample}/FASTQ_DNA/
        mkdir -p ${params.resultsDir}/${condition}/${sample}/nanoid/

        Rscript ${params.scripts}/nanoID.fastq.extraction.R FAST5paths=${params.resultsDir}/${condition}/${sample}/FAST5/ genomeFasta=${params.genome_fasta} mc.cores=${task.cpus} scriptsFolder=${params.scripts} resultsFolder=${params.resultsDir}/${condition}/${sample}/

        cat ${params.resultsDir}/${condition}/${sample}/FASTQ_DNA/*.fastq > ${params.resultsDir}/${condition}/${sample}/FASTQ_DNA_ALL.fastq
    """
    else
    """
        cat ${params.resultsDir}/${condition}/${sample}/FASTQ_DNA/*.fastq > ${params.resultsDir}/${condition}/${sample}/FASTQ_DNA_ALL.fastq
    """
}

// Fastq filtering.
process fastqFiltering {
    input:
        tuple val(sample), val(condition) from fastqExtraction_fastqFiltering

    output:
        tuple val(sample), val(condition) into fastqFiltering_genomicAlignment
        tuple val(sample), val(condition) into fastqFiltering_matureTranscriptomeAlignment
        tuple val(sample), val(condition) into fastqFiltering_prematureTranscriptomeAlignment

    script:
    if(params.fastqExtraction)
    """
        cat ${params.resultsDir}/${condition}/${sample}/FASTQ_DNA_ALL.fastq | NanoFilt -q ${params.qvalue} > ${params.resultsDir}/${condition}/${sample}/FASTQ_DNA.fastq
    """
    else
    """
        cp ${params.resultsDir}/${condition}/${sample}/FASTQ_DNA_ALL.fastq ${params.resultsDir}/${condition}/${sample}/FASTQ_DNA.fastq
    """
}

// Genomic alignment.
process genomicAlignment {
    input:
        tuple val(sample), val(condition) from fastqFiltering_genomicAlignment

    output:
        tuple val(sample), val(condition) into genomicAlignment_NID_alignmentExtraction
        tuple val(sample), val(condition) into genomicAlignment_NID_readExtraction

    script:
    if(params.genomicAlignment)
    """
        mkdir -p ${params.resultsDir}/${condition}/${sample}/genomeAlignment/
        
        minimap2 -ax splice -k14 --seed 1 -t ${task.cpus} ${params.genome_fasta} ${params.resultsDir}/${condition}/${sample}/FASTQ_DNA.fastq | samtools view -hSb | samtools sort -@ ${task.cpus} -o ${params.resultsDir}/${condition}/${sample}/genomeAlignment/minimap.G.bam
        samtools view ${params.resultsDir}/${condition}/${sample}/genomeAlignment/minimap.G.bam -bh -t genome.fa.fai -F 2308 -q ${params.alignmentqValue} | samtools sort -@ ${task.cpus} -o ${params.resultsDir}/${condition}/${sample}/genomeAlignment/minimap.filt.G.bam
        samtools index -@ ${task.cpus} ${params.resultsDir}/${condition}/${sample}/genomeAlignment/minimap.filt.G.bam

        ln -s ${params.resultsDir}/${condition}/${sample}/genomeAlignment/minimap.filt.G.bam ./minimap.filt.G.bam
        ln -s ${params.resultsDir}/${condition}/${sample}/genomeAlignment/minimap.filt.G.bam.bai ./minimap.filt.G.bam.bai
    """
    else
    """
        ln -s ${params.resultsDir}/${condition}/${sample}/genomeAlignment/minimap.filt.G.bam ./minimap.filt.G.bam
        ln -s ${params.resultsDir}/${condition}/${sample}/genomeAlignment/minimap.filt.G.bam.bai ./minimap.filt.G.bam.bai
    """
}

// nanoID - alignments extraction.
process NID_alignmentExtraction {
    input:
        tuple val(sample), val(condition) from genomicAlignment_NID_alignmentExtraction

    output:
        tuple val(sample), val(condition) into NID_alignmentExtraction_NID_readExtraction

    script:
    if(params.NID_alignmentExtraction)
    """
        Rscript ${params.scripts}/nanoID.alignment.extraction.R genomeFasta=${params.genome_fasta} mc.cores=${task.cpus} scriptsFolder=${params.scripts} resultsFolder=${params.resultsDir}/${condition}/${sample}/ bamPath=${params.resultsDir}/${condition}/${sample}/genomeAlignment/minimap.filt.G.bam
    """
    else
    """
        echo "Skipped"
    """
}

// nanoID - reads extraction.
process NID_readExtraction {
    input:
        tuple val(sample), val(condition) from NID_alignmentExtraction_NID_readExtraction

    output:
        tuple val(sample), val(condition) into NID_readExtraction_NID_alignmentReconstruction
        tuple val(sample), val(condition) into NID_readExtraction_NID_fivemerAlignmentReconstruction
        tuple val(sample), val(condition) into NID_readExtraction_NID_traceModel
        tuple val(sample), val(condition) into NID_readExtraction_NID_traceModelAdd
        tuple val(sample), val(condition) into NID_readExtraction_NID_sequencingSummary         
    
    script:
    if(params.NID_readExtraction)
    """
       Rscript ${params.scripts}/nanoID.read.extraction.R FAST5paths=${params.resultsDir}/${condition}/${sample}/FAST5/ genomeFasta=${params.genome_fasta} mc.cores=${task.cpus} scriptsFolder=${params.scripts} resultsFolder=${params.resultsDir}/${condition}/${sample}/ bamPath=${params.resultsDir}/${condition}/${sample}/genomeAlignment/minimap.filt.G.bam
    """
    else
    """
        echo "Skipped"
    """
}

// nanoID - alignment reconstruction.
process NID_alignmentReconstruction {
    input:
        tuple val(sample), val(condition) from NID_readExtraction_NID_alignmentReconstruction

    output:
        tuple val(sample), val(condition) into NID_alignmentReconstruction_NID_fivemerAlignmentReconstruction
	tuple val(condition), val(sample) into NID_alignmentReconstruction_NID_dataFormatting

    script:
    if(params.NID_alignmentReconstruction)
    """
        Rscript ${params.scripts}/nanoID.alignment.reconstruction.R scriptsFolder=${params.scripts} resultsFolder=${params.resultsDir}/${condition}/${sample}/ mc.cores=${task.cpus}
    """
    else
    """
        echo "Skipped"
    """
}

// nanoID - fivemers alignment reconstruction.
process NID_fivemerAlignmentReconstruction {
    input:
        tuple val(sample), val(condition) from NID_alignmentReconstruction_NID_fivemerAlignmentReconstruction

    output:
        tuple val(sample), val(condition) into NID_fivemerAlignmentReconstruction_NID_mismatchReadIdentification

    script:
    if(params.NID_fivemerAlignmentReconstruction)
    """
        Rscript ${params.scripts}/nanoID.five.mer.alignment.reconstruction.R scriptsFolder=${params.scripts} resultsFolder=${params.resultsDir}/${condition}/${sample}/ mc.cores=${task.cpus} nanoIDSupplementals=${params.nanoIDSupplementals}
    """
    else
    """
        echo "Skipped"
    """
}

// nanoID - mismatch read identification.
process NID_mismatchReadIdentification {
    input:
        tuple val(sample), val(condition) from NID_fivemerAlignmentReconstruction_NID_mismatchReadIdentification

    output:
        tuple val(condition), val(sample) into NID_mismatchReadIdentification_NID_dataFormatting

    script:
    if(params.NID_mismatchReadIdentification)
    """
        Rscript ${params.scripts}/nanoID.mismatch.read.identification.R scriptsFolder=${params.scripts} resultsFolder=${params.resultsDir}/${condition}/${sample}/ mc.cores=${task.cpus} nanoIDSupplementals=${params.nanoIDSupplementals}
    """
    else
    """
        echo "Skipped"
    """
}

// nanoID - trace model.
process NID_traceModel {
    input:
        tuple val(sample), val(condition) from NID_readExtraction_NID_traceModel    

    output:
        tuple val(condition), val(sample) into NID_traceModel_NID_dataFormatting

    script:
    if(params.NID_traceModel)
    """
        Rscript ${params.scripts}/nanoID.trace.model.R scriptsFolder=${params.scripts} resultsFolder=${params.resultsDir}/${condition}/${sample}/ mc.cores=${task.cpus} nanoIDSupplementals=${params.nanoIDSupplementals} moveSlot=${params.moveSlot}
    """
    else
    """
        echo "Skipped"
    """
}

// nanoID - trace model add.
process NID_traceModelAdd {
    input:
        tuple val(sample), val(condition) from NID_readExtraction_NID_traceModelAdd    

    output:
        tuple val(condition), val(sample) into NID_traceModelAdd_NID_dataFormatting

    script:
    if(params.NID_traceModelAdd)
    """
        Rscript ${params.scripts}/nanoID.trace.model.add.R scriptsFolder=${params.scripts} resultsFolder=${params.resultsDir}/${condition}/${sample}/ mc.cores=${task.cpus} nanoIDSupplementals=${params.nanoIDSupplementals} moveSlot=${params.moveSlot}
    """
    else
    """
        echo "Skipped"
    """
}

// nanoID - sequencing summary.
process NID_sequencingSummary {
    input:
        tuple val(sample), val(condition) from NID_readExtraction_NID_sequencingSummary

    output:
        tuple val(sample), val(condition) into NID_sequencingSummary_NID_rawSignal
        tuple val(sample), val(condition) into NID_sequencingSummary_NID_rawSignalFive
        tuple val(sample), val(condition) into NID_sequencingSummary_NID_rawSignalFiveAdd

    script:
    if(params.NID_sequencingSummary)
    """        
        fast5=\$(grep ${sample} ${params.samples} | cut -f3)
        FAST5paths=\$fast5

        Rscript ${params.scripts}/nanoID.sequencing.summary.R FAST5paths=\$fast5 scriptsFolder=${params.scripts} resultsFolder=${params.resultsDir}/${condition}/${sample}/ mc.cores=${task.cpus}
    """
    else
    """
        echo "Skipped"
    """
}

// nanoID - raw signal.
process NID_rawSignal {
    input:
        tuple val(sample), val(condition) from NID_sequencingSummary_NID_rawSignal    

    output:
        tuple val(condition), val(sample) into NID_rawSignal_NID_dataFormatting

    script:
    if(params.NID_rawSignal)
    """
        Rscript ${params.scripts}/nanoID.raw.signal.R scriptsFolder=${params.scripts} resultsFolder=${params.resultsDir}/${condition}/${sample}/ mc.cores=${task.cpus} nanoIDSupplementals=${params.nanoIDSupplementals} moveSlot=${params.moveSlot}
    """
    else
    """
        echo "Skipped"
    """
}

// nanoID - raw signal five.
process NID_rawSignalFive {
    input:
        tuple val(sample), val(condition) from NID_sequencingSummary_NID_rawSignalFive    

    output:
        tuple val(condition), val(sample) into NID_rawSignalFive_NID_dataFormatting

    script:
    if(params.NID_rawSignalFive)
    """
        Rscript ${params.scripts}/nanoID.raw.signal.five.R scriptsFolder=${params.scripts} resultsFolder=${params.resultsDir}/${condition}/${sample}/ mc.cores=${task.cpus} nanoIDSupplementals=${params.nanoIDSupplementals} moveSlot=${params.moveSlot}
    """
    else
    """
        echo "Skipped"
    """
}

// nanoID - raw signal five add.
process NID_rawSignalFiveAdd {
    input:
        tuple val(sample), val(condition) from NID_sequencingSummary_NID_rawSignalFiveAdd    

    output:
        tuple val(condition), val(sample) into NID_rawSignalFiveAdd_NID_dataFormatting
        tuple val(condition), val(sample) into NID_rawSignalFiveAdd_NID_dataFormatting2

    script:
    if(params.NID_rawSignalFiveAdd)
    """
        Rscript ${params.scripts}/nanoID.raw.signal.five.add.R scriptsFolder=${params.scripts} resultsFolder=${params.resultsDir}/${condition}/${sample}/ mc.cores=${task.cpus} nanoIDSupplementals=${params.nanoIDSupplementals} moveSlot=${params.moveSlot}
    """
    else
    """
        echo "Skipped"
    """
}


// nanoID - data formatting.
process NID_dataFormatting {
    input:
        tuple val(condition), val(sample) from NID_alignmentReconstruction_NID_dataFormatting.collect()
        tuple val(condition), val(sample) from NID_mismatchReadIdentification_NID_dataFormatting.collect()
        tuple val(condition), val(sample) from NID_traceModel_NID_dataFormatting.collect()
        tuple val(condition), val(sample) from NID_traceModelAdd_NID_dataFormatting.collect()
        tuple val(condition), val(sample) from NID_rawSignal_NID_dataFormatting.collect()
        tuple val(condition), val(sample) from NID_rawSignalFive_NID_dataFormatting.collect()
        tuple val(condition), val(sample) from NID_rawSignalFiveAdd_NID_dataFormatting.collect()
    
    output:
        tuple val(condition), val(sample) into NID_dataFormatting_NID_training

    script:
    if(params.NID_dataFormatting)
    """
        cat ${params.samples} | tail -n+2 | while read line || [ -n "\$line" ]
        do
           sample=\$(echo \$line | cut -f1 -d ' ')
           condition=\$(echo \$line | cut -f2 -d ' ')
           str=${params.resultsDir}"/"\$condition"/"\$sample/
           echo \$str

           rm -f \$str/nanoid/read.based.parameters*

           Rscript ${params.scripts}/nanoID.data.formatting.R scriptsFolder=${params.scripts} resultsFolder=\$str mc.cores=${task.cpus}
        done
    """
    else
    """
        echo "Skipped"
    """
}

// nanoID - training
process NID_training {
    input:
        tuple val(condition), val(sample) from NID_dataFormatting_NID_training
    
    output:
        tuple val(condition), val(sample) into NID_training_NID_classification

    script:
    if(params.NID_training)
    """
        mkdir -p ${params.resultsDir}/models/
        Rscript ${params.scripts}/nanoID.training.R scriptsFolder=${params.scripts} resultsDir=${params.resultsDir} unlabeled_time=${params.unlabeled_time} fullylabeled_time=${params.fullylabeled_time} trainedModelFolder=${params.resultsDir}/models/ Nsamplings=${params.Nsamplings} Nseeds=${params.Nseeds} trainingSizes=${params.trainingSizes} genesBed=${params.genesBed}
    """
    else
    """
        if [ -f ${params.nanoIDtrainedModel} ]; then
            mkdir -p ${params.resultsDir}/models/
            cp ${params.nanoIDtrainedModel} ${params.resultsDir}/models/
        fi
    """
}

// nanoID - classification.
process NID_classification {
    input:
        tuple val(condition), val(sample) from NID_training_NID_classification
    
    output:
        tuple val(condition), val(sample) into NID_classification_kineticRatesEstimation

    script:
    if(params.NID_classification)
    """
        cat ${params.samples} | tail -n+2 | while read line || [ -n "\$line" ]
        do
           sample=\$(echo \$line | cut -f1 -d ' ')
           condition=\$(echo \$line | cut -f2 -d ' ')
           str=${params.resultsDir}"/"\$condition"/"\$sample/
           echo \$str
           Rscript ${params.scripts}/nanoID.classification.R scriptsFolder=${params.scripts} resultsFolder=\$str mc.cores=${task.cpus} trainedModelFolder=${params.resultsDir}/models/
        done
    """
    else
    """
        echo "Skipped"
    """
}

// Kinetic rates estimation.
process kineticRatesEstimation {
    input:
        tuple val(sample), val(condition) from NID_classification_kineticRatesEstimation.collect()
    
    output:

    script:
    if(params.kineticRatesEstimation)
    """
        mkdir -p ${params.resultsDir}/kineticRates/
        Rscript ${params.scripts}/kineticRatesEstimation.R scriptsFolder=${params.scripts} resultsDir=${params.resultsDir} CCL=${params.CCL} genesBed=${params.genesBed} unlabeled_time=${params.unlabeled_time} fullylabeled_time=${params.fullylabeled_time} trainedModelFolder=${params.resultsDir}/models/

    """
    else
    """
        echo "Skipped."
    """
}
