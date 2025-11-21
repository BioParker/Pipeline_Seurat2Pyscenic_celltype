#!/usr/bin/env nextflow

/*
 * Pipeline parameters
 */


//Output directories

params.loomdir="loom.dir"
params.grndir="grn.dir"
params.corrdir="add_corr.dir"
params.ctxdir="ctx.dir"
params.aucdir="auc.dir"
params.outfilesdir="outfiles.dir"


/*
 * Processes
 */


//split seurat object into cell type loom objects
process mk_loom_files {

        conda params.rconda

        publishDir   params.loomdir, mode: 'symlink'

	input:
            val celltype
            path seurat
            val anncol

	output:
	    tuple path("${celltype}.loom"), val("${celltype}")

	script:
	"""
        Rscript '$projectDir'/scripts/mk_loom_files.R --seurat '$seurat' --celltype '$celltype' --cellcol '$anncol'
	"""
}


//Derive co-expression modules from expression matrix
process grn {

        conda params.psconda

        publishDir  params.grndir, mode: 'symlink'

        input:
            tuple path(loomfile), val(celltype)
            path tffile

        output:
            tuple path("${celltype}_adj.csv"), path("${loomfile}"), val("${celltype}")

        script:
        """
        pyscenic grn \
        --num_workers 10 \
        -o '$celltype'_adj.csv \
        '$loomfile' \
        '$tffile'
        """
}

//Add Pearson correlations based on TF-gene expression 
//to the network adjacencies output from the GRN step, 
//and output these to a new adjacencies file
process add_cor {

        conda params.psconda

        publishDir params.corrdir, mode: 'symlink'
        
        input:
            tuple path(adj_file), path(loomfile), val(celltype)

        output:
            tuple val("${celltype}"), path("${loomfile}"), path("${celltype}_adj_wcor.csv")

        script:
        """
        pyscenic add_cor \
        '$adj_file' \
        '$loomfile' \
        -o '$celltype'_adj_wcor.csv \
        --mask_dropouts \
        --sparse \
        """
}

//Find enriched motifs for a gene signature and
//optionally prune targets from this signature 
//based on cis-regulatory cues
process ctx {

        conda params.psconda

        publishDir params.ctxdir, mode: 'symlink'

        input:
            tuple val(celltype), path(loomfile), path(adj_file)
            path fthrfiles
            path motif2tf

        output:
            tuple val("${celltype}"), path("$loomfile"), path("${celltype}_regulons.csv")

        script:
        """
        pyscenic ctx \
        '$adj_file' \
	${fthrfiles.collect{"'$it'"}.join(' ')} \
        --annotations_fname '$motif2tf' \
        --expression_mtx_fname '$loomfile' \
        --mode "custom_multiprocessing" \
        --output '$celltype'_regulons.csv \
	--num_workers 10
        """

}

//
process auc {

        conda params.psconda

        publishDir params.aucdir, mode: 'symlink'

        input:
            tuple val(celltype), path(loomfile), path(regfile)

        output:
            tuple val("${celltype}"), path("${celltype}_auc_mtx.csv")

        script:
        """
        pyscenic aucell \
        $loomfile \
        '$regfile' \
        -o '$celltype'_auc_mtx.csv \
        --num_workers 10
        """

}

process consolidate_files {

        conda params.psconda

        publishDir params.outfilesdir

        input:
            tuple val(celltype), path(aucfile), path(loomfile), path(regfile), path(corrfile)

        output:
            tuple path("${celltype}/${loomfile}"), path("${celltype}/${corrfile}"), path("${celltype}/${regfile}"), path("${celltype}/${aucfile}")
            
        script:
        """
        mkdir $celltype;
        cp '$loomfile' '$celltype'/'$loomfile';
        cp '$corrfile' '$celltype'/'$corrfile';
        cp '$regfile' '$celltype'/'$regfile';
        cp '$aucfile' '$celltype'/'$aucfile'
        """

}

/*
 * Workflow
 */


workflow {

    //create input channel
    celltype = Channel.fromPath(params.ctfile)
	       .splitCsv()
	       .flatten()
    
    //outdirs

    //main files
    seurat = params.seuratobject

    //script path

    //accessory files
    tffile = file(params.tffile)
    fthrfiles = file(params.fthrfiles)
    motif2tf = file(params.motif2tf)

    //script options
    anncol = params.ctcol
    
    mk_loom_files(celltype,
	          seurat,
                  anncol)

    grn(mk_loom_files.out,
        tffile)

    add_cor(grn.out)

    ctx(add_cor.out,
        fthrfiles,
        motif2tf)

    auc(ctx.out)

    def publSets = auc.out.collect(flat: false)
			  .flatten()
                          .collate(2)
			  .concat(ctx.out, add_cor.out.map { it -> [ it[0], it[2] ] } )
			  .groupTuple()
			  .map { it -> it.flatten() }
			  

    consolidate_files(publSets)

}


