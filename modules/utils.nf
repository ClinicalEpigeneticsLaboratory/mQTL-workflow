def loadMynorm ( dir ) {
    if ( file("$dir/flow_two/mynorm_corrected.parquet").exists() )  {
        return file( "$dir/flow_two/mynorm_corrected.parquet" )
    }
    else if ( file("$dir/flow_two/mynorm.parquet").exists() ){
        return file( "$dir/flow_two/mynorm.parquet" )
    }
    else {
        error("Mynorm file does not exists in $dir/flow_two directory! Run flow_two.nf to generate necessary files.")
    }
    
}

def loadGenotypeFrame ( dir ) {
    if ( file("$dir/flow_one/genotype_table.parquet").exists() ) {
        return file("$dir/flow_one/genotype_table.parquet")
    }
    else {
        error("genotype table does not exists in $dir/flow_one/ directory!")
    }
}

def loadVCF ( dir ) {
    if ( file("$dir/flow_one/genotypes.vcf.gz").exists() ) {
        return file("$dir/flow_one/genotypes.vcf.gz")
    }
    else {
        error("VCF file does not exists in $dir/flow_one/ directory!")
    }
 
}

def adjustCPUs ( CPUs ) {
    if ( CPUs > 1 ){
        return Math.round( CPUs / 2 )
    } 
    else {
        return 1
    }
}