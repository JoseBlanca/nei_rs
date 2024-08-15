use nei_rs::read_vcf_file;

fn main() {
    use std::path::PathBuf;
    let mut test_vcf = PathBuf::new();
    let path = "/home/jose/analyses/g2psol/variants/as.haplotypes.vcf.gz.all";
    test_vcf.push(path);
    match read_vcf_file(&test_vcf) {
        Ok(variants) => {
            for var in variants.vars_iter {
                //println!("{:?}", var);
            }
        }
        Err(e) => {
            eprintln!("Error reading VCF file: {}", e);
        }
    };
}
