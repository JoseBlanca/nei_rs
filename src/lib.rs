use std::{
    env::var_os,
    io::{BufRead, BufReader, Empty, Read},
};

#[derive(Debug)]
struct Variant {}

pub struct Variants<'a> {
    samples: Option<Vec<String>>,
    vars_iter: Option<Box<dyn Iterator<Item = Result<Variant, VCFParseError>> + 'a>>,
}

impl<'a> Default for Variants<'a> {
    fn default() -> Variants<'a> {
        Variants {
            samples: None,
            vars_iter: None,
        }
    }
}

#[derive(thiserror::Error, Debug)]
pub enum VCFParseError {
    #[error("The line failed to define the fields and include the samples: `{0}`")]
    InvalidSampleLine(String),
    #[error("Error reading VCF line number: {0}")]
    ReadLineError(u64),
    #[error("The file is empty")]
    EmptyFile,
}

fn read_sample_line(line: &str) -> Result<Vec<String>, VCFParseError> {
    if !line.starts_with("#CHROM") {
        return Err(VCFParseError::InvalidSampleLine(line.to_string()));
    }

    let samples = line.split("\t").skip(9).map(|s| s.to_string()).collect();
    Ok(samples)
}

fn parse_variant_line(line: String) -> Result<Variant, VCFParseError> {
    println!("Parsing line: {}", line);
    Ok(Variant {})
}

pub fn parse_vcf<'a, T: Read + 'a>(mut file: BufReader<T>) -> Result<Variants<'a>, VCFParseError> {
    let mut vars = Variants::default();

    loop {
        let mut line = String::new();
        match file.read_line(&mut line) {
            Ok(0) => return Err(VCFParseError::EmptyFile),
            Ok(_) => (),
            Err(_) => return Err(VCFParseError::ReadLineError(0)),
        }
        if line.starts_with("##") {
        } else if line.starts_with("#CHROM") {
            let samples = read_sample_line(&line)?;
            vars.samples = Some(samples);
            break;
        } else {
            return Err(VCFParseError::InvalidSampleLine(line));
        }
    }

    let vars_iter = file.lines().map(|line_res| {
        let line = match line_res {
            Ok(line) => line,
            Err(_) => return Err(VCFParseError::ReadLineError(0)),
        };
        parse_variant_line(line)
    });
    vars.vars_iter = Some(Box::new(vars_iter));

    return Ok(vars);
}

#[cfg(test)]
mod tests {
    use self::super::*;

    // VCF example taken from https://samtools.github.io/hts-specs/VCFv4.5.pdf
    const VCF_45: &str = "##fileformat=VCFv4.5
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species=\"Homo sapiens\",taxonomy=x>
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">
##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">
##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">
##INFO=<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">
##FILTER=<ID=q10,Description=\"Quality below 10\">
##FILTER=<ID=s50,Description=\"Less than 50% of samples have data\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003
20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,.
20\t17330\t.\tT\tA\t3\tq10\tNS=3;DP=11;AF=0.017\tGT:GQ:DP:HQ\t0|0:49:3:58,50\t0|1:3:5:65,3\t0/0:41:3
20\t1110696\trs6040355\tA\tG,T\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4
20\t1230237\t.\tT\t.\t47\tPASS\tNS=3;DP=13;AA=T\tGT:GQ:DP:HQ\t0|0:54:7:56,60\t0|0:48:4:51,51\t0/0:61:2
20\t1234567\tmicrosat1\tGTC\tG,GTCT\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t0/1:35:4\t0/2:17:2\t1/1:40:3";

    #[test]
    fn it_works() {
        let mock_file = BufReader::new(VCF_45.as_bytes());
        let vars = parse_vcf(mock_file).expect("Error");
        let vars_iter = match vars.vars_iter {
            Some(iter) => iter,
            None => panic!("Vars iter is None"),
        };
        for var_res in vars_iter {
            let var = var_res.expect("Error reading variant");
            println!("{:?}", var);
        }
    }
}
