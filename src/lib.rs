use flate2::read::MultiGzDecoder;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::PathBuf;

const GT_FIELD_ID: &str = "GT";
const MISSING_ALLELE: i16 = -1;

#[derive(thiserror::Error, Debug)]
pub enum VCFParseError {
    #[error("The line failed to define the fields and include the samples: `{0}`")]
    InvalidSampleLine(String),
    #[error("Error reading VCF line number: {0}")]
    ReadLineError(u64),
    #[error("The file is empty")]
    EmptyFile,
    #[error("No variants found in the file")]
    NoVariantsError,
    #[error("Position not a valid integer in line: `{0}`")]
    PosNotInt(String),
    #[error("Qual is not a valid float in line: `{0}`")]
    QualNotFloat(String),
    #[error("Error parsing genotype format in line: `{0}`")]
    NoGenotypeFormatDefinition(String),
    #[error("There's no GT field in the GT format definition in line: `{0}`")]
    GenotypeNotFoundInFormatDefinition(String),
    #[error("GT format definition and GT do not match in line: `{0}`")]
    GtOutsideBounds(String),
    #[error("Incorrect allele `{0}` in line: `{1}`")]
    IncorrectAllele(String, String),
    #[error("Different ploidies found in line: `{0}`")]
    DifferentPloidiesError(String),
    #[error("Error parsing GTs in line: `{0}`")]
    GtParseError(String),
    #[error("File is not gzip and does not start with ##: `{0}`")]
    InvalidVCFFile(String),
    #[error("File is gzip, but does not start with ##: `{0}`")]
    InvalidGzipVCFFile(String),
    #[error("First GT `{0}` does not define ploidy in first variant line: `{1}`")]
    FirstGtDoesNotDefinePloidy(String, String),
}

#[derive(Debug)]
pub struct Variant {
    chrom: String,
    pos: u64,
    id: String,
    alleles: Vec<String>,
    qual: f64,
    filters: Vec<String>,
    gts: Vec<Vec<i16>>,
    ploidy: u8,
}

struct GtFormatCache {
    gt_string: String,
    gt_format_idxs: HashMap<String, usize>,
    gt_field_idx: usize,
    num_samples: usize,
    ploidy: u8,
}

fn get_ploidy_form_first_gt(
    gt: &str,
    gt_format_cache: &mut GtFormatCache,
) -> Result<u8, VCFParseError> {
    let gt = get_gt_item_from_gt_string(gt, gt_format_cache)?;
    let alleles: Vec<&str> = gt.split(|c| c == '/' || c == '|').collect();
    let ploidy = alleles.len();
    Ok(ploidy as u8)
}

fn parse_gt<'a>(
    gt: &'a str,
    sample_idx: usize,
    parsed_gts: &mut Vec<Vec<i16>>,
    line: &String,
) -> Result<u8, VCFParseError> {
    if gt == "0/0" {
        return Ok(2);
    } else if gt == "1/1" {
        parsed_gts[sample_idx][0] = 1;
        parsed_gts[sample_idx][1] = 1;
        return Ok(2);
    }

    let mut allele = 0;
    let mut ploidy_idx = 0;
    let mut allele_was_missing = false;
    for chr in gt.bytes() {
        allele *= 10;
        let digit = chr & 0b0000_1111;
        if digit < 10 {
            allele += digit as i16;
        } else if digit == 12 || digit == 15 {
            // chr is / or |
            parsed_gts[sample_idx][ploidy_idx] = allele;
            allele = 0;
            ploidy_idx += 1;
            allele_was_missing = false;
        } else if digit == 14 && allele == 0 {
            // chr is .
            parsed_gts[sample_idx][ploidy_idx] = MISSING_ALLELE;
            allele_was_missing = true;
            ploidy_idx += 1;
        } else {
            return Err(VCFParseError::IncorrectAllele(
                chr.to_string(),
                line.to_string(),
            ));
        }
    }
    if !allele_was_missing {
        parsed_gts[sample_idx][ploidy_idx] = allele
    };
    Ok((ploidy_idx + 1) as u8)
}

fn get_gt_item_from_gt_string<'a>(
    gt_str: &'a str,
    gt_format_cache: &mut GtFormatCache,
) -> Result<&'a str, VCFParseError> {
    let desired_field_idx = gt_format_cache.gt_field_idx;
    let mut idx = 0 as usize;
    for gt_item in gt_str.split(":") {
        if idx == desired_field_idx {
            return Ok(gt_item);
        }
        idx += 1;
    }
    Err(VCFParseError::NoGenotypeFormatDefinition(
        gt_str.to_string(),
    ))
}

fn parse_gts(
    gts: std::slice::Iter<&str>,
    gt_format_cache: &mut GtFormatCache,
    line: &String,
) -> Result<Vec<Vec<i16>>, VCFParseError> {
    let mut parsed_gts =
        vec![vec![0; gt_format_cache.ploidy as usize]; gt_format_cache.num_samples];

    let mut sample_idx = 0;
    for gt_str in gts {
        let gt = get_gt_item_from_gt_string(gt_str, gt_format_cache)?;

        let this_ploidy = match parse_gt(gt, sample_idx, &mut parsed_gts, line) {
            Ok(alleles) => alleles,
            Err(e) => return Err(e),
        };

        if gt_format_cache.ploidy != this_ploidy as u8 {
            return Err(VCFParseError::DifferentPloidiesError(line.to_string()));
        }
        sample_idx += 1;
    }
    Ok(parsed_gts)
}

fn parse_variant_line(
    line: String,
    gt_format_cache: &mut GtFormatCache,
) -> Result<Variant, VCFParseError> {
    let fields = line.split("\t").collect::<Vec<&str>>();

    let pos = match fields[1].parse::<u64>() {
        Ok(pos) => pos,
        Err(_) => return Err(VCFParseError::PosNotInt(line)),
    };

    let mut alleles = Vec::new();
    alleles.push(fields[3].to_string());
    alleles.extend(fields[4].split(",").map(|s| s.to_string()));

    let qual: f64;
    if fields[5] == "." {
        qual = 0.0;
    } else {
        qual = match fields[5].parse::<f64>() {
            Ok(pos) => pos,
            Err(_) => return Err(VCFParseError::QualNotFloat(line)),
        }
    };

    let mut filters = Vec::new();
    if fields[6] != "PASS" {
        filters.extend(fields[6].split(";").map(|s| s.to_string()));
    }

    let gt_format_str = fields[8].to_string();
    if gt_format_str != gt_format_cache.gt_string {
        let iter = fields[8]
            .split(":")
            .enumerate()
            .map(|(i, s)| (s.to_string(), i));
        let gt_format_idxs = HashMap::<String, usize>::from_iter(iter);
        gt_format_cache.gt_string = gt_format_str;
        gt_format_cache.gt_format_idxs = gt_format_idxs;
        gt_format_cache.gt_field_idx = match gt_format_cache.gt_format_idxs.get(GT_FIELD_ID) {
            Some(idx) => *idx,
            None => return Err(VCFParseError::GenotypeNotFoundInFormatDefinition(line)),
        };
    }

    if gt_format_cache.ploidy == 0 {
        gt_format_cache.ploidy = match get_ploidy_form_first_gt(&fields[9], gt_format_cache) {
            Ok(ploidy) => ploidy,
            Err(_) => {
                return Err(VCFParseError::FirstGtDoesNotDefinePloidy(
                    fields[9].to_string(),
                    line.to_string(),
                ))
            }
        };
    }

    let gts = parse_gts(fields[9..].iter(), gt_format_cache, &line)?;

    let ploidy = gts[0].len() as u8;

    let var = Variant {
        chrom: fields[0].to_string(),
        pos,
        id: fields[2].to_string(),
        alleles,
        qual,
        filters,
        gts,
        ploidy: ploidy,
    };
    Ok(var)
}

pub struct Variants<'a> {
    pub samples: Vec<String>,
    pub vars_iter: Box<dyn Iterator<Item = Result<Variant, VCFParseError>> + 'a>,
    pub ploidy: u8,
}

fn read_sample_line(line: &str) -> Result<Vec<String>, VCFParseError> {
    if !line.starts_with("#CHROM") {
        return Err(VCFParseError::InvalidSampleLine(line.to_string()));
    }

    let samples = line.split("\t").skip(9).map(|s| s.to_string()).collect();
    Ok(samples)
}

fn parse_vcf_buffer<'a, T: Read + 'a>(
    mut file: BufReader<T>,
) -> Result<Variants<'a>, VCFParseError> {
    let samples;
    loop {
        let mut line = String::new();
        match file.read_line(&mut line) {
            Ok(0) => return Err(VCFParseError::EmptyFile),
            Ok(_) => (),
            Err(_) => return Err(VCFParseError::ReadLineError(0)),
        }
        if line.starts_with("##") {
        } else if line.starts_with("#CHROM") {
            samples = read_sample_line(&line)?;
            break;
        } else {
            return Err(VCFParseError::InvalidSampleLine(line));
        }
    }

    let mut gt_format_cache = GtFormatCache {
        gt_string: "".to_string(),
        gt_format_idxs: HashMap::new(),
        gt_field_idx: 0,
        num_samples: samples.len(),
        ploidy: 0,
    };

    let mut vars_iter = file
        .lines()
        .map(move |line_res| {
            let line = match line_res {
                Ok(line) => line,
                Err(_) => return Err(VCFParseError::ReadLineError(0)),
            };
            parse_variant_line(line, &mut gt_format_cache)
        })
        .peekable();

    let first_var = match vars_iter.peek() {
        Some(Ok(var)) => Variant {
            chrom: var.chrom.clone(),
            pos: var.pos,
            id: var.id.clone(),
            alleles: var.alleles.clone(),
            qual: var.qual,
            filters: var.filters.clone(),
            gts: var.gts.clone(),
            ploidy: var.ploidy,
        },
        Some(Err(_)) => return Err(VCFParseError::NoVariantsError),
        None => return Err(VCFParseError::EmptyFile),
    };

    let vars = Variants {
        samples: samples,
        vars_iter: Box::new(vars_iter),
        ploidy: first_var.ploidy,
    };

    return Ok(vars);
}

#[derive(PartialEq)]
pub enum VcfFileKind {
    PlainTextVcf,
    GzippedVcf,
}

pub fn guess_vcf_file_kind(fpath: &PathBuf) -> Result<VcfFileKind, Box<dyn std::error::Error>> {
    let mut file = match File::open(fpath.clone()) {
        Ok(file) => file,
        Err(e) => return Err(Box::new(e)),
    };
    let mut buffer = vec![0; 2];
    file.read_exact(&mut buffer)?;

    if buffer == [0x23, 0x23] {
        return Ok(VcfFileKind::PlainTextVcf);
    }
    if buffer != [0x1f, 0x8b] {
        return Err(Box::new(VCFParseError::InvalidVCFFile(
            fpath.to_string_lossy().to_string(),
        )));
    }

    let file = File::open(fpath)?;
    let mut file = MultiGzDecoder::new(file);
    file.read_exact(&mut buffer)?;
    if buffer == [0x23, 0x23] {
        return Ok(VcfFileKind::GzippedVcf);
    }
    Err(Box::new(VCFParseError::InvalidGzipVCFFile(
        fpath.to_string_lossy().to_string(),
    )))
}

pub fn read_vcf_file(fpath: &PathBuf) -> Result<Variants, Box<dyn std::error::Error>> {
    let kind = match guess_vcf_file_kind(fpath) {
        Ok(kind) => kind,
        Err(e) => return Err(e),
    };

    let file = File::open(fpath)?;

    if kind == VcfFileKind::PlainTextVcf {
        let file = BufReader::new(file);
        match parse_vcf_buffer(file) {
            Ok(vars) => return Ok(vars),
            Err(e) => return Err(Box::new(e)),
        }
    } else if kind == VcfFileKind::GzippedVcf {
        let file = MultiGzDecoder::new(file);
        let file = BufReader::new(file);
        match parse_vcf_buffer(file) {
            Ok(vars) => return Ok(vars),
            Err(e) => return Err(Box::new(e)),
        }
    }
    Err(Box::new(VCFParseError::InvalidVCFFile(
        fpath.to_string_lossy().to_string(),
    )))
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
20\t1234567\tmicrosat1\tGTC\tG,GTCT\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t0/1:35:4\t0/2:17:2\t1/1:40:3
20\t1234567\tmicrosat1\tGTC\tG,GTCT\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t0/1:35:4\t0/2:17:2\t1/1:40:3";

    #[test]
    fn it_works() {
        let mock_file = BufReader::new(VCF_45.as_bytes());
        let vars = parse_vcf_buffer(mock_file).expect("Error");
        for var_res in vars.vars_iter {
            let _var = var_res.expect("Error reading variant");
        }
    }
}
