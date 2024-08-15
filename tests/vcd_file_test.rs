use nei_rs::{guess_vcf_file_kind, VcfFileKind};
use std::path::Path;

#[test]
fn vcf_file_type() {
    println!("Hello from vcd_file_test.rs");
    let data_dir = Path::new(file!()).parent().unwrap().join("data");
    let vcf_text_fpath = data_dir.join("format_example_4_5.vcf");
    println!("dir: {:?}", vcf_text_fpath);
    let file_type = guess_vcf_file_kind(&vcf_text_fpath).unwrap();
    assert!(file_type == VcfFileKind::PlainTextVcf);

    let vcf_gz_fpath = data_dir.join("format_example_4_5.vcf.gz");
    let file_type = guess_vcf_file_kind(&vcf_gz_fpath).unwrap();
    assert!(file_type == VcfFileKind::GzippedVcf);
}
