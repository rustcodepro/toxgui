use crate::structtox::ProteinCompareExtractSeq;
use crate::structtox::ProteinEqual;
use crate::structtox::ProteinEqualSeq;
use crate::structtox::ProteinEqualSeqCompare;
use crate::tensor::read_fasta;
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader};

/*
Gaurav Sablok
codeprog@icloud.com
*/

impl ProteinCompareExtractSeq {
    pub fn proteincompareseq(&self) -> Result<String, Box<dyn Error>> {
        let file1open = File::open(self.pathfile1.clone()).expect("file not present");
        let file2open = File::open(self.pathfile2.clone()).expect("file not present");
        let file1read = BufReader::new(file1open);
        let file2read = BufReader::new(file2open);
        let file1fasta = read_fasta(self.fastafile1.clone()).unwrap();
        let file2fasta = read_fasta(self.fastafile2.clone()).unwrap();
        let mut filecompare_1: Vec<ProteinEqual> = Vec::new();
        let mut filecompare_2: Vec<ProteinEqual> = Vec::new();
        for i in file1read.lines() {
            let line = i.expect("file not present");
            if !line.starts_with("#") {
                let linevec = line.split("\t").collect::<Vec<_>>();
                if linevec[2] == "protein_coding_gene" {
                    filecompare_1.push(ProteinEqual {
                        name: linevec[8].split(";").collect::<Vec<_>>()[0].replace("ID=", ""),
                        start: linevec[3].parse::<usize>().unwrap(),
                        stop: linevec[4].parse::<usize>().unwrap(),
                        strand: linevec[6].to_string(),
                    })
                }
            }
        }

        let mut filecompareseq_1: Vec<ProteinEqualSeq> = Vec::new();
        let mut filecompareseq_2: Vec<ProteinEqualSeq> = Vec::new();
        for i in file2read.lines() {
            let line = i.expect("file not present");
            if !line.starts_with("#") {
                let linevec = line.split("\t").collect::<Vec<_>>();
                if linevec[2] == "protein_coding_gene" {
                    filecompare_2.push(ProteinEqual {
                        name: linevec[8].split(";").collect::<Vec<_>>()[0].replace("ID=", ""),
                        start: linevec[3].parse::<usize>().unwrap(),
                        stop: linevec[4].parse::<usize>().unwrap(),
                        strand: linevec[6].to_string(),
                    });
                }
            }
        }

        for i in filecompare_1.iter() {
            for (val, keys) in file1fasta.iter() {
                if *i.name == val.to_string() {
                    filecompareseq_1.push(ProteinEqualSeq {
                        name: i.name.clone(),
                        start: i.start,
                        stop: i.stop,
                        strand: i.strand.clone(),
                        seq: keys.seq.clone(),
                    })
                }
            }
        }

        for i in filecompare_2.iter() {
            for (val, keys) in file2fasta.iter() {
                if *i.name == val.to_string() {
                    filecompareseq_2.push(ProteinEqualSeq {
                        name: i.name.clone(),
                        start: i.start,
                        stop: i.stop,
                        strand: i.strand.clone(),
                        seq: keys.seq.clone(),
                    });
                }
            }
        }

        let mut compareevery: Vec<ProteinEqualSeqCompare> = Vec::new();
        for i in filecompareseq_1.iter() {
            for val in filecompareseq_2.iter() {
                if i.name == val.name
                    && i.strand == val.strand
                    && i.start == val.start
                    && i.stop == val.stop
                    && i.seq == val.seq
                {
                    compareevery.push(ProteinEqualSeqCompare {
                        name1: i.name.clone(),
                        name2: val.name.clone(),
                        start1: i.start,
                        start2: val.start,
                        stop1: i.stop,
                        stop2: val.stop,
                        strand1: i.strand.parse::<usize>().unwrap(),
                        strand2: val.strand.parse::<usize>().unwrap(),
                        seq1: i.seq.clone(),
                        seq2: val.seq.clone(),
                    })
                }
            }
        }

        let mut compareevery_except_seq: Vec<ProteinEqualSeqCompare> = Vec::new();
        for i in filecompareseq_1.iter() {
            for val in filecompareseq_2.iter() {
                if i.name == val.name
                    && i.strand == val.strand
                    && i.start == val.start
                    && i.stop == val.stop
                    && i.seq != val.seq
                {
                    compareevery_except_seq.push(ProteinEqualSeqCompare {
                        name1: i.name.clone(),
                        name2: val.name.clone(),
                        start1: i.start,
                        start2: val.start,
                        stop1: i.stop,
                        stop2: val.stop,
                        strand1: i.strand.parse::<usize>().unwrap(),
                        strand2: val.strand.parse::<usize>().unwrap(),
                        seq1: i.seq.clone(),
                        seq2: val.seq.clone(),
                    })
                }
            }
        }

        let mut compareevery_strand: Vec<ProteinEqualSeqCompare> = Vec::new();
        for i in filecompareseq_1.iter() {
            for val in filecompareseq_2.iter() {
                if i.name == val.name
                    && i.strand != val.strand
                    && i.start == val.start
                    && i.stop == val.stop
                    && i.seq == val.seq
                {
                    compareevery_strand.push(ProteinEqualSeqCompare {
                        name1: i.name.clone(),
                        name2: val.name.clone(),
                        start1: i.start,
                        start2: val.start,
                        stop1: i.stop,
                        stop2: val.stop,
                        strand1: i.strand.parse::<usize>().unwrap(),
                        strand2: val.strand.parse::<usize>().unwrap(),
                        seq1: i.seq.clone(),
                        seq2: val.seq.clone(),
                    })
                }
            }
        }

        let mut compareevery_strand_seq: Vec<ProteinEqualSeqCompare> = Vec::new();
        for i in filecompareseq_1.iter() {
            for val in filecompareseq_2.iter() {
                if i.name == val.name
                    && i.strand != val.strand
                    && i.start == val.start
                    && i.stop == val.stop
                    && i.seq != val.seq
                {
                    compareevery_strand_seq.push(ProteinEqualSeqCompare {
                        name1: i.name.clone(),
                        name2: val.name.clone(),
                        start1: i.start,
                        start2: val.start,
                        stop1: i.stop,
                        stop2: val.stop,
                        strand1: i.strand.parse::<usize>().unwrap(),
                        strand2: val.strand.parse::<usize>().unwrap(),
                        seq1: i.seq.clone(),
                        seq2: val.seq.clone(),
                    })
                }
            }
        }

        let mut compareevery_write = File::create("compareevery.txt").expect("file not present");
        for i in compareevery.iter() {
            writeln!(
                compareevery_write,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                i.name1.clone(),
                i.name2.clone(),
                i.start1,
                i.start2,
                i.stop1,
                i.stop2,
                i.strand1,
                i.strand2,
                i.seq1,
                i.seq2
            )
            .expect("file not present");
        }

        let mut compareevery_except_seq_write =
            File::create("compareevery.txt").expect("file not present");
        for i in compareevery_except_seq.iter() {
            writeln!(
                compareevery_except_seq_write,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                i.name1.clone(),
                i.name2.clone(),
                i.start1,
                i.start2,
                i.stop1,
                i.stop2,
                i.strand1,
                i.strand2,
                i.seq1,
                i.seq2
            )
            .expect("file not present");
        }

        let mut compareevery_strand_write =
            File::create("compareevery.txt").expect("file not present");
        for i in compareevery_strand.iter() {
            writeln!(
                compareevery_strand_write,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                i.name1.clone(),
                i.name2.clone(),
                i.start1,
                i.start2,
                i.stop1,
                i.stop2,
                i.strand1,
                i.strand2,
                i.seq1,
                i.seq2
            )
            .expect("file not present");
        }

        let mut compareevery_strand_seq_write =
            File::create("compareevery.txt").expect("file not present");
        for i in compareevery_strand_seq.iter() {
            writeln!(
                compareevery_strand_seq_write,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                i.name1.clone(),
                i.name2.clone(),
                i.start1,
                i.start2,
                i.stop1,
                i.stop2,
                i.strand1,
                i.strand2,
                i.seq1,
                i.seq2
            )
            .expect("file not present");
        }

        Ok("The comparative file has been written".to_string())
    }
}
