use crate::structtox::CommonSeq;
use crate::structtox::FastaStruct;
use crate::structtox::ProteinEqual;
use crate::structtox::SeqCompare;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader};
use std::path::Path;

/*
 Gaurav Sablok
 codeprog@icloud.com
*/

impl SeqCompare {
    pub fn seqcompare(&self) -> Result<String, Box<dyn Error>> {
        let file1fasta = read_fasta(self.fasta_1.clone()).unwrap();
        let file2fasta = read_fasta(self.fasta_2.clone()).unwrap();
        let gff_1 = BufReader::new(File::open(self.gff_1.clone()).expect("file not present"));
        let gff_2 = BufReader::new(File::open(self.gff_2.clone()).expect("file not present"));
        let mut protein_1: Vec<ProteinEqual> = Vec::new();
        let mut protein_2: Vec<ProteinEqual> = Vec::new();
        for i in gff_1.lines() {
            let line = i.expect("file not present");
            if !line.starts_with("#") {
                let linevec = line.split("\t").collect::<Vec<_>>();
                if linevec[2] == "protein_coding_gene" {
                    protein_1.push(ProteinEqual {
                        name: linevec[8].split(";").collect::<Vec<_>>()[0].replace("ID=", ""),
                        start: linevec[3].parse::<usize>().unwrap(),
                        stop: linevec[4].parse::<usize>().unwrap(),
                        strand: linevec[6].to_string(),
                    })
                }
            }
        }
        for i in gff_2.lines() {
            let line = i.expect("file not present");
            if !line.starts_with("#") {
                let linevec = line.split("\t").collect::<Vec<_>>();
                if linevec[2] == "protein_coding_gene" {
                    protein_2.push(ProteinEqual {
                        name: linevec[8].split(";").collect::<Vec<_>>()[0].replace("ID=", ""),
                        start: linevec[3].parse::<usize>().unwrap(),
                        stop: linevec[4].parse::<usize>().unwrap(),
                        strand: linevec[6].to_string(),
                    })
                }
            }
        }

        let mut commonseq_1: Vec<(String, String)> = Vec::new();
        let mut commonseq_2: Vec<(String, String)> = Vec::new();
        for i in protein_1.iter() {
            for (_val, keys) in file1fasta.iter() {
                if i.name == keys.id {
                    let vecinsert: (String, String) = (i.name.clone(), keys.seq.clone());
                    commonseq_1.push(vecinsert);
                }
            }
        }
        for i in protein_2.iter() {
            for (_val, keys) in file2fasta.iter() {
                if i.name == keys.id {
                    let vecinsert: (String, String) = (i.name.clone(), keys.seq.clone());
                    commonseq_2.push(vecinsert);
                }
            }
        }

        let mut commonseq_sameid_same_seq: Vec<CommonSeq> = Vec::new();
        let mut commonseq_sameid_different_seq: Vec<CommonSeq> = Vec::new();
        let mut differentid_sameseq: Vec<CommonSeq> = Vec::new();

        for i in commonseq_1.iter() {
            for val in commonseq_2.iter() {
                if i.0.clone() == val.0.clone() && i.1.clone() == val.1.clone() {
                    commonseq_sameid_same_seq.push(CommonSeq {
                        name1: i.0.clone(),
                        name2: val.0.clone(),
                        seq1: i.1.clone(),
                        seq2: val.1.clone(),
                    })
                }
                if i.0.clone() == val.0.clone() && i.1.clone() != val.1.clone() {
                    commonseq_sameid_different_seq.push(CommonSeq {
                        name1: i.0.clone(),
                        name2: val.0.clone(),
                        seq1: i.1.clone(),
                        seq2: val.1.clone(),
                    })
                }
                if i.0.clone() != val.0.clone() && i.1.clone() == val.1.clone() {
                    differentid_sameseq.push(CommonSeq {
                        name1: i.0.clone(),
                        name2: val.0.clone(),
                        seq1: i.1.clone(),
                        seq2: val.1.clone(),
                    })
                }
            }
        }

        let mut commonseq_sameid_same_seq_write =
            File::create("commonseq_sameid_same_seq.fasta").expect("file not present");
        let mut commonseq_sameid_different_seq_write =
            File::create("commonseq_sameid_different_seq.fasta").expect("file not present");
        let mut differentid_sameseq_write =
            File::create("differentid_sameseq.fasta").expect("file not present");
        for i in commonseq_sameid_same_seq.iter() {
            writeln!(
                commonseq_sameid_same_seq_write,
                ">{}-1\n{}\n>{}-2\n{}",
                i.name1, i.name2, i.seq1, i.seq2
            )
            .expect("file not present");
        }
        for i in commonseq_sameid_different_seq.iter() {
            writeln!(
                commonseq_sameid_different_seq_write,
                ">{}-1\n{}\n>{}-2\n{}",
                i.name1, i.name2, i.seq1, i.seq2
            )
            .expect("file not present");
        }
        for i in differentid_sameseq.iter() {
            writeln!(
                differentid_sameseq_write,
                ">{}-1\n{}\n{}-2\n{}",
                i.name1, i.name2, i.seq1, i.seq2
            )
            .expect("file not present");
        }

        Ok("The file has been written".to_string())
    }
}

pub fn read_fasta<P: AsRef<Path>>(path: P) -> std::io::Result<HashMap<String, FastaStruct>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut records = HashMap::new();
    let mut current_id = String::new();
    let mut current_sequence = String::new();
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if !current_id.is_empty() {
                records.insert(
                    current_id.clone().replace(">", "").to_string(),
                    FastaStruct {
                        id: current_id.clone().replace(">", "").to_string(),
                        seq: current_sequence.clone(),
                        tag: current_id.clone().replace(">", "").to_string(),
                    },
                );
                current_sequence.clear();
            }
            current_id = line[1..].to_string();
        } else {
            current_sequence.push_str(&line);
        }
    }

    if !current_id.is_empty() {
        records.insert(
            current_id.clone().replace(">", "").to_string(),
            FastaStruct {
                id: current_id.clone().replace(">", "").to_string(),
                seq: current_sequence,
                tag: current_id.clone().replace(">", "").to_string(),
            },
        );
    }

    Ok(records)
}
