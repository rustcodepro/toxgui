use crate::structtox::GenomeTable;
use crate::structtox::ProteinCompareExtract;
use crate::structtox::ProteinEqual;
use crate::structtox::ProteomeRest;
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader};

/*
Gaurav Sablok
codeprog@icloud.com
*/

impl ProteinCompareExtract {
    pub fn proteincompare(&self) -> Result<String, Box<dyn Error>> {
        let file1open = File::open(self.pathfile1.clone()).expect("file not present");
        let file2open = File::open(self.pathfile2.clone()).expect("file not present");
        let file1read = BufReader::new(file1open);
        let file2read = BufReader::new(file2open);
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

        let mut commongenes: Vec<ProteomeRest> = Vec::new();
        let mut dissimilargenes: Vec<ProteomeRest> = Vec::new();
        for i in filecompare_1.iter() {
            for val in filecompare_2.iter() {
                if i.name == val.name {
                    commongenes.push(ProteomeRest {
                        name1: i.name.clone(),
                        name2: val.name.clone(),
                        start1: i.start,
                        start2: val.start,
                        stop1: i.stop,
                        stop2: val.stop,
                        strand1: i.strand.clone(),
                        strand2: val.strand.clone(),
                    });
                }
            }
        }
        for i in filecompare_1.iter() {
            for val in filecompare_2.iter() {
                if i.name != val.name {
                    dissimilargenes.push(ProteomeRest {
                        name1: i.name.clone(),
                        name2: val.name.clone(),
                        start1: i.start,
                        start2: val.start,
                        stop1: i.stop,
                        stop2: val.stop,
                        strand1: i.strand.clone(),
                        strand2: val.strand.clone(),
                    });
                }
            }
        }

        let mut common_genometablevec: Vec<GenomeTable> = Vec::new();
        for i in commongenes.iter() {
            common_genometablevec.push(GenomeTable {
                name: i.name1.clone(),
                start1: i.start1,
                start2: i.start2,
                end1: i.stop1,
                end2: i.stop2,
                strand1: i.strand1.clone(),
                strand2: i.strand2.clone(),
                startdifference: if i.start2 >= i.start1 {
                    (i.start2 - i.start1).to_string()
                } else {
                    (i.start1 - i.start2).to_string()
                },
                enddifference: if i.stop1 >= i.stop2 {
                    (i.stop1 - i.stop2).to_string()
                } else {
                    (i.stop2 - i.stop1).to_string()
                },
            })
        }

        let mut commongenes_same_strand: Vec<ProteomeRest> = Vec::new();
        let mut commonggenes_difference_strand: Vec<ProteomeRest> = Vec::new();
        for i in filecompare_1.iter() {
            for val in filecompare_2.iter() {
                if i.name == val.name && i.strand == val.strand {
                    commongenes_same_strand.push(ProteomeRest {
                        name1: i.name.clone(),
                        name2: val.name.clone(),
                        start1: i.start,
                        start2: val.start,
                        stop1: i.stop,
                        stop2: val.stop,
                        strand1: i.strand.clone(),
                        strand2: val.strand.clone(),
                    });
                }
            }
        }
        for i in filecompare_1.iter() {
            for val in filecompare_2.iter() {
                if i.name == val.name && i.strand != val.strand {
                    commonggenes_difference_strand.push(ProteomeRest {
                        name1: i.name.clone(),
                        name2: val.name.clone(),
                        start1: i.start,
                        start2: val.start,
                        stop1: i.stop,
                        stop2: val.stop,
                        strand1: i.strand.clone(),
                        strand2: val.strand.clone(),
                    });
                }
            }
        }

        let mut filewrite =
            File::create("comparison-result-common-genes.txt").expect("file not present");
        writeln!(
            filewrite,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            "GeneName",
            "Start1",
            "Start2",
            "End1",
            "End2",
            "Strand1",
            "Strand2",
            "Start Difference",
            "End Difference"
        )
        .expect("line not present");
        for i in common_genometablevec.iter() {
            writeln!(
                filewrite,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                i.name,
                i.start1,
                i.start2,
                i.end1,
                i.end2,
                i.strand1,
                i.strand2,
                i.startdifference,
                i.enddifference
            )
            .expect("line not present");
        }
        let mut dissimilar_genes = File::create("dissimilar-genes.txt").expect("file not present");
        for i in dissimilargenes.iter() {
            writeln!(
                dissimilar_genes,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                i.name1, i.name2, i.start1, i.start2, i.stop1, i.stop2, i.strand1, i.strand2
            )
            .expect("file not present");
        }

        let mut commongene_same_strand_write =
            File::create("commonggenes-same-strand.txt").expect("file not present");
        let mut commongene_different_strand_write =
            File::create("commongenes-different-strand.txt").expect("file not present");
        for i in commongenes_same_strand.iter() {
            writeln!(
                commongene_same_strand_write,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                i.name1.clone(),
                i.name2.clone(),
                i.start1,
                i.start2,
                i.stop1,
                i.stop2,
                i.strand1,
                i.strand2
            )
            .expect("line not present");
        }

        for i in commonggenes_difference_strand.iter() {
            writeln!(
                commongene_different_strand_write,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                i.name1.clone(),
                i.name2.clone(),
                i.start1,
                i.start2,
                i.stop1,
                i.stop2,
                i.strand1,
                i.strand2
            )
            .expect("file not present");
        }
        Ok("The comparative file has been written".to_string())
    }
}
