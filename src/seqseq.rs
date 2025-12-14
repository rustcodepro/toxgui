use crate::structtox::Extractplot;
use crate::structtox::SeqInfo;
use std::collections::HashSet;
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader};

/*
 Gaurav Sablok
 codeprog@icloud.com
*/

impl Extractplot {
    pub fn seqhash(&self) -> Result<HashSet<String>, Box<dyn Error>> {
        let mut hashid: HashSet<String> = HashSet::new();
        let fileopen = File::open(&self.pathfile1).expect("File not present");
        let fileread = BufReader::new(fileopen);
        for i in fileread.lines() {
            let line = i.expect("line not present");
            if !line.starts_with("#") {
                let linevec = line.split("\t").collect::<Vec<_>>();
                if linevec[2] == "protein_coding_gene" {
                    let name = linevec[8].split(";").collect::<Vec<_>>()[0].replace("ID=", "");
                    hashid.insert(name);
                }
            }
        }
        Ok(hashid)
    }

    pub fn seqhashadd(&self) -> Result<Vec<String>, Box<dyn Error>> {
        let mut seqstore: Vec<String> = Vec::new();
        let fileopen = File::open(&self.pathfile1).expect("File not present");
        let fileread = BufReader::new(fileopen);
        for i in fileread.lines() {
            let line = i.expect("line not present");
            if !line.starts_with("#") {
                seqstore.push(line.clone());
            }
        }
        Ok(seqstore)
    }

    pub fn extractseq(&self, newfile: &str) -> Result<String, Box<dyn Error>> {
        let fileopen = File::open(newfile).expect("file not present");
        let fileread = BufReader::new(fileopen);
        let mut filevec: Vec<_> = Vec::new();
        for i in fileread.lines() {
            let line = i.expect("file not present");
            filevec.push(line);
        }
        let hashidcloned = self.seqhash().unwrap();
        let hashidsseq = self.seqhashadd().unwrap();
        let mut seqvector: Vec<SeqInfo> = Vec::new();
        for i in hashidcloned.iter() {
            let mut exonvec: Vec<(usize, usize)> = Vec::new();
            let mut cdsvec: Vec<(usize, usize)> = Vec::new();
            let mut proteincodingvec: Vec<(usize, usize)> = Vec::new();
            let mut three_prime_utrvec: Vec<(usize, usize)> = Vec::new();
            let mut five_prime_utrvec: Vec<(usize, usize)> = Vec::new();
            for val in hashidcloned.iter() {
                for seq in hashidsseq.iter() {
                    let valinter = seq.split("\t").collect::<Vec<_>>();
                    if valinter[2] == "protein_coding_gene"
                        && *val == valinter[8].split(";").collect::<Vec<_>>()[0].replace("ID=", "")
                    {
                        let value: (usize, usize) = (
                            valinter[3].parse::<usize>().unwrap(),
                            valinter[4].parse::<usize>().unwrap(),
                        );
                        proteincodingvec.push(value);
                    }
                    if valinter[2] == "exon"
                        && *val
                            == valinter[8].split(";").collect::<Vec<_>>()[2].replace("gene_id=", "")
                    {
                        let exonpush: (usize, usize) = (
                            valinter[3].parse::<usize>().unwrap(),
                            valinter[4].parse::<usize>().unwrap(),
                        );
                        exonvec.push(exonpush);
                    }
                    if valinter[2] == "CDS"
                        && *val
                            == valinter[8].split(";").collect::<Vec<_>>()[2].replace("gene_id=", "")
                    {
                        let cdspush: (usize, usize) = (
                            valinter[3].parse::<usize>().unwrap(),
                            valinter[4].parse::<usize>().unwrap(),
                        );
                        cdsvec.push(cdspush);
                    }
                    if valinter[2] == "three_prime_UTR"
                        && *val
                            == valinter[8].split(";").collect::<Vec<_>>()[1]
                                .split("-")
                                .collect::<Vec<_>>()[0]
                                .replace("Parent=", "")
                    {
                        let threeutr: (usize, usize) = (
                            valinter[3].parse::<usize>().unwrap(),
                            valinter[4].parse::<usize>().unwrap(),
                        );
                        three_prime_utrvec.push(threeutr);
                    }
                    if valinter[2] == "five_prime_UTR"
                        && *val
                            == valinter[8].split(";").collect::<Vec<_>>()[1]
                                .split("-")
                                .collect::<Vec<_>>()[0]
                                .replace("Parent=", "")
                    {
                        let fiveutr: (usize, usize) = (
                            valinter[3].parse::<usize>().unwrap(),
                            valinter[4].parse::<usize>().unwrap(),
                        );
                        five_prime_utrvec.push(fiveutr);
                    }
                }
            }
            seqvector.push(SeqInfo {
                name: i.clone(),
                protein_coding: proteincodingvec,
                exon: exonvec,
                cds: cdsvec,
                three_prime: three_prime_utrvec,
                five_prime: five_prime_utrvec,
            });
        }
        let mut selectedones: Vec<SeqInfo> = Vec::new();
        for i in filevec.iter() {
            for val in seqvector.iter() {
                if *i == val.name {
                    selectedones.push(SeqInfo {
                        name: i.clone(),
                        protein_coding: val.protein_coding.clone(),
                        exon: val.exon.clone(),
                        cds: val.cds.clone(),
                        three_prime: val.three_prime.clone(),
                        five_prime: val.five_prime.clone(),
                    });
                }
            }
        }
        let mut filewrite = File::create("extractseq.txt").expect("file not present");
        for i in selectedones.iter() {
            writeln!(
                filewrite,
                "{}\t{:?}\t{:?}\t{:?}\t{:?}\t{:?}",
                i.name, i.protein_coding, i.exon, i.cds, i.three_prime, i.five_prime,
            )
            .expect("file not present");
        }
        Ok("The file has been written".to_string())
    }
}
