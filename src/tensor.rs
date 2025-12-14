use crate::structtox::FastaStruct;
use crate::structtox::PathFile;
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

impl PathFile {
    pub fn tensor(&self) -> Result<(), Box<dyn Error>> {
        let fastafile = read_fasta(&self.inputpath).expect("file not present");
        let mut inputtensor: Vec<String> = Vec::new();
        for (_val, seq) in fastafile.iter() {
            inputtensor.push(seq.seq.clone());
        }
        let mut tensor: Vec<Vec<Vec<f32>>> = Vec::new();
        for i in inputtensor.iter() {
            let pushtensor: Vec<Vec<f32>> = proteinencoder(i);
            tensor.push(pushtensor);
        }
        let mut filetensor = File::open("tensorprotein.txt").expect("file not present");
        writeln!(filetensor, "{:?}", tensor).expect("line not present");
        Ok(())
    }

    pub fn padded_tensor(&self) -> Result<(), Box<dyn Error>> {
        let fastafile = read_fasta(&self.inputpath).expect("file not present");
        let mut inputtensor: Vec<String> = Vec::new();
        for (_val, seq) in fastafile.iter() {
            inputtensor.push(seq.seq.clone());
        }
        let lengthtensor: Vec<usize> = inputtensor.iter().map(|x| x.len()).collect::<Vec<_>>();
        let maxlength = lengthtensor.iter().max().unwrap();
        let mut tensor: Vec<Vec<Vec<f32>>> = Vec::new();
        for i in inputtensor.iter() {
            let pushtensor: Vec<Vec<f32>> = proteinencoder(i);
            tensor.push(pushtensor);
        }
        let mut lengthtensor = Vec::new();
        for i in tensor.iter() {
            if i.len() == *maxlength {
                lengthtensor.push(i.clone());
            } else if i.len() != *maxlength {
                let mut newarray = i.clone();
                let valueinput = *maxlength as f32 - i.len() as f32;
                let iconst = 0f32;
                while iconst < valueinput {
                    newarray.push(
                        [
                            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                        ]
                        .to_vec(),
                    );
                }
            }
        }
        let mut filetensor = File::open("tensorprotein_padded.txt").expect("file not present");
        writeln!(filetensor, "{:?}", tensor).expect("line not present");
        Ok(())
    }
}

pub fn proteinencoder(seq: &str) -> Vec<Vec<f32>> {
    let value = seq.to_ascii_uppercase().chars().collect::<Vec<_>>();
    let mut vectorchar: Vec<Vec<f32>> = Vec::new();
    for i in value.iter() {
        if *i == 'A' {
            vectorchar.push(
                [
                    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ]
                .to_vec(),
            )
        }
        if *i == 'R' {
            vectorchar.push(
                [
                    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ]
                .to_vec(),
            )
        }
        if *i == 'N' {
            vectorchar.push(
                [
                    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ]
                .to_vec(),
            )
        }
        if *i == 'D' {
            vectorchar.push(
                [
                    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ]
                .to_vec(),
            )
        }
        if *i == 'C' {
            vectorchar.push(
                [
                    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ]
                .to_vec(),
            )
        }
        if *i == 'E' {
            vectorchar.push(
                [
                    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ]
                .to_vec(),
            )
        }
        if *i == 'Q' {
            vectorchar.push(
                [
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ]
                .to_vec(),
            )
        }
        if *i == 'G' {
            vectorchar.push(
                [
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ]
                .to_vec(),
            )
        }
        if *i == 'H' {
            vectorchar.push(
                [
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ]
                .to_vec(),
            )
        }
        if *i == 'I' {
            vectorchar.push(
                [
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ]
                .to_vec(),
            )
        }
        if *i == 'L' {
            vectorchar.push(
                [
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ]
                .to_vec(),
            )
        }
        if *i == 'K' {
            vectorchar.push(
                [
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ]
                .to_vec(),
            )
        }
        if *i == 'M' {
            vectorchar.push(
                [
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ]
                .to_vec(),
            )
        }
        if *i == 'F' {
            vectorchar.push(
                [
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ]
                .to_vec(),
            )
        }
        if *i == 'P' {
            vectorchar.push(
                [
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ]
                .to_vec(),
            )
        }
        if *i == 'S' {
            vectorchar.push(
                [
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ]
                .to_vec(),
            )
        }
        if *i == 'T' {
            vectorchar.push(
                [
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ]
                .to_vec(),
            )
        }
        if *i == 'W' {
            vectorchar.push(
                [
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
                ]
                .to_vec(),
            )
        }
        if *i == 'Y' {
            vectorchar.push(
                [
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
                ]
                .to_vec(),
            )
        }
        if *i == 'V' {
            vectorchar.push(
                [
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                ]
                .to_vec(),
            )
        }
    }
    vectorchar
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
                    current_id.clone().clone().split("|").collect::<Vec<_>>()[2]
                        .replace("gene=", "")
                        .replace(" ", "")
                        .to_string(),
                    FastaStruct {
                        id: current_id.clone().split("|").collect::<Vec<_>>()[2]
                            .replace("gene=", "")
                            .replace(" ", "")
                            .to_string(),
                        seq: current_sequence.clone(),
                        tag: current_id.clone().split("|").collect::<Vec<_>>()[2]
                            .replace("gene=", "")
                            .replace(" ", "")
                            .to_string(),
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
            current_id.clone().split("|").collect::<Vec<_>>()[2]
                .replace("gene=", "")
                .replace(" ", "")
                .to_string(),
            FastaStruct {
                id: current_id.clone().split("|").collect::<Vec<_>>()[2]
                    .replace("gene=", "")
                    .replace(" ", "")
                    .to_string(),
                seq: current_sequence,
                tag: current_id.clone().split("|").collect::<Vec<_>>()[2]
                    .replace("gene=", "")
                    .replace(" ", "")
                    .to_string(),
            },
        );
    }

    Ok(records)
}
