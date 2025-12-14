use crate::structtox::ComparativemRNA;
use crate::structtox::MRNA;
use crate::structtox::PathmRNA;
use std::collections::HashSet;
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader};

/*
Gaurav Sablok
codeprog@icloud.com
*/

impl PathmRNA {
    pub fn mrna_analysis(&self) -> Result<String, Box<dyn Error>> {
        let fileopen_1 = File::open(self.pathfile1.clone()).expect("File not present");
        let fileread_1 = BufReader::new(fileopen_1);
        let fileopen_2 = File::open(self.pathfile2.clone()).expect("file not present");
        let fileread_2 = BufReader::new(fileopen_2);
        let mut mrna_ids_1: HashSet<String> = HashSet::new();
        let mut mrnaread_1: Vec<MRNA> = Vec::new();
        let mut mrna_ids_2: HashSet<String> = HashSet::new();
        let mut mrnaread_2: Vec<MRNA> = Vec::new();
        for i in fileread_1.lines() {
            let line = i.expect("file not present");
            if !line.starts_with("#") {
                let linevec = line.split("\t").collect::<Vec<_>>();
                if linevec[2] == "mRNA" {
                    let nameinsert = linevec[8].split(";").collect::<Vec<_>>()[2]
                        .to_string()
                        .replace("Parent=", "");
                    let hashinsert = nameinsert.clone();
                    let startinsert: usize = linevec[3].parse::<usize>().unwrap();
                    let stopinsert: usize = linevec[4].parse::<usize>().unwrap();
                    let strandinsert: String = linevec[6].to_string();
                    mrnaread_1.push(MRNA {
                        name: nameinsert,
                        start: startinsert,
                        stop: stopinsert,
                        strand: strandinsert,
                    });
                    mrna_ids_1.insert(hashinsert);
                }
            }
        }

        for i in fileread_2.lines() {
            let line = i.expect("file not present");
            let linevec = line.split("\t").collect::<Vec<_>>();
            if !line.starts_with("#") {
                if linevec[2] == "mRNA" {
                    let nameinsert = linevec[8].split(";").collect::<Vec<_>>()[2]
                        .to_string()
                        .replace("Parent=", "");
                    let hashinsert = nameinsert.clone();
                    let startinsert: usize = linevec[3].parse::<usize>().unwrap();
                    let stopinsert: usize = linevec[4].parse::<usize>().unwrap();
                    let strandinsert: String = linevec[6].to_string();
                    mrnaread_2.push(MRNA {
                        name: nameinsert,
                        start: startinsert,
                        stop: stopinsert,
                        strand: strandinsert,
                    });
                    mrna_ids_2.insert(hashinsert);
                }
            }
        }

        let mut common_mrnas: Vec<ComparativemRNA> = Vec::new();
        for i in mrnaread_1.iter() {
            for val in mrnaread_2.iter() {
                if i.name == val.name {
                    common_mrnas.push(ComparativemRNA {
                        name1: i.name.clone(),
                        name2: val.name.clone(),
                        start1: i.start,
                        start2: val.start,
                        stop1: i.stop,
                        stop2: val.stop,
                        strand1: i.strand.clone(),
                        strand2: val.strand.clone(),
                    })
                }
            }
        }

        let mut comparative_mrna_positive_strand: Vec<ComparativemRNA> = Vec::new();
        let mut comparative_mrna_negative_strand: Vec<ComparativemRNA> = Vec::new();

        for i in mrnaread_1.iter() {
            for val in mrnaread_2.iter() {
                if i.name == val.name
                    && i.strand == "+".to_string()
                    && val.strand == "+".to_string()
                {
                    comparative_mrna_positive_strand.push(ComparativemRNA {
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
                if i.name == val.name
                    && i.strand == "-".to_string()
                    && val.strand == "-".to_string()
                {
                    comparative_mrna_negative_strand.push(ComparativemRNA {
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

        let mut positive_strand_start_difference: Vec<(String, usize)> = Vec::new();
        for i in comparative_mrna_positive_strand.iter() {
            if i.start1 != i.start2 {
                let valueinsert: (String, usize) = (i.name1.clone(), (i.start1 - i.start2));
                positive_strand_start_difference.push(valueinsert);
            }
        }

        let mut positive_strand_end_difference: Vec<(String, usize)> = Vec::new();
        for i in comparative_mrna_positive_strand.iter() {
            if i.stop1 != i.stop2 {
                let valueinsert: (String, usize) = (i.name1.clone(), (i.stop1 - i.stop2));
                positive_strand_end_difference.push(valueinsert);
            }
        }

        let mut negative_strand_start_difference: Vec<(String, usize)> = Vec::new();
        for i in comparative_mrna_negative_strand.iter() {
            if i.start1 != i.start2 {
                let valueinsert: (String, usize) = (i.name1.clone(), (i.start1 - i.start2));
                negative_strand_start_difference.push(valueinsert);
            }
        }

        let mut negative_strand_end_difference: Vec<(String, usize)> = Vec::new();
        for i in comparative_mrna_negative_strand.iter() {
            if i.stop1 != i.stop2 {
                let valueinsert: (String, usize) = (i.name1.clone(), (i.stop1 - i.stop2));
                negative_strand_end_difference.push(valueinsert);
            }
        }

        let mut filewrite1 =
            File::create("positve_strand_start_difference.txt").expect("file not present");
        writeln!(filewrite1, "{}\t{}", "Gene Name", "Difference").expect("line not present");
        for i in positive_strand_start_difference.iter() {
            writeln!(filewrite1, "{}\t{}", i.0, i.1).expect("line not present");
        }

        let mut filewrite2 =
            File::create("positve_strand_end_difference.txt").expect("file not present");
        writeln!(filewrite2, "{}\t{}", "Gene Name", "Difference").expect("line not present");
        for i in positive_strand_end_difference.iter() {
            writeln!(filewrite2, "{}\t{}", i.0, i.1).expect("line not present");
        }

        let mut filewrite3 =
            File::create("negative_strand_start_difference.txt").expect("file not present");
        writeln!(filewrite3, "{}\t{}", "Gene Name", "Difference").expect("line not present");
        for i in negative_strand_start_difference.iter() {
            writeln!(filewrite3, "{}\t{}", i.0, i.1).expect("line not present");
        }

        let mut filewrite4 =
            File::create("positve_strand_end_difference.txt").expect("file not present");
        writeln!(filewrite4, "{}\t{}", "Gene Name", "Difference").expect("line not present");
        for i in positive_strand_end_difference.iter() {
            writeln!(filewrite4, "{}\t{}", i.0, i.1).expect("line not present");
        }

        let mut positive_strand_start_bins_0_200: Vec<usize> = Vec::new();
        let mut positive_strand_start_bins_200_500: Vec<usize> = Vec::new();
        let mut positive_strand_start_bins_500_800: Vec<usize> = Vec::new();
        let mut positive_strand_start_bins_800_1000: Vec<usize> = Vec::new();
        let mut positive_stand_start_bins_1000_1500: Vec<usize> = Vec::new();
        let mut positive_strand_start_bins_greater_than_1500: Vec<usize> = Vec::new();

        for i in positive_strand_start_difference.iter() {
            if i.1 <= 200usize {
                positive_strand_start_bins_0_200.push(i.1);
            }
            if i.1 > 200usize && i.1 <= 500usize {
                positive_strand_start_bins_200_500.push(i.1)
            }
            if i.1 > 500 && i.1 <= 800usize {
                positive_strand_start_bins_500_800.push(i.1)
            }
            if i.1 > 800usize && i.1 <= 1000usize {
                positive_strand_start_bins_800_1000.push(i.1)
            }
            if i.1 > 1000usize && i.1 <= 1500usize {
                positive_stand_start_bins_1000_1500.push(i.1)
            }
            if i.1 > 1500usize {
                positive_strand_start_bins_greater_than_1500.push(i.1)
            } else {
                continue;
            }
        }

        let mut negative_strand_start_bins_0_200: Vec<usize> = Vec::new();
        let mut negative_strand_start_bins_200_500: Vec<usize> = Vec::new();
        let mut negative_strand_start_bins_500_800: Vec<usize> = Vec::new();
        let mut negative_strand_start_bins_800_1000: Vec<usize> = Vec::new();
        let mut negative_stand_start_bins_1000_1500: Vec<usize> = Vec::new();
        let mut negative_strand_start_bins_greater_than_1500: Vec<usize> = Vec::new();

        for i in negative_strand_start_difference.iter() {
            if i.1 <= 200usize {
                negative_strand_start_bins_0_200.push(i.1);
            }
            if i.1 > 200usize && i.1 <= 500usize {
                negative_strand_start_bins_200_500.push(i.1)
            }
            if i.1 > 500 && i.1 <= 800usize {
                negative_strand_start_bins_500_800.push(i.1)
            }
            if i.1 > 800usize && i.1 <= 1000usize {
                negative_strand_start_bins_800_1000.push(i.1)
            }
            if i.1 > 1000usize && i.1 <= 1500usize {
                negative_stand_start_bins_1000_1500.push(i.1)
            }
            if i.1 > 1500usize {
                negative_strand_start_bins_greater_than_1500.push(i.1)
            } else {
                continue;
            }
        }

        let mut positive_strand_end_bins_0_200: Vec<usize> = Vec::new();
        let mut positive_strand_end_bins_200_500: Vec<usize> = Vec::new();
        let mut positive_strand_end_bins_500_800: Vec<usize> = Vec::new();
        let mut positive_strand_end_bins_800_1000: Vec<usize> = Vec::new();
        let mut positive_stand_end_bins_1000_1500: Vec<usize> = Vec::new();
        let mut positive_strand_end_bins_greater_than_1500: Vec<usize> = Vec::new();

        for i in positive_strand_end_difference.iter() {
            if i.1 <= 200usize {
                positive_strand_end_bins_0_200.push(i.1);
            }
            if i.1 > 200usize && i.1 <= 500usize {
                positive_strand_end_bins_200_500.push(i.1)
            }
            if i.1 > 500 && i.1 <= 800usize {
                positive_strand_end_bins_500_800.push(i.1)
            }
            if i.1 > 800usize && i.1 <= 1000usize {
                positive_strand_end_bins_800_1000.push(i.1)
            }
            if i.1 > 1000usize && i.1 <= 1500usize {
                positive_stand_end_bins_1000_1500.push(i.1)
            }
            if i.1 > 1500usize {
                positive_strand_end_bins_greater_than_1500.push(i.1)
            } else {
                continue;
            }
        }

        let mut negative_strand_end_bins_0_200: Vec<usize> = Vec::new();
        let mut negative_strand_end_bins_200_500: Vec<usize> = Vec::new();
        let mut negative_strand_end_bins_500_800: Vec<usize> = Vec::new();
        let mut negative_strand_end_bins_800_1000: Vec<usize> = Vec::new();
        let mut negative_stand_end_bins_1000_1500: Vec<usize> = Vec::new();
        let mut negative_strand_end_bins_greater_than_1500: Vec<usize> = Vec::new();

        for i in negative_strand_end_difference.iter() {
            if i.1 <= 200usize {
                negative_strand_end_bins_0_200.push(i.1);
            }
            if i.1 > 200usize && i.1 <= 500usize {
                negative_strand_end_bins_200_500.push(i.1)
            }
            if i.1 > 500 && i.1 <= 800usize {
                negative_strand_end_bins_500_800.push(i.1)
            }
            if i.1 > 800usize && i.1 <= 1000usize {
                negative_strand_end_bins_800_1000.push(i.1)
            }
            if i.1 > 1000usize && i.1 <= 1500usize {
                negative_stand_end_bins_1000_1500.push(i.1)
            }
            if i.1 > 1500usize {
                negative_strand_end_bins_greater_than_1500.push(i.1)
            } else {
                continue;
            }
        }

        let mut positive_strand_start_1 =
            File::create("positive_strand_start_0_200.txt").expect("file not present");
        for i in positive_strand_start_bins_0_200.iter() {
            writeln!(positive_strand_start_1, "{}\n", i).expect("line not present");
        }
        let mut positive_strand_start_2 =
            File::create("positive_strand_start_200_500.txt").expect("file not present");
        for i in positive_strand_start_bins_200_500.iter() {
            writeln!(positive_strand_start_2, "{}\n", i).expect("line not present");
        }
        let mut positive_strand_start_3 =
            File::create("positive_strand_start_500_800.txt").expect("file not present");
        for i in positive_strand_start_bins_500_800.iter() {
            writeln!(positive_strand_start_3, "{}\n", i).expect("line not present");
        }
        let mut positive_strand_start_4 =
            File::create("positive_strand_start_800_100.txt").expect("file not present");
        for i in positive_strand_start_bins_800_1000.iter() {
            writeln!(positive_strand_start_4, "{}\n", i).expect("line not present");
        }
        let mut positive_strand_start_5 =
            File::create("positive_strand_start_1000_1500.txt").expect("file not present");
        for i in positive_strand_start_bins_800_1000.iter() {
            writeln!(positive_strand_start_5, "{}\n", i).expect("line not present");
        }
        let mut positive_strand_start_6 =
            File::create("positive_strand_start_1500.txt").expect("file not present");
        for i in positive_strand_start_bins_greater_than_1500.iter() {
            writeln!(positive_strand_start_6, "{}\n", i).expect("line not present");
        }

        let mut positive_strand_end_1 =
            File::create("positive_strand_end_0_200.txt").expect("file not present");
        for i in positive_strand_end_bins_0_200.iter() {
            writeln!(positive_strand_end_1, "{}\n", i).expect("line not present");
        }
        let mut positive_strand_end_2 =
            File::create("positive_strand_end_200_500.txt").expect("file not present");
        for i in positive_strand_end_bins_200_500.iter() {
            writeln!(positive_strand_end_2, "{}\n", i).expect("line not present");
        }
        let mut positive_strand_end_3 =
            File::create("positive_strand_end_500_800.txt").expect("file not present");
        for i in positive_strand_end_bins_500_800.iter() {
            writeln!(positive_strand_end_3, "{}\n", i).expect("line not present");
        }
        let mut positive_strand_end_4 =
            File::create("positive_strand_end_800_100.txt").expect("file not present");
        for i in positive_strand_end_bins_800_1000.iter() {
            writeln!(positive_strand_end_4, "{}\n", i).expect("line not present");
        }
        let mut positive_strand_end_5 =
            File::create("positive_strand_end_1000_1500.txt").expect("file not present");
        for i in positive_strand_end_bins_800_1000.iter() {
            writeln!(positive_strand_end_5, "{}\n", i).expect("line not present");
        }
        let mut positive_strand_end_6 =
            File::create("positive_strand_end_1500.txt").expect("file not present");
        for i in positive_strand_end_bins_greater_than_1500.iter() {
            writeln!(positive_strand_end_6, "{}\n", i).expect("line not present");
        }

        let mut negative_strand_end_1 =
            File::create("negative_strand_end_0_200.txt").expect("file not present");
        for i in negative_strand_end_bins_0_200.iter() {
            writeln!(negative_strand_end_1, "{}\n", i).expect("line not present");
        }
        let mut negative_strand_end_2 =
            File::create("negative_strand_end_200_500.txt").expect("file not present");
        for i in negative_strand_end_bins_200_500.iter() {
            writeln!(negative_strand_end_2, "{}\n", i).expect("line not present");
        }
        let mut negative_strand_end_3 =
            File::create("negative_strand_end_500_800.txt").expect("file not present");
        for i in negative_strand_end_bins_500_800.iter() {
            writeln!(negative_strand_end_3, "{}\n", i).expect("line not present");
        }
        let mut negative_strand_end_4 =
            File::create("negative_strand_end_800_100.txt").expect("file not present");
        for i in negative_strand_end_bins_800_1000.iter() {
            writeln!(negative_strand_end_4, "{}\n", i).expect("line not present");
        }
        let mut negative_strand_end_5 =
            File::create("negative_strand_end_1000_1500.txt").expect("file not present");
        for i in negative_strand_end_bins_800_1000.iter() {
            writeln!(negative_strand_end_5, "{}\n", i).expect("line not present");
        }
        let mut negative_strand_end_6 =
            File::create("negative_strand_end_1500.txt").expect("file not present");
        for i in negative_strand_end_bins_greater_than_1500.iter() {
            writeln!(negative_strand_end_6, "{}\n", i).expect("line not present");
        }

        let mut negative_strand_start_1 =
            File::create("negative_strand_start_0_200.txt").expect("file not present");
        for i in negative_strand_start_bins_0_200.iter() {
            writeln!(negative_strand_start_1, "{}\n", i).expect("line not present");
        }
        let mut negative_strand_start_2 =
            File::create("negative_strand_start_200_500.txt").expect("file not present");
        for i in negative_strand_start_bins_200_500.iter() {
            writeln!(negative_strand_start_2, "{}\n", i).expect("line not present");
        }
        let mut negative_strand_start_3 =
            File::create("negative_strand_start_500_800.txt").expect("file not present");
        for i in negative_strand_start_bins_500_800.iter() {
            writeln!(negative_strand_start_3, "{}\n", i).expect("line not present");
        }
        let mut negative_strand_start_4 =
            File::create("negative_strand_end_800_100.txt").expect("file not present");
        for i in negative_strand_start_bins_800_1000.iter() {
            writeln!(negative_strand_start_4, "{}\n", i).expect("line not present");
        }
        let mut negative_strand_start_5 =
            File::create("negative_strand_start_1000_1500.txt").expect("file not present");
        for i in negative_strand_start_bins_800_1000.iter() {
            writeln!(negative_strand_start_5, "{}\n", i).expect("line not present");
        }
        let mut negative_strand_start_6 =
            File::create("negative_strand_start_1500.txt").expect("file not present");
        for i in negative_strand_start_bins_greater_than_1500.iter() {
            writeln!(negative_strand_start_6, "{}\n", i).expect("line not present");
        }

        Ok("The analysis has been written".to_string())
    }
}
