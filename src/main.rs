mod args;
mod mrna;
mod proteincompare;
mod proteincompareseq;
mod proteinplotter;
mod seqcompare;
mod seqseq;
mod structtox;
mod tensor;
use self::structtox::Extractplot;
use self::structtox::PathFile;
use self::structtox::PathmRNA;
use self::structtox::ProteinCompareExtract;
use self::structtox::ProteinCompareExtractSeq;
use self::structtox::SeqCompare;
use self::structtox::ToxPath;
use crate::args::CommandParse;
use crate::args::Commands;
use clap::Parser;
use figlet_rs::FIGfont;

/*
Gaurav Sablok
codeprog@icloud.com
*/

fn main() {
    let standard_font = FIGfont::standard().unwrap();
    let figure = standard_font.convert("Toxannotator");
    assert!(figure.is_some());
    println!("{}", figure.unwrap());

    let argparse = CommandParse::parse();
    match &argparse.command {
        Commands::ProteinPlotter {
            inputfile1,
            inputfile2,
            threads,
        } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let command = ToxPath {
                    filepath1: inputfile1.clone(),
                    filepath2: inputfile2.clone(),
                };
                let proteinplot = command.proteinplotter().unwrap();
                println!("The command has finished:{:?}", proteinplot);
            });
        }
        Commands::ProteinTensor { inputfile, threads } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let command = PathFile {
                    inputpath: inputfile.clone(),
                };
                let _ = command.tensor();
                let _ = command.padded_tensor();
                println!("The padded tensor have been written:{:?}", command);
            });
        }
        Commands::ExtractSeq {
            annotationfile,
            idsfile,
            threads,
        } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let filepathread = Extractplot {
                    pathfile1: annotationfile.clone(),
                    pathfile2: idsfile.clone(),
                };
                let command = filepathread.extractseq(&filepathread.pathfile2).unwrap();
                println!("The command has finished:{:?}", command);
            });
        }
        Commands::ProteinCompare {
            gfffile1,
            gfffile2,
            threads,
        } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let filepath = ProteinCompareExtract {
                    pathfile1: gfffile1.clone(),
                    pathfile2: gfffile2.clone(),
                };
                let filepathrun = filepath.proteincompare().unwrap();
                println!(
                    "The file for the comparison for the proteome has been written: {:?}",
                    filepathrun
                );
            });
        }
        Commands::ProteinCompareSeqCommand {
            gfffile1,
            gfffile2,
            fastafile_1,
            fastafile_2,
            threads,
        } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let filepath = ProteinCompareExtractSeq {
                    pathfile1: gfffile1.clone(),
                    pathfile2: gfffile2.clone(),
                    fastafile1: fastafile_1.clone(),
                    fastafile2: fastafile_2.clone(),
                };
                let filepathrun = filepath.proteincompareseq().unwrap();
                println!(
                    "The file for the comparison for the proteome has been written: {:?}",
                    filepathrun
                );
            });
        }
        Commands::ComparemRNA {
            gff_1,
            gff_2,
            threads,
        } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let filepath = PathmRNA {
                    pathfile1: gff_1.clone(),
                    pathfile2: gff_2.clone(),
                };
                let command = filepath.mrna_analysis().unwrap();
                print!("The command has finished:{:?}", command);
            });
        }
        Commands::CompareFastaGFF {
            gff_1_input,
            gff_2_input,
            fasta_1_input,
            fasta_2_input,
            threads,
        } => {
            let command = rayon::ThreadPoolBuilder::new()
                .num_threads(threads.parse::<usize>().unwrap())
                .build()
                .unwrap();
            command.install(|| {
                let filepath = SeqCompare {
                    gff_1: gff_1_input.clone(),
                    gff_2: gff_2_input.clone(),
                    fasta_1: fasta_1_input.clone(),
                    fasta_2: fasta_2_input.clone(),
                };
                let newcommand = filepath.seqcompare().unwrap();
                print!("The command has finished:{:?}", newcommand);
            });
        }
    }
}
