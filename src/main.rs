use clap::{Parser, Subcommand};
use helicase::*;
use helicase::input::*;
use xorf::{Filter, BinaryFuse8};
use bitvec::prelude::*;
use serde::{Serialize, Deserialize};
use rmp_serde::Serializer;
use std::fs;


const CONFIG: Config = ParserOptions::default().config();

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Args {
    ///Kmer length. Caution, only k=31 is usable.
    #[arg(short)]
    k: usize,

    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand)]
enum Commands {
    Build {
        #[arg(short = 'f', long = "input")]
        input: String,

        #[arg(short = 'o', long = "output")]
        output: String,
    },
    Query {
        #[arg(short = 'i', long = "index")]
        index: String,

        #[arg(short = 'q', long = "query")]
        query: String,
    },
}

#[derive(Serialize, Deserialize, Debug)]
struct FuzFil {
    x: BinaryFuse8
}

fn main() {
    let args = Args::parse();

    match &args.command {
        Some(Commands::Build { input, output}) => {
            let mut reader = FastaParser::<CONFIG, _>::from_file(&input).expect("Error during fasta reading");

            let mut kmers:Vec<u64> = vec![];
            while let Some(_event) = reader.next(){
                let seq = reader.get_dna_string_owned();
                for i in 0..seq.len()-args.k+1{
                    let kmer = (&seq[i..i+args.k]).to_vec();
                    kmers.push(vectransformer(&kmer));
                }
            }
            let filter = BinaryFuse8::try_from(&kmers).unwrap();
            let packfil = FuzFil{x:filter};
            let mut serialized = Vec::new();
            packfil.serialize(&mut Serializer::new(&mut serialized)).unwrap();
            fs::write(output, serialized).expect("Error during writing");
        }

        Some(Commands::Query { index, query}) => {
            let mut qfile = FastaParser::<CONFIG, _>::from_file(query).expect("Error during fasta reading");
            let mpack = fs::read(index).expect("Error while reading index");
            let deserialized: FuzFil = rmp_serde::from_slice(&mpack).unwrap();
            let kmers: BinaryFuse8 = deserialized.x;

            let mut count_pos = 0;
            let mut count_neg = 0;
            while let Some(_event) = qfile.next(){
                let seq = qfile.get_dna_string_owned();
                for i in 0..seq.len()-args.k+1{
                    let kmer = (&seq[i..i+args.k]).to_vec();
                    let rc = revcomp(&kmer);
                    if kmers.contains(&vectransformer(&kmer)) || kmers.contains(&vectransformer(&rc)) {
                        count_pos += 1;
                    }
                    else {
                        count_neg += 1;
                    }
                }
            }
            println!("{:?} positive kmers", count_pos);
            println!("{:?} negative kmers", count_neg);
        }

        None => {
            println!("Please tell me what to do!")
        }
    }
}

fn vectransformer(kmer:&Vec<u8>) -> u64{
    let mut v = bitvec![];
    for nuc in kmer{
        let b = nuc.view_bits::<Msb0>();
        v.push(b[5]);
        v.push(b[6]);
    }
    v.push(false);
    v.push(false);
    return v.as_raw_slice()[0].try_into().unwrap();
}

fn revcomp(kmer:&Vec<u8>)->Vec<u8> {
    let mut rc:Vec<u8> = Vec::new();
    for nuc in kmer{
        match nuc{
            65 => rc.push(84),
            84 => rc.push(65),
            67 => rc.push(71),
            71 => rc.push(67),
            _ => println!("Wrong character in a query kmer, correct it and start again"),
        }
    }
    rc.reverse();
    return rc;
}
