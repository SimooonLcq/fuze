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
    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand)]
enum Commands {
    Build {
        ///Kmer length
        #[arg(short)]
        k: usize,

        ///Input fasta file
        #[arg(short = 'f', long = "input")]
        input: String,

        ///Output index file
        #[arg(short = 'o', long = "output")]
        output: String,
    },
    Query {
        ///Kmer length
        #[arg(short)]
        k: usize,

        ///z (=k-s, with k the lenght of query k-mers and s the length of indexed s-mers)
        #[arg(short)]
        z: usize,

        ///Input index file
        #[arg(short = 'i', long = "index")]
        index: String,

        ///Input query fasta file
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
        Some(Commands::Build { k, input, output}) => {
            assert!(*k <= 32, "k must be inferior or equal to 32!");
            let mut reader = FastaParser::<CONFIG, _>::from_file(&input).expect("Error during fasta reading");

            let mut kmers:Vec<u64> = vec![];
            while let Some(_event) = reader.next(){
                let seq = reader.get_dna_string_owned();
                for i in 0..seq.len()-k+1{
                    let kmer = (&seq[i..i+k]).to_vec();
                    kmers.push(vectransformer(&canonical(&kmer), k));
                }
            }
            kmers.sort();
            kmers.dedup();

            let filter = BinaryFuse8::try_from(&kmers).unwrap();
            let packfil = FuzFil{x:filter};
            let mut serialized = Vec::new();
            packfil.serialize(&mut Serializer::new(&mut serialized)).unwrap();
            fs::write(output, serialized).expect("Error during writing");
        }

        Some(Commands::Query { k, z, index, query}) => {
            assert!((k-z) <= 32, "s must be inferior or equal to 32!");
            let mut qfile = FastaParser::<CONFIG, _>::from_file(query).expect("Error during fasta reading");
            let mpack = fs::read(index).expect("Error while reading index");
            let deserialized: FuzFil = rmp_serde::from_slice(&mpack).unwrap();
            let kmers: BinaryFuse8 = deserialized.x;

            let mut count_pos = 0;
            let mut count_neg = 0;
            while let Some(_event) = qfile.next(){
                let seq = qfile.get_dna_string_owned();
                for i in 0..seq.len()-k+1{
                    let kmer = (&seq[i..i+k]).to_vec();
                    let smerx: Vec<Vec<u8>> = smers(&kmer, z, k);
                    
                    if start_checking(&kmers, &smerx) {
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

fn canonical(kmer:&Vec<u8>) -> Vec<u8>{
    let rc = revcomp(kmer);
    if &rc >= kmer{
        return kmer.to_vec();
    }
    else{
        return rc;
    }
}

fn vectransformer(kmer:&Vec<u8>, k:&usize) -> u64{
    let mut v = bitvec![];
    for nuc in kmer{
        let b = nuc.view_bits::<Msb0>();
        v.push(b[5]);
        v.push(b[6]);
    }
    for _i in 0..((32-k)*2){
        v.push(false);
    }
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

fn smers(kmer:&Vec<u8>, z:&usize, k:&usize) ->Vec<Vec<u8>> {
    let mut smers:Vec<Vec<u8>> = vec![];
    for i in 0..*z+1 {
        smers.push(kmer[i..i+k-z].to_vec());
    }
    return smers;
}

fn start_checking(filter:&BinaryFuse8, smers:&Vec<Vec<u8>>) -> bool{
    let mut i = 0;
    let mut positive = true;
    while (i < smers.len()) && (positive){
        positive = filter.contains(&vectransformer(&canonical(&smers[i]), &(smers[i].len())));
        i += 1;
    }
    return positive;
}
/*
fn smers_to_str(smers:&Vec<Vec<u8>>) -> Vec<String>{
    let mut vs:Vec<String> = vec![];
    for s in smers{
        vs.push(String::from_utf8(s.to_vec()).unwrap())
    }
    return vs;
}
    */