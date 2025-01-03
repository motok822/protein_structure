use plotters::prelude::*;
use rand::Rng;
use std::collections::HashMap;
// mod aco;
mod animation;
mod beam;
mod lib;
// use aco::ACO;
use animation::animation;
use beam::next_permutation;
use beam::Beam;
use clap::Parser;
use lib::{rotate_left, rotate_right, Amino, AminoAcid, Direction, Protein};
use std::cmp::Ordering;
use std::hash::{Hash, Hasher};

#[derive(Debug, PartialEq, PartialOrd)]
struct OrdF32(f32);

impl Eq for OrdF32 {}

impl Ord for OrdF32 {
    fn cmp(&self, other: &Self) -> Ordering {
        // NaNの扱いに注意しつつ、PartialOrdの結果を利用
        self.partial_cmp(other).unwrap_or(Ordering::Equal)
    }
}
#[derive(Debug, PartialEq, PartialOrd)]
struct FloatKey(f32);

impl Eq for FloatKey {}

impl Hash for FloatKey {
    fn hash<H: Hasher>(&self, state: &mut H) {
        // ビットパターンをハッシュ値に変換
        self.0.to_bits().hash(state);
    }
}
static PROTEIN_DATA: [&str; 22] = [
    "H4",
    "(HP)2PH2PHP2HPH2P2HPH",
    "H2(P2H)7H",
    "P2HP2(H2P4)3H2",
    "P3H2P2H2P5H7P2H2P4H2P2HP2",
    "P2H(P2H2)2P5H10P6(H2P2)2HP2H5",
    "H2(PH)3PH4PH(P3H)2P4H(P3H)2PHPH4(HP)3H2",
    "P2H3PH8P3H10PHP3H12P4H6PH2PHP",
    "H12(PH)2(P2H2)2P2HP2H2PPH2P2HP2(H2P2)2(HP)2H12",
    "H4P4H12P6(H12P3)3HP2(H2P2)2HPH",
    "P3H2P2H4P2H3(PH2)2PH4P8H6P2H6P9HPH2PH11P2H3PH2PHP2HPH3P6H3",
    "P6HPH2P5H3PH5PH2P4H2P2H2PH5PH10PH2PH7p11H7P2HPH3P6HPH2",
    "HPH2P2H4PH3P2H2P2HPH2PHPH2P2H2P3HP8H2",
    "H4PH2PH5P2HP2H2P2HP6HP2HP3HP2H2P2H3PH",
    "PHPH2PH6P2HPHP2HPH2(PH)2P3H(P2H2)2P2HPHP2HP",
    "PHPH2P2HPH3P2H2PH2P3H5P2HPH2(PH)2P4HP2(HP)2",
    "P2HP3HPH4P2H4PH2PH3P2(HP)2HP2HP6H2PH2PH",
    "H3P3H2PH(PH2)3PHP7HPHP2HP3HP2H6PH",
    "PHP4HPH3PHPH4PH2PH2P3HPHP3H3(P2H2)2P3H",
    "PH2PH3PH4P2H3P6HPH2P2H2PHP3H2(PH)2PH2P3",
    "(PH)2P4(HP)2HP2HPH6P2H3PHP2HPH2P2HPH3P4H",
    "PH2P6H2P3H3PHP2HPH2(P2H)2P2H2P2H7P2H2",
];

fn parse_amino_str(input: &str) -> Vec<Amino> {
    let mut result = Vec::new();
    let mut i = 0;
    let chars: Vec<char> = input.chars().collect();
    while i < chars.len() {
        match chars[i] {
            'H' => result.push(Amino::H),
            'P' => result.push(Amino::P),
            '(' => {
                let mut j = i + 1;
                while j < chars.len() && chars[j] != ')' {
                    j += 1;
                }
                let sequence = &chars[i + 1..j];
                let mut inner_result = parse_amino_str(&sequence.iter().collect::<String>());
                i = j + 1;
                let mut multiplier = 0;
                while i < chars.len() && chars[i].is_numeric() {
                    multiplier = multiplier * 10 + chars[i].to_digit(10).unwrap() as i32;
                    i += 1;
                }
                for _ in 0..multiplier {
                    result.append(&mut inner_result.clone());
                }
                continue;
            }
            c if c.is_ascii_digit() => {
                let mut multiplier = 0;
                while i < chars.len() && chars[i].is_numeric() {
                    multiplier = multiplier * 10 + chars[i].to_digit(10).unwrap() as i32;
                    i += 1;
                }
                if let Some(&last) = result.last() {
                    for _ in 0..multiplier - 1 {
                        result.push(last.clone());
                    }
                }
                continue;
            }
            _ => {}
        }
        i += 1;
    }
    result
}

static SAMPLE_PROTEIN_POINTS: [i32; 22] = [
    4, 9, 9, 8, 14, 23, 21, 36, 42, 53, 50, 48, 32, 34, 34, 33, 32, 32, 32, 31, 34, 33,
];
fn setup() -> Vec<Protein> {
    let mut sample_proteins = Vec::new();
    for i in 0..PROTEIN_DATA.len() {
        let amino_str = PROTEIN_DATA[i];
        let aminos = parse_amino_str(amino_str);
        let mut amino_acids = Vec::new();
        for j in 0..aminos.len() {
            amino_acids.push(AminoAcid {
                amino: aminos[j],
                pos: (0, 0, 0),
            });
        }
        sample_proteins.push(Protein {
            size: aminos.len() as i32,
            aminos: amino_acids,
            direct: Vec::new(),
            predict: 0,
        });
    }
    for i in 0..sample_proteins.len() {
        sample_proteins[i].aminos[0].pos = (0, 0, 0);
        sample_proteins[i].aminos[1].pos = (1, 0, 0);
    }
    sample_proteins[0].direct = vec![Direction::L, Direction::L];
    sample_proteins[1].direct = vec![
        Direction::L,
        Direction::S,
        Direction::L,
        Direction::L,
        Direction::R,
        Direction::R,
        Direction::L,
        Direction::R,
        Direction::L,
        Direction::L,
        Direction::S,
        Direction::L,
        Direction::R,
        Direction::R,
        Direction::L,
        Direction::L,
        Direction::S,
        Direction::L,
    ];
    sample_proteins
}

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    #[arg(short, long)]
    vis: bool,

    #[arg(short, long, default_value_t = 5)]
    step: u8,

    #[arg(short, long, default_value_t = 2)]
    dim: u8,

    #[arg(short, long, default_value_t = 10)]
    id: u8,
}

fn main() {
    let args = Args::parse();
    print!("{:?}", args);
    let mut sample_proteins = setup();
    let mut protein = &mut sample_proteins[args.id as usize];

    if args.vis {
        let mut beam = Beam {
            beam_width: 200,
            nodes: vec![protein.clone()],
            best_score: 0,
            best_ans: protein.clone(),
            num_direct: if args.dim == 2 { 3 } else { 5 },
        };
        beam.vis_one_step(args.step, args.dim);
    } else {
        let mut best_ans = protein.clone();
        let mut best_score = 0;
        for _ in 0..4 {
            let mut beam = Beam {
                beam_width: 200,
                nodes: vec![protein.clone()],
                best_score: 0,
                best_ans: protein.clone(),
                num_direct: if args.dim == 2 { 3 } else { 5 },
            };

            beam.first_step();
            println!("first step completed");

            for i in 0..10 {
                beam.one_step();
                println!("one step completed");
                println!("beam {}: {}", i, beam.best_score);
                beam.local_one_step();
                println!("local step completed");
                println!("beam {}: {}", i, beam.best_score);
            }
            if beam.best_score > best_score {
                best_ans = beam.best_ans.clone();
                best_score = beam.best_score;
            }
            animation(&mut beam.best_ans, args.dim);
        }
        animation(&mut best_ans, args.dim);
    }
}
