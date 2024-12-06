use plotters::prelude::*;
use rand::Rng;
use std::collections::HashMap;
mod anneal;
mod beam;
mod lib;
use anneal::Annealing;
use beam::Beam;
use lib::{rotate_left, rotate_right, Amino, AminoAcid, Direction, Heuristics, Protein};

static PROTEIN_DATA: [&str; 12] = [
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

static SAMPLE_PROTEIN_POINTS: [i32; 12] = [4, 9, 9, 8, 14, 23, 21, 36, 42, 53, 50, 48];
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
            point: SAMPLE_PROTEIN_POINTS[i],
        });
    }
    for i in 0..sample_proteins.len() {
        sample_proteins[i].aminos[0].pos = (0, 0, 0);
        sample_proteins[i].aminos[1].pos = (1, 0, 0);
    }
    sample_proteins[0].direct = vec![Direction::L, Direction::L];
    println!("{}", sample_proteins[0].calc_predict()); //test
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
    println!("{}", sample_proteins[1].calc_predict()); // test
    sample_proteins
}
fn main() {
    let mut sample_proteins = setup();
    let mut protein = &mut sample_proteins[11];
    //焼きなまし法
    let mut annealing = Annealing {
        temperature: 0.6,
        max_iter: 10000,
        now_ans: Vec::new(),
        now_score: 0,
        max_distance: 0.0,
        best_ans: Vec::new(),
        best_score: 0,
    };
    annealing.first_step(&mut protein);
    for i in 0..annealing.max_iter {
        annealing.one_step(&mut protein);
        println!("{}: {}", i, annealing.now_score);
    }
    println!("best score: {}", annealing.best_score);

    //ビームサーチ
    protein.direct = annealing.best_ans.clone();
    let mut beam = Beam {
        beam_width: 40,
        nodes: vec![protein.clone()],
        best_score: 0,
        best_ans: protein.clone(),
    };
    // beam.first_step(protein);
    for j in 0..5 {
        for i in 0..protein.direct.len() {
            beam.one_step(i);
            println!("now best {}: {}", i, beam.best_score);
        }
    }
    println!("best score: {}", beam.best_score);
}
