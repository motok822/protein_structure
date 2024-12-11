use plotters::prelude::*;
use rand::Rng;
use std::collections::HashMap;
// mod aco;
mod anneal;
mod beam;
mod lib;
// use aco::ACO;
use anneal::Annealing;
use beam::Beam;
use lib::{rotate_left, rotate_right, Amino, AminoAcid, Direction, Heuristics, Protein};
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::hash::{Hash, Hasher};
extern crate piston_window;
use piston_window::*;

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
    // let mut annealing = Annealing {
    //     temperature: 0.9,
    //     max_iter: 1000000,
    //     now_ans: protein.clone(),
    //     now_score: 0,
    //     best_ans: protein.clone(),
    //     best_score: 0,
    //     num_direct: 3,
    // };
    // //ビームサーチ
    // let mut beam = Beam {
    //     beam_width: 200,
    //     nodes: vec![protein.clone()],
    //     best_score: 0,
    //     best_ans: protein.clone(),
    //     num_direct: 3, // 2d
    // };
    // annealing.first_step();

    // for _ in 0..annealing.max_iter {
    //     annealing.one_step();
    //     println!("{}", annealing.best_score);
    // }

    // beam.first_step();

    // for _ in 0..10 {
    //     for i in 0..5 {
    //         beam.one_step();
    //         println!("beam {}: {}", i, beam.best_score);
    //     }
    //     annealing.best_ans = beam.best_ans.clone();
    //     annealing.now_ans = beam.best_ans.clone();
    //     annealing.best_score = beam.best_score;
    //     let mut heap = BinaryHeap::new();
    //     let mut map = HashMap::new();
    //     for i in 0..10000 {
    //         annealing.one_step();
    //         println!("anneal {}: {}", i, annealing.now_score);
    //         heap.push(OrdF32(annealing.now_ans.get_value()));
    //         map.insert(
    //             FloatKey(annealing.now_ans.get_value()),
    //             annealing.now_ans.clone(),
    //         );
    //     }
    //     let mut new_nodes = Vec::new();
    //     for _ in 0..beam.beam_width {
    //         if let Some(value) = heap.pop() {
    //             new_nodes.push(map.get(&FloatKey(value.0)).unwrap().clone());
    //         }
    //     }
    //     beam.best_ans = annealing.best_ans.clone();
    //     beam.best_score = annealing.best_score;
    //     beam.nodes = new_nodes.clone();
    // }

    // let mut heap = BinaryHeap::new();
    // let mut map = HashMap::new();
    // for _ in 0..10 {
    //     let mut beam = Beam {
    //         beam_width: 200,
    //         nodes: vec![protein.clone()],
    //         best_score: 0,
    //         best_ans: protein.clone(),
    //         num_direct: 3, // 2d
    //     };
    //     beam.first_step();
    //     for i in 0..4 {
    //         beam.one_step();
    //         println!("beam {}: {}", i, beam.best_score);
    //     }
    //     for i in 0..beam.nodes.len() {
    //         let value = beam.nodes[i].get_value();
    //         if !map.contains_key(&FloatKey(value)) {
    //             map.insert(FloatKey(value), beam.nodes[i].clone());
    //             heap.push(OrdF32(value));
    //         }
    //     }
    // }

    // let mut new_nodes = Vec::new();
    // for _ in 0..200 {
    //     if let Some(value) = heap.pop() {
    //         new_nodes.push(map.get(&FloatKey(value.0)).unwrap().clone());
    //     }
    // }

    //蟻コロニー最適化
    // let mut aco = ACO {
    //     pheromone: HashMap::new(),
    //     best_score: 0,
    //     protein: protein.clone(),
    //     alpha: 1.0,
    //     beta: 1.0,
    //     evaporation: 0.9,
    //     gamma: 1.0,
    //     num_of_ants: 30,
    // };
    // let max_iter = 10000;
    // aco.first_step(&mut protein);
    // for i in 0..max_iter {
    //     aco.one_step();
    //     println!("{}: {}", i, aco.best_score);
    // }
    let WINDOW_TITLE = "Protein Folding".to_string();
    let WINDOW_SIZE = [800, 800];
    let mut window: PistonWindow = WindowSettings::new(WINDOW_TITLE, WINDOW_SIZE)
        .exit_on_esc(true) //Escキーを押したら終了する
        .vsync(true) //垂直同期を有効にする
        .resizable(false) //ウィンドウのリサイズをさせない
        .samples(4) //アンチエイリアスのサンプル数
        .build()
        .unwrap_or_else(|e| panic!("Failed to build PistonWindow: {}", e));
    window.events.set_max_fps(60); //描画の最大FPS
    window.events.set_ups(60); //1秒間に何回アップデートするか
    let points = vec![
        (0.0, 0.0, 0.0), // 原点
        (100.0, 50.0, 50.0),
        (-50.0, -50.0, 100.0),
        (200.0, 100.0, -50.0),
    ];
    while let Some(e) = window.next() {
        window.draw_3d(&e, |window| {
            let args = e.render_args().unwrap();

            window
                .encoder
                .clear(&window.output_color, [0.3, 0.3, 0.3, 1.0]);
            window.encoder.clear_depth(&window.output_stencil, 1.0);

            data.u_model_view_proj = model_view_projection(
                model,
                first_person.camera(args.ext_dt).orthogonal(),
                projection,
            );
            window.encoder.draw(slice, pipeline, user_data);
            window.encoder.draw(&slice, &pso, &data);
        });
    }
}
