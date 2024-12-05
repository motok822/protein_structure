use plotters::prelude::*;
use std::collections::HashMap;

#[derive(PartialEq, Debug, Clone, Copy)]
enum Amino {
    H = 1,
    P = 2,
}
enum Direction {
    S = 1,
    L = 2,
    R = 3,
}
struct AminoAcid {
    amino: Amino,
    pos: (i32, i32),
}
struct Protein {
    size: i32,
    aminos: Vec<AminoAcid>,
    direct: Vec<Direction>,
    predict: i32,
    point: i32,
}
fn rotate_left((x, y, z): (f64, f64, f64)) -> (f64, f64, f64) {
    (y, -x, z)
}
fn rotate_right((x, y, z): (f64, f64, f64)) -> (f64, f64, f64) {
    (-y, x, z)
}

impl Protein {
    fn plot_2d(&mut self) {}
    fn calc_predict(&mut self) -> i32 {
        let mut map: HashMap<(i32, i32), Amino> = HashMap::new();
        map.insert((0, 0), self.aminos[0].amino);
        map.insert((1, 0), self.aminos[1].amino);
        let mut previous_direct = (1, 0);
        let mut last_pos = (1, 0);
        let mut result = 0;
        for i in 0..self.direct.len() {
            let (x, y) = match self.direct[i] {
                Direction::S => (
                    last_pos.0 + previous_direct.0,
                    last_pos.1 + previous_direct.1,
                ),
                Direction::L => (
                    last_pos.0 - previous_direct.1,
                    last_pos.1 + previous_direct.0,
                ),
                Direction::R => (
                    last_pos.0 + previous_direct.1,
                    last_pos.1 - previous_direct.0,
                ),
            };
            if map.contains_key(&(x, y)) {
                self.predict = 0;
                return 0;
            }
            let now_amino = self.aminos[i + 2].amino;
            let mut count = 0;
            let dxdy = vec![(1, 0), (0, 1), (-1, 0), (0, -1)];
            previous_direct = (x - last_pos.0, y - last_pos.1);
            if (now_amino == Amino::H) {
                for j in 0..4 {
                    let (dx, dy) = dxdy[j];
                    let (nx, ny) = (x + dx, y + dy);
                    if (map.contains_key(&(nx, ny))
                        && map[&(nx, ny)] == Amino::H
                        && (-dx, -dy) != previous_direct)
                    {
                        count += 1;
                    }
                }
            }
            last_pos = (x, y);
            map.insert((x, y), now_amino);
            result += count;
        }
        self.predict = result;
        result
    }

    fn animation(&mut self) {
        let area = BitMapBackend::gif(
            "./animated.gif", // アニメーションファイルの名前。この名前で保存される
            (1200, 800),      //  グラフのサイズ（幅x高さ)
            100,              //  1フレームの時間。単位は [ms]
        )
        .unwrap()
        .into_drawing_area();

        let mut data = Vec::new();
        data.push((0.0, 0.0, 0.0, self.aminos[0].amino));
        data.push((5.0, 0.0, 0.0, self.aminos[1].amino));
        for i in 0..self.direct.len() {
            let mut prev_direction = (
                data[i + 1].0 - data[i].0,
                data[i + 1].1 - data[i].1,
                data[i + 1].2 - data[i].2,
            );
            let (x, y, z, a) = match self.direct[i] {
                Direction::S => (
                    data[i + 1].0 + prev_direction.0,
                    data[i + 1].1 + prev_direction.1,
                    data[i + 1].2 + prev_direction.2,
                    self.aminos[i + 2].amino,
                ),
                Direction::L => (
                    data[i + 1].0 + rotate_left(prev_direction).0,
                    data[i + 1].1 + rotate_left(prev_direction).1,
                    data[i + 1].2 + rotate_left(prev_direction).2,
                    self.aminos[i + 2].amino,
                ),
                Direction::R => (
                    data[i + 1].0 + rotate_right(prev_direction).0,
                    data[i + 1].1 + rotate_right(prev_direction).1,
                    data[i + 1].2 + rotate_right(prev_direction).2,
                    self.aminos[i + 2].amino,
                ),
            };
            data.push((x, y, z, a));
        }

        for _ in 0..=20 {
            area.fill(&WHITE).unwrap();

            let mut chart = ChartBuilder::on(&area)
                .margin(20)
                .caption("protein structure", ("sans-serif", 40))
                .build_cartesian_3d(-20.0..20.0, -20.0..20.0, -20.0..20.0)
                .unwrap();

            chart.configure_axes().draw().unwrap();

            for i in 0..data.len() {
                data[i].0 += 1.0;
            }

            chart
                .draw_series(data.iter().map(|&(x, y, z, a)| {
                    if (a == Amino::H) {
                        TriangleMarker::new((x, y, z), 5, &RED)
                    } else {
                        TriangleMarker::new((x, y, z), 5, &BLUE)
                    }
                }))
                .unwrap();

            chart
                .draw_series(LineSeries::new(
                    data.iter().map(|&(x, y, z, _)| (x, y, z)),
                    &BLACK.mix(0.3),
                ))
                .unwrap();

            // グラフを更新する
            area.present().unwrap();
        }
    }
}

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
                pos: (0, 0),
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
        sample_proteins[i].aminos[0].pos = (0, 0);
        sample_proteins[i].aminos[1].pos = (1, 0);
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
    sample_proteins[1].animation();
}
