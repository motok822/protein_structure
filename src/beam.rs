use crate::lib::{rotate_down, rotate_left, rotate_right, rotate_up, Amino, Direction, Protein};
use core::num;
use plotters::prelude::DrawingArea;
use plotters::prelude::*;
use rand::seq::SliceRandom;
use rand::thread_rng;
use rand::Rng;
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::collections::HashMap;
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

pub fn next_permutation<T: Ord>(a: &mut [T]) -> bool {
    let Some(i) = a.windows(2).rposition(|w| w[0] < w[1]) else {
        return false;
    };
    let j = a.iter().rposition(|x| x > &a[i]).unwrap();
    a.swap(i, j);
    a[i + 1..].reverse();
    true
}

pub struct Beam {
    pub beam_width: i32,
    pub nodes: Vec<Protein>,
    pub best_score: i32,
    pub best_ans: Protein,
    pub num_direct: i32,
}

impl Beam {
    pub fn first_step(&mut self) {
        while true {
            let mut rng = rand::thread_rng();
            let mut direct = Vec::new();
            let mut protein = self.best_ans.clone();
            for _ in 0..protein.size - 2 {
                let r = rng.gen_range(0..self.num_direct);
                match r {
                    0 => direct.push(Direction::S),
                    1 => direct.push(Direction::L),
                    2 => direct.push(Direction::R),
                    3 => direct.push(Direction::U),
                    4 => direct.push(Direction::D),
                    _ => {}
                }
            }
            protein.direct = direct.clone();
            let mut score = protein.calc_predict();
            if score != -1 {
                self.best_score = score;
                self.best_ans = protein.clone();
                break;
            }
        }
        self.nodes = vec![self.best_ans.clone()];
    }
    pub fn local_one_step(&mut self) {
        let step_size = 5;
        let mut heap = BinaryHeap::new();
        let mut map = HashMap::new();
        for i in (0..self.nodes.len()).step_by(3) {
            let mut node = self.nodes[i].clone();
            for c in 0..self.best_ans.direct.len() - step_size {
                let mut directions = vec![Direction::S, Direction::L, Direction::R];

                if self.num_direct == 5 {
                    directions.push(Direction::U);
                    directions.push(Direction::D);
                }
                let first_pos = node.clone().aminos[c].pos;
                let target_pos = node.clone().aminos[c + step_size].pos;

                let x_diff = target_pos.0 - first_pos.0;
                let y_diff = target_pos.1 - first_pos.1;
                let z_diff = target_pos.2 - first_pos.2;
                let mut residue =
                    step_size as u32 - (x_diff.abs() + y_diff.abs() + z_diff.abs()) as u32;

                let mut d = Vec::new();
                for _ in 0..x_diff.abs() {
                    if x_diff > 0 {
                        d.push((1, 0, 0));
                    } else {
                        d.push((-1, 0, 0));
                    }
                }
                for _ in 0..y_diff.abs() {
                    if y_diff > 0 {
                        d.push((0, 1, 0));
                    } else {
                        d.push((0, -1, 0));
                    }
                }
                for _ in 0..z_diff.abs() {
                    if z_diff > 0 {
                        d.push((0, 0, 1));
                    } else {
                        d.push((0, 0, -1));
                    }
                }
                let mut ds = Vec::new();
                let mut base = 2 as u32;
                if self.num_direct == 5 {
                    base = 3 as u32;
                }
                for i in 0..(base.pow(residue / 2)) {
                    let mut _d = d.clone();
                    let mut cnt = i;
                    let mut _residue = residue;
                    while _residue > 0 {
                        let r = cnt % base;
                        cnt /= base;
                        _residue -= 2;
                        if r == 0 {
                            _d.push((1, 0, 0));
                            _d.push((-1, 0, 0));
                        } else if r == 1 {
                            _d.push((0, 1, 0));
                            _d.push((0, -1, 0));
                        } else if self.num_direct == 5 && r == 2 {
                            _d.push((0, 0, 1));
                            _d.push((0, 0, -1));
                        }
                    }
                    ds.push(_d);
                }
                // println!("{:?}", ds);
                let mut directions = vec![Direction::S, Direction::L, Direction::R];
                if self.num_direct == 5 {
                    directions.push(Direction::U);
                    directions.push(Direction::D);
                }

                for i in 0..ds.len() {
                    let mut d = ds[i].clone();
                    while true {
                        let mut new_node = node.clone();
                        for j in c..c + step_size {
                            let now_direct = (d[j - c].0, d[j - c].1, d[j - c].2);
                            new_node.aminos[j + 2].pos = (
                                new_node.aminos[j + 1].pos.0 + now_direct.0,
                                new_node.aminos[j + 1].pos.1 + now_direct.1,
                                new_node.aminos[j + 1].pos.2 + now_direct.2,
                            );
                            let prev_direct = (
                                new_node.aminos[j + 1].pos.0 - new_node.aminos[j].pos.0,
                                new_node.aminos[j + 1].pos.1 - new_node.aminos[j].pos.1,
                                new_node.aminos[j + 1].pos.2 - new_node.aminos[j].pos.2,
                            );
                            if rotate_left(prev_direct) == now_direct {
                                new_node.direct[j] = Direction::L;
                            } else if rotate_right(prev_direct) == now_direct {
                                new_node.direct[j] = Direction::R;
                            } else if self.num_direct == 5 && rotate_up(prev_direct) == now_direct {
                                new_node.direct[j] = Direction::U;
                            } else if self.num_direct == 5 && rotate_down(prev_direct) == now_direct
                            {
                                new_node.direct[j] = Direction::D;
                            } else {
                                new_node.direct[j] = Direction::S;
                            }
                        }
                        let score = new_node.calc_predict();
                        if score != -1 {
                            let mut value = new_node.get_value();
                            if !map.contains_key(&value) {
                                map.insert(value, (new_node, score));
                                heap.push(value);
                            }
                        }
                        if next_permutation(&mut d) == false {
                            break;
                        }
                    }
                }
            }
        }
        let mut new_nodes: Vec<Protein> = Vec::new();
        for i in 0..self.beam_width {
            if let Some(value) = heap.pop() {
                let (node, score) = map.remove(&value).unwrap();
                if i == 0 {
                    if self.best_score < score {
                        self.best_score = score;
                        self.best_ans = node.clone();
                    }
                }
                new_nodes.push(node.clone());
            }
        }
        self.nodes = new_nodes;
    }
    pub fn one_step(&mut self) {
        for c in 0..self.best_ans.direct.len() {
            let mut heap = BinaryHeap::new();
            let mut map = HashMap::new();
            let mut directions = vec![Direction::S, Direction::L, Direction::R];

            if self.num_direct == 5 {
                directions.push(Direction::U);
                directions.push(Direction::D);
            }

            for i in 0..self.nodes.len() {
                let mut node = self.nodes[i].clone();
                for j in 0..directions.len() {
                    let direct = directions[j];
                    let mut new_node = node.clone();
                    new_node.direct[c] = direct;
                    let score = new_node.calc_predict();
                    if score != -1 {
                        let mut value = new_node.get_value();
                        if !map.contains_key(&value) {
                            map.insert(value, (new_node, score));
                            heap.push(value);
                        }
                    }
                }
            }
            let mut new_nodes: Vec<Protein> = Vec::new();
            let mut x = 0;
            while heap.len() > 0 {
                if let Some(value) = heap.pop() {
                    let (node, score) = map.remove(&value).unwrap();
                    if self.best_score < score {
                        self.best_score = score;
                        self.best_ans = node.clone();
                    }
                    let exp = (-8.0 * (x * x) as f32
                        / ((self.beam_width * self.beam_width) as f32))
                        .exp();
                    let mut rng = rand::thread_rng();
                    let prob = rng.gen_range(0.0..1.0) as f32;
                    if prob < exp {
                        x += 1;
                        new_nodes.push(node.clone());
                    }
                }
                if new_nodes.len() >= self.beam_width as usize {
                    break;
                }
                x += 1;
            }
            self.nodes = new_nodes;
        }
    }
    pub fn vis_one_step(&mut self, vis_step: u8, dim: u8) {
        let area = BitMapBackend::gif(
            "./animated.gif", // アニメーションファイルの名前。この名前で保存される
            (1200, 800),      //  グラフのサイズ（幅x高さ)
            1,                //  1フレームの時間。単位は [ms]
        )
        .unwrap()
        .into_drawing_area();
        for _ in 0..4 {
            self.first_step();
            println!("first step completed");

            for i in 0..10 {
                for c in 0..self.best_ans.direct.len() {
                    let mut heap = BinaryHeap::new();
                    let mut map = HashMap::new();
                    let mut directions = vec![Direction::S, Direction::L, Direction::R];

                    if self.num_direct == 5 {
                        directions.push(Direction::U);
                        directions.push(Direction::D);
                    }

                    for i in 0..self.nodes.len() {
                        let mut node = self.nodes[i].clone();
                        for j in 0..directions.len() {
                            let direct = directions[j];
                            let mut new_node = node.clone();
                            new_node.direct[c] = direct;
                            let score = new_node.calc_predict();
                            if score != -1 {
                                let mut value = new_node.get_value();
                                if !map.contains_key(&value) {
                                    map.insert(value, (new_node, score));
                                    heap.push(value);
                                }
                            }
                        }
                    }
                    let mut new_nodes: Vec<Protein> = Vec::new();
                    let mut x = 0;
                    let mut count = 0;
                    while heap.len() > 0 {
                        if let Some(value) = heap.pop() {
                            let (node, score) = map.remove(&value).unwrap();
                            if self.best_score < score {
                                self.best_score = score;
                                self.best_ans = node.clone();
                            }
                            let exp = (-8.0 * (x * x) as f32
                                / ((self.beam_width * self.beam_width) as f32))
                                .exp();
                            let mut rng = rand::thread_rng();
                            let prob = rng.gen_range(0.0..1.0) as f32;
                            if prob < exp {
                                x += 1;
                                new_nodes.push(node.clone());
                                count += 1;
                                if count == vis_step {
                                    count = 0;
                                    //visualize の処理
                                    let mut data = Vec::new();
                                    data.push((0, 0, 0, node.aminos[0].amino));
                                    if dim == 2 {
                                        data.push((5, 0, 0, node.aminos[1].amino));
                                    } else {
                                        data.push((10, 0, 0, node.aminos[1].amino));
                                    }
                                    for i in 0..node.direct.len() {
                                        let mut prev_direction = (
                                            data[i + 1].0 - data[i].0,
                                            data[i + 1].1 - data[i].1,
                                            data[i + 1].2 - data[i].2,
                                        );
                                        let (x, y, z, a) = match node.direct[i] {
                                            Direction::S => (
                                                data[i + 1].0 + prev_direction.0,
                                                data[i + 1].1 + prev_direction.1,
                                                data[i + 1].2 + prev_direction.2,
                                                node.aminos[i + 2].amino,
                                            ),
                                            Direction::L => (
                                                data[i + 1].0 + rotate_left(prev_direction).0,
                                                data[i + 1].1 + rotate_left(prev_direction).1,
                                                data[i + 1].2 + rotate_left(prev_direction).2,
                                                node.aminos[i + 2].amino,
                                            ),
                                            Direction::R => (
                                                data[i + 1].0 + rotate_right(prev_direction).0,
                                                data[i + 1].1 + rotate_right(prev_direction).1,
                                                data[i + 1].2 + rotate_right(prev_direction).2,
                                                node.aminos[i + 2].amino,
                                            ),
                                            Direction::U => (
                                                data[i + 1].0 + rotate_up(prev_direction).0,
                                                data[i + 1].1 + rotate_up(prev_direction).1,
                                                data[i + 1].2 + rotate_up(prev_direction).2,
                                                node.aminos[i + 2].amino,
                                            ),
                                            Direction::D => (
                                                data[i + 1].0 + rotate_down(prev_direction).0,
                                                data[i + 1].1 + rotate_down(prev_direction).1,
                                                data[i + 1].2 + rotate_down(prev_direction).2,
                                                node.aminos[i + 2].amino,
                                            ),
                                        };
                                        data.push((x, y, z, a));
                                    }

                                    area.fill(&WHITE).unwrap();

                                    let mut chart = ChartBuilder::on(&area)
                                        .margin(20)
                                        .caption("protein structure", ("sans-serif", 40))
                                        .build_cartesian_3d(-100..100, -100..100, -100..100)
                                        .unwrap();

                                    chart.configure_axes().draw().unwrap();

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
                        if new_nodes.len() >= self.beam_width as usize {
                            break;
                        }
                        x += 1;
                    }
                    self.nodes = new_nodes;
                }

                println!("one step completed");
                println!("beam {}: {}", i, self.best_score);
            }
        }
    }
}
