use crate::lib::{rotate_left, rotate_right, Direction, Heuristics, Protein};
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
                self.nodes[0] = protein.clone();
                break;
            }
        }
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
                    let direct = directions[j].clone();
                    let mut new_node = node.clone();
                    new_node.direct[c] = direct.clone();
                    let score = new_node.calc_predict();
                    if score != -1 {
                        let mut value = new_node.get_value();
                        if !map.contains_key(&FloatKey(value)) {
                            map.insert(FloatKey(value), (new_node, score));
                            heap.push(OrdF32(value));
                        }
                    }
                }
            }
            let mut new_nodes = Vec::new();
            for i in 0..self.beam_width {
                if let Some(value) = heap.pop() {
                    let (node, score) = map.get(&FloatKey(value.0)).unwrap().clone();
                    if i == 0 {
                        if self.best_score < score {
                            self.best_score = score;
                            self.best_ans = node;
                        }
                    }
                    new_nodes.push(map.get(&FloatKey(value.0)).unwrap().0.clone());
                }
            }
            self.nodes = new_nodes;
        }
    }
}
