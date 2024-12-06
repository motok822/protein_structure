use crate::lib::{rotate_left, rotate_right, Direction, Heuristics, Protein};
use rand::Rng;
use std::collections::BinaryHeap;
use std::collections::HashMap;

pub struct Beam {
    pub beam_width: i32,
    pub nodes: Vec<Protein>,
    pub best_score: i32,
    pub best_ans: Protein,
}

impl Beam {
    pub fn first_step(&mut self, protein: &mut Protein) {
        while true {
            let mut rng = rand::thread_rng();
            let mut direct = Vec::new();
            for _ in 0..protein.size - 2 {
                let r = rng.gen_range(0..3);
                match r {
                    0 => direct.push(Direction::S),
                    1 => direct.push(Direction::L),
                    2 => direct.push(Direction::R),
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
    pub fn one_step(&mut self, count: usize) {
        let mut heap = BinaryHeap::new();
        let mut map: HashMap<i32, (Protein, i32)> = HashMap::new();
        let directions = vec![Direction::S, Direction::L, Direction::R];
        for i in 0..self.nodes.len() {
            let mut node = self.nodes[i].clone();
            for j in 0..directions.len() {
                let direct = directions[j].clone();
                let mut new_node = node.clone();
                new_node.direct[count] = direct.clone();
                let score = new_node.calc_predict();
                if score != -1 {
                    let max_distance = new_node.calc_max_distance().round() as i32;
                    let mut value = score - (max_distance / 3); //評価関数
                    heap.push(value);
                    if !map.contains_key(&value) {
                        map.insert(value, (new_node, score));
                    }
                }
            }
        }
        let mut new_nodes = Vec::new();
        for i in 0..self.beam_width {
            if let Some(value) = heap.pop() {
                let (node, score) = map.get(&value).unwrap().clone();
                if i == 0 {
                    if self.best_score < score {
                        self.best_score = score;
                        self.best_ans = node;
                    }
                }
                new_nodes.push(map.get(&value).unwrap().clone().0);
            }
        }
        self.nodes = new_nodes;
    }
}
