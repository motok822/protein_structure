use crate::lib::{rotate_left, rotate_right, Direction, Heuristics, Protein};
use rand::Rng;
use std::collections::BinaryHeap;
use std::collections::HashMap;
use std::ops::Not;

pub struct ACO {
    pub pheromone: HashMap<(i32, i32), f64>,
    pub best_score: i32,
    pub protein: Protein,
    pub alpha: f64,
    pub beta: f64,
    pub evaporation: f64,
    pub gamma: f64,
    pub num_of_ants: i32,
}

impl ACO {
    pub fn first_step(&mut self, protein: &mut Protein) {
        let mut cnt = 0;
        while true {
            let mut rng = rand::thread_rng();
            let mut direct = Vec::new();
            for i in 0..protein.size - 2 {
                let r = rng.gen_range(0..3);
                match r {
                    0 => direct.push(Direction::S),
                    1 => direct.push(Direction::L),
                    2 => direct.push(Direction::R),
                    _ => {}
                }
                // self.pheromone.insert((i as i32, r as i32), 1.0);
            }
            protein.direct = direct.clone();
            let mut score = protein.calc_predict();
            if score != -1 {
                self.best_score = score;
                self.protein = protein.clone();
                break;
            }
            cnt += 1;
        }
    }
    pub fn one_step(&mut self) {
        let mut directions = vec![Direction::S, Direction::L, Direction::R];
        let mut protein = self.protein.clone();
        let mut new_proteins = Vec::new();
        let mut all_routes = Vec::new();
        let mut all_scores = Vec::new();
        let mut cnt = 0;
        for _ in 0..self.num_of_ants {
            let mut route = Vec::new();
            for i in 0..protein.direct.len() {
                let mut probs = Vec::new();
                let mut sum = 0.0;
                for j in 0..directions.len() {
                    let direct = directions[j].clone();
                    if (!self.pheromone.contains_key(&(i as i32, direct as i32))) {
                        self.pheromone.insert((i as i32, direct as i32), 1.0);
                    }
                    let p = self.pheromone.get(&(i as i32, direct as i32)).unwrap();
                    let mut new_protein = protein.clone();
                    new_protein.direct[i] = direct.clone();
                    let mut score = new_protein.calc_predict();
                    if score == -1 {
                        score = 0;
                    }
                    let dist = (-(score as f64) / self.gamma).exp();
                    let prob = p.powf(self.alpha) * dist.powf(self.beta);
                    probs.push(prob);
                    sum += prob;
                }
                for j in 0..probs.len() {
                    probs[j] /= sum;
                }
                let mut rng = rand::thread_rng();
                let mut r = rng.gen_range(0.0..1.0);
                for j in 0..probs.len() {
                    r -= probs[j];
                    if r <= 0.0 {
                        route.push(directions[j].clone());
                        break;
                    }
                }
            }
            let mut new_protein = protein.clone();
            new_protein.direct = route.clone();
            let mut score = new_protein.calc_predict();
            if score == -1 {
                score = 0;
            }
            new_proteins.push(new_protein.clone());
            all_routes.push(route.clone());
            all_scores.push(score);
        }
        // println!("score cnt: {}", cnt);
        self.update_pheromone(all_routes.clone(), all_scores.clone());

        for i in 0..new_proteins.len() {
            if self.best_score < all_scores[i] {
                self.best_score = all_scores[i];
                self.protein = new_proteins[i].clone();
            }
        }
    }
    pub fn update_pheromone(&mut self, route: Vec<Vec<Direction>>, score: Vec<i32>) {
        let mut update_map: HashMap<(i32, i32), f64> = HashMap::new();
        for i in 0..route.len() {
            let route = route[i].clone();
            let score = score[i];
            for j in 0..route.len() {
                let direct = route[j].clone();
                if (!update_map.contains_key(&(j as i32, direct as i32))) {
                    update_map.insert((j as i32, direct as i32), 0.0);
                }
                let value = update_map.get(&(j as i32, direct as i32)).unwrap();
                update_map.insert((j as i32, direct as i32), value + score as f64);
                j
            }
        }
        for i in 0..route.len() {
            for j in 0..route[i].len() {
                let direct = route[i][j].clone();
                let value = update_map.get(&(j as i32, direct as i32)).unwrap();
                if !self.pheromone.contains_key(&(j as i32, direct as i32)) {
                    self.pheromone.insert((j as i32, direct as i32), 1.0);
                }
                let p = self.pheromone.get(&(j as i32, direct as i32)).unwrap();
                let mut next_p = p * self.evaporation + value;
                if next_p < 0.0 {
                    next_p = 0.0;
                }
                self.pheromone.insert((j as i32, direct as i32), next_p);
            }
        }
    }
}
