use crate::lib::{Direction, Heuristics, Protein};
use rand::Rng;
pub struct Annealing {
    pub temperature: f64,
    pub max_iter: i32,
    pub now_ans: Vec<Direction>,
    pub now_score: i32,
    pub max_distance: f32,
    pub best_ans: Vec<Direction>,
    pub best_score: i32,
}

impl Heuristics for Annealing {
    fn first_step(&mut self, protein: &mut Protein) {
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
                self.now_score = score;
                self.now_ans = direct.clone();
                self.max_distance = protein.calc_max_distance();
                self.best_ans = direct.clone();
                self.best_score = score;
                break;
            }
        }
    }
    fn get_value(&self) -> f32 {
        self.now_score as f32 - self.max_distance / 3.0
    }
    fn one_step(&mut self, protein: &mut Protein) {
        let mut rng = rand::thread_rng();
        let mut direct = self.now_ans.clone();
        let mut score = self.now_score;
        let mut max_distance = self.max_distance;
        let mut new_direct = direct.clone();
        let mut r = rng.gen_range(0..direct.len());
        match direct[r] {
            Direction::S => {
                let mut new_r = rng.gen_range(0..3);
                while new_r == 0 {
                    new_r = rng.gen_range(0..3);
                }
                match new_r {
                    1 => new_direct[r] = Direction::L,
                    2 => new_direct[r] = Direction::R,
                    _ => {}
                }
            }
            Direction::L => {
                let mut new_r = rng.gen_range(0..3);
                while new_r == 1 {
                    new_r = rng.gen_range(0..3);
                }
                match new_r {
                    0 => new_direct[r] = Direction::S,
                    2 => new_direct[r] = Direction::R,
                    _ => {}
                }
            }
            Direction::R => {
                let mut new_r = rng.gen_range(0..3);
                while new_r == 2 {
                    new_r = rng.gen_range(0..3);
                }
                match new_r {
                    0 => new_direct[r] = Direction::S,
                    1 => new_direct[r] = Direction::L,
                    _ => {}
                }
            }
        }
        protein.direct = new_direct.clone();
        let new_score = protein.calc_predict();
        let new_distance = protein.calc_max_distance();
        let new_value = new_score as f32 - new_distance / 3.0;
        println!("{}: {}", new_score, new_distance);

        if new_score != -1 {
            if new_value > self.get_value() {
                self.now_score = new_score;
                self.max_distance = max_distance;
                self.now_ans = new_direct.clone();
            } else {
                let diff = (new_value - self.get_value()) as f64;
                let prob = (diff / self.temperature).exp();
                let p = rng.gen_range(0.0..1.0);
                if p < prob {
                    self.now_score = new_score;
                    self.max_distance = max_distance;
                    self.now_ans = new_direct.clone();
                }
            }
        }

        if (self.now_score > self.best_score) {
            self.best_score = self.now_score;
            self.best_ans = self.now_ans.clone();
        }
    }
}
