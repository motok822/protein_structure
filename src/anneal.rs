use crate::lib::{Direction, Heuristics, Protein};
use rand::Rng;
pub struct Annealing {
    pub temperature: f64,
    pub max_iter: i32,
    pub now_ans: Protein,
    pub now_score: i32,
    pub best_ans: Protein,
    pub best_score: i32,
}

impl Heuristics for Annealing {
    fn first_step(&mut self) {
        while true {
            let mut rng = rand::thread_rng();
            let mut direct = Vec::new();
            let mut protein = self.now_ans.clone();
            for _ in 0..protein.size - 2 {
                let r = rng.gen_range(0..3);
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
                self.now_score = score;
                self.now_ans = protein.clone();
                self.best_ans = protein.clone();
                self.best_score = score;
                break;
            }
        }
    }
    fn one_step(&mut self) {
        let mut rng = rand::thread_rng();
        let mut direct = self.now_ans.direct.clone();
        let mut score = self.now_score;
        let mut new_direct = direct.clone();
        let mut protein = self.now_ans.clone();
        for _ in 0..5 {
            let mut r = rng.gen_range(0..direct.len());
            let mut new_r = rng.gen_range(0..5);
            if direct[r] as i32 == new_r {
                new_r += 1;
                new_r %= 5;
            }
            match new_r {
                0 => new_direct[r] = Direction::S,
                1 => new_direct[r] = Direction::L,
                2 => new_direct[r] = Direction::R,
                3 => new_direct[r] = Direction::U,
                4 => new_direct[r] = Direction::D,
                _ => {}
            }
        }
        protein.direct = new_direct.clone();
        let new_score = protein.calc_predict();
        let new_value = protein.get_value();

        if new_score != -1 {
            if new_value > self.best_ans.get_value() {
                self.now_score = new_score;
                self.now_ans = protein.clone();
                self.now_ans.direct = new_direct.clone();
            } else {
                let diff = (new_value - self.best_ans.get_value()) as f64;
                let prob = (diff / self.temperature).exp();
                let p = rng.gen_range(0.0..1.0);
                if p < prob {
                    self.now_score = new_score;
                    self.now_ans = protein.clone();
                    self.now_ans.direct = new_direct.clone();
                }
            }
        }

        if (self.now_score > self.best_score) {
            self.best_score = self.now_score;
            self.best_ans = self.now_ans.clone();
        }
    }
}
