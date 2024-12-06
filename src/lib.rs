use std::collections::HashMap;
#[derive(PartialEq, Debug, Clone, Copy)]
pub enum Amino {
    H = 1,
    P = 2,
}
#[derive(PartialEq, Debug, Clone, Copy)]
pub enum Direction {
    S = 1,
    L = 2,
    R = 3,
}
#[derive(PartialEq, Debug, Clone, Copy)]
pub struct AminoAcid {
    pub amino: Amino,
    pub pos: (i32, i32, i32),
}
#[derive(PartialEq, Debug, Clone)]
pub struct Protein {
    pub size: i32,
    pub aminos: Vec<AminoAcid>,
    pub direct: Vec<Direction>,
    pub predict: i32,
    pub point: i32,
}
pub fn rotate_left((x, y, z): (i32, i32, i32)) -> (i32, i32, i32) {
    (y, -x, z)
}
pub fn rotate_right((x, y, z): (i32, i32, i32)) -> (i32, i32, i32) {
    (-y, x, z)
}

impl Protein {
    pub fn calc_predict(&mut self) -> i32 {
        let mut map: HashMap<(i32, i32, i32), Amino> = HashMap::new();
        map.insert((0, 0, 0), self.aminos[0].amino);
        map.insert((1, 0, 0), self.aminos[1].amino);
        let mut previous_direct = (1, 0, 0);
        let mut last_pos = (1, 0, 0);
        let mut result = 0;
        self.aminos[0].pos = (0, 0, 0);
        self.aminos[1].pos = (1, 0, 0);
        for i in 0..self.direct.len() {
            let (x, y, z) = match self.direct[i] {
                Direction::S => (
                    last_pos.0 + previous_direct.0,
                    last_pos.1 + previous_direct.1,
                    last_pos.2 + previous_direct.2,
                ),
                Direction::L => (
                    last_pos.0 + rotate_left(previous_direct).0,
                    last_pos.1 + rotate_left(previous_direct).1,
                    last_pos.2 + rotate_left(previous_direct).2,
                ),
                Direction::R => (
                    last_pos.0 + rotate_right(previous_direct).0,
                    last_pos.1 + rotate_right(previous_direct).1,
                    last_pos.2 + rotate_right(previous_direct).2,
                ),
            };
            if map.contains_key(&(x, y, z)) {
                self.predict = -1;
                return -1;
            }
            let now_amino = self.aminos[i + 2].amino;
            let mut count = 0;
            let dxdy = vec![(1, 0, 0), (0, 1, 0), (-1, 0, 0), (0, -1, 0)];
            previous_direct = (x - last_pos.0, y - last_pos.1, z - last_pos.2);
            if (now_amino == Amino::H) {
                for j in 0..4 {
                    let (dx, dy, dz) = dxdy[j];
                    let (nx, ny, nz) = (x + dx, y + dy, z + dz);
                    if (map.contains_key(&(nx, ny, nz))
                        && map[&(nx, ny, nz)] == Amino::H
                        && (-dx, -dy, -dz) != previous_direct)
                    {
                        count += 1;
                    }
                }
            }
            last_pos = (x, y, z);
            map.insert((x, y, z), now_amino);
            self.aminos[i + 2].pos = (x, y, z);
            result += count;
        }
        self.predict = result;
        result
    }
    pub fn calc_max_distance(&mut self) -> f32 {
        let mut max_distance = 0.0;
        let mut all_pos = Vec::new();
        all_pos.push((0, 0, 0));
        all_pos.push((1, 0, 0));
        let mut last_pos = (1, 0, 0);
        let mut last_direct = (1, 0, 0);
        for i in 0..self.direct.len() {
            let (x, y, z) = match self.direct[i] {
                Direction::S => (
                    last_pos.0 + last_direct.0,
                    last_pos.1 + last_direct.1,
                    last_pos.2 + last_direct.2,
                ),
                Direction::L => (
                    last_pos.0 + rotate_left(last_direct).0,
                    last_pos.1 + rotate_left(last_direct).1,
                    last_pos.2 + rotate_left(last_direct).2,
                ),
                Direction::R => (
                    last_pos.0 + rotate_right(last_direct).0,
                    last_pos.1 + rotate_right(last_direct).1,
                    last_pos.2 + rotate_right(last_direct).2,
                ),
            };
            all_pos.push((x, y, z));
            last_direct = (x - last_pos.0, y - last_pos.1, z - last_pos.2);
            last_pos = (x, y, z);
        }
        for i in 0..all_pos.len() {
            for j in 0..all_pos.len() {
                let (x1, y1, z1) = all_pos[i];
                let (x2, y2, z2) = all_pos[j];
                let distance = (((x1 - x2) * (x1 - x2)
                    + (y1 - y2) * (y1 - y2)
                    + (z1 - z2) * (z1 - z2)) as f32)
                    .sqrt();
                if distance > max_distance {
                    max_distance = distance;
                }
            }
        }
        max_distance
    }
}

pub trait Heuristics {
    fn first_step(&mut self, protein: &mut Protein);
    fn one_step(&mut self, protein: &mut Protein);
    fn get_value(&self) -> f32;
}
