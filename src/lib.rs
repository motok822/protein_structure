use std::cmp::{max, min};
use std::collections::HashMap;
use std::path::absolute;
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
    U = 4,
    D = 5,
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
}
pub fn rotate_left((x, y, z): (i32, i32, i32)) -> (i32, i32, i32) {
    if z == 0 {
        return (y, -x, z);
    } else {
        return (x, -z, y);
    }
}
pub fn rotate_right((x, y, z): (i32, i32, i32)) -> (i32, i32, i32) {
    if z == 0 {
        return (-y, x, z);
    } else {
        return (x, z, -y);
    }
}
pub fn rotate_up((x, y, z): (i32, i32, i32)) -> (i32, i32, i32) {
    if y == 0 {
        return (-z, y, x);
    } else {
        return (-y, x, z);
    }
}
pub fn rotate_down((x, y, z): (i32, i32, i32)) -> (i32, i32, i32) {
    if y == 0 {
        return (z, y, -x);
    } else {
        return (y, -x, z);
    }
}

impl Protein {
    pub fn get_value(&mut self) -> i32 {
        let mut score = self.calc_predict() as f32;
        let mut max_distance = self.calc_max_distance();
        // let (above_cube, area) = self.get_above_cube();
        // score as f32 - area / 5.0 - above_cube / 5.0
        ((score - max_distance / 3.0) * 10.0).round() as i32
    }

    pub fn calc_predict(&mut self) -> i32 {
        let mut map: HashMap<(i32, i32, i32), Amino> = HashMap::new();
        map.insert((0, 0, 0), self.aminos[0].amino);
        map.insert((1, 0, 0), self.aminos[1].amino);
        let mut previous_direct = (1, 0, 0);
        let mut last_pos = (1, 0, 0);
        let mut result = 0;
        self.aminos[0].pos = (0, 0, 0);
        self.aminos[1].pos = (1, 0, 0);
        assert!(self.aminos.len() == self.direct.len() + 2);
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
                Direction::U => (
                    last_pos.0 + rotate_up(previous_direct).0,
                    last_pos.1 + rotate_up(previous_direct).1,
                    last_pos.2 + rotate_up(previous_direct).2,
                ),
                Direction::D => (
                    last_pos.0 + rotate_down(previous_direct).0,
                    last_pos.1 + rotate_down(previous_direct).1,
                    last_pos.2 + rotate_down(previous_direct).2,
                ),
            };
            if map.contains_key(&(x, y, z)) {
                self.predict = -1;
                return -1;
            }
            let now_amino = self.aminos[i + 2].amino;
            let mut count = 0;
            let dxdy = vec![
                (1, 0, 0),
                (-1, 0, 0),
                (0, 1, 0),
                (0, -1, 0),
                (0, 0, 1),
                (0, 0, -1),
            ];
            previous_direct = (x - last_pos.0, y - last_pos.1, z - last_pos.2);
            if (now_amino == Amino::H) {
                for j in 0..6 {
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
    pub fn get_above_cube(&mut self) -> (f32, f32) {
        let mut all_pos = Vec::new();
        all_pos.push((0, 0, 0));
        all_pos.push((1, 0, 0));
        let mut last_pos = (1, 0, 0);
        let mut last_direct = (1, 0, 0);
        let mut min_x = 1000;
        let mut max_x = -1000;
        let mut min_y = 1000;
        let mut max_y = -1000;
        let mut min_z = 1000;
        let mut max_z = -1000;
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
                Direction::U => (
                    last_pos.0 + rotate_up(last_direct).0,
                    last_pos.1 + rotate_up(last_direct).1,
                    last_pos.2 + rotate_up(last_direct).2,
                ),
                Direction::D => (
                    last_pos.0 + rotate_down(last_direct).0,
                    last_pos.1 + rotate_down(last_direct).1,
                    last_pos.2 + rotate_down(last_direct).2,
                ),
            };
            all_pos.push((x, y, z));
            last_direct = (x - last_pos.0, y - last_pos.1, z - last_pos.2);
            last_pos = (x, y, z);
        }
        for i in 0..all_pos.len() {
            if self.aminos[i].amino != Amino::H {
                continue;
            }
            let (x1, y1, z1) = all_pos[i];
            min_x = min(min_x, x1);
            max_x = max(max_x, x1);
            min_y = min(min_y, y1);
            max_y = max(max_y, y1);
            min_z = min(min_z, z1);
            max_z = max(max_z, z1);
        }
        if min_z == max_z {
            let average = (((max_x - min_x) as f32) + ((max_y - min_y) as f32)) / 2.0;
            return (
                average * average - (max_x - min_x) as f32 * (max_y - min_y) as f32,
                (max_x - min_x) as f32 * (max_y - min_y) as f32,
            );
        }

        let average = ((max_x as f32 - min_x as f32)
            + (max_y as f32 - min_y as f32)
            + (max_z as f32 - min_z as f32))
            / 3.0;
        return (
            average * average * average
                - (max_x - min_x) as f32 * (max_y - min_y) as f32 * (max_z - min_z) as f32,
            (max_x - min_x) as f32 * (max_y - min_y) as f32 * (max_z - min_z) as f32,
        );
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
                Direction::U => (
                    last_pos.0 + rotate_up(last_direct).0,
                    last_pos.1 + rotate_up(last_direct).1,
                    last_pos.2 + rotate_up(last_direct).2,
                ),
                Direction::D => (
                    last_pos.0 + rotate_down(last_direct).0,
                    last_pos.1 + rotate_down(last_direct).1,
                    last_pos.2 + rotate_down(last_direct).2,
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
