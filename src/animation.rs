use crate::lib::{
    rotate_down, rotate_left, rotate_right, rotate_up, Amino, AminoAcid, Direction, Protein,
};
use plotters::prelude::*;
pub fn animation(protein: &mut Protein, dim: u8) {
    let area = BitMapBackend::gif(
        "./animated.gif", // アニメーションファイルの名前。この名前で保存される
        (1200, 800),      //  グラフのサイズ（幅x高さ)
        100,              //  1フレームの時間。単位は [ms]
    )
    .unwrap()
    .into_drawing_area();

    let mut data = Vec::new();
    data.push((0, 0, 0, protein.aminos[0].amino));
    if dim == 2 {
        data.push((5, 0, 0, protein.aminos[1].amino));
    } else {
        data.push((10, 0, 0, protein.aminos[1].amino));
    }
    for i in 0..protein.direct.len() {
        let mut prev_direction = (
            data[i + 1].0 - data[i].0,
            data[i + 1].1 - data[i].1,
            data[i + 1].2 - data[i].2,
        );
        let (x, y, z, a) = match protein.direct[i] {
            Direction::S => (
                data[i + 1].0 + prev_direction.0,
                data[i + 1].1 + prev_direction.1,
                data[i + 1].2 + prev_direction.2,
                protein.aminos[i + 2].amino,
            ),
            Direction::L => (
                data[i + 1].0 + rotate_left(prev_direction).0,
                data[i + 1].1 + rotate_left(prev_direction).1,
                data[i + 1].2 + rotate_left(prev_direction).2,
                protein.aminos[i + 2].amino,
            ),
            Direction::R => (
                data[i + 1].0 + rotate_right(prev_direction).0,
                data[i + 1].1 + rotate_right(prev_direction).1,
                data[i + 1].2 + rotate_right(prev_direction).2,
                protein.aminos[i + 2].amino,
            ),
            Direction::U => (
                data[i + 1].0 + rotate_up(prev_direction).0,
                data[i + 1].1 + rotate_up(prev_direction).1,
                data[i + 1].2 + rotate_up(prev_direction).2,
                protein.aminos[i + 2].amino,
            ),
            Direction::D => (
                data[i + 1].0 + rotate_down(prev_direction).0,
                data[i + 1].1 + rotate_down(prev_direction).1,
                data[i + 1].2 + rotate_down(prev_direction).2,
                protein.aminos[i + 2].amino,
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
