use lib::{rotate_left, rotate_right, Amino, AminoAcid, Direction, Heuristics, Protein};
use plotters::prelude::*;

fn animation(protein: &mut Protein) {
    let area = BitMapBackend::gif(
        "./animated.gif", // アニメーションファイルの名前。この名前で保存される
        (1200, 800),      //  グラフのサイズ（幅x高さ)
        100,              //  1フレームの時間。単位は [ms]
    )
    .unwrap()
    .into_drawing_area();

    let mut data = Vec::new();
    data.push((0, 0, 0, protein.aminos[0].amino));
    data.push((5, 0, 0, protein.aminos[1].amino));
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
        };
        data.push((x, y, z, a));
    }

    for _ in 0..=20 {
        area.fill(&WHITE).unwrap();

        let mut chart = ChartBuilder::on(&area)
            .margin(20)
            .caption("protein structure", ("sans-serif", 40))
            .build_cartesian_3d(-20..20, -20..20, -20..20)
            .unwrap();

        chart.configure_axes().draw().unwrap();

        for i in 0..data.len() {
            data[i].0 += 1;
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
