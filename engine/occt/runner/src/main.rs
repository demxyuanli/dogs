//! Pure Rust runner — calls occt-core crate directly.
use std::io::{self, BufRead, Write};
use occt_core::*;
// Method name aliases: GpXyz::multiplied = multiply_scalar, GpMat::m = data

fn main() {
    let stdin = io::stdin();
    let mut out = io::stdout().lock();
    for line in stdin.lock().lines() {
        let line = line.unwrap();
        let mut it = line.split_whitespace();
        let op = it.next().unwrap_or("");
        let args: Vec<f64> = it.filter_map(|t| t.parse().ok()).collect();
        writeln!(out, "{}", dispatch(op, &args)).unwrap();
    }
}

fn dispatch(op: &str, a: &[f64]) -> String {
    if a.is_empty() && op != "echo" { return "err no args".into(); }
    match op {
        "echo" => format!("ok {}", a[0] as i32),

        "gp_XYZ::Dot" => {let v=GpXyz::new(a[0],a[1],a[2]); format!("ok {}", v.dot(&GpXyz::new(a[3],a[4],a[5])))}
        "gp_XYZ::Crossed" => {let v=GpXyz::new(a[0],a[1],a[2]).crossed(&GpXyz::new(a[3],a[4],a[5])); format!("ok {} {} {}", v.x, v.y, v.z)}
        "gp_XYZ::Modulus" => format!("ok {}", GpXyz::new(a[0],a[1],a[2]).modulus()),
        "gp_XYZ::Added" => {let v=GpXyz::new(a[0],a[1],a[2]).added(&GpXyz::new(a[3],a[4],a[5])); format!("ok {} {} {}", v.x, v.y, v.z)}
        "gp_XYZ::MultipliedScalar" => {let v=GpXyz::new(a[0],a[1],a[2]).multiplied(a[3]); format!("ok {} {} {}", v.x, v.y, v.z)}
        "gp_XYZ::CrossSquareMagnitude" => {let v=GpXyz::new(a[0],a[1],a[2]); format!("ok {}", v.cross_square_magnitude(&GpXyz::new(a[3],a[4],a[5])))}
        "gp_XYZ::DotCross" => {let v=GpXyz::new(a[0],a[1],a[2]); format!("ok {}", v.dot_cross(&GpXyz::new(a[3],a[4],a[5]), &GpXyz::new(a[6],a[7],a[8])))}
        "gp_XYZ::MultipliedMat" => {
            let v = GpXyz::new(a[0],a[1],a[2]);
            let m = GpMat::new(a[3],a[4],a[5],a[6],a[7],a[8],a[9],a[10],a[11]);
            let r = v.multiplied_mat(&m);
            format!("ok {} {} {}", r.x, r.y, r.z)
        }

        "gp_Mat::Multiply" => {
            let ma = GpMat::new(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8]);
            let mb = GpMat::new(a[9],a[10],a[11],a[12],a[13],a[14],a[15],a[16],a[17]);
            let r = ma.multiply(&mb);
            format!("ok {} {} {} {} {} {} {} {} {}", r.m[0][0],r.m[0][1],r.m[0][2], r.m[1][0],r.m[1][1],r.m[1][2], r.m[2][0],r.m[2][1],r.m[2][2])
        }
        "gp_Mat::Determinant" => format!("ok {}", GpMat::new(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8]).determinant()),
        "gp_Mat::Transposed" => {
            let r = GpMat::new(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8]).transpose();
            format!("ok {} {} {} {} {} {} {} {} {}", r.m[0][0],r.m[0][1],r.m[0][2], r.m[1][0],r.m[1][1],r.m[1][2], r.m[2][0],r.m[2][1],r.m[2][2])
        }

        "gp_Vec::Dot" => format!("ok {}", GpVec::new(a[0],a[1],a[2]).dot(&GpVec::new(a[3],a[4],a[5]))),
        "gp_Vec::Crossed" => {let v=GpVec::new(a[0],a[1],a[2]).crossed(&GpVec::new(a[3],a[4],a[5])); format!("ok {} {} {}", v.x(), v.y(), v.z())}
        "gp_Vec::Magnitude" => format!("ok {}", GpVec::new(a[0],a[1],a[2]).magnitude()),

        "gp_Dir::Crossed" => {
            match GpDir::new(a[0],a[1],a[2]).and_then(|d1| {
                GpDir::new(a[3],a[4],a[5]).and_then(|d2| d1.crossed(&d2))
            }) {
                Ok(r) => format!("ok {} {} {}", r.x(), r.y(), r.z()),
                Err(_) => "err zero cross".into(),
            }
        }

        // ---- 2D ops ----
        "gp_XY::Dot" => format!("ok {}", GpXY::new(a[0],a[1]).dot(&GpXY::new(a[2],a[3]))),
        "gp_XY::Crossed" => format!("ok {}", GpXY::new(a[0],a[1]).crossed(&GpXY::new(a[2],a[3]))),
        "gp_XY::Modulus" => format!("ok {}", GpXY::new(a[0],a[1]).modulus()),
        "gp_XY::Added" => {let v=GpXY::new(a[0],a[1]).added(&GpXY::new(a[2],a[3])); format!("ok {} {}", v.x, v.y)},
        "gp_XY::MultipliedScalar" => {let v=GpXY::new(a[0],a[1]).multiplied(a[2]); format!("ok {} {}", v.x, v.y)},
        "gp_XY::Normalized" => match GpXY::new(a[0],a[1]).normalized() { Ok(v)=>format!("ok {} {}",v.x,v.y), Err(_)=>"err zero".into() },
        "gp_Vec2d::Dot" => format!("ok {}", GpVec2d::new(a[0],a[1]).dot(&GpVec2d::new(a[2],a[3]))),
        "gp_Vec2d::Magnitude" => format!("ok {}", GpVec2d::new(a[0],a[1]).magnitude()),
        "gp_Vec2d::Crossed" => format!("ok {}", GpVec2d::new(a[0],a[1]).crossed(&GpVec2d::new(a[2],a[3]))),

        "gp_Dir2d::new" => match GpDir2d::new(a[0],a[1]) { Ok(d)=>format!("ok {} {}",d.x,d.y), Err(_)=>"err zero".into() },
        "gp_Dir2d::Dot" => format!("ok {}", GpDir2d::new(a[0],a[1]).unwrap_or_default().dot(&GpDir2d::new(a[2],a[3]).unwrap_or_default())),

        _ => "err unknown op".into(),
    }
}
