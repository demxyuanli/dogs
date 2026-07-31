//! OCCT driver — calls C++ reference via FFI. Inline-only ops.
use std::io::{self, BufRead, Write};

extern "C" {
    fn echo(x: i32) -> i32;

    fn gp_XYZ_Dot(ax:f64,ay:f64,az:f64,bx:f64,by:f64,bz:f64, out:*mut f64) -> i32;
    fn gp_XYZ_Crossed(ax:f64,ay:f64,az:f64,bx:f64,by:f64,bz:f64,ox:*mut f64,oy:*mut f64,oz:*mut f64) -> i32;
    fn gp_XYZ_Modulus(x:f64,y:f64,z:f64, out:*mut f64) -> i32;
    fn gp_XYZ_Added(ax:f64,ay:f64,az:f64,bx:f64,by:f64,bz:f64,ox:*mut f64,oy:*mut f64,oz:*mut f64) -> i32;
    fn gp_XYZ_MultipliedScalar(x:f64,y:f64,z:f64,s:f64,ox:*mut f64,oy:*mut f64,oz:*mut f64) -> i32;
    fn gp_XYZ_CrossSquareMagnitude(ax:f64,ay:f64,az:f64,bx:f64,by:f64,bz:f64,out:*mut f64) -> i32;
    fn gp_XYZ_DotCross(ax:f64,ay:f64,az:f64,bx:f64,by:f64,bz:f64,cx:f64,cy:f64,cz:f64,out:*mut f64) -> i32;
    fn gp_XYZ_MultipliedMat(x:f64,y:f64,z:f64,m00:f64,m01:f64,m02:f64,m10:f64,m11:f64,m12:f64,m20:f64,m21:f64,m22:f64,ox:*mut f64,oy:*mut f64,oz:*mut f64) -> i32;

    fn gp_Mat_Multiply(a11:f64,a12:f64,a13:f64,a21:f64,a22:f64,a23:f64,a31:f64,a32:f64,a33:f64,
        b11:f64,b12:f64,b13:f64,b21:f64,b22:f64,b23:f64,b31:f64,b32:f64,b33:f64,
        r00:*mut f64,r01:*mut f64,r02:*mut f64,r10:*mut f64,r11:*mut f64,r12:*mut f64,r20:*mut f64,r21:*mut f64,r22:*mut f64) -> i32;
    fn gp_Mat_Determinant(a11:f64,a12:f64,a13:f64,a21:f64,a22:f64,a23:f64,a31:f64,a32:f64,a33:f64,out:*mut f64) -> i32;
    fn gp_Mat_Transposed(a11:f64,a12:f64,a13:f64,a21:f64,a22:f64,a23:f64,a31:f64,a32:f64,a33:f64,
        r00:*mut f64,r01:*mut f64,r02:*mut f64,r10:*mut f64,r11:*mut f64,r12:*mut f64,r20:*mut f64,r21:*mut f64,r22:*mut f64) -> i32;

    fn gp_Vec_Dot(ax:f64,ay:f64,az:f64,bx:f64,by:f64,bz:f64,out:*mut f64) -> i32;
    fn gp_Vec_Crossed(ax:f64,ay:f64,az:f64,bx:f64,by:f64,bz:f64,ox:*mut f64,oy:*mut f64,oz:*mut f64) -> i32;
    fn gp_Vec_Magnitude(x:f64,y:f64,z:f64,out:*mut f64) -> i32;

    fn gp_Dir_Crossed(ax:f64,ay:f64,az:f64,bx:f64,by:f64,bz:f64,ox:*mut f64,oy:*mut f64,oz:*mut f64) -> i32;

    // 2D ops
    fn gp_XY_Dot(ax:f64,ay:f64,bx:f64,by:f64,out:*mut f64) -> i32;
    fn gp_XY_Crossed(ax:f64,ay:f64,bx:f64,by:f64,out:*mut f64) -> i32;
    fn gp_XY_Modulus(x:f64,y:f64,out:*mut f64) -> i32;
    fn gp_XY_Added(ax:f64,ay:f64,bx:f64,by:f64,ox:*mut f64,oy:*mut f64) -> i32;
    fn gp_XY_MultipliedScalar(x:f64,y:f64,s:f64,ox:*mut f64,oy:*mut f64) -> i32;
    fn gp_XY_Normalized(x:f64,y:f64,ox:*mut f64,oy:*mut f64,err:*mut i32) -> i32;
    fn gp_Vec2d_Dot(ax:f64,ay:f64,bx:f64,by:f64,out:*mut f64) -> i32;
    fn gp_Vec2d_Magnitude(x:f64,y:f64,out:*mut f64) -> i32;
    fn gp_Vec2d_Crossed(ax:f64,ay:f64,bx:f64,by:f64,out:*mut f64) -> i32;

    fn gp_Dir2d_New(x:f64,y:f64,ox:*mut f64,oy:*mut f64,err:*mut i32) -> i32;
    fn gp_Dir2d_Dot(ax:f64,ay:f64,bx:f64,by:f64,out:*mut f64) -> i32;
}

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
        "echo" => format!("ok {}", unsafe { echo(a[0] as i32) }),

        "gp_XYZ::Dot" => {let mut r=0.0;unsafe{gp_XYZ_Dot(a[0],a[1],a[2],a[3],a[4],a[5],&mut r)};format!("ok {r}")}
        "gp_XYZ::Crossed" => {let(mut x,mut y,mut z)=(0.0,0.0,0.0);unsafe{gp_XYZ_Crossed(a[0],a[1],a[2],a[3],a[4],a[5],&mut x,&mut y,&mut z)};format!("ok {x} {y} {z}")}
        "gp_XYZ::Modulus" => {let mut r=0.0;unsafe{gp_XYZ_Modulus(a[0],a[1],a[2],&mut r)};format!("ok {r}")}
        "gp_XYZ::Added" => {let(mut x,mut y,mut z)=(0.0,0.0,0.0);unsafe{gp_XYZ_Added(a[0],a[1],a[2],a[3],a[4],a[5],&mut x,&mut y,&mut z)};format!("ok {x} {y} {z}")}
        "gp_XYZ::MultipliedScalar" => {let(mut x,mut y,mut z)=(0.0,0.0,0.0);unsafe{gp_XYZ_MultipliedScalar(a[0],a[1],a[2],a[3],&mut x,&mut y,&mut z)};format!("ok {x} {y} {z}")}
        "gp_XYZ::CrossSquareMagnitude" => {let mut r=0.0;unsafe{gp_XYZ_CrossSquareMagnitude(a[0],a[1],a[2],a[3],a[4],a[5],&mut r)};format!("ok {r}")}
        "gp_XYZ::DotCross" => {let mut r=0.0;unsafe{gp_XYZ_DotCross(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],&mut r)};format!("ok {r}")}
        "gp_XYZ::MultipliedMat" => {let(mut x,mut y,mut z)=(0.0,0.0,0.0);unsafe{gp_XYZ_MultipliedMat(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],a[10],a[11],&mut x,&mut y,&mut z)};format!("ok {x} {y} {z}")}

        "gp_Mat::Multiply" => {
            let(mut r00,mut r01,mut r02,mut r10,mut r11,mut r12,mut r20,mut r21,mut r22)=(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
            unsafe{gp_Mat_Multiply(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],a[10],a[11],a[12],a[13],a[14],a[15],a[16],a[17],
                &mut r00,&mut r01,&mut r02,&mut r10,&mut r11,&mut r12,&mut r20,&mut r21,&mut r22)};
            format!("ok {r00} {r01} {r02} {r10} {r11} {r12} {r20} {r21} {r22}")
        }
        "gp_Mat::Determinant" => {let mut r=0.0;unsafe{gp_Mat_Determinant(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],&mut r)};format!("ok {r}")}
        "gp_Mat::Transposed" => {
            let(mut r00,mut r01,mut r02,mut r10,mut r11,mut r12,mut r20,mut r21,mut r22)=(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
            unsafe{gp_Mat_Transposed(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],
                &mut r00,&mut r01,&mut r02,&mut r10,&mut r11,&mut r12,&mut r20,&mut r21,&mut r22)};
            format!("ok {r00} {r01} {r02} {r10} {r11} {r12} {r20} {r21} {r22}")
        }

        "gp_Vec::Dot" => {let mut r=0.0;unsafe{gp_Vec_Dot(a[0],a[1],a[2],a[3],a[4],a[5],&mut r)};format!("ok {r}")}
        "gp_Vec::Crossed" => {let(mut x,mut y,mut z)=(0.0,0.0,0.0);unsafe{gp_Vec_Crossed(a[0],a[1],a[2],a[3],a[4],a[5],&mut x,&mut y,&mut z)};format!("ok {x} {y} {z}")}
        "gp_Vec::Magnitude" => {let mut r=0.0;unsafe{gp_Vec_Magnitude(a[0],a[1],a[2],&mut r)};format!("ok {r}")}

        "gp_Dir::Crossed" => {let(mut x,mut y,mut z)=(0.0,0.0,0.0);unsafe{gp_Dir_Crossed(a[0],a[1],a[2],a[3],a[4],a[5],&mut x,&mut y,&mut z)};format!("ok {x} {y} {z}")}

        // ---- 2D ops ----
        "gp_XY::Dot" => {let mut r=0.0;unsafe{gp_XY_Dot(a[0],a[1],a[2],a[3],&mut r)};format!("ok {r}")}
        "gp_XY::Crossed" => {let mut r=0.0;unsafe{gp_XY_Crossed(a[0],a[1],a[2],a[3],&mut r)};format!("ok {r}")}
        "gp_XY::Modulus" => {let mut r=0.0;unsafe{gp_XY_Modulus(a[0],a[1],&mut r)};format!("ok {r}")}
        "gp_XY::Added" => {let(mut x,mut y)=(0.0,0.0);unsafe{gp_XY_Added(a[0],a[1],a[2],a[3],&mut x,&mut y)};format!("ok {x} {y}")}
        "gp_XY::MultipliedScalar" => {let(mut x,mut y)=(0.0,0.0);unsafe{gp_XY_MultipliedScalar(a[0],a[1],a[2],&mut x,&mut y)};format!("ok {x} {y}")}
        "gp_XY::Normalized" => {let(mut x,mut y,mut e)=(0.0,0.0,0i32);unsafe{gp_XY_Normalized(a[0],a[1],&mut x,&mut y,&mut e)};if e!=0{"err zero".into()}else{format!("ok {x} {y}")}}
        "gp_Vec2d::Dot" => {let mut r=0.0;unsafe{gp_Vec2d_Dot(a[0],a[1],a[2],a[3],&mut r)};format!("ok {r}")}
        "gp_Vec2d::Magnitude" => {let mut r=0.0;unsafe{gp_Vec2d_Magnitude(a[0],a[1],&mut r)};format!("ok {r}")}
        "gp_Vec2d::Crossed" => {let mut r=0.0;unsafe{gp_Vec2d_Crossed(a[0],a[1],a[2],a[3],&mut r)};format!("ok {r}")}

        "gp_Dir2d::new" => {let(mut x,mut y,mut e)=(0.0,0.0,0i32);unsafe{gp_Dir2d_New(a[0],a[1],&mut x,&mut y,&mut e)};if e!=0{"err zero".into()}else{format!("ok {x} {y}")}}
        "gp_Dir2d::Dot" => {let mut r=0.0;unsafe{gp_Dir2d_Dot(a[0],a[1],a[2],a[3],&mut r)};format!("ok {r}")}

        _ => "err unknown op".into(),
    }
}
