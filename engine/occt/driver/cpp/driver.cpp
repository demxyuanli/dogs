// OCCT reference driver — golden values for differential testing.
// Inline OCCT ops + linked Standard_Failure.cxx + Standard.cxx for exception support.
#include <cmath>
#include <gp_XYZ.hxx>
#include <gp_Mat.hxx>
#include <gp_Vec.hxx>
#include <gp_Dir.hxx>
#include <gp_XY.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Mat2d.hxx>
#include <gp_Vec2d.hxx>
#include <gp_Dir2d.hxx>
#include <gp_Trsf2d.hxx>
#include <Precision.hxx>

extern "C" {

int echo(int x) { return x; }

// ---- gp_XYZ (all inline) ----

int gp_XYZ_Dot(double ax, double ay, double az, double bx, double by, double bz, double* out) {
    gp_XYZ a(ax, ay, az), b(bx, by, bz);
    *out = a.Dot(b);
    return 0;
}

int gp_XYZ_Crossed(double ax, double ay, double az, double bx, double by, double bz,
                   double* ox, double* oy, double* oz) {
    gp_XYZ a(ax, ay, az), b(bx, by, bz);
    gp_XYZ c = a.Crossed(b);
    *ox = c.X(); *oy = c.Y(); *oz = c.Z();
    return 0;
}

int gp_XYZ_Modulus(double x, double y, double z, double* out) {
    gp_XYZ a(x, y, z);
    *out = a.Modulus();
    return 0;
}

int gp_XYZ_Added(double ax, double ay, double az, double bx, double by, double bz,
                 double* ox, double* oy, double* oz) {
    gp_XYZ a(ax, ay, az), b(bx, by, bz);
    gp_XYZ c = a.Added(b);
    *ox = c.X(); *oy = c.Y(); *oz = c.Z();
    return 0;
}

int gp_XYZ_MultipliedScalar(double x, double y, double z, double s,
                             double* ox, double* oy, double* oz) {
    gp_XYZ a(x, y, z);
    gp_XYZ c = a.Multiplied(s);
    *ox = c.X(); *oy = c.Y(); *oz = c.Z();
    return 0;
}

int gp_XYZ_CrossSquareMagnitude(double ax, double ay, double az,
                                 double bx, double by, double bz, double* out) {
    gp_XYZ a(ax, ay, az), b(bx, by, bz);
    *out = a.CrossSquareMagnitude(b);
    return 0;
}

int gp_XYZ_DotCross(double ax, double ay, double az,
                     double bx, double by, double bz,
                     double cx, double cy, double cz, double* out) {
    gp_XYZ a(ax, ay, az), b(bx, by, bz), c(cx, cy, cz);
    *out = a.DotCross(b, c);
    return 0;
}

int gp_XYZ_MultipliedMat(double x, double y, double z,
                          double m00,double m01,double m02,
                          double m10,double m11,double m12,
                          double m20,double m21,double m22,
                          double* ox, double* oy, double* oz) {
    gp_Mat m(m00,m01,m02,m10,m11,m12,m20,m21,m22);
    gp_XYZ a(x,y,z);
    gp_XYZ r = a.Multiplied(m);
    *ox = r.X(); *oy = r.Y(); *oz = r.Z();
    return 0;
}

// ---- gp_Mat (inline only) ----

int gp_Mat_Multiply(double a11,double a12,double a13,double a21,double a22,double a23,double a31,double a32,double a33,
                    double b11,double b12,double b13,double b21,double b22,double b23,double b31,double b32,double b33,
                    double* r00,double* r01,double* r02,
                    double* r10,double* r11,double* r12,
                    double* r20,double* r21,double* r22) {
    gp_Mat A(a11,a12,a13,a21,a22,a23,a31,a32,a33);
    gp_Mat B(b11,b12,b13,b21,b22,b23,b31,b32,b33);
    gp_Mat C = A.Multiplied(B);
    *r00=C(1,1); *r01=C(1,2); *r02=C(1,3);
    *r10=C(2,1); *r11=C(2,2); *r12=C(2,3);
    *r20=C(3,1); *r21=C(3,2); *r22=C(3,3);
    return 0;
}

int gp_Mat_Determinant(double a11,double a12,double a13,
                        double a21,double a22,double a23,
                        double a31,double a32,double a33, double* out) {
    gp_Mat m(a11,a12,a13,a21,a22,a23,a31,a32,a33);
    *out = m.Determinant();
    return 0;
}

int gp_Mat_Transposed(double a11,double a12,double a13,
                       double a21,double a22,double a23,
                       double a31,double a32,double a33,
                       double* r00,double* r01,double* r02,
                       double* r10,double* r11,double* r12,
                       double* r20,double* r21,double* r22) {
    gp_Mat m(a11,a12,a13,a21,a22,a23,a31,a32,a33);
    gp_Mat t = m.Transposed();
    *r00=t(1,1); *r01=t(1,2); *r02=t(1,3);
    *r10=t(2,1); *r11=t(2,2); *r12=t(2,3);
    *r20=t(3,1); *r21=t(3,2); *r22=t(3,3);
    return 0;
}

// ---- gp_Vec (inline, delegates to gp_XYZ) ----

int gp_Vec_Dot(double ax, double ay, double az, double bx, double by, double bz, double* out) {
    gp_Vec a(ax, ay, az), b(bx, by, bz);
    *out = a.Dot(b);
    return 0;
}

int gp_Vec_Crossed(double ax, double ay, double az, double bx, double by, double bz,
                   double* ox, double* oy, double* oz) {
    gp_Vec a(ax, ay, az), b(bx, by, bz);
    gp_Vec c = a.Crossed(b);
    *ox = c.X(); *oy = c.Y(); *oz = c.Z();
    return 0;
}

int gp_Vec_Magnitude(double x, double y, double z, double* out) {
    gp_Vec a(x, y, z);
    *out = a.Magnitude();
    return 0;
}

// ---- gp_Dir (inline Crossed only — constructors may throw which needs Standard_Failure) ----

int gp_Dir_Crossed(double ax, double ay, double az, double bx, double by, double bz,
                   double* ox, double* oy, double* oz) {
    gp_Dir a(ax, ay, az), b(bx, by, bz);
    gp_Dir c = a.Crossed(b);
    *ox = c.X(); *oy = c.Y(); *oz = c.Z();
    return 0;
}

// ---- gp_XY (all inline) ----

int gp_XY_Dot(double ax, double ay, double bx, double by, double* out) {
    gp_XY a(ax, ay), b(bx, by); *out = a.Dot(b); return 0;
}

int gp_XY_Crossed(double ax, double ay, double bx, double by, double* out) {
    gp_XY a(ax, ay), b(bx, by); *out = a.Crossed(b); return 0;
}

int gp_XY_Modulus(double x, double y, double* out) {
    *out = gp_XY(x, y).Modulus(); return 0;
}

int gp_XY_Added(double ax, double ay, double bx, double by, double* ox, double* oy) {
    gp_XY r = gp_XY(ax, ay).Added(gp_XY(bx, by));
    *ox = r.X(); *oy = r.Y(); return 0;
}

int gp_XY_MultipliedScalar(double x, double y, double s, double* ox, double* oy) {
    gp_XY r = gp_XY(x, y).Multiplied(s);
    *ox = r.X(); *oy = r.Y(); return 0;
}

int gp_XY_Normalized(double x, double y, double* ox, double* oy, int* err) {
    gp_XY a(x, y);
    double d = a.Modulus();
    if (d <= gp::Resolution()) { *err = 1; return 0; }
    gp_XY n = a.Normalized();
    *ox = n.X(); *oy = n.Y(); *err = 0; return 0;
}

// gp_XY_MultipliedMat2d deferred — gp_Mat2d constructor non-inline

// ---- gp_Mat2d (deferred — constructors are non-inline, need OCCT .lib) ----

// ---- gp_Vec2d (inline, delegates to gp_XY) ----

int gp_Vec2d_Dot(double ax, double ay, double bx, double by, double* out) {
    gp_Vec2d a(ax, ay), b(bx, by); *out = a.Dot(b); return 0;
}

int gp_Vec2d_Magnitude(double x, double y, double* out) {
    *out = gp_Vec2d(x, y).Magnitude(); return 0;
}

int gp_Vec2d_Crossed(double ax, double ay, double bx, double by, double* out) {
    *out = gp_Vec2d(ax, ay).Crossed(gp_Vec2d(bx, by)); return 0;
}

// ---- gp_Dir2d ----

int gp_Dir2d_New(double x, double y, double* ox, double* oy, int* err) {
    double sq = x*x + y*y;
    if (sq <= gp::Resolution() * gp::Resolution()) { *err = 1; return 0; }
    gp_Dir2d d(x, y);
    *ox = d.X(); *oy = d.Y(); *err = 0; return 0;
}

int gp_Dir2d_Dot(double ax, double ay, double bx, double by, double* out) {
    gp_Dir2d a(ax, ay), b(bx, by); *out = a.Dot(b); return 0;
}

// ---- gp_Trsf2d (deferred — SetRotation needs gp_Pnt2d non-inline) ----

} // extern "C"
