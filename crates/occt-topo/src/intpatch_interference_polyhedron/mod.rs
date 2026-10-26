//! `IntPatch_InterferencePolyhedron`.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/IntPatch/IntPatch_InterferencePolyhedron.hxx/.cxx`.
//!
//! Computes the interference between two `IntPatch_Polyhedron`s (or the
//! self-interference of one): section points, section lines and tangent zones.
//! The class derives from `Intf_Interference`; the port composes [`IntfInterference`]
//! as `base`, the same choice as the other `Intf_*` ports in `occt-core`.
//!
//! Ported members (OCCT 8.0.0):
//! * both ctors (`cxx:34-96`), `Perform` x2 (`cxx:98-124`)
//! * `Interference(Obje1, Obje2)` (`cxx:128-146`), the BVH candidate walk
//! * `Intersect` (`cxx:151-838`) - see `intersect.rs`
//! * `TangentZoneValue` (`cxx:842-1117`) - see `tangent_zone_value.rs`
//! * `CoupleCharacteristics` (`cxx:1121-1178`)
//!
//! Translation notes:
//! * `Interference(const IntPatch_Polyhedron&)` (`cxx:127`, empty body) is not
//!   translated; the self-interference ctor calls the two-polyhedron overload
//!   with the same polyhedron twice, exactly as `cxx:90`.
//! * `memset` zeroing of `OI` / `TI` / `dpOeT` / `dpOpT` / `deOpT`
//!   (`cxx:40-45`) becomes the zeroed arrays of [`Default`].
//! * `TangentZoneValue` calls `TheTZ.Insert(...)`, which `Intf_TangentZone`
//!   implements as a hard `false` (`Intf_TangentZone.cxx:79-118`); the port
//!   keeps that call so the `piToInsert` bookkeeping below it stays
//!   transliterated rather than deleted.
#![allow(dead_code)]

mod intersect;
mod tangent_zone_value;

use occt_core::gp::GpXyz;
use occt_core::intf::{IntfInterference, IntfTangentZone};
use occt_core::precision::epsilon;

use crate::intpatch_polyhedron::IntPatchPolyhedron;
use crate::intpatch_polyhedron_bvh::{IntPatchBVHTraversal, IntPatchPolyhedronBVH};

/// `static const int Pourcent3[9]` (`cxx:33`).
pub(crate) const POURCENT3: [usize; 9] = [0, 1, 2, 0, 1, 2, 0, 1, 2];

/// `IntPatch_InterferencePolyhedron` (`hxx:31-91`).
#[derive(Clone, Debug)]
pub struct IntPatchInterferencePolyhedron {
    /// The `Intf_Interference` base sub-object (`mySPoins` / `mySLines` /
    /// `myTZones` / `SelfIntf` / `Tolerance`).
    pub base: IntfInterference,
    /// `OI[3]` (`hxx:82`).
    oi: [i32; 3],
    /// `TI[3]` (`hxx:83`).
    ti: [i32; 3],
    /// `dpOpT[3][3]` (`hxx:84`).
    dp_op_t: [[f64; 3]; 3],
    /// `dpOeT[3][3]` (`hxx:85`).
    dp_oe_t: [[f64; 3]; 3],
    /// `deOpT[3][3]` (`hxx:86`).
    de_op_t: [[f64; 3]; 3],
    /// `voo[3]` (`hxx:87`).
    voo: [GpXyz; 3],
    /// `vtt[3]` (`hxx:88`).
    vtt: [GpXyz; 3],
    /// `Incidence` (`hxx:89`).
    incidence: f64,
}

impl Default for IntPatchInterferencePolyhedron {
    /// `IntPatch_InterferencePolyhedron()` (`cxx:34-47`).
    fn default() -> Self {
        Self {
            base: IntfInterference::new(false),
            oi: [0; 3],
            ti: [0; 3],
            dp_op_t: [[0.0; 3]; 3],
            dp_oe_t: [[0.0; 3]; 3],
            de_op_t: [[0.0; 3]; 3],
            voo: [GpXyz::new(0.0, 0.0, 0.0); 3],
            vtt: [GpXyz::new(0.0, 0.0, 0.0); 3],
            incidence: 0.0,
        }
    }
}

impl IntPatchInterferencePolyhedron {
    /// `IntPatch_InterferencePolyhedron()` (`cxx:34-47`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `IntPatch_InterferencePolyhedron(FirstPol, SeconPol)` (`cxx:51-73`).
    pub fn with_polyhedra(
        first_pol: &IntPatchPolyhedron,
        secon_pol: &IntPatchPolyhedron,
    ) -> Self {
        let mut me = Self::default();
        if !first_pol
            .bounding()
            .is_out_box(secon_pol.bounding())
        {
            me.base.tolerance = first_pol.deflection_over_estimation()
                + secon_pol.deflection_over_estimation();
            if me.base.tolerance == 0.0 {
                me.base.tolerance = epsilon(1000.0);
            }
            me.interference(first_pol, secon_pol);
        }
        me
    }

    /// `IntPatch_InterferencePolyhedron(Objet)` (`cxx:77-91`).
    pub fn with_self(objet: &IntPatchPolyhedron) -> Self {
        let mut me = Self {
            base: IntfInterference::new(true),
            ..Self::default()
        };
        me.base.tolerance = objet.deflection_over_estimation() * 2.0;
        if me.base.tolerance == 0.0 {
            me.base.tolerance = epsilon(1000.0);
        }
        // `Interference(Objet, Objet)` (`cxx:90`).
        me.interference(objet, objet);
        me
    }

    /// `Perform(FirstPol, SeconPol)` (`cxx:98-110`).
    pub fn perform_polyhedra(
        &mut self,
        first_pol: &IntPatchPolyhedron,
        secon_pol: &IntPatchPolyhedron,
    ) {
        self.base.self_interference(false);
        if !first_pol
            .bounding()
            .is_out_box(secon_pol.bounding())
        {
            self.base.tolerance = first_pol.deflection_over_estimation()
                + secon_pol.deflection_over_estimation();
            if self.base.tolerance == 0.0 {
                self.base.tolerance = epsilon(1000.0);
            }
            self.interference(first_pol, secon_pol);
        }
    }

    /// `Perform(Objet)` (`cxx:114-124`).
    pub fn perform_self(&mut self, objet: &IntPatchPolyhedron) {
        self.base.self_interference(true);
        self.base.tolerance = objet.deflection_over_estimation() * 2.0;
        if self.base.tolerance == 0.0 {
            self.base.tolerance = epsilon(1000.0);
        }
        self.interference(objet, objet);
    }

    /// `Interference(theFirstPol, theSecondPol)` (`cxx:128-146`).
    fn interference(&mut self, first_pol: &IntPatchPolyhedron, second_pol: &IntPatchPolyhedron) {
        // Build BVH sets for both polyhedra
        let mut set1 = IntPatchPolyhedronBVH::with_poly(first_pol);
        let mut set2 = IntPatchPolyhedronBVH::with_poly(second_pol);

        // Find candidate triangle pairs via BVH traversal
        let mut traversal = IntPatchBVHTraversal::new();
        traversal.perform(&mut set1, &mut set2, self.base.self_intf);

        // Process each candidate pair
        for pair in traversal.pairs() {
            self.intersect(pair.first, first_pol, pair.second, second_pol);
        }
    }

    /// `CoupleCharacteristics(FirstPol, SeconPol)` (`cxx:1121-1178`).
    fn couple_characteristics(
        &mut self,
        first_pol: &IntPatchPolyhedron,
        secon_pol: &IntPatchPolyhedron,
    ) {
        for n1 in 0..3 {
            let n2 = POURCENT3[n1 + 1];
            self.voo[n1] = first_pol
                .point(self.oi[n2])
                .xyz()
                .subtracted(&first_pol.point(self.oi[n1]).xyz());
            self.vtt[n1] = secon_pol
                .point(self.ti[n2])
                .xyz()
                .subtracted(&secon_pol.point(self.ti[n1]).xyz());
        }

        let mut vvec = self.voo[0]
            .crossed(&self.voo[1])
            .added(&self.voo[1].crossed(&self.voo[2]))
            .added(&self.voo[2].crossed(&self.voo[0]));
        let vnor_t = self.vtt[0]
            .crossed(&self.vtt[1])
            .added(&self.vtt[1].crossed(&self.vtt[2]))
            .added(&self.vtt[2].crossed(&self.vtt[0]));
        if vnor_t.modulus() > vvec.modulus() {
            vvec = vnor_t;
        }

        for n1 in 0..3 {
            for n2 in 0..3 {
                let vto = first_pol
                    .point(self.oi[n1])
                    .xyz()
                    .subtracted(&secon_pol.point(self.ti[n2]).xyz());
                self.dp_op_t[n1][n2] = vto.modulus();

                let mut lg = self.vtt[n2].modulus();
                if lg > 1e-16 {
                    //-- RealEpsilon()
                    let vv = vto.crossed(&self.vtt[n2]);
                    lg = if vvec.dot(&vv) > 0.0 { lg } else { -lg };
                    self.dp_oe_t[n1][n2] = vv.modulus() / lg;
                } else {
                    self.dp_oe_t[n1][n2] = self.dp_op_t[n1][n2];
                }

                lg = self.voo[n1].modulus();
                if lg > 1e-16 {
                    //-- RealEpsilon()
                    let vv = vto.crossed(&self.voo[n1]);
                    lg = if vvec.dot(&vv) > 0.0 { -lg } else { lg };
                    self.de_op_t[n1][n2] = vv.modulus() / lg;
                } else {
                    self.de_op_t[n1][n2] = self.dp_op_t[n1][n2];
                }
            }
        }
    }

    /// One entry point for the [`IntfTangentZone`] returned by
    /// [`Self::tangent_zone_value`], so the caller stays a transliteration of
    /// `cxx:239-242`.
    fn store_tangent_zone(&mut self, the_tz: &IntfTangentZone) {
        if !self.base.insert_zone(the_tz) {
            self.base.t_zones.push(the_tz.clone());
        }
    }
}
