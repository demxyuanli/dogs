//! Candidate search and distance evaluation of `BRepBuilderAPI_Sewing`:
//! `EvaluateAngulars` (`cxx:1189-1280`), `EvaluateDistances`
//! (`cxx:1287-1522`), `AnalysisNearestEdges` (`cxx:1618-1782`),
//! `FindCandidates` (`cxx:1785-2086`) and `MergedNearestEdges`
//! (`cxx:4242-4431`).

mod abscissa;
mod angulars;
mod distances;
mod find_candidates;
mod helpers;
mod nearest_edges;
