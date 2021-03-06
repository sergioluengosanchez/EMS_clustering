# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

c_hybrid_SEM <- function(data_r, is_directional_r, num_clusters, max_parents, directional_parents) {
    .Call('_SomaSimulation_c_hybrid_SEM', PACKAGE = 'SomaSimulation', data_r, is_directional_r, num_clusters, max_parents, directional_parents)
}

geodesic_distance <- function(input_path, vertex_idx) {
    .Call('_SomaSimulation_geodesic_distance', PACKAGE = 'SomaSimulation', input_path, vertex_idx)
}

poly_volume <- function(x) {
    .Call('_SomaSimulation_poly_volume', PACKAGE = 'SomaSimulation', x)
}

