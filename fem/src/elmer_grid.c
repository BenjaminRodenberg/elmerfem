/*  
   Elmer, A Finite Element Software for Multiphysical Problems
  
   Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
   
   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.
  
   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.
   
   You should have received a copy of the GNU Lesser General Public
   License along with this library (in file ../../LGPL-2.1); if not, write 
   to the Free Software Foundation, Inc., 51 Franklin Street, 
   Fifth Floor, Boston, MA  02110-1301  USA
*/

/* Authors: Moritz Hanke (DKRZ), adapted by Thomas Zwinger
    original code provided by Moritz Hanke (DKRZ) in the frame of the TerraDT project
  Email:   thomas.zwinger@csc.fi
  Web:     http://www.csc.fi/elmer
  Address: CSC - IT Center for Science Ltd.
           Keilaranta 14
           02101 Espoo, Finland 

  Original Date: 10.6.2025
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <assert.h>

#include <proj.h>

#include "elmer_grid.h"

struct glb2loc {
  int global_id, local_id;
};

static void _compute_cell_centers(
  int nbr_cells, int * cell_to_vertex, const int * num_vertices_per_cell,
  const double * x_vertices, const double * y_vertices,
  double * x_cells, double * y_cells);

static int compare_glb2loc_glb (const void * a, const void * b) {
  return ((struct glb2loc*)a)->global_id - ((struct glb2loc*)b)->global_id;
}

void convert2rad(
  double * x_vertices, double * y_vertices, int nbr_vertices) {
  // define transformation
  PJ * P =
    proj_create_crs_to_crs(
      PJ_DEFAULT_CTX, "EPSG:3413", "+proj=longlat +datum=WGS84", NULL);

  if (!P) {
    fputs("failed to create transformation", stderr);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  // transform all vertices
  for (int i = 0; i < nbr_vertices; ++i) {
    PJ_COORD src_coord = proj_coord(x_vertices[i], y_vertices[i], 0, 0);
    PJ_COORD tgt_coord = proj_trans(P, PJ_FWD, src_coord);
    x_vertices[i] = proj_torad(tgt_coord.lp.lam);
    y_vertices[i] = proj_torad(tgt_coord.lp.phi);
  }

  // clean up
  proj_destroy(P);
}

static inline void LLtoXYZ(double lon, double lat, double p_out[]) {

   while (lon < -M_PI) lon += 2.0 * M_PI;
   while (lon >= M_PI) lon -= 2.0 * M_PI;

   double cos_lat = cos(lat);
   p_out[0] = cos_lat * cos(lon);
   p_out[1] = cos_lat * sin(lon);
   p_out[2] = sin(lat);
}

static inline void XYZtoLL (double const p_in[], double * lon, double * lat) {

   *lon = atan2(p_in[1] , p_in[0]);
   *lat = M_PI_2 - acos(p_in[2]);
}

static inline void normalise_vector(double v[]) {

   double norm = 1.0 / sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);

   v[0] *= norm;
   v[1] *= norm;
   v[2] *= norm;
}

void compute_cell_centers(
  int nbr_cells, int * cell_to_vertex, const int * num_vertices_per_cell,
  const double * x_vertices, const double * y_vertices,
  double * x_cells, double * y_cells) {

  _compute_cell_centers(
    nbr_cells, cell_to_vertex, num_vertices_per_cell,
    x_vertices, y_vertices, x_cells, y_cells);
}

static void _compute_cell_centers(
  int nbr_cells, int * cell_to_vertex, const int * num_vertices_per_cell,
  const double * x_vertices, const double * y_vertices,
  double * x_cells, double * y_cells) {

  for (int i = 0; i < nbr_cells; ++i) {

    int num_vertices = num_vertices_per_cell[i];
    double cell_center[3] = {0.0, 0.0, 0.0};

    for (int j = 0; j < num_vertices; ++j) {
      double vertex[3];
      LLtoXYZ(
        x_vertices[cell_to_vertex[j]], y_vertices[cell_to_vertex[j]], vertex);
      cell_center[0] += vertex[0];
      cell_center[1] += vertex[1];
      cell_center[2] += vertex[2];
    }

    normalise_vector(cell_center);
    XYZtoLL(cell_center, x_cells + i, y_cells + i);

    cell_to_vertex += num_vertices;
  }
}

