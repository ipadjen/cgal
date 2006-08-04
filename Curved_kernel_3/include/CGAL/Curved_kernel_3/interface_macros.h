// Copyright (c) 2000-2004  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source: /CVSROOT/CGAL/Packages/Curved_kernel/include/CGAL/Curved_kernel/interface_macros.h,v $
// $Revision$ $Date$
// $Name:  $
//
// Author(s)     : Herve Bronnimann, Sylvain Pion, Susan Hert

// This file is intentionally not protected against re-inclusion.
// It's aimed at being included from within a kernel traits class, this
// way we share more code.

// It is the responsability of the including file to correctly set the 2
// macros CGAL_Kernel_pred and CGAL_Kernel_cons.
// And they are #undefed at the end of this file.

  CGAL_Spherical_Kernel_cons(Get_equation, get_equation_object) 
  CGAL_Spherical_Kernel_cons(Construct_circular_arc_point_3, construct_circular_arc_point_3_object)
  CGAL_Spherical_Kernel_cons(Construct_sphere_3, construct_sphere_3_object)
  CGAL_Spherical_Kernel_cons(Construct_plane_3, construct_plane_3_object)
  CGAL_Spherical_Kernel_cons(Construct_line_3, construct_line_3_object)
  CGAL_Spherical_Kernel_cons(Construct_circle_3, construct_circle_3_object)
  CGAL_Spherical_Kernel_cons(Construct_diametral_sphere_3, construct_diametral_sphere_3_object)
  CGAL_Spherical_Kernel_cons(Construct_supporting_plane_3, construct_supporting_plane_3_object)
  CGAL_Spherical_Kernel_cons(Compute_circular_x_3, compute_circular_x_3_object)
  CGAL_Spherical_Kernel_cons(Compute_circular_y_3, compute_circular_y_3_object)
  CGAL_Spherical_Kernel_cons(Compute_circular_z_3, compute_circular_z_3_object)
  CGAL_Spherical_Kernel_cons(Intersect_3, intersect_3_object)

  CGAL_Spherical_Kernel_pred(Compare_x_3,  compare_x_3_object)
  CGAL_Spherical_Kernel_pred(Compare_y_3,  compare_y_3_object)
  CGAL_Spherical_Kernel_pred(Compare_z_3,  compare_z_3_object)
  CGAL_Spherical_Kernel_pred(Compare_xy_3,  compare_xy_3_object)
  CGAL_Spherical_Kernel_pred(Compare_xyz_3, compare_xyz_3_object)
  CGAL_Spherical_Kernel_pred(Equal_3,  equal_3_object)
  CGAL_Spherical_Kernel_pred(Has_on_3,  has_on_3_object)

#undef CGAL_Spherical_Kernel_pred
#undef CGAL_Spherical_Kernel_cons
