/* Copyright (c) 2012
   Julian Pfeifle, Universitat Politecnica de Catalunya
   julian.pfeifle@upc.edu

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version: http://www.gnu.org/licenses/gpl.txt.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
--------------------------------------------------------------------------------
*/

#include "polymake/client.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/linalg.h"
#include "polymake/Array.h"
#include "polymake/Set.h"
#include "polymake/polytope/sympol_interface.h"
#include "polymake/group/permlib.h"
#include "permlib/permlib_api.h"
#include <sstream>
#include <iostream>
#include <fstream>

namespace polymake { namespace polytope {

typedef Set< Set<int> > SimplexSet;
typedef Array< Set<int> > SimplexArray;
typedef Array< Set<int> > Orbits;

void write_set(const Set<int>& ridge, const SimplexSet& simplex_set, const std::string& dir)
{
   std::ostringstream filename;
   filename << dir << "/dim" << ridge.size();
   for (Entire<Set<int> >::const_iterator rit = entire(ridge); !rit.at_end(); ++rit)
      filename << "_" << *rit;
   std::ofstream outfile(filename.str().c_str(), std::ios_base::trunc);
   for (Entire<SimplexSet>::const_iterator sit = entire(simplex_set); !sit.at_end(); ++sit) {
      outfile << "{";
      int ct(0);
      for (Entire<Set<int> >::const_iterator s = entire(*sit); !s.at_end(); ++s, ++ct) {
         if (ct) outfile << " ";
         outfile << *s;
      }
      outfile << "}" << endl;
   }
}

void write_simplex_reps_from_ridge(perl::Object p, Set<int> ridge, perl::OptionSet options)
{
   const Matrix<Rational> V = p.give("VERTICES");
   const Array<Array<int> > gens = p.give(options["group_generators"]);
   const group::PermlibGroup sym_group(gens);
   Set<int> simplex;
   SimplexSet simplex_set;
   const int k = ridge.size();
   Orbits orbits = sym_group.setwise_stabilizer(ridge).orbits();
   for (Entire<Orbits>::const_iterator oit = entire(orbits); !oit.at_end(); ++oit) {
      simplex = ridge;
      simplex += *(oit->begin());
      if (simplex.size() == k+1 and rank(V.minor(simplex, All)) == k+1) 
         simplex_set += sym_group.lex_min_representative(simplex);
   }
   const std::string dir = options["directory"];
   write_set(ridge, simplex_set, dir);
}

UserFunction4perl("# @category Producing from scratch"
                  "# Find the equivalence classes of //k//-dimensional simplices"
                  "# contained in P, modulo its symmetry group"
                  "# @param Polytope P"
                  "# @param Set<Int> facet the facet to take as starting point"
                  "# @param string the directory where the output should be written",
                  &write_simplex_reps_from_ridge, "write_simplex_reps_from_ridge(Polytope Set<Int> { group_generators => 'QUOTIENT_SPACE.SYMMETRY_GROUP.GENERATORS', directory => '.' } )");

} }

// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End:
