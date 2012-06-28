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
#include "polymake/Array.h"
#include "polymake/Set.h"
#include "polymake/Matrix.h"
#include "polymake/polytope/sympol_interface.h"
#include "polymake/group/permlib.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <bits/error_constants.h>

namespace polymake { namespace polytope {

typedef Array<Set<int> > RepArray;
typedef Array<Set<int> > Orbits;

std::string set2string(const Set<int>& s, char glue = '_') 
{
    int ct(0);
    std::ostringstream os;
    for (Entire<Set<int> >::const_iterator sit = entire(s); !sit.at_end(); ++sit) {
	if (ct++)
	    os << glue;
	os << *sit;
    }
    return os.str();
}

std::string point_to_arrow(const std::string sp)
{
    std::ostringstream os;
    for (Entire<std::string>::const_iterator sit = entire(sp); !sit.at_end(); ++sit)
        if (*sit == '.') os << "->"; 
        else os << *sit;
    return os.str();
}

template <typename Scalar>
int make_facet_makefile(perl::Object p, perl::Object ref, std::string refname, perl::OptionSet options)
{
    const std::string group_str = options["group"];
    std::ostringstream gen_str;
    gen_str << group_str << ".GENERATORS";
    const Array<Array<int> > gens = ref.give(gen_str.str().c_str());
    const group::PermlibGroup sym_group(gens);
    RepArray ridges;
    if (!(p.lookup("MAXIMAL_SIMPLEX_REPRESENTATIVES") >> ridges)) {
        Set<Set<int> > vertices; 
        for (Entire<Orbits>::const_iterator oit = entire(sym_group.orbits()); !oit.at_end(); ++oit)
            vertices += scalar2set(*(oit->begin()));
        ridges = RepArray(vertices.size(), entire(vertices));
    }
    const int this_dim = ridges[0].size();
    std::ostringstream polyname, prefix, repfilename;
    polyname << "step" << this_dim << ".poly";
    prefix << "dim" << this_dim << "_";
    repfilename << prefix.str() << "reps.txt";
    const std::string dir = options["out_dir"];
    int status = mkdir(dir.c_str(), ACCESSPERMS);
    if (status != 0 && errno != EEXIST) { 
        cerr << "couldn't create directory " << dir << endl;
        throw std::runtime_error("Halt");
    }
    const std::string subdir = dir + "/" + prefix.str();
    status = mkdir(subdir.c_str(), ACCESSPERMS);
    if (status != 0 && errno != EEXIST) {
        cerr << "couldn't create directory " << subdir << endl;
        throw std::runtime_error("Halt");
    }
    std::ofstream makefile((subdir + std::string("/Makefile")).c_str(), std::ios_base::trunc);

    makefile << "all: " << polyname.str() << endl << endl;

    makefile << polyname.str() << ": " << repfilename.str()
	     << "\n\t polymake 'my $$p=load(\"../../" << refname << "\"); "
	     << "my $$a=read_sets_to_array(\"" << repfilename.str() << "\"); "
	     << "my $$q=new " << p.type().name() << "(VERTICES=>$$p->VERTICES, "
	     << "MAXIMAL_SIMPLEX_REPRESENTATIVES=>$$a); "
	     << "save($$q, \"" << polyname.str() << "\");'" << endl
             << "\t ln -sf " << prefix.str() << "/" << polyname.str() << " ../" << polyname.str()
	     << endl << endl;

    makefile << repfilename.str() << ": unsorted\n\t cat " << prefix.str() 
	     << "* | sort | uniq > " << repfilename.str() << endl << endl;

    makefile << "objects = ";

    for (Entire<RepArray>::const_iterator sit = entire(ridges); !sit.at_end(); ++sit) 
        makefile << prefix.str() << set2string(*sit) << " ";
    makefile << endl << endl;

    makefile << "comma:=," << endl << endl;

    makefile << "unsorted: $(objects)" << endl << endl;

    makefile << "$(objects): " << prefix.str() << "%:" << endl
             << "\tpolymake 'my $$p=load(\"../../" << refname << "\"); "
             << "my $$s=new Set<Int>($(subst _,$(comma),$(subst " << prefix.str() << ",,$@))); "
             << "write_simplex_reps_from_ridge($$p, $$s);'" 
             << endl << endl;
    return 1;
}

UserFunctionTemplate4perl("# @category Producing from scratch"
                          "# Find the equivalence classes of //k//-dimensional simplices"
                          "# contained in P, modulo its symmetry group"
                          "# @param Polytope P"
                          "# @param refname a filename that contains just the VERTICES and the @option group of the polytope"
                          "# @param directory the name of a directory where to write the makefile",
                          "make_facet_makefile<Scalar>(Polytope<Scalar> Polytope<Scalar> $ { out_dir => '.', group => 'QUOTIENT_SPACE.SYMMETRY_GROUP' } )");

} }

// Local Variables:
// mode:C++
// c-basic-offset:4
// indent-tabs-mode:nil
// End:
