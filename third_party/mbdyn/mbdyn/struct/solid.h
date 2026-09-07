/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2023
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
 * Paolo Mantegazza     <mantegazza@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation (version 2 of the License).
 *
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
        Copyright (C) 2022(-2023) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

#ifndef ___SOLID_H__INCLUDED___
#define ___SOLID_H__INCLUDED___

#include "dataman.h"
#include "elem.h"
#include "gravity.h"

// 3D elements
class Hexahedron8;
class Hexahedron8f;
class Hexahedron8upc;
class Hexahedron20;
class Hexahedron20f;
class Hexahedron20upc;
class Hexahedron20r;
class Hexahedron20fr;
class Hexahedron20upcr;
class Hexahedron27;
class Hexahedron27f;
class Pentahedron6;
class Pentahedron15;
class Pentahedron15f;
class Pentahedron15upc;
class Pentahedron18;
class Pentahedron18f;
class Tetrahedron10;
class Tetrahedron10f;
class Tetrahedron10upc;
class Tetrahedron20;
class Tetrahedron20f;

// 3D collocation rules
class Gauss2x2x2;
class Gauss3x3x3;
class GaussH20r;
class CollocPenta15;
class CollocPenta18;
class CollocTet10h;
class CollocTet20;

// 3D base class
class SolidElem: public InitialAssemblyElem, public GravityOwner {
public:
     SolidElem(unsigned uLabel,
               flag fOut);
     virtual ~SolidElem();

     virtual Elem::Type GetElemType() const override;

     virtual void
     SetValue(DataManager *pDM,
              VectorHandler& X, VectorHandler& XP,
              SimulationEntity::Hints *ph) override;

     virtual std::ostream& Restart(std::ostream& out) const override;

     virtual unsigned int iGetInitialNumDof() const override;

     virtual bool bIsDeformable() const override;

     virtual unsigned int iGetNumPrivData() const override;

     virtual unsigned int iGetPrivDataIdx(const char *s) const override;

     virtual doublereal dGetPrivData(unsigned int i) const override;

     virtual doublereal dGetE() const;

     virtual doublereal dGetU() const = 0;

     virtual doublereal dGetDissipatedEnergy() const = 0;

     virtual doublereal dGetRequestedTimeStep() const = 0;
};

template <typename ElementType, typename CollocationType>
SolidElem*
ReadSolid(DataManager* pDM, MBDynParser& HP, unsigned int uLabel);

#endif
