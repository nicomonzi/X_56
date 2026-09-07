/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2023
 *
 * Pierangelo Masarati  <pierangelo.masarati@polimi.it>
 * Paolo Mantegazza     <paolo.mantegazza@polimi.it>
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
  Shared Tpetra scalar/ordinal/node type aliases.
  Every translation unit that touches the Tpetra wrapper should
  include this header instead of repeating the typedefs.

  Tpetra port: converted from Epetra by the MBDyn project.
*/

#ifndef ___TPETRA_TYPES__INCLUDED___
#define ___TPETRA_TYPES__INCLUDED___

#ifdef USE_TRILINOS

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <Tpetra_Core.hpp>
#pragma GCC diagnostic pop

#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Import.hpp>

#pragma GCC diagnostic pop

/* -----------------------------------------------------------------------
 * Scalar / ordinal / node selection
 * MBDyn always uses double precision and sequential (or MPI) execution.
 * ----------------------------------------------------------------------- */
using TpetraSC   = double;
using TpetraLO   = int;
using TpetraGO   = long long;   /* 64-bit global indices */
using TpetraNode = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType;

/* Convenience type aliases */
using TpetraComm   = Teuchos::Comm<int>;
using TpetraMap    = Tpetra::Map<TpetraLO, TpetraGO, TpetraNode>;
using TpetraCrs    = Tpetra::CrsMatrix<TpetraSC, TpetraLO, TpetraGO, TpetraNode>;
using TpetraVector = Tpetra::Vector<TpetraSC, TpetraLO, TpetraGO, TpetraNode>;
using TpetraMV     = Tpetra::MultiVector<TpetraSC, TpetraLO, TpetraGO, TpetraNode>;
using TpetraImport = Tpetra::Import<TpetraLO, TpetraGO, TpetraNode>;
using TpetraOp     = Tpetra::Operator<TpetraSC, TpetraLO, TpetraGO, TpetraNode>;

/* Helper: build the serial (non-MPI) communicator */
inline Teuchos::RCP<const TpetraComm> tpetraSerialComm()
{
     return Teuchos::rcp(new Teuchos::SerialComm<int>());
}

#ifdef USE_MPI
#include <Teuchos_DefaultMpiComm.hpp>
/* Helper: wrap an MPI communicator */
inline Teuchos::RCP<const TpetraComm> tpetraMpiComm(MPI_Comm mpiComm)
{
     return Teuchos::rcp(new Teuchos::MpiComm<int>(mpiComm));
}
#endif  /* USE_MPI */

#endif  /* USE_TRILINOS */
#endif  /* ___TPETRA_TYPES__INCLUDED___ */
