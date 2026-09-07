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
  AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
  Copyright (C) 2022(-2023) all rights reserved.

  The copyright of this code is transferred
  to Pierangelo Masarati and Paolo Mantegazza
  for use in the software MBDyn as described
  in the GNU Public License version 2.1

  Tpetra port: converted from Epetra by the MBDyn project.
*/

#ifndef __TPETRA_SPARSE_MATRIX_HANDLER_H__INCLUDED__
#define __TPETRA_SPARSE_MATRIX_HANDLER_H__INCLUDED__

#ifdef USE_TRILINOS
#include "myassert.h"
#include "spmh.h"
#include "cscmhtpl.h"
#include "tpetra_types.h"

class TpetraSparseMatrixHandler: public SparseMatrixHandler {
public:
     TpetraSparseMatrixHandler(const integer& iNumRows,
                               const integer& iNumCols,
                               integer iNumColsAlloc,
                               const Teuchos::RCP<const TpetraComm>& pComm);
     virtual ~TpetraSparseMatrixHandler();

#ifdef DEBUG
     virtual void IsValid() const override;
#endif

     virtual void Resize(integer, integer) override;

     virtual void ResizeReset(integer, integer) override;

     virtual void Reset() override;

     virtual const doublereal&
     operator()(integer iRow, integer iCol) const override;

     virtual doublereal&
     operator()(integer iRow, integer iCol) override;

     using SparseMatrixHandler::MakeCompressedColumnForm;
     using SparseMatrixHandler::MakeCompressedRowForm;

     virtual
     int32_t MakeCompressedRowForm(std::vector<doublereal>& Ax,
                                   std::vector<int32_t>& Ai,
                                   std::vector<int32_t>& Ap,
                                   int offset = 0) const override;

     virtual
     int64_t MakeCompressedRowForm(std::vector<doublereal>& Ax,
                                   std::vector<int64_t>& Ai,
                                   std::vector<int64_t>& Ap,
                                   int offset = 0) const override;

     virtual
     int32_t MakeCompressedRowForm(doublereal *const Ax,
                                   int32_t *const Ai,
                                   int32_t *const Ap,
                                   int offset = 0) const override;

     virtual
     int64_t MakeCompressedRowForm(doublereal *const Ax,
                                   int64_t *const Ai,
                                   int64_t *const Ap,
                                   int offset = 0) const override;

     virtual
     int32_t MakeCompressedColumnForm(std::vector<doublereal>& Ax,
                                      std::vector<int32_t>& Ai,
                                      std::vector<int32_t>& Ap,
                                      int offset = 0) const override;

     virtual
     int64_t MakeCompressedColumnForm(std::vector<doublereal>& Ax,
                                      std::vector<int64_t>& Ai,
                                      std::vector<int64_t>& Ap,
                                      int offset = 0) const override;

     virtual
     int32_t MakeCompressedColumnForm(doublereal *const Ax,
                                      int32_t *const Ai,
                                      int32_t *const Ap,
                                      int offset = 0) const override;

     virtual
     int64_t MakeCompressedColumnForm(doublereal *const Ax,
                                      int64_t *const Ai,
                                      int64_t *const Ap,
                                      int offset = 0) const override;

     virtual
     VectorHandler& GetCol(integer icol,
                           VectorHandler& out) const override;

     virtual void Scale(const std::vector<doublereal>& oRowScale,
                        const std::vector<doublereal>& oColScale) override;

     virtual bool AddItem(integer iRow,
                          const sp_grad::SpGradient& oItem) override;

     virtual void
     IncCoef(integer iRow, integer iCol,
             const doublereal& dCoef) override final;

     virtual void
     DecCoef(integer iRow, integer iCol,
             const doublereal& dCoef) override final;

     virtual void EnumerateNz(
          const std::function<EnumerateNzCallback>& func) const override;

     virtual doublereal Norm(Norm_t eNorm = NORM_1) const override;

     virtual integer Nz() const override;

     virtual TpetraSparseMatrixHandler* Copy() const override;

     virtual integer PacMat() override;

     /* ---------------------------------------------------------------
      * Iterator over non-zeros (same interface as the Epetra version)
      * --------------------------------------------------------------- */
     class const_iterator {
          friend class TpetraSparseMatrixHandler;

          const integer* const rowptr;
          const integer* const colind;
          const doublereal* const values;
#ifdef DEBUG
          const integer NRows;
#endif
          const integer NZ;
          integer iIdx;
          SparseMatrixHandler::SparseMatrixElement elem;

#ifdef DEBUG
          bool bInvariant() const {
               ASSERT(iIdx >= 0);
               ASSERT(iIdx <= NZ);
               ASSERT(iIdx == NZ || (elem.iRow >= 0 && elem.iRow < NRows));
               ASSERT(iIdx == NZ || rowptr[elem.iRow] <= iIdx);
               ASSERT(iIdx == NZ || rowptr[elem.iRow + 1] > iIdx);
               return true;
          }
#endif
          void UpdateElem() {
               if (iIdx < NZ) {
                    elem.iCol  = colind[iIdx];
                    elem.dCoef = values[iIdx];
                    if (iIdx >= rowptr[elem.iRow + 1]) {
                         ++elem.iRow;
                    }
               } else {
                    ASSERT(iIdx == NZ);
#ifdef DEBUG
                    elem.iRow  = std::numeric_limits<decltype(elem.iRow)>::min();
                    elem.iCol  = std::numeric_limits<decltype(elem.iCol)>::min();
                    elem.dCoef = -std::numeric_limits<decltype(elem.dCoef)>::max();
#endif
               }
          }

          const_iterator(const integer* rowptr_a,
                         const integer* colind_a,
                         const doublereal* values_a,
                         integer NRows_a,
                         integer iIdx_a,
                         integer iRow)
               :rowptr(rowptr_a),
                colind(colind_a),
                values(values_a),
#ifdef DEBUG
                NRows(NRows_a),
#endif
                NZ(rowptr[NRows_a] - rowptr_a[0]),
                iIdx(iIdx_a) {
               elem.iRow = iRow;
               UpdateElem();
               ASSERT(bInvariant());
          }

     public:
          ~const_iterator() { ASSERT(bInvariant()); }

          const const_iterator& operator++() {
               ASSERT(bInvariant());
               ++iIdx;
               UpdateElem();
               ASSERT(bInvariant());
               return *this;
          }

          const SparseMatrixHandler::SparseMatrixElement* operator->() const {
               ASSERT(bInvariant());
               return &elem;
          }

          const SparseMatrixHandler::SparseMatrixElement& operator*() const {
               ASSERT(bInvariant());
               return elem;
          }

          bool operator==(const const_iterator& op) const {
               ASSERT(bInvariant());
               ASSERT(rowptr == op.rowptr);
               ASSERT(colind == op.colind);
               ASSERT(values == op.values);
               return iIdx == op.iIdx;
          }

          bool operator!=(const const_iterator& op) const {
               ASSERT(bInvariant());
               ASSERT(rowptr == op.rowptr);
               ASSERT(colind == op.colind);
               ASSERT(values == op.values);
               return iIdx != op.iIdx;
          }
     };

     const_iterator begin() const;
     const_iterator end() const;

     Teuchos::RCP<const TpetraCrs> pGetTpetraCrsMatrixConst() const { return pMat; }
     Teuchos::RCP<const TpetraCrs> pGetTpetraCrsMatrix() const { return pMat; }
     Teuchos::RCP<TpetraCrs>       pGetTpetraCrsMatrix()       { return pMat; }

     using MatrixHandler::operator=;

protected:
     virtual VectorHandler&
     MatVecMul_base(void (VectorHandler::*op)(integer iRow,
                                              const doublereal& dCoef),
                    VectorHandler& out,
                    const VectorHandler& in) const override;

     virtual VectorHandler&
     MatTVecMul_base(void (VectorHandler::*op)(integer iRow,
                                               const doublereal& dCoef),
                     VectorHandler& out,
                     const VectorHandler& in) const override;

private:
     /*
      * After fillComplete() the matrix data lives in a Kokkos device view.
      * We mirror it to host once (in PacMat / EnsureFilled) and cache
      * raw pointers for the iterator and CRS-form extraction.
      */
     void EnsureFilled() const;
     void InsertOrSumValues(TpetraGO globalRow,
                            integer nEntries,
                            const doublereal* vals,
                            const TpetraGO* cols);

     CSCMatrixHandlerTpl<doublereal, integer, 0>& GetTransposedCSC() const;

     template <typename idx_type>
     idx_type MakeCompressedRowFormTpl(doublereal *const Ax,
                                       idx_type *const Ai,
                                       idx_type *const Ap,
                                       int offset) const;

     template <typename idx_type>
     idx_type MakeCompressedColumnFormTpl(doublereal *const Ax,
                                          idx_type *const Ai,
                                          idx_type *const Ap,
                                          int offset) const;

     inline integer PacMat() const {
          return const_cast<TpetraSparseMatrixHandler*>(this)->PacMat();
     }

     Teuchos::RCP<const TpetraComm>  pComm;
     Teuchos::RCP<const TpetraMap>   pRowMap;
     Teuchos::RCP<const TpetraMap>   pColMap;
     mutable Teuchos::RCP<TpetraCrs> pMat;

     const integer iNumColsAlloc;
     mutable bool bFilled;
     mutable bool bHostCacheDirty;

     /*
      * Host-side CRS arrays populated once after fillComplete().
      * These back the const_iterator and the MakeCompressed*Form methods
      * without further device round-trips.
      * After fillComplete() the CRS structure is fixed; only values change.
      * bHostCacheDirty is set when sumIntoGlobalValues modifies the Tpetra
      * matrix, so that the next EnsureFilled() call re-extracts values.
      */
     mutable std::vector<integer>    oRowPtr;   /* length NRows+1 */
     mutable std::vector<integer>    oColInd;   /* length Nz      */
     mutable std::vector<doublereal> oValues;   /* length Nz      */

     mutable CSCMatrixHandlerTpl<doublereal, integer, 0> oCscT;

     /*
      * Pre-fillComplete entry buffer.  Tpetra::CrsMatrix allocates a
      * fixed number of entries per row at construction (iNumColsAlloc)
      * and insertGlobalValues cannot grow beyond that.  We accumulate
      * entries in this map and flush them into the CrsMatrix at
      * fillComplete time in EnsureFilled().
      */
     std::vector<std::map<TpetraGO, TpetraSC>> oEntryBuffer;
};

#endif  /* USE_TRILINOS */
#endif  /* __TPETRA_SPARSE_MATRIX_HANDLER_H__INCLUDED__ */
