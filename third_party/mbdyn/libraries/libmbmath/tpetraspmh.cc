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

#include "mbconfig.h"

#ifdef USE_TRILINOS
#include <algorithm>
#include <iomanip>
#include <new>

#include "sp_gradient.h"
#undef HAVE_BLAS
#include "tpetraspmh.h"
#include "tpetravh.h"
#include "cscmhtpl.h"

/* -----------------------------------------------------------------------
 * Construction
 * ----------------------------------------------------------------------- */
TpetraSparseMatrixHandler::TpetraSparseMatrixHandler(
     const integer& iNumRows_a,
     const integer& iNumCols_a,
     integer iNumColsAlloc_a,
     const Teuchos::RCP<const TpetraComm>& pComm_a)
     :SparseMatrixHandler(iNumRows_a, iNumCols_a),
      pComm(pComm_a),
      iNumColsAlloc(iNumColsAlloc_a),
      bFilled(false),
      bHostCacheDirty(false),
      oEntryBuffer(iNumRows_a)
{
     if (iNumRows_a != iNumCols_a) {
          silent_cerr("TpetraSparseMatrixHandler: matrix must be square!\n");
          throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
     }

     const TpetraGO nGlobal   = static_cast<TpetraGO>(iNumRows_a);
     const TpetraGO indexBase = 0;  /* 0-based global indices */

     pRowMap = Teuchos::rcp(new TpetraMap(nGlobal, indexBase, pComm));
     pColMap = pRowMap;             /* square matrix – same map for rows and cols */

     /* CrsMatrix is NOT created here; it is built from the entry buffer
      * in EnsureFilled() once the full sparsity pattern is known. */

#ifdef DEBUG
     IsValid();
#endif
}

TpetraSparseMatrixHandler::~TpetraSparseMatrixHandler()
{
#ifdef DEBUG
     IsValid();
#endif
}

/* -----------------------------------------------------------------------
 * Debug
 * ----------------------------------------------------------------------- */
#ifdef DEBUG
void TpetraSparseMatrixHandler::IsValid() const
{
     if (bFilled) {
          ASSERT(pMat.get() != nullptr);
          ASSERT(static_cast<integer>(pMat->getGlobalNumRows()) == NRows);
          ASSERT(static_cast<integer>(pMat->getGlobalNumCols()) == NCols);
     }
}
#endif

/* -----------------------------------------------------------------------
 * Internal helpers
 * ----------------------------------------------------------------------- */

/*
 * Call fillComplete() and cache the host-side CRS arrays once so that
 * subsequent calls to iterator, MakeCompressed*Form, Norm, Nz etc. are
 * free of additional Kokkos syncs.
 */
void TpetraSparseMatrixHandler::EnsureFilled() const
{
     if (bFilled && !bHostCacheDirty) {
          return;
     }

     if (!bFilled) {
          /* Build the CrsMatrix from the entry buffer.
           * First build a CrsGraph with the exact sparsity pattern, then
           * create the CrsMatrix from it so that isStaticGraph()==true.
           * This ensures sumIntoGlobalValues works after fillComplete. */
          using TpetraGraph = Tpetra::CrsGraph<TpetraLO, TpetraGO, TpetraNode>;

          const integer nRows = NRows;
          Teuchos::Array<size_t> numEntriesPerRow(nRows);
          for (integer i = 0; i < nRows; ++i) {
               numEntriesPerRow[i] = oEntryBuffer[i].size();
          }

          auto pGraph = Teuchos::rcp(new TpetraGraph(
               pRowMap, pColMap, numEntriesPerRow()));

          /* Insert column indices row by row. */
          for (integer i = 0; i < nRows; ++i) {
               const auto& rowEntries = oEntryBuffer[i];
               if (rowEntries.empty()) {
                    continue;
               }
               std::vector<TpetraGO> cols;
               cols.reserve(rowEntries.size());
               for (const auto& entry : rowEntries) {
                    cols.push_back(entry.first);
               }
               Teuchos::ArrayView<const TpetraGO> colView(cols.data(), cols.size());
               pGraph->insertGlobalIndices(static_cast<TpetraGO>(i), colView);
          }

          pGraph->fillComplete(pColMap, pRowMap);

          /* Create the matrix from the static graph and fill values. */
          pMat = Teuchos::rcp(new TpetraCrs(pGraph));

          for (integer i = 0; i < nRows; ++i) {
               const auto& rowEntries = oEntryBuffer[i];
               if (rowEntries.empty()) {
                    continue;
               }
               std::vector<TpetraGO> cols;
               std::vector<TpetraSC> vals;
               cols.reserve(rowEntries.size());
               vals.reserve(rowEntries.size());
               for (const auto& entry : rowEntries) {
                    cols.push_back(entry.first);
                    vals.push_back(entry.second);
               }
               Teuchos::ArrayView<const TpetraGO> colView(cols.data(), cols.size());
               Teuchos::ArrayView<const TpetraSC> valView(vals.data(), vals.size());
               pMat->replaceGlobalValues(static_cast<TpetraGO>(i), colView, valView);
          }

          pMat->fillComplete(pColMap, pRowMap);
          bFilled = true;
     }

     /* Ensure the matrix is fill-complete before extracting host arrays.
      * Reset() calls resumeFill() to allow value modification;
      * we must re-finalize here. */
     if (pMat->isFillActive()) {
          pMat->fillComplete(pColMap, pRowMap);
     }

     /* Extract into host std::vectors for cheap random access later. */
     const auto& rowPtrView = pMat->getLocalRowPtrsHost();
     const auto& colIndView = pMat->getLocalIndicesHost();
     const auto& valuesView = pMat->getLocalValuesHost(Tpetra::Access::ReadOnly);

     const integer nRows = NRows;
     const integer nNz   = static_cast<integer>(pMat->getLocalNumEntries());

     oRowPtr.resize(nRows + 1);
     oColInd.resize(nNz);
     oValues.resize(nNz);

     for (integer i = 0; i <= nRows; ++i) {
          oRowPtr[i] = static_cast<integer>(rowPtrView[i]);
     }
     for (integer k = 0; k < nNz; ++k) {
          oColInd[k] = static_cast<integer>(colIndView[k]);
          oValues[k] = static_cast<doublereal>(valuesView[k]);
     }

     /* Build the transposed CSC for MatVecMul/MatTVecMul and operator(). */
     oCscT = CSCMatrixHandlerTpl<doublereal, integer, 0>(
          oValues.data(), oColInd.data(), oRowPtr.data(), nRows, nNz);

     bHostCacheDirty = false;
}

void TpetraSparseMatrixHandler::InsertOrSumValues(TpetraGO globalRow,
                                                   integer nEntries,
                                                   const doublereal* vals,
                                                   const TpetraGO* cols)
{
     if (nEntries == 0) {
          return;
     }

     if (bFilled) {
          /* After fillComplete: use sumIntoGlobalValues */
          Teuchos::ArrayView<const TpetraGO> colView(cols,  nEntries);
          Teuchos::ArrayView<const TpetraSC> valView(vals,  nEntries);

          const TpetraLO err =
               pMat->sumIntoGlobalValues(globalRow, colView, valView);
          if (err < 0) {
               ASSERT(0);
               throw ErrGeneric(MBDYN_EXCEPT_ARGS);
          }
          if (err != nEntries) {
               /* Some entries were not in the sparsity pattern;
                * the matrix structure must be rebuilt from scratch. */
               throw MatrixHandler::ErrRebuildMatrix(MBDYN_EXCEPT_ARGS);
          }
          bHostCacheDirty = true;
     } else {
          /* Before fillComplete: accumulate in entry buffer.
           * Duplicate column indices are summed (matching IncCoef semantics). */
          ASSERT(globalRow >= 0 && globalRow < static_cast<TpetraGO>(NRows));
          auto& rowMap = oEntryBuffer[globalRow];
          for (integer k = 0; k < nEntries; ++k) {
               rowMap[cols[k]] += vals[k];
          }
     }
}

CSCMatrixHandlerTpl<doublereal, integer, 0>&
TpetraSparseMatrixHandler::GetTransposedCSC() const
{
     ASSERT(bFilled);
#ifdef DEBUG
     IsValid();
     oCscT.IsValid();
#endif
     return oCscT;
}

/* -----------------------------------------------------------------------
 * SparseMatrixHandler interface
 * ----------------------------------------------------------------------- */
void TpetraSparseMatrixHandler::Resize(integer, integer)
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

void TpetraSparseMatrixHandler::ResizeReset(integer, integer)
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

void TpetraSparseMatrixHandler::Reset()
{
#ifdef DEBUG
     IsValid();
#endif

     if (bFilled) {
          /*
           * The matrix structure is fixed after fillComplete; only values
           * need to be zeroed.  We must resumeFill() first because
           * sumIntoGlobalValues (used later in IncCoef) requires
           * isFillActive()==true.  fillComplete() will be called in
           * PacMat()/EnsureFilled() before the matrix is used by a solver.
           */
          if (!pMat->isFillActive()) {
               pMat->resumeFill();
          }
          pMat->setAllToScalar(0.);
          std::fill(oValues.begin(), oValues.end(), 0.);

          /* oCscT shares the same values buffer; it stays valid. */
          bHostCacheDirty = false;
     } else {
          /* Haven't filled yet – clear the entry buffer for fresh assembly. */
          for (auto& rowMap : oEntryBuffer) {
               rowMap.clear();
          }
          pMat = Teuchos::null;
          oRowPtr.clear();
          oColInd.clear();
          oValues.clear();
          bHostCacheDirty = false;
     }

#ifdef DEBUG
     IsValid();
#endif
}

integer TpetraSparseMatrixHandler::PacMat()
{
#ifdef DEBUG
     IsValid();
#endif
     EnsureFilled();
#ifdef DEBUG
     IsValid();
#endif
     return Nz();
}

integer TpetraSparseMatrixHandler::Nz() const
{
#ifdef DEBUG
     IsValid();
#endif
     ASSERT(bFilled);
     return static_cast<integer>(pMat->getGlobalNumEntries());
}

/* -----------------------------------------------------------------------
 * Element access (1-based MBDyn convention)
 * ----------------------------------------------------------------------- */
const doublereal&
TpetraSparseMatrixHandler::operator()(integer iRow, integer iCol) const
{
     if (!bFilled) {
          throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
     }
     return GetTransposedCSC()(iCol, iRow);
}

doublereal&
TpetraSparseMatrixHandler::operator()(integer iRow, integer iCol)
{
     if (!bFilled) {
          throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
     }
     return GetTransposedCSC()(iCol, iRow);
}

void
TpetraSparseMatrixHandler::IncCoef(integer iRow, integer iCol,
                                    const doublereal& dCoef)
{
     const TpetraGO globalRow = static_cast<TpetraGO>(iRow - 1); /* 0-based */
     const TpetraGO globalCol = static_cast<TpetraGO>(iCol - 1);

     InsertOrSumValues(globalRow, 1, &dCoef, &globalCol);
}

void
TpetraSparseMatrixHandler::DecCoef(integer iRow, integer iCol,
                                    const doublereal& dCoef)
{
     IncCoef(iRow, iCol, -dCoef);
}

bool TpetraSparseMatrixHandler::AddItem(integer iRow,
                                         const sp_grad::SpGradient& oItem)
{
#ifdef DEBUG
     IsValid();
#endif
     ASSERT(iRow >= 1);
     ASSERT(iRow <= iGetNumRows());

     std::vector<doublereal>         rgValues;
     std::vector<TpetraGO>           rgColIndex;

     rgValues.reserve(oItem.iGetSize());
     rgColIndex.reserve(oItem.iGetSize());

     for (const auto& oDer: oItem) {
          rgValues.push_back(oDer.dDer);
          rgColIndex.push_back(static_cast<TpetraGO>(oDer.iDof - 1)); /* 0-based */
     }

     const TpetraGO globalRow = static_cast<TpetraGO>(iRow - 1);
     InsertOrSumValues(globalRow,
                       static_cast<integer>(rgValues.size()),
                       rgValues.data(),
                       rgColIndex.data());
     return true;
}

/* -----------------------------------------------------------------------
 * Matrix-vector products via the transposed CSC cache
 * ----------------------------------------------------------------------- */
VectorHandler&
TpetraSparseMatrixHandler::MatVecMul_base(
     void (VectorHandler::*op)(integer, const doublereal&),
     VectorHandler& out, const VectorHandler& in) const
{
     EnsureFilled();
     return GetTransposedCSC().MatTVecMul_base(op, out, in);
}

VectorHandler&
TpetraSparseMatrixHandler::MatTVecMul_base(
     void (VectorHandler::*op)(integer, const doublereal&),
     VectorHandler& out, const VectorHandler& in) const
{
     EnsureFilled();
     return GetTransposedCSC().MatVecMul_base(op, out, in);
}

/* -----------------------------------------------------------------------
 * Norm
 * ----------------------------------------------------------------------- */
doublereal TpetraSparseMatrixHandler::Norm(Norm_t eNorm) const
{
#ifdef DEBUG
     IsValid();
#endif
     EnsureFilled();

     const integer nNz   = Nz();
     const integer nRows = NRows;
     const integer nCols = NCols;

     switch (eNorm) {
     case NORM_1: {
          /* Column sum norm: max over columns of sum of |a_ij| */
          std::vector<doublereal> colSum(nCols, 0.);
          for (integer k = 0; k < nNz; ++k) {
               colSum[oColInd[k]] += std::fabs(oValues[k]);
          }
          return *std::max_element(colSum.begin(), colSum.end());
     }
     case NORM_INF: {
          /* Row sum norm: max over rows of sum of |a_ij| */
          std::vector<doublereal> rowSum(nRows, 0.);
          for (integer i = 0; i < nRows; ++i) {
               for (integer k = oRowPtr[i]; k < oRowPtr[i + 1]; ++k) {
                    rowSum[i] += std::fabs(oValues[k]);
               }
          }
          return *std::max_element(rowSum.begin(), rowSum.end());
     }
     default:
          throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
     }
}

/* -----------------------------------------------------------------------
 * Scale
 * ----------------------------------------------------------------------- */
void TpetraSparseMatrixHandler::Scale(
     const std::vector<doublereal>& oRowScale,
     const std::vector<doublereal>& oColScale)
{
#ifdef DEBUG
     IsValid();
#endif
     EnsureFilled();

     CSCMatrixHandlerTpl<doublereal, integer, 0> A_T(
          oValues.data(), oColInd.data(), oRowPtr.data(), NRows, Nz());

     A_T.Scale(oColScale, oRowScale); /* exchange rows ↔ cols */

     /* Push scaled values back into the Tpetra matrix */
     auto valHost = pMat->getLocalValuesHost(Tpetra::Access::ReadWrite);
     ASSERT(static_cast<size_t>(Nz()) == valHost.size());
     for (integer k = 0; k < Nz(); ++k) {
          valHost[k] = oValues[k];
     }

#ifdef DEBUG
     IsValid();
#endif
}

/* -----------------------------------------------------------------------
 * EnumerateNz
 * ----------------------------------------------------------------------- */
void TpetraSparseMatrixHandler::EnumerateNz(
     const std::function<EnumerateNzCallback>& func) const
{
#ifdef DEBUG
     IsValid();
#endif
     EnsureFilled();

     const integer nRows = iGetNumRows();
     for (integer iRow = 0; iRow < nRows; ++iRow) {
          for (integer k = oRowPtr[iRow]; k < oRowPtr[iRow + 1]; ++k) {
               func(iRow + 1, oColInd[k] + 1, oValues[k]); /* 1-based */
          }
     }
}

/* -----------------------------------------------------------------------
 * GetCol (not implemented – same as Epetra version)
 * ----------------------------------------------------------------------- */
VectorHandler& TpetraSparseMatrixHandler::GetCol(integer, VectorHandler&) const
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

/* -----------------------------------------------------------------------
 * Copy
 * ----------------------------------------------------------------------- */
TpetraSparseMatrixHandler* TpetraSparseMatrixHandler::Copy() const
{
#ifdef DEBUG
     IsValid();
#endif
     TpetraSparseMatrixHandler* pMH = nullptr;

     SAFENEWWITHCONSTRUCTOR(pMH,
                            TpetraSparseMatrixHandler,
                            TpetraSparseMatrixHandler(iGetNumRows(),
                                                      iGetNumCols(),
                                                      iNumColsAlloc,
                                                      pComm));
     return pMH;
}

/* -----------------------------------------------------------------------
 * Iterators
 * ----------------------------------------------------------------------- */
TpetraSparseMatrixHandler::const_iterator
TpetraSparseMatrixHandler::begin() const
{
     EnsureFilled();
     return const_iterator(oRowPtr.data(), oColInd.data(), oValues.data(),
                           NRows, 0, 0);
}

TpetraSparseMatrixHandler::const_iterator
TpetraSparseMatrixHandler::end() const
{
     EnsureFilled();
     return const_iterator(oRowPtr.data(), oColInd.data(), oValues.data(),
                           NRows, Nz(), NRows - 1);
}

/* -----------------------------------------------------------------------
 * MakeCompressedRowForm  (template)
 * ----------------------------------------------------------------------- */
template <typename idx_type>
idx_type TpetraSparseMatrixHandler::MakeCompressedRowFormTpl(
     doublereal *const Ax,
     idx_type *const Ai,
     idx_type *const Ap,
     int offset) const
{
#ifdef DEBUG
     IsValid();
#endif
     EnsureFilled();

     const integer nRows = iGetNumRows();
     const integer nNz   = Nz();

     for (integer i = 0; i <= nRows; ++i) {
          Ap[i] = static_cast<idx_type>(oRowPtr[i]) + static_cast<idx_type>(offset);
     }
     for (integer k = 0; k < nNz; ++k) {
          Ax[k] = oValues[k];
          Ai[k] = static_cast<idx_type>(oColInd[k]) + static_cast<idx_type>(offset);
     }
     return static_cast<idx_type>(nNz);
}

int32_t TpetraSparseMatrixHandler::MakeCompressedRowForm(
     std::vector<doublereal>& Ax, std::vector<int32_t>& Ai,
     std::vector<int32_t>& Ap, int offset) const
{
     EnsureFilled();
     Ai.resize(Nz()); Ax.resize(Nz()); Ap.resize(iGetNumCols() + 1);
     return MakeCompressedRowForm(Ax.data(), Ai.data(), Ap.data(), offset);
}

int64_t TpetraSparseMatrixHandler::MakeCompressedRowForm(
     std::vector<doublereal>& Ax, std::vector<int64_t>& Ai,
     std::vector<int64_t>& Ap, int offset) const
{
     EnsureFilled();
     Ai.resize(Nz()); Ax.resize(Nz()); Ap.resize(iGetNumCols() + 1);
     return MakeCompressedRowForm(Ax.data(), Ai.data(), Ap.data(), offset);
}

int32_t TpetraSparseMatrixHandler::MakeCompressedRowForm(
     doublereal *const Ax, int32_t *const Ai, int32_t *const Ap, int offset) const
{
     return MakeCompressedRowFormTpl(Ax, Ai, Ap, offset);
}

int64_t TpetraSparseMatrixHandler::MakeCompressedRowForm(
     doublereal *const Ax, int64_t *const Ai, int64_t *const Ap, int offset) const
{
     return MakeCompressedRowFormTpl(Ax, Ai, Ap, offset);
}

/* -----------------------------------------------------------------------
 * MakeCompressedColumnForm  (template)
 * ----------------------------------------------------------------------- */
template <typename idx_type>
idx_type TpetraSparseMatrixHandler::MakeCompressedColumnFormTpl(
     doublereal *const Ax,
     idx_type *const Ai,
     idx_type *const Ap,
     int offset) const
{
#ifdef DEBUG
     IsValid();
#endif
     EnsureFilled();

     const integer nRows = iGetNumRows();
     const integer nCols = iGetNumCols();
     const integer nNz   = Nz();

     std::vector<idx_type> colSize(nCols, 0);

     for (integer k = 0; k < nNz; ++k) {
          ASSERT(oColInd[k] >= 0);
          ASSERT(oColInd[k] < nCols);
          ++colSize[oColInd[k]];
     }

     idx_type iPtr = static_cast<idx_type>(offset);
     for (integer iCol = 0; iCol < nCols; ++iCol) {
          Ap[iCol] = iPtr;
          iPtr    += colSize[iCol];
     }
     Ap[nCols] = iPtr;
     ASSERT(iPtr - offset == static_cast<idx_type>(nNz));

     std::fill(colSize.begin(), colSize.end(), 0);

#ifdef DEBUG
     constexpr idx_type iInvalidIndex = -1;
     constexpr doublereal dInvalidValue = -123456789.;
     std::fill(Ai, Ai + nNz, iInvalidIndex);
     std::fill(Ax, Ax + nNz, dInvalidValue);
#endif

     for (integer iRow = 0; iRow < nRows; ++iRow) {
          for (integer k = oRowPtr[iRow]; k < oRowPtr[iRow + 1]; ++k) {
               const integer iColIdx = oColInd[k];
               const idx_type iSlot  = Ap[iColIdx] - offset + colSize[iColIdx]++;

               ASSERT(iSlot >= 0);
               ASSERT(iSlot < static_cast<idx_type>(nNz));

               Ai[iSlot] = static_cast<idx_type>(iRow) + static_cast<idx_type>(offset);
               Ax[iSlot] = oValues[k];
          }
     }

     return iPtr - static_cast<idx_type>(offset);
}

int32_t TpetraSparseMatrixHandler::MakeCompressedColumnForm(
     std::vector<doublereal>& Ax, std::vector<int32_t>& Ai,
     std::vector<int32_t>& Ap, int offset) const
{
     EnsureFilled();
     Ai.resize(Nz()); Ax.resize(Nz()); Ap.resize(iGetNumCols() + 1);
     return MakeCompressedColumnForm(Ax.data(), Ai.data(), Ap.data(), offset);
}

int64_t TpetraSparseMatrixHandler::MakeCompressedColumnForm(
     std::vector<doublereal>& Ax, std::vector<int64_t>& Ai,
     std::vector<int64_t>& Ap, int offset) const
{
     EnsureFilled();
     Ai.resize(Nz()); Ax.resize(Nz()); Ap.resize(iGetNumCols() + 1);
     return MakeCompressedColumnForm(Ax.data(), Ai.data(), Ap.data(), offset);
}

int32_t TpetraSparseMatrixHandler::MakeCompressedColumnForm(
     doublereal *const Ax, int32_t *const Ai, int32_t *const Ap, int offset) const
{
     return MakeCompressedColumnFormTpl(Ax, Ai, Ap, offset);
}

int64_t TpetraSparseMatrixHandler::MakeCompressedColumnForm(
     doublereal *const Ax, int64_t *const Ai, int64_t *const Ap, int offset) const
{
     return MakeCompressedColumnFormTpl(Ax, Ai, Ap, offset);
}

#endif  /* USE_TRILINOS */
