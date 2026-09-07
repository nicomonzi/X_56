/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2025
 *
 * Pierangelo Masarati	<pierangelo.masarati@polimi.it>
 * Paolo Mantegazza	<paolo.mantegazza@polimi.it>
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

#ifndef ___RESTART_DATA_H__INCLUDED___
#define ___RESTART_DATA_H__INCLUDED___

#include <map>
#include <iostream>
#include <sstream>

#include "ac/f2c.h"
#include "except.h"
#include "myassert.h"
#include "binary_conversion.h"
#include "matvec3.h"
#include "matvec6.h"
#include "matvec3n.h"
#include "sp_matrix_base.h"

class RestartData {
public:
     enum RestartEntity: unsigned {
          DATA_MANAGER,
          INITVAL_TIME,
          NODES_ABSTRACT,
          NODES_STRUCT,
          NODES_ELECTRIC,
          NODES_THERMAL,
          NODES_HYDRAULIC,
          NODES_PARAMETER,
          ELEM_AUTOSTRUCT,
          ELEM_BODIES,
          ELEM_JOINTS,
          ELEM_BEAMS,
          ELEM_SOLIDS,
          ELEM_GENELS,
          ELEM_LOADABLE,
          ELEM_INERTIA,
          ELEM_HYDRAULIC,
          ELEM_ELECTRIC,
          ELEM_THERMAL
     };

     enum RestartAction {
          RESTART_SAVE,
          RESTART_RESTORE
     };

     class Exception: public MBDynErrBase {
     public:
          Exception(MBDYN_EXCEPT_ARGS_DECL): MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {
          }
     };

     class ExceptionRead: public Exception {
     public:
          ExceptionRead(MBDYN_EXCEPT_ARGS_DECL): Exception(MBDYN_EXCEPT_ARGS_PASSTHRU) {
          }
     };

     class ExceptionWrite: public Exception {
     public:
          ExceptionWrite(MBDYN_EXCEPT_ARGS_DECL): Exception(MBDYN_EXCEPT_ARGS_PASSTHRU) {
          }
     };

     void WriteFile(const std::string& filename);
     void ReadFile(const std::string& filename);

     void WriteFile(std::ostream& os);
     void ReadFile(std::istream& is);

     // To be used by nodes and elements
     template <typename T>
     void Sync(RestartEntity eEntity, unsigned int uLabel, const std::string& strName, T& oData, RestartAction eAction) {
          Sync(eEntity, uLabel, -1, strName, oData, eAction);
     }

     // To be used by constitutive laws
     template <typename T>
     void Sync(RestartEntity eEntity, unsigned int uLabel, integer iIndex, const std::string& strName, T& oData, RestartAction eAction) {
          std::string strRawData;

          if (eAction == RESTART_SAVE) {
               strRawData = ExtractRawData(oData);
          }

          const Key oKey{eEntity, uLabel, iIndex, strName};
          
          if (Sync(oKey, strRawData, eAction)) {
               RestoreRawData(oData, strRawData);
          }
     }

     class ErrDuplicateEntry: public ErrGeneric {
     public:
          ErrDuplicateEntry(MBDYN_EXCEPT_ARGS_DECL): ErrGeneric(MBDYN_EXCEPT_ARGS_PASSTHRU) {
          }
     };

     class ErrInvalidFileFormat: public ErrGeneric {
     public:
          ErrInvalidFileFormat(MBDYN_EXCEPT_ARGS_DECL): ErrGeneric(MBDYN_EXCEPT_ARGS_PASSTHRU) {
          }
     };
public:
     struct Key {
          RestartEntity eEntity;
          unsigned uLabel;
          integer iIndex;
          std::string strName;

          bool operator<(const Key& oKey) const {
               if (eEntity < oKey.eEntity) {
                    return true;
               } else if (eEntity > oKey.eEntity) {
                    return false;
               }

               if (uLabel < oKey.uLabel) {
                    return true;
               } else if (uLabel > oKey.uLabel) {
                    return false;
               }

               if (iIndex < oKey.iIndex) {
                    return true;
               } else if (iIndex > oKey.iIndex) {
                    return false;
               }

               return strName < oKey.strName;
          }
     };
private:
     template <typename T>
     static std::string ExtractRawData(const T& oData) {
          std::ostringstream os(std::ios::binary | std::ios::out);

          BinaryConversion::WriteBinary(os, oData);

          return os.str();
     }

     template <typename T>
     static void RestoreRawData(T& oData, const std::string& strRawData) {
          std::istringstream is(strRawData, std::ios::binary | std::ios::in);

          BinaryConversion::ReadBinary(is, oData);
     }
     
     bool Sync(const Key& oKey, std::string& strData, RestartAction eAction);

     std::map<Key, std::string> m_Map;

     static constexpr char szHeaderString[] = "MBDyn binary restart file version 0.0.1\n";
};

#endif
