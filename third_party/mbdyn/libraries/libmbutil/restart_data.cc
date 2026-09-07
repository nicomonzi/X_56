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

#include "mbconfig.h"

#include <fstream>

#include "myassert.h"
#include "restart_data.h"

void RestartData::WriteFile(const std::string& filename)
{
     std::ofstream os(filename, std::ios::binary);

     WriteFile(os);
}

void RestartData::ReadFile(const std::string& filename)
{
     std::ifstream is(filename, std::ios::binary);

     ReadFile(is);
}

void RestartData::WriteFile(std::ostream& os)
{
     os.exceptions(std::ios::badbit | std::ios::failbit);

     using namespace BinaryConversion;

     os << szHeaderString;

     const size_t uNumItems = m_Map.size();

     DEBUGCERR("Number of items to write: " << uNumItems << "\n");

     WriteBinary(os, uNumItems);

#ifdef DEBUG
     size_t i = 0;
#endif
     for (const auto& oItem: m_Map) {
          DEBUGCERR("Writing item " << i << "\n");
          WriteBinary(os, static_cast<unsigned>(oItem.first.eEntity));
          WriteBinary(os, oItem.first.uLabel);
          WriteBinary(os, oItem.first.iIndex);
          WriteBinary(os, oItem.first.strName);
          WriteBinary(os, oItem.second);
          DEBUGCERR("Item written " << i++ << " " << oItem.first.eEntity << " " << oItem.first.uLabel << " " << oItem.first.iIndex << " " << oItem.first.strName << "\n");
     }
}

void RestartData::ReadFile(std::istream& is)
{
     is.exceptions(std::ios::badbit | std::ios::failbit);

     using namespace BinaryConversion;

     std::string strHeader;

     std::getline(is, strHeader, '\n');

     constexpr size_t nHeaderLength = sizeof(szHeaderString) / sizeof(szHeaderString[0]) - 1;

     if (0 == strHeader.compare(0, nHeaderLength, szHeaderString)) {
          throw ErrInvalidFileFormat(MBDYN_EXCEPT_ARGS);
     }

     size_t uNumItems = 0;

     ReadBinary(is, uNumItems);

     DEBUGCERR("Number of items to read: " << uNumItems << "\n");

     Key oKey;
     std::string strData;

     m_Map.clear();

     for (size_t i = 0; i < uNumItems; ++i) {
          DEBUGCERR("Reading item " << i << "\n");

          unsigned uEntity;
          ReadBinary(is, uEntity);
          oKey.eEntity = static_cast<RestartEntity>(uEntity);
          ReadBinary(is, oKey.uLabel);
          ReadBinary(is, oKey.iIndex);
          ReadBinary(is, oKey.strName);
          ReadBinary(is, strData);
          m_Map.insert(std::make_pair(oKey, strData));

          DEBUGCERR("Insert item " << i << " " << uEntity << " " << oKey.uLabel << " " << oKey.iIndex << " " << oKey.strName << "\n");
     }
}

bool RestartData::Sync(const Key& oKey, std::string& strData, RestartAction eAction)
{
     auto it = m_Map.find(oKey);

     if (eAction == RESTART_SAVE) {
          if (it != m_Map.end()) {
               throw ErrDuplicateEntry(MBDYN_EXCEPT_ARGS);
          }

          m_Map.insert(std::make_pair(oKey, strData));
     } else {
          ASSERT(eAction == RESTART_RESTORE);

          if (it != m_Map.end()) {
               strData = it->second;
               return true;
          } else {
               // Do not alter strData
          }
     }

     return false;
}
