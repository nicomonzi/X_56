/* $Header$ */
/*
 * This library comes with MBDyn (C), a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2023
 *
 * Pierangelo Masarati  <pierangelo.masarati@polimi.it>
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

/* Codice relativo alla classe File_name, per l'handling di nomi di files
 * in ambiente DOS e UNIX */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "filename.h"
#include "myassert.h"
#include "mynewmem.h"


FileName::FileName(const std::string sFName, int iExtSepNum)
{
   	if (sFName != "") {
      		iInit(sFName, iExtSepNum);
   	}
}

FileName::~FileName(void)
{
	NO_OP;
}

int
FileName::iInit(const  std::string sFName, int iExtSepNum)
{

   	if (sFName == "") {
      		return 0;
   	}

	std::string s(sFName);
	std::string::size_type pos = std::string::npos;

	/* the extension separator is looked for
	 * after the last directory separator, if any */
	std::string::size_type lim = s.find_last_of(DIR_SEP);
	if (lim == std::string::npos) {
		lim = 0;
	}

	if (iExtSepNum > 0) {
		/* iExtSepNum-th extension separator, from the beginning */
		int iCnt = 0;
		for (std::string::size_type i = lim; i < s.size(); i++) {
			if (s[i] == EXT_SEP && ++iCnt == iExtSepNum) {
				pos = i;
				break;
			}
		}

	} else if (iExtSepNum < 0) {
		/* |iExtSepNum|-th extension separator, from the end;
		 * a leading extension separator (e.g. ".profile")
		 * is not treated as such */
		int iCnt = 0;
		for (std::string::size_type i = s.size(); i-- > lim + 1; ) {
			if (s[i] == EXT_SEP && ++iCnt == -iExtSepNum) {
				pos = i;
				break;
			}
		}
	}

	if (pos == std::string::npos) {
		sBase = s;
		sExt.clear();

	} else {
		sBase = s.substr(0, pos);
		sExt = s.substr(pos);
	}

	return int(sBase.size());
}

const std::string
FileName::_sPutExt(std::string sEName)
{
   	if (sEName == "") {
      		sEName = sExt;
   	}

	sName = sBase;
   	if (sEName[0] != '\0') {
      		if (sEName[0] != EXT_SEP) {
			sName += EXT_SEP;
      		}
      		sName += sEName;
   	}

   	return sName;
}

const std::string
FileName::sGet(void) const
{
	std::string tmps("");
   	return const_cast<FileName *>(this)->_sPutExt(tmps);
}

int
is_abs_path(const char *const p)
{
	ASSERT(p != 0);

#ifdef _WIN32
	if ((strstr(p, ":\\") != 0) || (strstr(p, ":/") != 0)) {
		return 1;
	}
#else // ! _WIN32
	if (p[0] == '/') {
		return 1;
	}
#endif // ! _WIN32

	return 0;
}

