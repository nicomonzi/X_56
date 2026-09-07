/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2023
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cstdio>
#include <cstring>
#include "mbstrbuf.h"
//#include <stdio.h>

std::ostream&
mbstrbuf::stats(std::ostream& out)
{
	return out << "len=" << buf.capacity() << "; cursor=" << buf.size() << std::endl;
}


void
mbstrbuf::make_room(unsigned newlen)
{
	buf.reserve(buf.size() + newlen);
}

void
mbstrbuf::return_cursor(unsigned newcursor)
{
	buf.resize(newcursor);
}


void
mbstrbuf::print_str(const char *str)
{
	buf += str;
}

void
mbstrbuf::print_double(const char *fmt, double d)
{
	int dlen = snprintf(NULL, 0, fmt, d);
	if (dlen <= 0) {
		return;
	}

	std::string::size_type oldlen = buf.size();
	buf.resize(oldlen + dlen + 1);
	snprintf(&buf[oldlen], dlen + 1, fmt, d);
	buf.resize(oldlen + dlen);
}

const char *
mbstrbuf::get_buf(void) const
{
	return buf.c_str();
}

unsigned
mbstrbuf::get_len(void) const
{
	return unsigned(buf.size());
}

std::ostream&
operator << (std::ostream& out, const mbstrbuf& buf)
{
	return out << buf.buf;
}

#ifdef MAIN

int
main(void)
{
	mbstrbuf buf(10);
	buf.stats(std::cout) << std::endl;

	buf.print_double("%e", 10.2);

	std::cout << buf << std::endl;
	buf.stats(std::cout) << std::endl;

	buf.print_str(",");
	buf.print_double("%e", 12345678.9);

	std::cout << buf << std::endl;
	buf.stats(std::cout) << std::endl;

	buf.print_str(",");
	buf.print_double("%e", 12345678.9);

	std::cout << buf << std::endl;
	buf.stats(std::cout) << std::endl;

	buf.print_str(",");
	buf.print_double("%e", 12345678.9);

	std::cout << buf << std::endl;
	buf.stats(std::cout) << std::endl;

	return 0;
}

#endif // MAIN
