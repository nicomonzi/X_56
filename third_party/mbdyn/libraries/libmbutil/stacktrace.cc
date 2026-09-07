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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "stacktrace.h"

#include <iostream>
#include <csignal>

#if defined(HAVE_STACKTRACE) && (__cplusplus >= 202302L) // we should be able to use stacktrace
#include <stacktrace>

void segfault_handler(int signal) {
    std::cerr << "Error: Signal " << signal << " caught. Backtrace:\n";
    std::cerr << std::stacktrace::current() << '\n';
    std::_Exit(signal); // Exit immediately to avoid looping if a segfault happens here
}

#elif defined(HAVE_EXECINFO_H) // use execinfo

#include <execinfo.h>
#include <unistd.h>
void segfault_handler(int signal) {
    void* array[400];
    size_t size = backtrace(array, 400);

    std::cerr << "Error: Signal " << signal << " caught. Backtrace:\n";
    backtrace_symbols_fd(array, size, STDERR_FILENO);
    _exit(signal);
}

#else // do nothing
void segfault_handler(int signal) {
     return;
}
#endif

void set_stacktrace_callback() {
    std::signal(SIGABRT, segfault_handler);
    // std::signal(SIGFPE, segfault_handler);
    std::signal(SIGSEGV, segfault_handler);
    std::signal(SIGTERM, segfault_handler);
}
