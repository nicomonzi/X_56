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

/*
 AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
        Copyright (C) 2022(-2024) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

    module load: "libmodule-mcp_test1.la";

    user defined: 0, mcp test1,
    m1, k1, q1_0, qdot1_0, m2, k2, q2_0, qdot2_0;

    user defined: 1, mcp test2,
    m1, k1, q1_0, qdot1_0, string, "-m1 * g * sin(omega * Time)",
    output, yes;

    user defined: 2, mcp test2,
    m1, k1, q1_0, qdot1_0, string, "-m1 * g * cos(omega * Time)",
    output, yes;

    user defined: 3, mcp test2,
    m1, k1, q1_0, qdot1_0, string, "m1 * g * sin(omega * Time)",
    output, yes;

    user defined: 4, mcp test2,
    m1, k1, q1_0, qdot1_0, string, "m1 * g * cos(omega * Time)",
    output, yes;

    user defined: 5, mcp test2,
    m1, k1, q1_0, qdot1_0, string, "m1 * g * cos(omega * Time)^2",
    output, yes;

    user defined: 6, mcp test2,
    m1, k1, q1_0, qdot1_0, string, "m1 * g * sin(omega * Time)^2",
    output, yes;

    user defined: 7, mcp test2,
    m1, k1, q1_0, qdot1_0, string, "-m1 * g * cos(omega * Time)^2",
    output, yes;

    user defined: 8, mcp test2,
    m1, k1, q1_0, qdot1_0, string, "-m1 * g * sin(omega * Time)^2",
    output, yes;

    user defined: 9, mcp test2,
    m2, k2, q2_0, qdot2_0, string, "-m2 * g * sin(omega * Time)",
    output, yes;

    user defined: 10, mcp test2,
    m2, k2, q2_0, qdot2_0, string, "-m2 * g * cos(omega * Time)",
    output, yes;

    user defined: 11, mcp test2,
    m2, k2, q2_0, qdot2_0, string, "m2 * g * sin(omega * Time)",
    output, yes;

    user defined: 12, mcp test2,
    m2, k2, q2_0, qdot2_0, string, "m2 * g * cos(omega * Time)",
    output, yes;

    user defined: 13, mcp test2,
    m2, k2, q2_0, qdot2_0, string, "m2 * g * cos(omega * Time)^2",
    output, yes;

    user defined: 14, mcp test2,
    m2, k2, q2_0, qdot2_0, string, "m2 * g * sin(omega * Time)^2",
    output, yes;

    user defined: 15, mcp test2,
    m2, k2, q2_0, qdot2_0, string, "-m2 * g * cos(omega * Time)^2",
    output, yes;

    user defined: 16, mcp test2,
    m2, k2, q2_0, qdot2_0, string, "-m2 * g * sin(omega * Time)^2",
    output, yes;

    genel: 2, clamp,
           1, abstract, algebraic, 0.;