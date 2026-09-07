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

#include "mbconfig.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>

extern char *optarg;
extern int optind;

#define BUFFER_SIZE 256 
#define MAX_NAME_LENGTH 256 

/*
@param input_name: Nastran output basename (for .op2 and .mat files)
@param output_name: Output .fem filename (if empty, uses input_name)
@param is_modal_displacement: initial modal displacements flag (1: enabled, 0: disabled)
@param is_modal_velocity: initial modal velocities flag (1: enabled, 0: disabled)
@param mass_direction: index of "mass" component to be used (-1: check consistency, 1: x, 2: y, 3: z)  
@return: exit code 
@brief: femgen main routine that converts Nastran .op2 and .mat files into an
        ASCII finite element model file with extension .fem suitable for MBDyn
        modal joint.  
*/
extern int femgen(const char input_name[MAX_NAME_LENGTH],
                  const char output_name[MAX_NAME_LENGTH],
                  int32_t is_modal_displacement, 
                  int32_t is_modal_velocity, 
                  int32_t mass_direction);

static void
usage(FILE *outf, int rc)
{
	fprintf(outf,
"usage: femgen <nastran_name> [-o <output_fem>] [-d] [-m {cxyz}] [-h]\n"
"\n"
"  <nastran_name>   Nastran output basename (reads <name>.op2 and <name>.mat)\n"
"  -o <output_fem>  Output .fem filename (default: <nastran_name>.fem)\n"
"  -d               Disable initial modal displacements/velocities\n"
"  -m <idx>         Index of 'mass' component ('x', 'y', 'z'; 'c' to check)\n"
"  -h               Show this help message\n"
"\n"
"Examples:\n"
"  femgen blade              # reads blade.op2, blade.mat -> blade.fem\n"
"  femgen blade -o modal     # reads blade.op2, blade.mat -> modal.fem\n"
"  femgen blade.op2          # reads blade.op2, blade.mat -> blade.fem\n"
		);
	exit(rc);
}

int main(int argc, char *argv[])
{
	char input_name[MAX_NAME_LENGTH] = { '\0' };
	char output_name[MAX_NAME_LENGTH] = { '\0' };

	int32_t is_modal_displacement = 1;
	int32_t is_modal_velocity = 1; 
	int32_t mass_direction = 0;
	
	for (;;) {
		int opt = getopt(argc, argv, "dhm:o:");
		if (opt == -1) {
			break;
		}

		switch (opt) {
		case 'd':
			is_modal_displacement = 0;
			is_modal_velocity = 0;
			break;

		case '?':
		case 'h':
			usage(stdout, EXIT_SUCCESS);

		case 'm':
			if (optarg[1] != '\0') {
				usage(stderr, EXIT_FAILURE);
			}

			switch (optarg[0]) {
			case 'c':
				mass_direction = -1;
				break;

			case 'x':
				mass_direction = 1;
				break;

			case 'y':
				mass_direction = 2;
				break;

			case 'z':
				mass_direction = 3;
				break;

			default:
				fprintf(stderr, "femgen: unhandled parameter '%c' for option '-m'\n", optarg[0]);
				usage(stderr, EXIT_FAILURE);
			}
			break;

		case 'o': {
			size_t len = strlen(optarg);
			if (len >= sizeof(output_name)) {
				fprintf(stderr, "femgen: output file name '%s' too long\n", optarg);
				exit(EXIT_FAILURE);
			}
			strcpy(output_name, optarg);
			} break;

		default:
			fprintf(stderr, "femgen: unhandled option '%c'\n", opt);
			usage(stderr, EXIT_FAILURE);
		}
	}

	/* Get the required positional argument: Nastran input name */
	if (optind < argc) {
		size_t len = strlen(argv[optind]);
		if (len >= sizeof(input_name)) {
			fprintf(stderr, "femgen: input file name '%s' too long\n", argv[optind]);
			exit(EXIT_FAILURE);
		}
		strcpy(input_name, argv[optind]);
		optind++;
	}

	/* Check if input name was provided */
	if (input_name[0] == '\0') {
		fprintf(stderr, "femgen: missing required argument <nastran_name>\n\n");
		usage(stderr, EXIT_FAILURE);
	}

	if (optind < argc) {
		fprintf(stderr, "femgen: extra arguments ignored\n");
	}

	return femgen(input_name, output_name, is_modal_displacement, is_modal_velocity, mass_direction);
}

