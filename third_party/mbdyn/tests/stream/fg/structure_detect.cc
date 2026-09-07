/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
 *
 * Pierangelo Masarati	<masarati@aero.polimi.it>
 * Paolo Mantegazza	<mantegazza@aero.polimi.it>
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

#include <iostream>
#include <cstddef>
#include <cstdlib>
#include <string>
#include <cstring>
#include <typeinfo>
#ifdef HAVE_CXXABI_H
#include <cxxabi.h>
#endif

// -I flightgear-$(FGVERSION)/src/ -I simgear-$(FGVERSION)/

// controls structure
#include "Network/net_ctrls.hxx"

// flight data model structure
#include "Network/net_fdm.hxx"

std::string
mbdyn_demangle(const char *name)
{
	std::string ret;

#ifdef HAVE_CXXABI_H
	const char* demangled_name = name;
	int status = -1;
	char* res = abi::__cxa_demangle(name, NULL, NULL, &status);
	if (status == 0) {
		demangled_name = res;
	}
	ret = demangled_name;
#else
	ret = name;
#endif

	return ret;
}

std::string
mbdyn_demangle(const std::type_info& t)
{
	return mbdyn_demangle(t.name());
}

#define MYFG_ENUM(type,name) \
        std::cout << #name " = " << type::name
#define MYFG_OFFSETOF(type,datum,name) \
        std::cout << #name ": type, size, offset = " << mbdyn_demangle(typeid(datum.name)) << ", " << sizeof(datum.name) << ", " << offsetof(type, name)

static int
detect_ctrls(void)
{
	FGNetCtrls ctrls;

	std::cout << "Controls System" << std::endl;
	std::cout << std::endl;
	std::cout << "sizeof(" << mbdyn_demangle(typeid(ctrls)) << ") = " << sizeof(ctrls) << std::endl;
	std::cout << std::endl;

	std::cout << "FG_NET_CTRLS_VERSION = " << FG_NET_CTRLS_VERSION << std::endl;
	std::cout << std::endl;

	int rc = 0;

#if FGVER == 020400

#if 0 // 2.4.0
// net_ctrls.hxx -- defines a common net I/O interface to the flight
//                  sim controls
//
// Written by Curtis Olson - http://www.flightgear.org/~curt
// Started July 2001.
//
// This file is in the Public Domain, and comes with no warranty.
//
// $Id: structure_detect.cc,v 1.1 2017/02/01 15:31:27 zanoni Exp $


#ifndef _NET_CTRLS_HXX
#define _NET_CTRLS_HXX

#include <simgear/misc/stdint.hxx>

// NOTE: this file defines an external interface structure.  Due to
// variability between platforms and architectures, we only used fixed
// length types here.  Specifically, integer types can vary in length.
// I am not aware of any platforms that don't use 4 bytes for float
// and 8 bytes for double.

//    !!! IMPORTANT !!!
/* There is some space reserved in the protocol for future use.
 * When adding a new type, add it just before the "reserved" definition
 * and subtract the size of this new type from the RESERVED_SPACE definition
 * (1 for (u)int32_t or float and 2 for double).
 *	
 * This way the protocol will be forward and backward compatible until
 * RESERVED_SPACE becomes zero.
 */

#define RESERVED_SPACE 25
const uint32_t FG_NET_CTRLS_VERSION = 27;


// Define a structure containing the control parameters

class FGNetCtrls {

public:

    enum {
        FG_MAX_ENGINES = 4,
        FG_MAX_WHEELS = 16,
        FG_MAX_TANKS = 8
    };

    uint32_t version;		         // increment when data values change

    // Aero controls
    double aileron;		         // -1 ... 1
    double elevator;		         // -1 ... 1
    double rudder;		         // -1 ... 1
    double aileron_trim;	         // -1 ... 1
    double elevator_trim;	         // -1 ... 1
    double rudder_trim;		         // -1 ... 1
    double flaps;		         //  0 ... 1
    double spoilers;
    double speedbrake;

    // Aero control faults
    uint32_t flaps_power;                 // true = power available
    uint32_t flap_motor_ok;

    // Engine controls
    uint32_t num_engines;		 // number of valid engines
    uint32_t master_bat[FG_MAX_ENGINES];
    uint32_t master_alt[FG_MAX_ENGINES];
    uint32_t magnetos[FG_MAX_ENGINES];
    uint32_t starter_power[FG_MAX_ENGINES];// true = starter power
    double throttle[FG_MAX_ENGINES];     //  0 ... 1
    double mixture[FG_MAX_ENGINES];      //  0 ... 1
    double condition[FG_MAX_ENGINES];    //  0 ... 1
    uint32_t fuel_pump_power[FG_MAX_ENGINES];// true = on
    double prop_advance[FG_MAX_ENGINES]; //  0 ... 1
    uint32_t feed_tank_to[4];
    uint32_t reverse[4];


    // Engine faults
    uint32_t engine_ok[FG_MAX_ENGINES];
    uint32_t mag_left_ok[FG_MAX_ENGINES];
    uint32_t mag_right_ok[FG_MAX_ENGINES];
    uint32_t spark_plugs_ok[FG_MAX_ENGINES];  // false = fouled plugs
    uint32_t oil_press_status[FG_MAX_ENGINES];// 0 = normal, 1 = low, 2 = full fail
    uint32_t fuel_pump_ok[FG_MAX_ENGINES];

    // Fuel management
    uint32_t num_tanks;                      // number of valid tanks
    uint32_t fuel_selector[FG_MAX_TANKS];    // false = off, true = on
    uint32_t xfer_pump[5];                   // specifies transfer from array
                                             // value tank to tank specified by
                                             // int value
    uint32_t cross_feed;                     // false = off, true = on

    // Brake controls
    double brake_left;
    double brake_right;
    double copilot_brake_left;
    double copilot_brake_right;
    double brake_parking;
    
    // Landing Gear
    uint32_t gear_handle; // true=gear handle down; false= gear handle up

    // Switches
    uint32_t master_avionics;
    
        // nav and Comm
    double	comm_1;
    double	comm_2;
    double	nav_1;
    double	nav_2;

    // wind and turbulance
    double wind_speed_kt;
    double wind_dir_deg;
    double turbulence_norm;

    // temp and pressure
    double temp_c;
    double press_inhg;

    // other information about environment
    double hground;		         // ground elevation (meters)
    double magvar;		         // local magnetic variation in degs.

    // hazards
    uint32_t icing;                      // icing status could me much
                                         // more complex but I'm
                                         // starting simple here.

    // simulation control
    uint32_t speedup;		         // integer speedup multiplier
    uint32_t freeze;		         // 0=normal
				         // 0x01=master
				         // 0x02=position
				         // 0x04=fuel

    // --- New since FlightGear 0.9.10 (FG_NET_CTRLS_VERSION = 27)

    // --- Add new variables just before this line.

    uint32_t reserved[RESERVED_SPACE];	 // 100 bytes reserved for future use.
};


#endif // _NET_CTRLS_HXX
#endif

/*
    enum {
        FG_MAX_ENGINES = 4,
        FG_MAX_WHEELS = 16,
        FG_MAX_TANKS = 8
    };
*/

	MYFG_ENUM(FGNetCtrls, FG_MAX_ENGINES) << std::endl;
	MYFG_ENUM(FGNetCtrls, FG_MAX_WHEELS) << std::endl;
	MYFG_ENUM(FGNetCtrls, FG_MAX_TANKS) << std::endl;
	std::cout << std::endl;

//    uint32_t version;		         // increment when data values change
	MYFG_OFFSETOF(FGNetCtrls, ctrls, version) << std::endl;
	std::cout << std::endl;

/*
    // Aero controls
    double aileron;		         // -1 ... 1
    double elevator;		         // -1 ... 1
    double rudder;		         // -1 ... 1
    double aileron_trim;	         // -1 ... 1
    double elevator_trim;	         // -1 ... 1
    double rudder_trim;		         // -1 ... 1
    double flaps;		         //  0 ... 1
    double spoilers;
    double speedbrake;
*/
	std::cout << "Aero controls" << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, aileron) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, elevator) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, rudder) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, aileron_trim) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, elevator_trim) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, rudder_trim) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, flaps) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, spoilers) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, speedbrake) << std::endl;
	std::cout << std::endl;

/*
    // Aero control faults
    uint32_t flaps_power;                 // true = power available
    uint32_t flap_motor_ok;
*/
	std::cout << "Aero control faults" << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, flaps_power) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, flap_motor_ok) << std::endl;
	std::cout << std::endl;

/*
    // Engine controls
    uint32_t num_engines;		 // number of valid engines
    uint32_t master_bat[FG_MAX_ENGINES];
    uint32_t master_alt[FG_MAX_ENGINES];
    uint32_t magnetos[FG_MAX_ENGINES];
    uint32_t starter_power[FG_MAX_ENGINES];// true = starter power
    double throttle[FG_MAX_ENGINES];     //  0 ... 1
    double mixture[FG_MAX_ENGINES];      //  0 ... 1
    double condition[FG_MAX_ENGINES];    //  0 ... 1
    uint32_t fuel_pump_power[FG_MAX_ENGINES];// true = on
    double prop_advance[FG_MAX_ENGINES]; //  0 ... 1
    uint32_t feed_tank_to[4];
    uint32_t reverse[4];
*/
	std::cout << "Engine controls" << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, num_engines) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, master_bat) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, master_alt) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, magnetos) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, starter_power) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, throttle) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, mixture) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, condition) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, fuel_pump_power) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, prop_advance) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, feed_tank_to) << " (note: 4 instead of FG_MAX_ENGINES)" << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, reverse) << " (note: 4 instead of FG_MAX_ENGINES)" << std::endl;
	std::cout << std::endl;

/*
    // Engine faults
    uint32_t engine_ok[FG_MAX_ENGINES];
    uint32_t mag_left_ok[FG_MAX_ENGINES];
    uint32_t mag_right_ok[FG_MAX_ENGINES];
    uint32_t spark_plugs_ok[FG_MAX_ENGINES];  // false = fouled plugs
    uint32_t oil_press_status[FG_MAX_ENGINES];// 0 = normal, 1 = low, 2 = full fail
    uint32_t fuel_pump_ok[FG_MAX_ENGINES];
*/
	std::cout << "Engine faults" << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, engine_ok) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, mag_left_ok) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, mag_right_ok) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, spark_plugs_ok) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, oil_press_status) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, fuel_pump_ok) << std::endl;
	std::cout << std::endl;

/*
    // Fuel management
    uint32_t num_tanks;                      // number of valid tanks
    uint32_t fuel_selector[FG_MAX_TANKS];    // false = off, true = on
    uint32_t xfer_pump[5];                   // specifies transfer from array
                                             // value tank to tank specified by
                                             // int value
    uint32_t cross_feed;                     // false = off, true = on
*/
	std::cout << "Engine faults" << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, num_tanks) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, fuel_selector) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, xfer_pump) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, cross_feed) << std::endl;
	std::cout << std::endl;

/*
    // Brake controls
    double brake_left;
    double brake_right;
    double copilot_brake_left;
    double copilot_brake_right;
    double brake_parking;
*/
	std::cout << "Brake controls" << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, brake_left) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, brake_right) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, copilot_brake_left) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, copilot_brake_right) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, brake_parking) << std::endl;
	std::cout << std::endl;

/*
    // Landing Gear
    uint32_t gear_handle; // true=gear handle down; false= gear handle up
*/
	std::cout << "Landing Gear" << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, gear_handle) << std::endl;
	std::cout << std::endl;

/*
    // Switches
    uint32_t master_avionics;
*/
	std::cout << "Switches" << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, master_avionics) << std::endl;
	std::cout << std::endl;

/*    
        // nav and Comm
    double	comm_1;
    double	comm_2;
    double	nav_1;
    double	nav_2;
*/
	std::cout << "nav and Comm" << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, comm_1) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, comm_2) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, nav_1) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, nav_2) << std::endl;
	std::cout << std::endl;

/*
    // wind and turbulance
    double wind_speed_kt;
    double wind_dir_deg;
    double turbulence_norm;
*/
	std::cout << "wind and turbulance" << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, wind_speed_kt) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, wind_dir_deg) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, turbulence_norm) << std::endl;
	std::cout << std::endl;

/*
    // temp and pressure
    double temp_c;
    double press_inhg;
*/
	std::cout << "temp and pressure" << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, temp_c) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, press_inhg) << std::endl;
	std::cout << std::endl;

/*
    // other information about environment
    double hground;		         // ground elevation (meters)
    double magvar;		         // local magnetic variation in degs.
*/
	std::cout << "other information about environment" << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, hground) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, magvar) << std::endl;
	std::cout << std::endl;

/*
    // hazards
    uint32_t icing;                      // icing status could me much
                                         // more complex but I'm
                                         // starting simple here.
*/
	std::cout << "hazards" << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, icing) << std::endl;
	std::cout << std::endl;

/*
    // simulation control
    uint32_t speedup;		         // integer speedup multiplier
    uint32_t freeze;		         // 0=normal
				         // 0x01=master
				         // 0x02=position
				         // 0x04=fuel
*/
	std::cout << "simulation control" << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, speedup) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, freeze) << std::endl;
	std::cout << std::endl;

/*
    // --- New since FlightGear 0.9.10 (FG_NET_CTRLS_VERSION = 27)

    // --- Add new variables just before this line.

    uint32_t reserved[RESERVED_SPACE];	 // 100 bytes reserved for future use.
*/
	std::cout << "reserved (New since FlightGear 0.9.10)" << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, reserved) << std::endl;
	std::cout << std::endl;

#elif (FGVER == 030400) || (FGVER == 030200)

#if 0 // 3.4.0
// net_ctrls.hxx -- defines a common net I/O interface to the flight
//                  sim controls
//
// Written by Curtis Olson - http://www.flightgear.org/~curt
// Started July 2001.
//
// This file is in the Public Domain, and comes with no warranty.
//
// $Id: structure_detect.cc,v 1.1 2017/02/01 15:31:27 zanoni Exp $


#ifndef _NET_CTRLS_HXX
#define _NET_CTRLS_HXX

#include <simgear/misc/stdint.hxx>

// NOTE: this file defines an external interface structure.  Due to
// variability between platforms and architectures, we only used fixed
// length types here.  Specifically, integer types can vary in length.
// I am not aware of any platforms that don't use 4 bytes for float
// and 8 bytes for double.

//    !!! IMPORTANT !!!
/* There is some space reserved in the protocol for future use.
 * When adding a new type, add it just before the "reserved" definition
 * and subtract the size of this new type from the RESERVED_SPACE definition
 * (1 for (u)int32_t or float and 2 for double).
 *	
 * This way the protocol will be forward and backward compatible until
 * RESERVED_SPACE becomes zero.
 */

#define RESERVED_SPACE 25
const uint32_t FG_NET_CTRLS_VERSION = 27;


// Define a structure containing the control parameters

class FGNetCtrls {

public:

    enum {
        FG_MAX_ENGINES = 4,
        FG_MAX_WHEELS = 16,
        FG_MAX_TANKS = 8
    };

    uint32_t version;		         // increment when data values change

    // Aero controls
    double aileron;		         // -1 ... 1
    double elevator;		         // -1 ... 1
    double rudder;		         // -1 ... 1
    double aileron_trim;	         // -1 ... 1
    double elevator_trim;	         // -1 ... 1
    double rudder_trim;		         // -1 ... 1
    double flaps;		         //  0 ... 1
    double spoilers;
    double speedbrake;

    // Aero control faults
    uint32_t flaps_power;                 // true = power available
    uint32_t flap_motor_ok;

    // Engine controls
    uint32_t num_engines;		 // number of valid engines
    uint32_t master_bat[FG_MAX_ENGINES];
    uint32_t master_alt[FG_MAX_ENGINES];
    uint32_t magnetos[FG_MAX_ENGINES];
    uint32_t starter_power[FG_MAX_ENGINES];// true = starter power
    double throttle[FG_MAX_ENGINES];     //  0 ... 1
    double mixture[FG_MAX_ENGINES];      //  0 ... 1
    double condition[FG_MAX_ENGINES];    //  0 ... 1
    uint32_t fuel_pump_power[FG_MAX_ENGINES];// true = on
    double prop_advance[FG_MAX_ENGINES]; //  0 ... 1
    uint32_t feed_tank_to[4];
    uint32_t reverse[4];


    // Engine faults
    uint32_t engine_ok[FG_MAX_ENGINES];
    uint32_t mag_left_ok[FG_MAX_ENGINES];
    uint32_t mag_right_ok[FG_MAX_ENGINES];
    uint32_t spark_plugs_ok[FG_MAX_ENGINES];  // false = fouled plugs
    uint32_t oil_press_status[FG_MAX_ENGINES];// 0 = normal, 1 = low, 2 = full fail
    uint32_t fuel_pump_ok[FG_MAX_ENGINES];

    // Fuel management
    uint32_t num_tanks;                      // number of valid tanks
    uint32_t fuel_selector[FG_MAX_TANKS];    // false = off, true = on
    uint32_t xfer_pump[5];                   // specifies transfer from array
                                             // value tank to tank specified by
                                             // int value
    uint32_t cross_feed;                     // false = off, true = on

    // Brake controls
    double brake_left;
    double brake_right;
    double copilot_brake_left;
    double copilot_brake_right;
    double brake_parking;
    
    // Landing Gear
    uint32_t gear_handle; // true=gear handle down; false= gear handle up

    // Switches
    uint32_t master_avionics;
    
        // nav and Comm
    double	comm_1;
    double	comm_2;
    double	nav_1;
    double	nav_2;

    // wind and turbulance
    double wind_speed_kt;
    double wind_dir_deg;
    double turbulence_norm;

    // temp and pressure
    double temp_c;
    double press_inhg;

    // other information about environment
    double hground;		         // ground elevation (meters)
    double magvar;		         // local magnetic variation in degs.

    // hazards
    uint32_t icing;                      // icing status could me much
                                         // more complex but I'm
                                         // starting simple here.

    // simulation control
    uint32_t speedup;		         // integer speedup multiplier
    uint32_t freeze;		         // 0=normal
				         // 0x01=master
				         // 0x02=position
				         // 0x04=fuel

    // --- New since FlightGear 0.9.10 (FG_NET_CTRLS_VERSION = 27)

    // --- Add new variables just before this line.

    uint32_t reserved[RESERVED_SPACE];	 // 100 bytes reserved for future use.
};


#endif // _NET_CTRLS_HXX


#endif

/*
    enum {
        FG_MAX_ENGINES = 4,
        FG_MAX_WHEELS = 16,
        FG_MAX_TANKS = 8
    };
*/

	MYFG_ENUM(FGNetCtrls, FG_MAX_ENGINES) << std::endl;
	MYFG_ENUM(FGNetCtrls, FG_MAX_WHEELS) << std::endl;
	MYFG_ENUM(FGNetCtrls, FG_MAX_TANKS) << std::endl;
	std::cout << std::endl;

//    uint32_t version;		         // increment when data values change
	MYFG_OFFSETOF(FGNetCtrls, ctrls, version) << std::endl;
	std::cout << std::endl;

/*
    // Aero controls
    double aileron;		         // -1 ... 1
    double elevator;		         // -1 ... 1
    double rudder;		         // -1 ... 1
    double aileron_trim;	         // -1 ... 1
    double elevator_trim;	         // -1 ... 1
    double rudder_trim;		         // -1 ... 1
    double flaps;		         //  0 ... 1
    double spoilers;
    double speedbrake;
*/
	std::cout << "Aero controls" << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, aileron) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, elevator) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, rudder) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, aileron_trim) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, elevator_trim) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, rudder_trim) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, flaps) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, spoilers) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, speedbrake) << std::endl;
	std::cout << std::endl;

/*
    // Aero control faults
    uint32_t flaps_power;                 // true = power available
    uint32_t flap_motor_ok;
*/
	std::cout << "Aero control faults" << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, flaps_power) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, flap_motor_ok) << std::endl;
	std::cout << std::endl;

/*
    // Engine controls
    uint32_t num_engines;		 // number of valid engines
    uint32_t master_bat[FG_MAX_ENGINES];
    uint32_t master_alt[FG_MAX_ENGINES];
    uint32_t magnetos[FG_MAX_ENGINES];
    uint32_t starter_power[FG_MAX_ENGINES];// true = starter power
    double throttle[FG_MAX_ENGINES];     //  0 ... 1
    double mixture[FG_MAX_ENGINES];      //  0 ... 1
    double condition[FG_MAX_ENGINES];    //  0 ... 1
    uint32_t fuel_pump_power[FG_MAX_ENGINES];// true = on
    double prop_advance[FG_MAX_ENGINES]; //  0 ... 1
    uint32_t feed_tank_to[4];
    uint32_t reverse[4];
*/
	std::cout << "Engine controls" << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, num_engines) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, master_bat) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, master_alt) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, magnetos) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, starter_power) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, throttle) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, mixture) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, condition) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, fuel_pump_power) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, prop_advance) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, feed_tank_to) << " (note: 4 instead of FG_MAX_ENGINES)" << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, reverse) << " (note: 4 instead of FG_MAX_ENGINES)" << std::endl;
	std::cout << std::endl;

/*
    // Engine faults
    uint32_t engine_ok[FG_MAX_ENGINES];
    uint32_t mag_left_ok[FG_MAX_ENGINES];
    uint32_t mag_right_ok[FG_MAX_ENGINES];
    uint32_t spark_plugs_ok[FG_MAX_ENGINES];  // false = fouled plugs
    uint32_t oil_press_status[FG_MAX_ENGINES];// 0 = normal, 1 = low, 2 = full fail
    uint32_t fuel_pump_ok[FG_MAX_ENGINES];
*/
	std::cout << "Engine faults" << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, engine_ok) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, mag_left_ok) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, mag_right_ok) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, spark_plugs_ok) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, oil_press_status) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, fuel_pump_ok) << std::endl;
	std::cout << std::endl;

/*
    // Fuel management
    uint32_t num_tanks;                      // number of valid tanks
    uint32_t fuel_selector[FG_MAX_TANKS];    // false = off, true = on
    uint32_t xfer_pump[5];                   // specifies transfer from array
                                             // value tank to tank specified by
                                             // int value
    uint32_t cross_feed;                     // false = off, true = on
*/
	std::cout << "Engine faults" << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, num_tanks) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, fuel_selector) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, xfer_pump) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, cross_feed) << std::endl;
	std::cout << std::endl;

/*
    // Brake controls
    double brake_left;
    double brake_right;
    double copilot_brake_left;
    double copilot_brake_right;
    double brake_parking;
*/
	std::cout << "Brake controls" << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, brake_left) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, brake_right) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, copilot_brake_left) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, copilot_brake_right) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, brake_parking) << std::endl;
	std::cout << std::endl;

/*
    // Landing Gear
    uint32_t gear_handle; // true=gear handle down; false= gear handle up
*/
	std::cout << "Landing Gear" << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, gear_handle) << std::endl;
	std::cout << std::endl;

/*
    // Switches
    uint32_t master_avionics;
*/
	std::cout << "Switches" << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, master_avionics) << std::endl;
	std::cout << std::endl;

/*    
        // nav and Comm
    double	comm_1;
    double	comm_2;
    double	nav_1;
    double	nav_2;
*/
	std::cout << "nav and Comm" << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, comm_1) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, comm_2) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, nav_1) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, nav_2) << std::endl;
	std::cout << std::endl;

/*
    // wind and turbulance
    double wind_speed_kt;
    double wind_dir_deg;
    double turbulence_norm;
*/
	std::cout << "wind and turbulance" << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, wind_speed_kt) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, wind_dir_deg) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, turbulence_norm) << std::endl;
	std::cout << std::endl;

/*
    // temp and pressure
    double temp_c;
    double press_inhg;
*/
	std::cout << "temp and pressure" << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, temp_c) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, press_inhg) << std::endl;
	std::cout << std::endl;

/*
    // other information about environment
    double hground;		         // ground elevation (meters)
    double magvar;		         // local magnetic variation in degs.
*/
	std::cout << "other information about environment" << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, hground) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, magvar) << std::endl;
	std::cout << std::endl;

/*
    // hazards
    uint32_t icing;                      // icing status could me much
                                         // more complex but I'm
                                         // starting simple here.
*/
	std::cout << "hazards" << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, icing) << std::endl;
	std::cout << std::endl;

/*
    // simulation control
    uint32_t speedup;		         // integer speedup multiplier
    uint32_t freeze;		         // 0=normal
				         // 0x01=master
				         // 0x02=position
				         // 0x04=fuel
*/
	std::cout << "simulation control" << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, speedup) << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, freeze) << std::endl;
	std::cout << std::endl;

/*
    // --- New since FlightGear 0.9.10 (FG_NET_CTRLS_VERSION = 27)

    // --- Add new variables just before this line.

    uint32_t reserved[RESERVED_SPACE];	 // 100 bytes reserved for future use.
*/
	std::cout << "reserved (New since FlightGear 0.9.10)" << std::endl;
	MYFG_OFFSETOF(FGNetCtrls, ctrls, reserved) << std::endl;
	std::cout << std::endl;

#else
	//std::cerr << "unsupported FGVERSION = " FGVERSION << std::endl;
	rc = 1;
#endif

#include "ctrls.hcc"

	return rc;
}

static int
detect_fdm(void)
{

	FGNetFDM fdm;

	std::cout << "Flight Dynamics Model" << std::endl;
	std::cout << std::endl;
	std::cout << "sizeof(" << mbdyn_demangle(typeid(fdm)) << ") = " << sizeof(fdm) << std::endl;
	std::cout << std::endl;

	std::cout << "FG_NET_FDM_VERSION = " << FG_NET_FDM_VERSION << std::endl;
	std::cout << std::endl;

	int rc = 0;

#if FGVER == 020400

#if 0 // 2.4.0
// net_fdm.hxx -- defines a common net I/O interface to the flight
//                dynamics model
//
// Written by Curtis Olson - http://www.flightgear.org/~curt
// Started September 2001.
//
// This file is in the Public Domain, and comes with no warranty.
//
// $Id: structure_detect.cc,v 1.1 2017/02/01 15:31:27 zanoni Exp $


#ifndef _NET_FDM_HXX
#define _NET_FDM_HXX


#include <time.h> // time_t
#include <simgear/misc/stdint.hxx>

// NOTE: this file defines an external interface structure.  Due to
// variability between platforms and architectures, we only used fixed
// length types here.  Specifically, integer types can vary in length.
// I am not aware of any platforms that don't use 4 bytes for float
// and 8 bytes for double.

const uint32_t FG_NET_FDM_VERSION = 24;


// Define a structure containing the top level flight dynamics model
// parameters

class FGNetFDM {

public:

    enum {
        FG_MAX_ENGINES = 4,
        FG_MAX_WHEELS = 3,
        FG_MAX_TANKS = 4
    };

    uint32_t version;		// increment when data values change
    uint32_t padding;		// padding

    // Positions
    double longitude;		// geodetic (radians)
    double latitude;		// geodetic (radians)
    double altitude;		// above sea level (meters)
    float agl;			// above ground level (meters)
    float phi;			// roll (radians)
    float theta;		// pitch (radians)
    float psi;			// yaw or true heading (radians)
    float alpha;                // angle of attack (radians)
    float beta;                 // side slip angle (radians)

    // Velocities
    float phidot;		// roll rate (radians/sec)
    float thetadot;		// pitch rate (radians/sec)
    float psidot;		// yaw rate (radians/sec)
    float vcas;		        // calibrated airspeed
    float climb_rate;		// feet per second
    float v_north;              // north velocity in local/body frame, fps
    float v_east;               // east velocity in local/body frame, fps
    float v_down;               // down/vertical velocity in local/body frame, fps
    float v_wind_body_north;    // north velocity in local/body frame
                                // relative to local airmass, fps
    float v_wind_body_east;     // east velocity in local/body frame
                                // relative to local airmass, fps
    float v_wind_body_down;     // down/vertical velocity in local/body
                                // frame relative to local airmass, fps

    // Accelerations
    float A_X_pilot;		// X accel in body frame ft/sec^2
    float A_Y_pilot;		// Y accel in body frame ft/sec^2
    float A_Z_pilot;		// Z accel in body frame ft/sec^2

    // Stall
    float stall_warning;        // 0.0 - 1.0 indicating the amount of stall
    float slip_deg;		// slip ball deflection

    // Pressure
    
    // Engine status
    uint32_t num_engines;	     // Number of valid engines
    uint32_t eng_state[FG_MAX_ENGINES];// Engine state (off, cranking, running)
    float rpm[FG_MAX_ENGINES];	     // Engine RPM rev/min
    float fuel_flow[FG_MAX_ENGINES]; // Fuel flow gallons/hr
    float fuel_px[FG_MAX_ENGINES];   // Fuel pressure psi
    float egt[FG_MAX_ENGINES];	     // Exhuast gas temp deg F
    float cht[FG_MAX_ENGINES];	     // Cylinder head temp deg F
    float mp_osi[FG_MAX_ENGINES];    // Manifold pressure
    float tit[FG_MAX_ENGINES];	     // Turbine Inlet Temperature
    float oil_temp[FG_MAX_ENGINES];  // Oil temp deg F
    float oil_px[FG_MAX_ENGINES];    // Oil pressure psi

    // Consumables
    uint32_t num_tanks;		// Max number of fuel tanks
    float fuel_quantity[FG_MAX_TANKS];

    // Gear status
    uint32_t num_wheels;
    uint32_t wow[FG_MAX_WHEELS];
    float gear_pos[FG_MAX_WHEELS];
    float gear_steer[FG_MAX_WHEELS];
    float gear_compression[FG_MAX_WHEELS];

    // Environment
    uint32_t cur_time;           // current unix time
                                 // FIXME: make this uint64_t before 2038
    int32_t warp;                // offset in seconds to unix time
    float visibility;            // visibility in meters (for env. effects)

    // Control surface positions (normalized values)
    float elevator;
    float elevator_trim_tab;
    float left_flap;
    float right_flap;
    float left_aileron;
    float right_aileron;
    float rudder;
    float nose_wheel;
    float speedbrake;
    float spoilers;
};


#endif // _NET_FDM_HXX
#endif

/*
    enum {
        FG_MAX_ENGINES = 4,
        FG_MAX_WHEELS = 3,
        FG_MAX_TANKS = 4
    };
*/
	MYFG_ENUM(FGNetFDM, FG_MAX_ENGINES) << std::endl;
	MYFG_ENUM(FGNetFDM, FG_MAX_WHEELS) << std::endl;
	MYFG_ENUM(FGNetFDM, FG_MAX_TANKS) << std::endl;
	std::cout << std::endl;

/*
    uint32_t version;		// increment when data values change
    uint32_t padding;		// padding
*/
	MYFG_OFFSETOF(FGNetFDM, fdm, version) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, padding) << " (unused)" << std::endl;
	std::cout << std::endl;

/*
    // Positions
    double longitude;		// geodetic (radians)
    double latitude;		// geodetic (radians)
    double altitude;		// above sea level (meters)
    float agl;			// above ground level (meters)
    float phi;			// roll (radians)
    float theta;		// pitch (radians)
    float psi;			// yaw or true heading (radians)
    float alpha;                // angle of attack (radians)
    float beta;                 // side slip angle (radians)
*/
	std::cout << "Positions" << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, longitude) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, latitude) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, altitude) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, agl) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, phi) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, theta) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, psi) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, alpha) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, beta) << std::endl;
	std::cout << std::endl;

/*
    // Velocities
    float phidot;		// roll rate (radians/sec)
    float thetadot;		// pitch rate (radians/sec)
    float psidot;		// yaw rate (radians/sec)
    float vcas;		        // calibrated airspeed
    float climb_rate;		// feet per second
    float v_north;              // north velocity in local/body frame, fps
    float v_east;               // east velocity in local/body frame, fps
    float v_down;               // down/vertical velocity in local/body frame, fps
    float v_wind_body_north;    // north velocity in local/body frame
                                // relative to local airmass, fps
    float v_wind_body_east;     // east velocity in local/body frame
                                // relative to local airmass, fps
    float v_wind_body_down;     // down/vertical velocity in local/body
                                // frame relative to local airmass, fps
*/

	std::cout << "Velocities" << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, phidot) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, thetadot) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, psidot) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, vcas) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, climb_rate) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, v_north) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, v_east) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, v_down) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, v_wind_body_north) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, v_wind_body_east) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, v_wind_body_down) << std::endl;
	std::cout << std::endl;

/*
    // Accelerations
    float A_X_pilot;		// X accel in body frame ft/sec^2
    float A_Y_pilot;		// Y accel in body frame ft/sec^2
    float A_Z_pilot;		// Z accel in body frame ft/sec^2
*/
	std::cout << "Accelerations" << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, A_X_pilot) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, A_Y_pilot) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, A_Z_pilot) << std::endl;
	std::cout << std::endl;

/*
    // Stall
    float stall_warning;        // 0.0 - 1.0 indicating the amount of stall
    float slip_deg;		// slip ball deflection
*/
	std::cout << "Stall" << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, stall_warning) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, slip_deg) << std::endl;
	std::cout << std::endl;

    // Pressure

/*    
    // Engine status
    uint32_t num_engines;	     // Number of valid engines
    uint32_t eng_state[FG_MAX_ENGINES];// Engine state (off, cranking, running)
    float rpm[FG_MAX_ENGINES];	     // Engine RPM rev/min
    float fuel_flow[FG_MAX_ENGINES]; // Fuel flow gallons/hr
    float fuel_px[FG_MAX_ENGINES];   // Fuel pressure psi
    float egt[FG_MAX_ENGINES];	     // Exhuast gas temp deg F
    float cht[FG_MAX_ENGINES];	     // Cylinder head temp deg F
    float mp_osi[FG_MAX_ENGINES];    // Manifold pressure
    float tit[FG_MAX_ENGINES];	     // Turbine Inlet Temperature
    float oil_temp[FG_MAX_ENGINES];  // Oil temp deg F
    float oil_px[FG_MAX_ENGINES];    // Oil pressure psi
*/
	std::cout << "Engine status" << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, num_engines) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, eng_state) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, rpm) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, fuel_flow) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, fuel_px) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, egt) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, cht) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, mp_osi) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, tit) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, oil_temp) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, oil_px) << std::endl;
	std::cout << std::endl;

/*
    // Consumables
    uint32_t num_tanks;		// Max number of fuel tanks
    float fuel_quantity[FG_MAX_TANKS];
*/
	std::cout << "Consumables" << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, num_tanks) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, fuel_quantity) << std::endl;
	std::cout << std::endl;

/*
    // Gear status
    uint32_t num_wheels;
    uint32_t wow[FG_MAX_WHEELS];
    float gear_pos[FG_MAX_WHEELS];
    float gear_steer[FG_MAX_WHEELS];
    float gear_compression[FG_MAX_WHEELS];
*/
	std::cout << "Gear status" << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, num_wheels) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, wow) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, gear_pos) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, gear_steer) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, gear_compression) << std::endl;
	std::cout << std::endl;

/*
    // Environment
    uint32_t cur_time;           // current unix time
                                 // FIXME: make this uint64_t before 2038
    int32_t warp;                // offset in seconds to unix time
    float visibility;            // visibility in meters (for env. effects)
*/
	std::cout << "Environment" << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, cur_time) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, warp) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, visibility) << std::endl;
	std::cout << std::endl;

/*
    // Control surface positions (normalized values)
    float elevator;
    float elevator_trim_tab;
    float left_flap;
    float right_flap;
    float left_aileron;
    float right_aileron;
    float rudder;
    float nose_wheel;
    float speedbrake;
    float spoilers;
*/
	std::cout << "Control surface positions" << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, elevator) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, elevator_trim_tab) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, left_flap) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, right_flap) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, left_aileron) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, right_aileron) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, rudder) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, nose_wheel) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, speedbrake) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, spoilers) << std::endl;
	std::cout << std::endl;

#elif (FGVER == 030400) || (FGVER == 030200)

#if 0 // 3.4.0
// net_fdm.hxx -- defines a common net I/O interface to the flight
//                dynamics model
//
// Written by Curtis Olson - http://www.flightgear.org/~curt
// Started September 2001.
//
// This file is in the Public Domain, and comes with no warranty.
//
// $Id: structure_detect.cc,v 1.1 2017/02/01 15:31:27 zanoni Exp $


#ifndef _NET_FDM_HXX
#define _NET_FDM_HXX


#include <time.h> // time_t
#include <simgear/misc/stdint.hxx>

// NOTE: this file defines an external interface structure.  Due to
// variability between platforms and architectures, we only used fixed
// length types here.  Specifically, integer types can vary in length.
// I am not aware of any platforms that don't use 4 bytes for float
// and 8 bytes for double.

const uint32_t FG_NET_FDM_VERSION = 24;


// Define a structure containing the top level flight dynamics model
// parameters

class FGNetFDM {

public:

    enum {
        FG_MAX_ENGINES = 4,
        FG_MAX_WHEELS = 3,
        FG_MAX_TANKS = 4
    };

    uint32_t version;		// increment when data values change
    uint32_t padding;		// padding

    // Positions
    double longitude;		// geodetic (radians)
    double latitude;		// geodetic (radians)
    double altitude;		// above sea level (meters)
    float agl;			// above ground level (meters)
    float phi;			// roll (radians)
    float theta;		// pitch (radians)
    float psi;			// yaw or true heading (radians)
    float alpha;                // angle of attack (radians)
    float beta;                 // side slip angle (radians)

    // Velocities
    float phidot;		// roll rate (radians/sec)
    float thetadot;		// pitch rate (radians/sec)
    float psidot;		// yaw rate (radians/sec)
    float vcas;		        // calibrated airspeed
    float climb_rate;		// feet per second
    float v_north;              // north velocity in local/body frame, fps
    float v_east;               // east velocity in local/body frame, fps
    float v_down;               // down/vertical velocity in local/body frame, fps
    float v_body_u;    // ECEF velocity in body frame
    float v_body_v;    // ECEF velocity in body frame 
    float v_body_w;    // ECEF velocity in body frame
    
    // Accelerations
    float A_X_pilot;		// X accel in body frame ft/sec^2
    float A_Y_pilot;		// Y accel in body frame ft/sec^2
    float A_Z_pilot;		// Z accel in body frame ft/sec^2

    // Stall
    float stall_warning;        // 0.0 - 1.0 indicating the amount of stall
    float slip_deg;		// slip ball deflection

    // Pressure
    
    // Engine status
    uint32_t num_engines;	     // Number of valid engines
    uint32_t eng_state[FG_MAX_ENGINES];// Engine state (off, cranking, running)
    float rpm[FG_MAX_ENGINES];	     // Engine RPM rev/min
    float fuel_flow[FG_MAX_ENGINES]; // Fuel flow gallons/hr
    float fuel_px[FG_MAX_ENGINES];   // Fuel pressure psi
    float egt[FG_MAX_ENGINES];	     // Exhuast gas temp deg F
    float cht[FG_MAX_ENGINES];	     // Cylinder head temp deg F
    float mp_osi[FG_MAX_ENGINES];    // Manifold pressure
    float tit[FG_MAX_ENGINES];	     // Turbine Inlet Temperature
    float oil_temp[FG_MAX_ENGINES];  // Oil temp deg F
    float oil_px[FG_MAX_ENGINES];    // Oil pressure psi

    // Consumables
    uint32_t num_tanks;		// Max number of fuel tanks
    float fuel_quantity[FG_MAX_TANKS];

    // Gear status
    uint32_t num_wheels;
    uint32_t wow[FG_MAX_WHEELS];
    float gear_pos[FG_MAX_WHEELS];
    float gear_steer[FG_MAX_WHEELS];
    float gear_compression[FG_MAX_WHEELS];

    // Environment
    uint32_t cur_time;           // current unix time
                                 // FIXME: make this uint64_t before 2038
    int32_t warp;                // offset in seconds to unix time
    float visibility;            // visibility in meters (for env. effects)

    // Control surface positions (normalized values)
    float elevator;
    float elevator_trim_tab;
    float left_flap;
    float right_flap;
    float left_aileron;
    float right_aileron;
    float rudder;
    float nose_wheel;
    float speedbrake;
    float spoilers;
};

#endif // _NET_FDM_HXX
#endif

/*
    enum {
        FG_MAX_ENGINES = 4,
        FG_MAX_WHEELS = 3,
        FG_MAX_TANKS = 4
    };
*/
	MYFG_ENUM(FGNetFDM, FG_MAX_ENGINES) << std::endl;
	MYFG_ENUM(FGNetFDM, FG_MAX_WHEELS) << std::endl;
	MYFG_ENUM(FGNetFDM, FG_MAX_TANKS) << std::endl;
	std::cout << std::endl;

/*
    uint32_t version;		// increment when data values change
    uint32_t padding;		// padding
*/
	MYFG_OFFSETOF(FGNetFDM, fdm, version) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, padding) << " (unused)" << std::endl;
	std::cout << std::endl;

/*
    // Positions
    double longitude;		// geodetic (radians)
    double latitude;		// geodetic (radians)
    double altitude;		// above sea level (meters)
    float agl;			// above ground level (meters)
    float phi;			// roll (radians)
    float theta;		// pitch (radians)
    float psi;			// yaw or true heading (radians)
    float alpha;                // angle of attack (radians)
    float beta;                 // side slip angle (radians)
*/
	std::cout << "Positions" << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, longitude) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, latitude) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, altitude) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, agl) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, phi) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, theta) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, psi) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, alpha) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, beta) << std::endl;
	std::cout << std::endl;

/*
    // Velocities
    float phidot;		// roll rate (radians/sec)
    float thetadot;		// pitch rate (radians/sec)
    float psidot;		// yaw rate (radians/sec)
    float vcas;		        // calibrated airspeed
    float climb_rate;		// feet per second
    float v_north;              // north velocity in local/body frame, fps
    float v_east;               // east velocity in local/body frame, fps
    float v_down;               // down/vertical velocity in local/body frame, fps
    float v_body_u;    // ECEF velocity in body frame
    float v_body_v;    // ECEF velocity in body frame 
    float v_body_w;    // ECEF velocity in body frame
*/

	std::cout << "Velocities" << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, phidot) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, thetadot) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, psidot) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, vcas) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, climb_rate) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, v_north) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, v_east) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, v_down) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, v_body_u) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, v_body_v) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, v_body_w) << std::endl;
	std::cout << std::endl;

/*
    // Accelerations
    float A_X_pilot;		// X accel in body frame ft/sec^2
    float A_Y_pilot;		// Y accel in body frame ft/sec^2
    float A_Z_pilot;		// Z accel in body frame ft/sec^2
*/
	std::cout << "Accelerations" << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, A_X_pilot) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, A_Y_pilot) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, A_Z_pilot) << std::endl;
	std::cout << std::endl;

/*
    // Stall
    float stall_warning;        // 0.0 - 1.0 indicating the amount of stall
    float slip_deg;		// slip ball deflection
*/
	std::cout << "Stall" << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, stall_warning) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, slip_deg) << std::endl;
	std::cout << std::endl;

    // Pressure

/*    
    // Engine status
    uint32_t num_engines;	     // Number of valid engines
    uint32_t eng_state[FG_MAX_ENGINES];// Engine state (off, cranking, running)
    float rpm[FG_MAX_ENGINES];	     // Engine RPM rev/min
    float fuel_flow[FG_MAX_ENGINES]; // Fuel flow gallons/hr
    float fuel_px[FG_MAX_ENGINES];   // Fuel pressure psi
    float egt[FG_MAX_ENGINES];	     // Exhuast gas temp deg F
    float cht[FG_MAX_ENGINES];	     // Cylinder head temp deg F
    float mp_osi[FG_MAX_ENGINES];    // Manifold pressure
    float tit[FG_MAX_ENGINES];	     // Turbine Inlet Temperature
    float oil_temp[FG_MAX_ENGINES];  // Oil temp deg F
    float oil_px[FG_MAX_ENGINES];    // Oil pressure psi
*/
	std::cout << "Engine status" << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, num_engines) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, eng_state) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, rpm) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, fuel_flow) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, fuel_px) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, egt) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, cht) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, mp_osi) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, tit) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, oil_temp) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, oil_px) << std::endl;
	std::cout << std::endl;

/*
    // Consumables
    uint32_t num_tanks;		// Max number of fuel tanks
    float fuel_quantity[FG_MAX_TANKS];
*/
	std::cout << "Consumables" << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, num_tanks) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, fuel_quantity) << std::endl;
	std::cout << std::endl;

/*
    // Gear status
    uint32_t num_wheels;
    uint32_t wow[FG_MAX_WHEELS];
    float gear_pos[FG_MAX_WHEELS];
    float gear_steer[FG_MAX_WHEELS];
    float gear_compression[FG_MAX_WHEELS];
*/
	std::cout << "Gear status" << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, num_wheels) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, wow) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, gear_pos) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, gear_steer) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, gear_compression) << std::endl;
	std::cout << std::endl;

/*
    // Environment
    uint32_t cur_time;           // current unix time
                                 // FIXME: make this uint64_t before 2038
    int32_t warp;                // offset in seconds to unix time
    float visibility;            // visibility in meters (for env. effects)
*/
	std::cout << "Environment" << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, cur_time) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, warp) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, visibility) << std::endl;
	std::cout << std::endl;

/*
    // Control surface positions (normalized values)
    float elevator;
    float elevator_trim_tab;
    float left_flap;
    float right_flap;
    float left_aileron;
    float right_aileron;
    float rudder;
    float nose_wheel;
    float speedbrake;
    float spoilers;
*/
	std::cout << "Control surface positions" << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, elevator) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, elevator_trim_tab) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, left_flap) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, right_flap) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, left_aileron) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, right_aileron) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, rudder) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, nose_wheel) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, speedbrake) << std::endl;
	MYFG_OFFSETOF(FGNetFDM, fdm, spoilers) << std::endl;
	std::cout << std::endl;

#else
	//std::cerr << "unsupported FGVERSION = " FGVERSION << std::endl;
	rc = 1;
#endif

#include "fdm.hcc"

	return rc;
}

int
main(int argc, char *argv[])
{
	std::cout << "flightgear-" FGVERSION << std::endl;

	std::cout << std::endl;
	detect_ctrls();

	std::cout << std::endl;
	detect_fdm();

	return 0;
}
