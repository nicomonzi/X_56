## $Header$
###############################################################################
#
# module-marsh: MBDyn as a MARSH (https://marsh-sim.github.io/) Flight
# Dynamics Module node.
#
###############################################################################

MODULE_DEPENDENCIES=mavlink_producers.lo mavlink_consumers.lo

# vendored MAVLink headers use packed structs -> GCC's
# -Waddress-of-packed-member fires constantly on generated code
MODULE_INCLUDE=-Wno-address-of-packed-member
