"""preCICE 2.x adapter between the MBDyn nodal socket and DUST."""

from __future__ import annotations

import numpy as np
import precice


WRITE_FIELDS = ("Position", "Velocity", "Rotation", "AngularVelocity")
READ_FIELDS = ("Force", "Moment")


class CouplingAdapter:
    """Drive the serial-implicit MBDyn/DUST coupling until final time."""

    def __init__(self, mbdyn, config_file: str) -> None:
        self.mbdyn = mbdyn
        self.interface = precice.Interface("MBDyn", config_file, 0, 1)
        self.mesh_id = self.interface.get_mesh_id("MBDynNodes")
        self.vertex_ids = self.interface.set_mesh_vertices(
            self.mesh_id, self.mbdyn.reference_nodes
        )
        self.data_ids = {
            name: self.interface.get_data_id(name, self.mesh_id)
            for name in (*WRITE_FIELDS, *READ_FIELDS)
        }
        self.precice_dt = self.interface.initialize()

        # Some coupling configurations request the initial kinematics before
        # the first time window. Write every structural field, not only one.
        write_initial = precice.action_write_initial_data()
        if self.interface.is_action_required(write_initial):
            self._write_kinematics()
            self.interface.mark_action_fulfilled(write_initial)
        self.interface.initialize_data()

    def _write_kinematics(self) -> None:
        for name in WRITE_FIELDS:
            self.interface.write_block_vector_data(
                self.data_ids[name], self.vertex_ids, self.mbdyn.data[name]
            )

    def _read_loads(self) -> None:
        for name in READ_FIELDS:
            self.mbdyn.data[name] = self.interface.read_block_vector_data(
                self.data_ids[name], self.vertex_ids
            )

    def run(self, structural_dt: float) -> None:
        write_checkpoint = precice.action_write_iteration_checkpoint()
        read_checkpoint = precice.action_read_iteration_checkpoint()
        saved = None
        simulation_time = 0.0
        iteration = 0

        try:
            while self.interface.is_coupling_ongoing():
                iteration += 1
                if self.interface.is_action_required(write_checkpoint):
                    saved = {
                        name: self.mbdyn.data[name].copy()
                        for name in WRITE_FIELDS
                    }
                    self.interface.mark_action_fulfilled(write_checkpoint)

                # Apply the last DUST loads to the structural trial iteration.
                if self.mbdyn.send_loads(converged=False):
                    break
                if self.mbdyn.receive_kinematics():
                    break

                self._write_kinematics()
                dt = min(structural_dt, self.precice_dt)
                self.precice_dt = self.interface.advance(dt)
                self._read_loads()

                if self.interface.is_action_required(read_checkpoint):
                    if saved is None:
                        raise RuntimeError("preCICE requested a missing checkpoint")
                    for name, values in saved.items():
                        self.mbdyn.data[name] = values
                    self.interface.mark_action_fulfilled(read_checkpoint)
                else:
                    # Tell MBDyn to accept the converged trial and start a new step.
                    if self.mbdyn.send_loads(converged=True):
                        break
                    if self.mbdyn.receive_kinematics():
                        break
                    simulation_time += dt
                    print(
                        f"Coupled time {simulation_time:8.3f} s "
                        f"(last implicit iteration {iteration})",
                        flush=True,
                    )
                    iteration = 0
        finally:
            self.interface.finalize()
            self.mbdyn.close()
