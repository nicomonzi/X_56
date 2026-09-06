"""Small MBDyn socket interface used by the local preCICE adapter."""

from __future__ import annotations

import numpy as np
from mbc_py_interface import mbcNodal


class MBDynInterface:
    """Exchange nodal kinematics and loads with an external structural force."""

    def __init__(self, socket_path: str, reference_file: str) -> None:
        self.reference_nodes = np.loadtxt(reference_file, dtype=float)
        if self.reference_nodes.ndim == 1:
            self.reference_nodes = self.reference_nodes.reshape(1, 3)
        self.node_count = int(self.reference_nodes.shape[0])

        # orientation-vector rotation and accelerations must match the MBDyn card.
        self.nodal = mbcNodal(
            socket_path, "", 0, -1, 0, 1, 0, self.node_count,
            0, 0x100, 1
        )
        self.nodal.negotiate()
        self.nodal.recv()

        self.data = {
            "Position": self._reshape(self.nodal.n_x),
            "Velocity": self._reshape(self.nodal.n_xp),
            "Rotation": self._reshape(self.nodal.n_theta),
            "AngularVelocity": self._reshape(self.nodal.n_omega),
            "Force": np.zeros((self.node_count, 3), dtype=float),
            "Moment": np.zeros((self.node_count, 3), dtype=float),
        }
        self.nodal.n_f[:] = 0.0
        self.nodal.n_m[:] = 0.0

    def _reshape(self, values) -> np.ndarray:
        return np.asarray(values, dtype=float).reshape(self.node_count, 3).copy()

    def receive_kinematics(self) -> bool:
        """Receive a trial MBDyn state; True means MBDyn ended the exchange."""
        if self.nodal.recv():
            return True
        self.data["Position"] = self._reshape(self.nodal.n_x)
        self.data["Velocity"] = self._reshape(self.nodal.n_xp)
        self.data["Rotation"] = self._reshape(self.nodal.n_theta)
        self.data["AngularVelocity"] = self._reshape(self.nodal.n_omega)
        return False

    def send_loads(self, converged: bool) -> bool:
        """Send DUST nodal loads and select trial or converged MBDyn advance."""
        self.nodal.n_f[:] = self.data["Force"].reshape(-1)
        self.nodal.n_m[:] = self.data["Moment"].reshape(-1)
        return bool(self.nodal.send(converged))

    def close(self) -> None:
        self.nodal.destroy()
