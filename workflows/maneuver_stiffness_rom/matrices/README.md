# Contratto delle matrici

Ogni sottocartella di caso deve contenere:

- `kh_fixed_basis.csv`: matrice 6 x 6, senza intestazione, separata da virgole;
- `metadata.json`: metadati che certificano base, modi, unita e fattore di carico.

Esempio di `metadata.json`:

```json
{
  "load_factor": 1.6,
  "basis_id": "x56_mbdyn_modal_60_phi0_modes_7_12",
  "modes": [7, 8, 9, 10, 11, 12],
  "fem_sha256": "63e23dc991d2f458d7665412abda8ad34ddbb1654822bfc0867e6512a1df252b",
  "units": "MBDyn FEM modal normalization; K_h in s^-2",
  "projection": "Phi0.T @ K_tangent(nz) @ Phi0"
}
```

`kh_fixed_basis.csv` deve essere realmente proiettata sulla stessa `Phi0` usata
dal `.fem` MBDyn. La diagonale `KHH` prodotta dall'ALTER SOL103 standard non e
intercambiabile: appartiene alla base propria prestressata del singolo run.
Anche scala e normalizzazione devono coincidere: nel `.fem` corrente la massa
modale e identita e, per esempio, `K_77 = (2*pi*3.217134 Hz)^2`.
