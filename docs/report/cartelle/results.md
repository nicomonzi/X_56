# Rapporto: `results`

## Ruolo e stato

Cartella globale molto piccola: otto file, circa 1 MB. Contiene cinque PNG e tre
CSV di confronto, fra cui grafici Nastran rigido vs MBDyn e risultati globali.
Non contiene codice né README che attribuisca con precisione ogni figura a un
workflow.

## Obsolescenza e azioni

I risultati non sono necessariamente obsoleti, ma **la collocazione è
obsoleta/ambigua**. Per ogni file va identificato lo script produttore tramite
nomi, cronologia Git e colonne CSV; poi spostarlo sotto
`results/processed/<studio>` o `GRAFICI/<studio>`, accompagnato da input hash,
commit e comando di generazione. Fino ad allora conservare tutto: lo spazio è
trascurabile e la provenienza non è ancora formalizzata.

