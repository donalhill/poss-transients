"""Palomar's POSS-I plate log (VizieR VI/25/nposs): for each numbered red plate, the observing night (evening date of
the double date) and the local start time. Writes data/palomar_log_vi25.csv."""
from pathlib import Path
from astroquery.vizier import Vizier

ROOT = Path(__file__).resolve().parents[1]
t = Vizier(row_limit=-1, columns=["**"]).get_catalogs("VI/25/nposs")[0].to_pandas()
t = t.rename(columns={"POSS": "plate_number", "Obs": "evening_date", "ObsE": "red_start_local", "EHA": "red_start_ha", "Eexp": "red_exposure_min"})
t["plate_number"] = t.plate_number.astype(str).str.strip()
t = t.rename(columns={"_RA.icrs": "ra_j2000", "_DE.icrs": "dec_j2000"})
t[["recno", "plate_number", "evening_date", "red_start_local", "red_start_ha", "red_exposure_min", "ra_j2000", "dec_j2000"]].to_csv(ROOT / "data/palomar_log_vi25.csv", index=False)
print(len(t), "rows written")
