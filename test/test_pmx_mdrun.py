import unittest
from pmxNCMC.pmx_mdrun import PMX_MDRUN_RE
from pathlib import Path
import os
import tempfile
import logging
import shutil

base=Path(__file__).parent

class MyTestCase(unittest.TestCase):
    def setUp(self):
        self.MDRUN="mpirun -np 2 --bind-to none gmx_mpi mdrun"
        self.GROMPP="gmx grompp"
        self.env=os.environ.copy()
        self.min_output=True

    def test_PMX_MDRUN_RE_init(self):
        print("Test New Start")
        os.chdir(base/'2-pentene/mdrun/test/01-new-sim')
        top = Path("../../../topol.top")
        csv = Path("./md.csv")
        if csv.is_file():
            csv.unlink()
        mdp_folder = Path("../../01-trans/rep_999/mdp/")
        folder_start = Path("000000")
        with tempfile.TemporaryDirectory(prefix="pmxRE_test_") as tmp_folder:
            tmp_folder=Path(tmp_folder)
            mdrun = PMX_MDRUN_RE(top, csv, mdp_folder, folder_start, self.MDRUN, self.GROMPP,
                                 tmp_folder, self.env, self.min_output)
            self.assertTrue(mdrun.safety_check())
            self.assertTrue(mdrun.safe_flag)
            mdrun.run_eq_mdrun()
            mdrun.run_ti_grompp()
            mdrun.run_ti_mdrun()
        csv.unlink()
        for i in range(2):
            for f in (folder_start/str(i)).iterdir():
                if f.name != "eq.tpr" and f.is_file():
                    f.unlink()

    def test_PMX_MDRUN_RE_restart(self):
        print("Test Restart")
        os.chdir(base / '2-pentene/mdrun/test/02-append')
        top = Path("../../../topol.top")
        shutil.copy("md_tmp.csv","md.csv")
        csv = Path("./md.csv")
        mdp_folder = Path("../../01-trans/rep_999/mdp/")
        folder_start = Path("000199")
        p_new = Path("000200")
        if p_new.is_dir():
            shutil.rmtree(p_new)
        with tempfile.TemporaryDirectory(prefix="pmxRE_test_") as tmp_folder:
            tmp_folder=Path(tmp_folder)
            mdrun = PMX_MDRUN_RE(top, csv, mdp_folder, folder_start, self.MDRUN, self.GROMPP,
                                 tmp_folder, self.env, self.min_output)
            self.assertTrue(mdrun.safety_check())
            self.assertTrue(mdrun.safe_flag)
            mdrun.run_eq_grompp()
            mdrun.run_eq_mdrun()
            mdrun.run_ti_grompp()
            mdrun.run_ti_mdrun()
        shutil.rmtree(p_new)

    def test_PMX_MDRUN_RE_new_cycle(self):
        print("Test New Cycle")
        os.chdir(base/'2-pentene/mdrun/test/01-new-sim')
        top = Path("../../../topol.top")
        csv = Path("./md.csv")
        if csv.is_file():
            csv.unlink()
        mdp_folder = Path("../../01-trans/rep_999/mdp/")
        folder_start = Path("000000")
        with tempfile.TemporaryDirectory(prefix="pmxRE_test_") as tmp_folder:
            tmp_folder = Path(tmp_folder)
            mdrun = PMX_MDRUN_RE(top, csv, mdp_folder, folder_start, self.MDRUN, self.GROMPP,
                                 tmp_folder, self.env, self.min_output)
            self.assertTrue(mdrun.safety_check())
            mdrun.log_sim_settings()
            self.assertTrue(mdrun.cycle(4, 0.5))
        csv.unlink()
        for i in range(2):
            for f in (folder_start/str(i)).iterdir():
                if f.name != "eq.tpr" and f.is_file():
                    f.unlink()
        shutil.rmtree("000001")
        shutil.rmtree("000002")
        shutil.rmtree("000003")



if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    unittest.main()
