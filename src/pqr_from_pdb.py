"""
make the data from the pdb and itp
"""

import sys
import pandas as pd

import logger
import itp_to_df
import pdb_to_df
import parse_charmm_data


class PdbToPqr:
    """
    preapre the file with positions and charges and radii
    """

    def __init__(self,
                 log: logger.logging.Logger
                 ) -> None:
        self.initiate(log)

    def initiate(self,
                 log: logger.logging.Logger
                 ) -> None:
        """get all the infos"""
        itp = itp_to_df.Itp('APT_COR.itp').atoms
        pdb = pdb_to_df.Pdb(sys.argv[1], log).pdb_df
        charmm = parse_charmm_data.ParseData('CHARMM.DAT', log).radius_df
        self.get_charges(pdb, itp)

    def get_charges(self,
                    pdb: pd.DataFrame,
                    itp: pd.DataFrame
                    ) -> None:
        """get charges of the atoms in the pdb file"""
        atoms: list[str] = pdb['atom_name']
        charges: dict[str, float] = {}
        self.set_aptes_charges(pdb, itp)
        for item in atoms:
            charges[item] = itp[itp['atomname'] == item]['charge'].values

    def set_aptes_charges(self,
                          pdb: pd.DataFrame,
                          itp: pd.DataFrame
                          ) -> None:
        """get charges for different aptes longs"""
        itp_df: pd.DataFrame = itp[itp['resname'] == 'APT']
        apt_df: pd.DataFrame = pdb[pdb['residue_name'] == 'APT']
        pro_charges: pd.DataFrame
        unpro_charges: pd.DataFrame
        pro_charges, unpro_charges = self.get_aptes_q(itp_df)
        dfs_by_residue: dict[int, pd.DataFrame] = {}
        for residue_number, group in apt_df.groupby('residue_number'):
            if len(group) == 13:
                print("IT'S PROTONATED")
            elif len(group) == 12:
                print("IT'S UNPROTONATED")
            dfs_by_residue[residue_number] = group

    @staticmethod
    def get_aptes_q(itp_df: pd.DataFrame
                    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """return charges for pro- and unprotonated aptes"""
        pro_found: bool = False
        unpro_found: bool = False
        for _, group in itp_df.groupby('resnr'):
            if len(group) == 13 and not pro_found:
                pro_charges: pd.DataFrame = \
                    group[['atomname', 'charge']].copy()
                pro_found = True
            elif len(group) == 12 and not unpro_found:
                unpro_charges: pd.DataFrame = \
                    group[['atomname', 'charge']].copy()
                unpro_found = True
            if pro_found and unpro_found:
                break
        return pro_charges, unpro_charges


if __name__ == '__main__':
    PdbToPqr(log=logger.setup_logger('pdb2pqr.log'))
