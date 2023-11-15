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
        itp: pd.DataFrame = itp_to_df.Itp('APT_COR.itp').atoms
        pdb: pd.DataFrame = pdb_to_df.Pdb(sys.argv[1], log).pdb_df
        charmm: pd.DataFrame = \
            parse_charmm_data.ParseData('CHARMM.DAT', log).radius_df
        self.get_charges(pdb, itp)

    def get_charges(self,
                    pdb: pd.DataFrame,
                    itp: pd.DataFrame
                    ) -> None:
        """get charges of the atoms in the pdb file"""
        atoms: list[str] = pdb['atom_name']
        charges: dict[str, float] = {}
        aptes_with_charges: pd.DataFrame = self.set_aptes_charges(pdb, itp)
        for item in atoms:
            charges[item] = itp[itp['atomname'] == item]['charge'].values

    def set_aptes_charges(self,
                          pdb: pd.DataFrame,
                          itp: pd.DataFrame
                          ) -> pd.DataFrame:
        """get charges for different aptes longs"""
        itp_df: pd.DataFrame = itp[itp['resname'] == 'APT']
        apt_df: pd.DataFrame = pdb[pdb['residue_name'] == 'APT']
        pro_charges: pd.DataFrame
        unpro_charges: pd.DataFrame
        pro_charges, unpro_charges = self.get_aptes_charges(itp_df)
        dfs_with_charges: list[pd.DataFrame] = []
        for _, group in apt_df.groupby('residue_number'):
            if len(group) == 13:
                dfs_with_charges.append(self.add_charge_column_to_df(
                    group, pro_charges))
            elif len(group) == 12:
                dfs_with_charges.append(self.add_charge_column_to_df(
                    group, unpro_charges))
        return pd.concat(dfs_with_charges, axis=0, ignore_index=True)

    @staticmethod
    def add_charge_column_to_df(main_df: pd.DataFrame,
                                charge_df: pd.DataFrame
                                ) -> pd.DataFrame:
        """add charge column to the dataframe"""
        tmp_df = main_df.copy()
        tmp_df = pd.merge(tmp_df,
                          charge_df[['atomname', 'charge']],
                          left_on='atom_name',
                          right_on='atomname',
                          how='left')
        return tmp_df.drop('atomname', axis=1)

    @staticmethod
    def get_aptes_charges(itp_df: pd.DataFrame
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
